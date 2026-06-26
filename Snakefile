# ==============================================================================
# guide_extraction — Main Snakefile Entry Point
# ==============================================================================
# scPerturb-seq Guide Extraction Workflow
#
# Pipeline (scprocess is an external pre-processing step; its H5 output is a
# direct input to this workflow):
#   guide_fasta → sgRNA_index ────────────────────────┐
#   gex_h5 (external) → extract_whitelist ────────────┤
#   sgRNA_fastq ─────────────────────────→ sgRNA_quant → merge → merged MEX
#
# Usage:
#   snakemake --cores 8          # local multi-core (default)
#   snakemake -n                 # dry-run preview
# ==============================================================================

import os
import sys
import glob as _pyglob
import yaml


# ---- Path constants ----
# Snakemake uses the Snakefile directory as the working directory
BASEDIR = os.getcwd()

# ---- Load global config ----
configfile: "config/config.yaml"

# ---- Load sample topology (GROUPS dict, available to all rules) ----
# v0.1.2: samples file path can be overridden via config key.
_samples_file = config.get("samples_file", os.path.join(BASEDIR, "config", "samples.yaml"))
if not os.path.isabs(_samples_file):
    _samples_file = os.path.join(BASEDIR, _samples_file)
with open(_samples_file, "r") as f:
    _samples_cfg = yaml.safe_load(f)

GROUPS = _samples_cfg.get("groups", {})
if not GROUPS:
    sys.exit(f"ERROR: No groups defined in {_samples_file}")


# ---- Resolve sgRNA FASTQ — expand glob patterns to concrete file lists ----
def _resolve_sgRNA_fastq(group_cfg: dict, read: str) -> list:
    """Expand sgRNA glob patterns from samples.yaml, returning sorted FASTQ file paths."""
    d = group_cfg["sgRNA_fastq_dir"]
    pat = group_cfg[f"sgRNA_{read}_pattern"]
    files = sorted(_pyglob.glob(os.path.join(d, pat)))
    if not files:
        print(f"WARNING: No sgRNA FASTQ found for {read} pattern '{pat}' in {d}",
              file=sys.stderr)
    return files


# Pre-resolve sgRNA FASTQ for all groups, inject into GROUPS dict
for gname, gcfg in GROUPS.items():
    gcfg["_sgRNA_r1"] = _resolve_sgRNA_fastq(gcfg, "r1")
    gcfg["_sgRNA_r2"] = _resolve_sgRNA_fastq(gcfg, "r2")
    # Pre-resolve pre-trimmed FASTQ only if preprocess.trimmed is enabled
    if config.get("preprocess", {}).get("trimmed", False):
        _trim_dir = os.path.join(gcfg["sgRNA_fastq_dir"], "trimmed")
        gcfg["_sgRNA_r1_trimmed"] = _resolve_sgRNA_fastq(
            {"sgRNA_fastq_dir": _trim_dir, "sgRNA_r1_pattern": gcfg["sgRNA_r1_pattern"]}, "r1"
        ) if os.path.isdir(_trim_dir) else []
        gcfg["_sgRNA_r2_trimmed"] = _resolve_sgRNA_fastq(
            {"sgRNA_fastq_dir": _trim_dir, "sgRNA_r2_pattern": gcfg["sgRNA_r2_pattern"]}, "r2"
        ) if os.path.isdir(_trim_dir) else []
    else:
        gcfg["_sgRNA_r1_trimmed"] = []
        gcfg["_sgRNA_r2_trimmed"] = []


# ---- Chemistry resolution (v0.1.2+) ----
# Load the single-source-of-truth chemistry specification and resolve
# tenx_chemistry → all downstream parameters (af_chemistry, whitelist,
# translation behaviour, geometry overrides, HAM chemistry, UMI length).
#
# Populates config["_chemistry"] for rules to consume.
# Backfills legacy config keys for backward compatibility with configs
# that lack a tenx_chemistry field (pre-v0.1.2).
_SPEC_PATH = os.path.join(BASEDIR, "config", "chemistry_spec.yaml")
if not os.path.exists(_SPEC_PATH):
    sys.exit(
        "ERROR: config/chemistry_spec.yaml not found. "
        "This file is required by sgprocess >= 0.1.2.")
with open(_SPEC_PATH, "r") as f:
    _CHEMISTRY_SPEC = yaml.safe_load(f)


def _resolve_chemistry(cfg):
    """Resolve tenx_chemistry → downstream parameters.

    Populates cfg["_chemistry"] with the resolved spec dict.
    Backfills legacy config keys so pre-v0.1.2 configs continue working.
    """
    chem = cfg.get("tenx_chemistry", None)
    if chem is None:
        # Legacy config — no tenx_chemistry field.
        # Rules fall back to their existing config paths unchanged.
        cfg["_chemistry"] = None
        return

    overrides = cfg.get("chemistry_overrides", {})

    if chem == "custom":
        spec = dict(cfg.get("custom_chemistry", {}))
        required = ["af_chemistry", "whitelist", "expected_ori",
                     "translation", "ham_chemistry", "umi_len"]
        for k in required:
            if k not in spec:
                raise ValueError(
                    f"custom_chemistry missing required key: {k}")
    else:
        spec = dict(_CHEMISTRY_SPEC.get(
            chem, _CHEMISTRY_SPEC.get("3v3", {})))
        if not spec:
            raise ValueError(f"Unknown tenx_chemistry: {chem}")
        spec.update(overrides)

    cfg["_chemistry"] = spec

    # ── Backfill legacy config paths ──

    # skip_translation
    if "skip_translation" not in cfg:
        cfg["skip_translation"] = not spec.get("translation", True)

    # simpleaf quant chemistry (use geometry_override if present)
    cfg.setdefault("simpleaf", {})
    cfg["simpleaf"].setdefault("quant", {})
    if "chemistry" not in cfg["simpleaf"]["quant"]:
        geom = spec.get("geometry_override")
        cfg["simpleaf"]["quant"]["chemistry"] = (
            geom if geom else spec["af_chemistry"])

    # hash_matcher chemistry
    cfg.setdefault("hash_matcher", {})
    if "chemistry" not in cfg["hash_matcher"]:
        cfg["hash_matcher"]["chemistry"] = spec.get(
            "ham_chemistry", "10xv3")


_resolve_chemistry(config)


# ---- Method selection ----
# Supported: "simpleaf" (default) | "hash_matcher"
METHOD = config.get("guide_extraction", {}).get("method", "simpleaf")
if METHOD not in ("simpleaf", "hash_matcher"):
    sys.exit(f"ERROR: Unknown guide_extraction.method='{METHOD}'. "
             f"Supported: simpleaf, hash_matcher")

# Per-method output directories (used by merge.smk)
if METHOD == "simpleaf":
    QUANT_OUT_SUBDIR = "simpleaf_quant/af_quant/alevin"
    QUANT_MTX_FILE   = "quants_mat.mtx"
    QUANT_ROWS_FILE  = "quants_mat_rows.txt"
    QUANT_COLS_FILE  = "quants_mat_cols.txt"
else:  # hash_matcher
    QUANT_OUT_SUBDIR = "guide_quant/matrix"
    QUANT_MTX_FILE   = "matrix.mtx.gz"
    QUANT_ROWS_FILE  = "barcodes.tsv.gz"
    QUANT_COLS_FILE  = "features.tsv.gz"


# ---- Final target ----
rule all:
    input:
        os.path.join(config["out_dir"], "merged", "merged_matrix.mtx.gz"),
        os.path.join(config["out_dir"], "merged", "merged_barcodes.tsv.gz"),
        os.path.join(config["out_dir"], "merged", "merged_features.tsv.gz"),
    

# ---- Import rule modules ----
include: "rules/reference.smk"
include: "rules/whitelist.smk"

if METHOD == "simpleaf":
    include: "rules/quant.smk"
else:
    include: "rules/guide_quant.smk"

include: "rules/merge.smk"
