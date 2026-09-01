# ==============================================================================
# scprocess-perturb — Snakefile
# ==============================================================================
# Perturb-seq guide extraction (+ optional assignment). GEX matrices come from
# an external pre-processing step and are a direct input.
#   guide_fasta → sgRNA_index ────────────────────────┐
#   gex_h5 (external) → extract_whitelist ────────────┤
#   sgRNA_fastq ──────────→ sgRNA_quant → merge → merged MEX → assignment
#
# Usage:
#   snakemake --cores 8          # run
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
_groups_file = config.get("groups_file", os.path.join(BASEDIR, "config", "groups.yaml"))
if not os.path.isabs(_groups_file):
    _groups_file = os.path.join(BASEDIR, _groups_file)
with open(_groups_file, "r") as f:
    _groups_cfg = yaml.safe_load(f)

GROUPS = _groups_cfg.get("groups", {})
if not GROUPS:
    sys.exit(f"ERROR: No groups defined in {_groups_file}")


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


# ---- Chemistry resolution ----
# Resolve tenx_chemistry → downstream parameters (af_chemistry, whitelist,
# translation, geometry override, HAM chemistry, UMI length) into
# config["_chemistry"] for rules to consume.
_SPEC_PATH = os.path.join(BASEDIR, "config", "chemistry_spec.yaml")
if not os.path.exists(_SPEC_PATH):
    sys.exit("ERROR: config/chemistry_spec.yaml not found.")
with open(_SPEC_PATH, "r") as f:
    _CHEMISTRY_SPEC = yaml.safe_load(f)


def _resolve_chemistry(cfg):
    """Resolve tenx_chemistry into cfg["_chemistry"] and backfill legacy keys."""
    chem = cfg.get("tenx_chemistry", None)
    if chem is None:
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

    # Backfill legacy config keys
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


# ---- Reference auto-derivation ----
# Reference paths default to {out_dir}/refs/; users need only provide guide_csv.
_refs = config.setdefault("references", {})
_ref_dir = os.path.join(config["out_dir"], "refs")
_refs.setdefault("guide_csv", config.get("guide_csv", ""))
_refs.setdefault("guide_fasta",       os.path.join(_ref_dir, "guides.fasta"))
_refs.setdefault("guide_t2g_2col",    os.path.join(_ref_dir, "t2g.tsv"))
_refs.setdefault("sgRNA_index_dir",   os.path.join(_ref_dir, "piscem_index"))
_refs.setdefault("guide_hash",        os.path.join(_ref_dir, "guide_hash.pkl"))
_refs.setdefault("whitelist_dir",     os.path.join(_ref_dir, "whitelist_cache"))

# ---- Translation table auto-derivation ----
# Class A (3' dual-oligo) chemistries need a RNA<->Feature translation table;
# derive its path into whitelist_dir (downloaded on demand). Class B (5') has
# translation_file: null. A user-set config["translation_table"] wins.
_chem_spec = config.get("_chemistry") or {}
if _chem_spec.get("translation_file") and "translation_table" not in config:
    config["translation_table"] = os.path.join(
        _refs["whitelist_dir"], _chem_spec["translation_file"])

# ---- Section defaults (any section may be omitted from config.yaml) ----
_wl = config.setdefault("whitelist", {})
_wl.setdefault("min_umi", 1000)
_wl.setdefault("min_genes", 500)

_res = config.setdefault("resources", {})
_res.setdefault("simpleaf_quant_threads", 12)
_res.setdefault("simpleaf_index_threads", 4)
_res.setdefault("hash_quant_threads", 4)

_ham = config.setdefault("hash_matcher", {})
_ham.setdefault("umi_threshold", 1)
_ham.setdefault("cb_max_hamming", 1)

_saf = config.setdefault("simpleaf", {})
_saf_idx = _saf.setdefault("index", {})
_saf_idx.setdefault("kmer_length", 15)
_saf_idx.setdefault("minimizer_length", 11)
_saf_quant = _saf.setdefault("quant", {})
_saf_quant.setdefault("resolution", "parsimony-gene")
_saf_quant.setdefault("use_knee", False)

# Assignment guide_csv falls back to the top-level guide_csv
_asgn = config.setdefault("assignment", {})
_asgn.setdefault("guide_csv", config.get("guide_csv", ""))
_asgn.setdefault("methods", [])


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
# Extend rule all inputs with assignment targets when methods configured
_assignment_methods = config.get("assignment", {}).get("methods", [])
_assignment_targets = []
for _m in _assignment_methods:
    _base = os.path.join(config["out_dir"], "assignment", _m)
    _assignment_targets.append(os.path.join(_base, "assignments.csv"))
    _assignment_targets.append(os.path.join(_base, "perturbation_obs.csv"))

rule all:
    input:
        os.path.join(config["out_dir"], "guide_matrix", "merged_matrix.mtx.gz"),
        os.path.join(config["out_dir"], "guide_matrix", "merged_barcodes.tsv.gz"),
        os.path.join(config["out_dir"], "guide_matrix", "merged_features.tsv.gz"),
        *_assignment_targets,


# ---- Import rule modules ----
include: "rules/reference.smk"
include: "rules/whitelist.smk"

if METHOD == "simpleaf":
    include: "rules/quant.smk"
else:
    include: "rules/guide_quant.smk"

include: "rules/merge.smk"

if _assignment_methods:
    include: "rules/assignment.smk"
