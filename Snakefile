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
with open(os.path.join(BASEDIR, "config", "samples.yaml"), "r") as f:
    _samples_cfg = yaml.safe_load(f)

GROUPS = _samples_cfg.get("groups", {})
if not GROUPS:
    sys.exit("ERROR: No groups defined in config/samples.yaml")


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
