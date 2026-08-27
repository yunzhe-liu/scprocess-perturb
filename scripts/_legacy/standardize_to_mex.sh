#!/usr/bin/env bash
# ==============================================================================
# standardize_to_mex.sh — simpleaf quant output → 10x MEX standard format
# ==============================================================================
# Usage: standardize_to_mex.sh <output_prefix> <alevin_dir> <mex_out_dir>
#
# Outputs:
#   {mex_out_dir}/{prefix}_matrix.mtx.gz     — count matrix
#   {mex_out_dir}/{prefix}_barcodes.tsv.gz   — cell barcodes
#   {mex_out_dir}/{prefix}_features.tsv.gz   — features (guide_id \t guide_id \t "CRISPR Guide Capture")
# ==============================================================================
set -euo pipefail

PREFIX="${1:?Usage: $0 <output_prefix> <alevin_dir> <mex_out_dir>}"
INDIR="${2:?}"
OUTDIR="${3:?}"

mkdir -p "${OUTDIR}"

echo "Converting: ${INDIR} → ${OUTDIR}/${PREFIX}_*"

# 1. Matrix — gzip compress
gzip -c "${INDIR}/quants_mat.mtx" > "${OUTDIR}/${PREFIX}_matrix.mtx.gz"
echo "  matrix:  $(du -h "${OUTDIR}/${PREFIX}_matrix.mtx.gz" | cut -f1)"

# 2. Barcodes — gzip compress
gzip -c "${INDIR}/quants_mat_rows.txt" > "${OUTDIR}/${PREFIX}_barcodes.tsv.gz"
echo "  barcodes: $(du -h "${OUTDIR}/${PREFIX}_barcodes.tsv.gz" | cut -f1)"

# 3. Features — 3 columns: id, name, type (CRISPR Guide Capture)
awk -F'\t' '{print $1 "\t" $1 "\t" "CRISPR Guide Capture"}' \
    "${INDIR}/quants_mat_cols.txt" | \
    gzip -c > "${OUTDIR}/${PREFIX}_features.tsv.gz"
echo "  features: $(du -h "${OUTDIR}/${PREFIX}_features.tsv.gz" | cut -f1)"

echo "Done: $(ls "${OUTDIR}/${PREFIX}_"*)"
