#!/usr/bin/env bash
# ==============================================================================
# preprocess_cutadapt.sh — Trim TSO+poly-G from sgRNA R2 reads (paired-end)
# ==============================================================================
# Usage:
#   bash preprocess_cutadapt.sh <samples.yaml> [--dry-run]
#
# Reads samples.yaml to discover per-group FASTQ directories and glob patterns.
# For each group, runs cutadapt in paired-end mode to trim the 5' TSO+poly-G
# adapter from R2 reads. R1 reads are passed through unchanged.
#
# Output is written to <sgRNA_fastq_dir>/trimmed/ for each group.
# The main workflow auto-detects these trimmed files and uses them in preference
# to raw FASTQ.
#
# Requires: cutadapt >= 5.2 (conda env: scp_analysis or pip install cutadapt)
#
# Adapter strategy:
#   Four poly-G length variants (3G, 4G, 5G, 6G) are used as concurrent 5'
#   anchors. The dynamic programming matrix in cutadapt calculates the
#   best-scoring alignment for each read, accommodating the natural 3-6bp
#   poly-G tail length variation from 10x template switching.
#
#   TSO (27bp):            AAGCAGTGGTATCAACGCAGAGTACAT
#   + poly-G 3-6bp:        GGG / GGGG / GGGGG / GGGGGG
#   = 4 anchors (30-33bp)
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJ_DIR="$(dirname "$SCRIPT_DIR")"

# --- Parse arguments ---
SAMPLES_YAML="${1:-$PROJ_DIR/config/samples.yaml}"
DRY_RUN=false
[[ "${2:-}" == "--dry-run" ]] && DRY_RUN=true

if [[ ! -f "$SAMPLES_YAML" ]]; then
    echo "ERROR: samples.yaml not found: $SAMPLES_YAML"
    exit 1
fi

# --- Conda setup ---
CONDA_BASE="${CONDA_BASE:-$HOME/software/miniconda3}"
source "$CONDA_BASE/etc/profile.d/conda.sh" 2>/dev/null || true
conda activate scp_analysis 2>/dev/null || true
if command -v cutadapt &>/dev/null; then
    echo "cutadapt found: $(cutadapt --version)"
else
    echo "ERROR: cutadapt not found. Install: conda activate scp_analysis && pip install cutadapt"
    exit 1
fi

# --- 5' anchor definitions (TSO + variable poly-G) ---
# These four anchors cover the natural poly-G length variation (3-6 G's).
ANCHOR_G3="AAGCAGTGGTATCAACGCAGAGTACATGGG"
ANCHOR_G4="AAGCAGTGGTATCAACGCAGAGTACATGGGG"
ANCHOR_G5="AAGCAGTGGTATCAACGCAGAGTACATGGGGG"
ANCHOR_G6="AAGCAGTGGTATCAACGCAGAGTACATGGGGGG"
ERR_TOL=0.13    # 13% error tolerance (~4 mismatches per 30bp)
MIN_OVERLAP=10  # Minimum 10bp overlap for adapter match
MIN_LEN=18      # Discard trimmed reads shorter than 18bp
THREADS=8

# --- Parse samples.yaml to discover groups ---
echo "Reading sample topology from: $SAMPLES_YAML"
GROUPS=$(python3 "$PROJ_DIR/scripts/parse_samples.py" "$SAMPLES_YAML")

if [[ -z "$GROUPS" ]]; then
    echo "ERROR: No groups parsed from samples.yaml"
    exit 1
fi

echo ""
echo "============================================"
echo "cutadapt preprocessing (5' TSO+poly-G trim)"
echo "============================================"
echo "Anchors: TSO+3G, TSO+4G, TSO+5G, TSO+6G"
echo "Error tolerance: ${ERR_TOL} (${MIN_OVERLAP}bp min overlap)"
echo "Min trimmed length: ${MIN_LEN}bp"
echo "Threads: ${THREADS}"
echo ""

# --- Process each group ---
TOTAL_GROUPS=0
TOTAL_READS=0
TOTAL_PASSED=0

while IFS=$'\t' read -r gname fastq_dir r1_pat r2_pat; do
    TOTAL_GROUPS=$((TOTAL_GROUPS + 1))
    TRIM_DIR="${fastq_dir}/trimmed"
    
    # Find FASTQ files by glob
    shopt -s nullglob
    r1_files=($fastq_dir/$r1_pat)
    r2_files=($fastq_dir/$r2_pat)
    shopt -u nullglob
    
    if [[ ${#r1_files[@]} -eq 0 ]]; then
        echo "[$gname] WARNING: No R1 files matching '$r1_pat' in $fastq_dir — skipping"
        continue
    fi
    if [[ ${#r1_files[@]} -ne ${#r2_files[@]} ]]; then
        echo "[$gname] WARNING: R1/R2 file count mismatch (${#r1_files[@]} vs ${#r2_files[@]}) — skipping"
        continue
    fi
    
    echo "[$gname] Found ${#r1_files[@]} FASTQ pair(s)"
    
    if $DRY_RUN; then
        echo "  (dry-run) Would trim ${#r1_files[@]} pairs → $TRIM_DIR/"
        continue
    fi
    
    mkdir -p "$TRIM_DIR"
    
    group_reads=0
    group_passed=0
    
    for i in "${!r1_files[@]}"; do
        r1_in="${r1_files[$i]}"
        r2_in="${r2_files[$i]}"
        
        # Derive output filenames — keep the same base name
        r1_base="$(basename "$r1_in")"
        r2_base="$(basename "$r2_in")"
        r1_out="$TRIM_DIR/$r1_base"
        r2_out="$TRIM_DIR/$r2_base"
        
        # Run cutadapt
        cutadapt \
            -G "$ANCHOR_G3" \
            -G "$ANCHOR_G4" \
            -G "$ANCHOR_G5" \
            -G "$ANCHOR_G6" \
            -e "$ERR_TOL" \
            -O "$MIN_OVERLAP" \
            --discard-untrimmed \
            -m "$MIN_LEN" \
            -j "$THREADS" \
            -o "$r1_out" \
            -p "$r2_out" \
            "$r1_in" "$r2_in" 2>&1 | grep -E "processed|adapter|pass|filter|short" || true
        
        # Count reads from cutadapt's stderr (hidden in the pipe above)
    done
    
    # Summary for this group
    n_trimmed=$(ls "$TRIM_DIR"/*.fastq.gz 2>/dev/null | wc -l)
    echo "  → $n_trimmed trimmed FASTQ files in $TRIM_DIR/"
    
done <<< "$GROUPS"

echo ""
echo "============================================"
echo "Preprocessing complete ($TOTAL_GROUPS groups)"
echo "Trimmed files are in <sgRNA_fastq_dir>/trimmed/"
echo "The workflow will auto-detect and use them."
echo "============================================"
