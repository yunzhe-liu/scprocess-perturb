#!/bin/bash
# ==============================================================================
# param_sweep_master.sh — simpleaf parameter sweep (k=13 baseline)
# ==============================================================================
# Sweeps:
#   1. UMI resolution (5 modes)   — k=13, m=9, whitelist mode
#   2. minimizer length (4 values) — k=13, parsimony-gene, whitelist mode
#
# Each run: build index (if needed) → extract_whitelist → simpleaf_quant ×48
#           → merge → barcode_translation
#
# Output:
#   /data/yunzliu/results/guide_extraction/param_sweep/resolution/{mode}/merged/
#   /data/yunzliu/results/guide_extraction/param_sweep/minimizer/{mXX}/merged/
# ==============================================================================

set -euo pipefail

SGPROCESS=/home/yunzliu/sgprocess
CONDA_BASE=/home/yunzliu/software/miniconda3
INDEX_BASE=/data/yunzliu/simpleaf_index_k13
REF_FASTA=/data/yunzliu/references/guides_fixed.fasta
T2G=/data/yunzliu/references/t2g_2col_guide.tsv

source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate scprocess

# ── Helper: generate snakemake config YAML ──
gen_config() {
    local out_dir=$1 log_dir=$2 index_dir=$3 resolution=$4 m_val=$5
    cat > /tmp/_param_sweep_config.yaml << EOF
proj_dir: ${SGPROCESS}
out_dir: ${out_dir}
log_dir: ${log_dir}
work_dir: /data/yunzliu/tmp/param_sweep_\$(basename ${out_dir})

guide_extraction:
  method: simpleaf

references:
  guide_csv: /data/yunzliu/references/raw_guides_k562_essential.csv
  guide_fasta: ${REF_FASTA}
  guide_hash: /data/yunzliu/references/guide_hash_v3.pkl
  guide_t2g_2col: ${T2G}
  sgRNA_index_dir: ${index_dir}

simpleaf:
  conda_env: simpleaf
  bin_dir: ${CONDA_BASE}/envs/simpleaf/bin
  af_home: ${CONDA_BASE}/envs/simpleaf/opt/alevin-fry
  index:
    kmer_length: 13
    minimizer_length: ${m_val}
  quant:
    chemistry: 10xv3
    resolution: ${resolution}

preprocess:
  trimmed: false
whitelist:
  min_umi: 1000
  min_genes: 500
resources:
  default_threads: 8
  simpleaf_index_threads: 4
  simpleaf_quant_threads: 12
  hash_quant_threads: 4
snakemake:
  latency_wait: 60
  rerun_triggers:
  - mtime
  - params
  - input
  - software-env
  - code
EOF
}

# ── Helper: run one sweep ──
run_sweep() {
    local out_dir=$1 log_dir=$2 index_dir=$3 resolution=$4 m_val=$5 label=$6

    echo ""
    echo "===================================================================="
    echo "  PARAM SWEEP: ${label}"
    echo "  resolution=${resolution}  m=${m_val}  index=${index_dir}"
    echo "  out_dir=${out_dir}"
    echo "===================================================================="

    gen_config "$out_dir" "$log_dir" "$index_dir" "$resolution" "$m_val"

    # Build index if not present
    if [ ! -f "${index_dir}/index/piscem_idx.ctab" ]; then
        echo "[1/4] Building piscem index (k=13, m=${m_val})..."
        source "$CONDA_BASE/etc/profile.d/conda.sh"
        conda activate simpleaf
        export ALEVIN_FRY_HOME="${CONDA_BASE}/envs/simpleaf/opt/alevin-fry"
        rm -rf "${index_dir}"
        simpleaf index \
            --output "${index_dir}" \
            --threads 4 \
            --ref-seq "${REF_FASTA}" \
            --kmer-length 13 \
            --minimizer-length "${m_val}"
        # Create parent-level symlinks (required by simpleaf quant)
        for f in "${index_dir}"/index/piscem_idx.*; do
            base=$(basename "$f")
            [ -e "${index_dir}/${base}" ] || ln -sf "index/${base}" "${index_dir}/${base}"
        done
        [ -f "${index_dir}/index/simpleaf_index.json" ] && \
            ln -sf "index/simpleaf_index.json" "${index_dir}/simpleaf_index.json"
        echo "  Index built: ${index_dir}"
        conda activate scprocess
    else
        echo "[1/4] Index already exists: ${index_dir}"
    fi

    # Unlock stale locks
    cd "$SGPROCESS"
    snakemake --configfile /tmp/_param_sweep_config.yaml --unlock 2>&1 | tail -1

    # Run snakemake
    echo "[2/4] Running snakemake (48 lanes quant + merge)..."
    snakemake --configfile /tmp/_param_sweep_config.yaml \
        --cores 48 --latency-wait 60 < /dev/null \
        > "${log_dir}/snakemake.log" 2>&1

    echo "[3/4] Snakemake done. Checking output..."

    # Verify
    N_MTX=$(find "${out_dir}" -name "quants_mat.mtx" 2>/dev/null | wc -l)
    echo "  Per-lane matrices: ${N_MTX}/48"
    if [ -f "${out_dir}/merged/merged_matrix.mtx.gz" ]; then
        echo "  Merged MEX: OK"
    else
        echo "  Merged MEX: MISSING!"
    fi

    echo "[4/4] ${label} — COMPLETE"
}

# =============================================================================
# PART 1: UMI Resolution Sweep (k=13, m=9, whitelist)
# =============================================================================
RES_INDEX=/data/yunzliu/simpleaf_index_k13
RES_LOG=/data/yunzliu/logs/guide_extraction/param_sweep/resolution

for mode in cr-like cr-like-em parsimony parsimony-em parsimony-gene-em; do
    OUT_DIR=/data/yunzliu/results/guide_extraction/param_sweep/resolution/${mode}
    run_sweep "$OUT_DIR" "$RES_LOG" "$RES_INDEX" "$mode" 9 "resolution=${mode}"
done

# =============================================================================
# PART 2: Minimizer Sweep (k=13, parsimony-gene, whitelist)
# =============================================================================
M_LOG=/data/yunzliu/logs/guide_extraction/param_sweep/minimizer
RES=parsimony-gene

for m in 3 5 7 11; do
    IDX=/data/yunzliu/simpleaf_index_k13_m${m}
    OUT_DIR=/data/yunzliu/results/guide_extraction/param_sweep/minimizer/m${m}
    run_sweep "$OUT_DIR" "$M_LOG" "$IDX" "$RES" "$m" "m=${m}"
done

echo ""
echo "===================================================================="
echo "  ALL PARAM SWEEPS COMPLETE"
echo "===================================================================="
echo ""
echo "Resolution sweep outputs:"
for mode in cr-like cr-like-em parsimony parsimony-em parsimony-gene-em; do
    echo "  /data/yunzliu/results/guide_extraction/param_sweep/resolution/${mode}/merged/"
done
echo ""
echo "Minimizer sweep outputs:"
for m in 3 5 7 11; do
    echo "  /data/yunzliu/results/guide_extraction/param_sweep/minimizer/m${m}/merged/"
done
