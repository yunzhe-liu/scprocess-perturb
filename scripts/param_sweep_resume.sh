#!/bin/bash
# ==============================================================================
# param_sweep_resume.sh — Resume parameter sweep after crash
# ==============================================================================
# Completed: cr-like
# Remaining: cr-like-em, parsimony, parsimony-em, parsimony-gene-em,
#            minimizer m=3,5,7,11

set -euo pipefail

SCPROCESS_PERTURB=/home/yunzliu/scprocess-perturb
CONDA_BASE=/home/yunzliu/software/miniconda3
INDEX_BASE=/data/yunzliu/simpleaf_index_k13
REF_FASTA=/data/yunzliu/references/guides_fixed.fasta
T2G=/data/yunzliu/references/t2g_2col_guide.tsv

source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate scprocess

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

run_sweep() {
    local out_dir=$1 log_dir=$2 index_dir=$3 resolution=$4 m_val=$5 label=$6
    echo ""
    echo "===================================================================="
    echo "  PARAM SWEEP: ${label}"
    echo "  resolution=${resolution}  m=${m_val}  index=${index_dir}"
    echo "===================================================================="
    gen_config "$out_dir" "$log_dir" "$index_dir" "$resolution" "$m_val"

    # Build index if needed
    if [ ! -f "${index_dir}/index/piscem_idx.ctab" ]; then
        echo "[1/4] Building index (k=13, m=${m_val})..."
        source "$CONDA_BASE/etc/profile.d/conda.sh"
        conda activate simpleaf
        export ALEVIN_FRY_HOME="${CONDA_BASE}/envs/simpleaf/opt/alevin-fry"
        rm -rf "${index_dir}"
        simpleaf index --output "${index_dir}" --threads 4 \
            --ref-seq "${REF_FASTA}" --kmer-length 13 --minimizer-length "${m_val}"
        for f in "${index_dir}"/index/piscem_idx.*; do
            base=$(basename "$f")
            [ -e "${index_dir}/${base}" ] || ln -sf "index/${base}" "${index_dir}/${base}"
        done
        [ -f "${index_dir}/index/simpleaf_index.json" ] && \
            ln -sf "index/simpleaf_index.json" "${index_dir}/simpleaf_index.json"
        conda activate scprocess
    else
        echo "[1/4] Index ready: ${index_dir}"
    fi

    cd "$SGPROCESS"
    snakemake --configfile /tmp/_param_sweep_config.yaml --unlock 2>&1 | tail -1

    echo "[2/4] Running snakemake..."
    snakemake --configfile /tmp/_param_sweep_config.yaml \
        --cores 48 --latency-wait 60 < /dev/null \
        > "${log_dir}/snakemake.log" 2>&1

    N_MTX=$(find "${out_dir}" -name "quants_mat.mtx" 2>/dev/null | wc -l)
    echo "[3/4] Per-lane: ${N_MTX}/48"
    [ -f "${out_dir}/merged/merged_matrix.mtx.gz" ] && echo "  Merged: OK" || echo "  Merged: MISSING!"
    echo "[4/4] ${label} — COMPLETE"
}

# ─── PART 1: Resolution sweep (remaining 4 modes) ───
RES_INDEX=/data/yunzliu/simpleaf_index_k13
RES_LOG=/data/yunzliu/logs/guide_extraction/param_sweep/resolution
M_LOG=/data/yunzliu/logs/guide_extraction/param_sweep/minimizer
RES=parsimony-gene

for mode in cr-like-em parsimony parsimony-em parsimony-gene-em; do
    OUT_DIR=/data/yunzliu/results/guide_extraction/param_sweep/resolution/${mode}
    run_sweep "$OUT_DIR" "$RES_LOG" "$RES_INDEX" "$mode" 9 "resolution=${mode}"
done

# ─── PART 2: Minimizer sweep (m=3,5,7,11) ───
for m in 3 5 7 11; do
    IDX=/data/yunzliu/simpleaf_index_k13_m${m}
    OUT_DIR=/data/yunzliu/results/guide_extraction/param_sweep/minimizer/m${m}
    run_sweep "$OUT_DIR" "$M_LOG" "$IDX" "$RES" "$m" "m=${m}"
done

echo ""
echo "ALL DONE."
echo "Resolution: $(ls -d /data/yunzliu/results/guide_extraction/param_sweep/resolution/*/merged/ 2>/dev/null | wc -l)/5"
echo "Minimizer:  $(ls -d /data/yunzliu/results/guide_extraction/param_sweep/minimizer/*/merged/ 2>/dev/null | wc -l)/4"
