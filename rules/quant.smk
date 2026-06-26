# ==============================================================================
# quant.smk — sgRNA Guide Quantification (simpleaf quant)
# ==============================================================================
# Supports two modes via config["simpleaf"]["quant"]["use_knee"]:
#   false (default) — GEX-whitelist-gated via --explicit-pl
#   true            — knee-calling via --knee (no whitelist required)
# ==============================================================================

import os


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _join_fastq(files: list) -> str:
    """Join FASTQ file paths with comma for simpleaf --reads1/--reads2."""
    return ",".join(files)


def _resolve_reads(wildcards, read: str) -> list:
    """Return trimmed FASTQ if available, otherwise raw FASTQ."""
    gcfg = GROUPS[wildcards.group]
    trimmed = gcfg.get(f"_sgRNA_{read}_trimmed", [])
    if trimmed:
        return trimmed
    return gcfg.get(f"_sgRNA_{read}", [])


# ---------------------------------------------------------------------------
# Rule: simpleaf_quant
# ---------------------------------------------------------------------------
USE_KNEE = config.get("simpleaf", {}).get("quant", {}).get("use_knee", False)

rule simpleaf_quant:
    input:
        r1_files = lambda wildcards: _resolve_reads(wildcards, "r1"),
        r2_files = lambda wildcards: _resolve_reads(wildcards, "r2"),
        wl       = (lambda wildcards: []) if USE_KNEE else os.path.join(config["out_dir"], "{group}", "barcode_whitelist_noheader.txt"),
        idx_ctab = os.path.join(config["references"]["sgRNA_index_dir"], "index", "piscem_idx.ctab"),
        t2g      = config["references"]["guide_t2g_2col"],
    output:
        mtx  = os.path.join(config["out_dir"], "{group}", "simpleaf_quant", "af_quant", "alevin", "quants_mat.mtx"),
        rows = os.path.join(config["out_dir"], "{group}", "simpleaf_quant", "af_quant", "alevin", "quants_mat_rows.txt"),
        cols = os.path.join(config["out_dir"], "{group}", "simpleaf_quant", "af_quant", "alevin", "quants_mat_cols.txt"),
    params:
        out_dir    = os.path.join(config["out_dir"], "{group}", "simpleaf_quant"),
        index_dir  = config["references"]["sgRNA_index_dir"],
        index_piscem_pref = config["references"]["sgRNA_index_dir"] + "/index/piscem_idx",
        # v0.1.2: chemistry resolved from _chemistry spec (if available),
        # with geometry_override taking precedence over the raw config value.
        chemistry  = (config.get("_chemistry") or {}).get("geometry_override") or config["simpleaf"]["quant"]["chemistry"],
        resolution = config["simpleaf"]["quant"]["resolution"],
        af_home    = config["simpleaf"]["af_home"],
        lock_yaml  = os.path.join(config["proj_dir"], "envs", "simpleaf.lock.yaml"),
        use_knee   = USE_KNEE,
        reads1 = lambda wildcards: _join_fastq(_resolve_reads(wildcards, "r1")),
        reads2 = lambda wildcards: _join_fastq(_resolve_reads(wildcards, "r2")),
    log:
        os.path.join(config["log_dir"], "quant", "{group}.log"),
    benchmark:
        os.path.join(config["log_dir"], "benchmark", "quant_{group}.tsv"),
    threads: config["resources"]["simpleaf_quant_threads"]
    conda:
        os.path.join(config["proj_dir"], "envs", "simpleaf.lock.yaml"),
    shell:"""
        set -euo pipefail
        exec &>> {log}
        CONDA_BASE="$HOME/software/miniconda3"
        ENV_NAME="simpleaf"
        source "$CONDA_BASE/etc/profile.d/conda.sh"
        if [ -d "$CONDA_BASE/envs/$ENV_NAME" ]; then
            conda activate "$ENV_NAME"
        else
            echo "Creating conda environment '$ENV_NAME' from lock file..."
            conda env create -f "{params.lock_yaml}"
            conda activate "$ENV_NAME"
        fi
        export ALEVIN_FRY_HOME="{params.af_home}"
        simpleaf set-paths

        echo "=== simpleaf quant: {wildcards.group} ==="
        echo "  Index:     {params.index_piscem_pref}"
        echo "  Chemistry: {params.chemistry}"
        echo "  Threads:   {threads}"

        rm -rf "{params.out_dir}"
        mkdir -p "{params.out_dir}"

        USE_KNEE="{params.use_knee}"
        if [ "$USE_KNEE" = "True" ]; then
            echo "  Mode:      knee (no whitelist)"
            simpleaf quant \
                --chemistry "{params.chemistry}" \
                --output "{params.out_dir}" \
                --threads {threads} \
                --index "{params.index_piscem_pref}" \
                --reads1 "{params.reads1}" \
                --reads2 "{params.reads2}" \
                --t2g-map "{input.t2g}" \
                --resolution "{params.resolution}" \
                --knee
        else
            echo "  Whitelist: {input.wl}"
            echo "  Mode:      whitelist"
            simpleaf quant \
                --chemistry "{params.chemistry}" \
                --output "{params.out_dir}" \
                --threads {threads} \
                --index "{params.index_piscem_pref}" \
                --reads1 "{params.reads1}" \
                --reads2 "{params.reads2}" \
                --t2g-map "{input.t2g}" \
                --resolution "{params.resolution}" \
                --explicit-pl "{input.wl}"
        fi

        echo "  Done."
        echo "  Output: {params.out_dir}/af_quant/alevin/"
    """
