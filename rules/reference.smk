# ==============================================================================
# reference.smk — Guide Reference Preparation & piscem Index Building
# ==============================================================================
# Rule 1: generate_guide_reference — one-time: CSV → FASTA + t2g
# Rule 2: build_sgRNA_index        — one-time: FASTA → piscem index
# ==============================================================================

# ---------------------------------------------------------------------------
# Rule: generate_guide_reference
# Decomposes wide-format dual-sgRNA library into single-guide FASTA + t2g map.
# Skips if outputs already exist (rule is idempotent).
# ---------------------------------------------------------------------------
rule generate_guide_reference:
    input:
        csv = config["references"]["guide_csv"],
    output:
        fasta = config["references"]["guide_fasta"],
        t2g   = config["references"]["guide_t2g_2col"],
    params:
        adapter  = os.path.join(config["proj_dir"], "scripts", "feature_reference_adapter.py"),
        lock_yaml = os.path.join(config["proj_dir"], "envs", "scp_analysis.lock.yaml"),
    conda:
        os.path.join(config["proj_dir"], "envs", "scp_analysis.lock.yaml"),
    shell:"""
        set -euo pipefail
        CONDA_BASE="$HOME/software/miniconda3"
        ENV_NAME="scp_analysis"
        source "$CONDA_BASE/etc/profile.d/conda.sh"
        if [ -d "$CONDA_BASE/envs/$ENV_NAME" ]; then
            conda activate "$ENV_NAME"
        else
            echo "Creating conda environment '$ENV_NAME' from lock file..."
            conda env create -f "{params.lock_yaml}"
            conda activate "$ENV_NAME"
        fi
        # Skip if outputs already exist (adapter is idempotent; avoid redundant runs)
        if [ -f "{output.fasta}" ] && [ -f "{output.t2g}" ]; then
            echo "Guide reference files exist, skipping generation."
            exit 0
        fi
        python {params.adapter} --csv {input.csv} \
            --out-fasta {output.fasta} \
            --out-t2g {output.t2g}
    """

# ---------------------------------------------------------------------------
# Rule: build_sgRNA_index
# Builds piscem index (k=15, m=11) from guides.fasta.
# The index is globally shared and reused by all sample groups.
# ---------------------------------------------------------------------------
rule build_sgRNA_index:
    input:
        fasta = config["references"]["guide_fasta"],
    output:
        # Sentinel file for the piscem index (pre-built; skipped if exists)
        ctab = os.path.join(config["references"]["sgRNA_index_dir"], "index", "piscem_idx.ctab"),
    params:
        out_dir   = config["references"]["sgRNA_index_dir"],
        kmer      = config["simpleaf"]["index"]["kmer_length"],
        minimizer = config["simpleaf"]["index"]["minimizer_length"],
        af_home   = config["simpleaf"]["af_home"],
        lock_yaml = os.path.join(config["proj_dir"], "envs", "simpleaf.lock.yaml"),
    threads: config["resources"]["simpleaf_index_threads"]
    conda:
        os.path.join(config["proj_dir"], "envs", "simpleaf.lock.yaml"),
    shell:"""
        set -euo pipefail
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
        mkdir -p "{params.out_dir}"

        simpleaf index \
            --output "{params.out_dir}" \
            --threads {threads} \
            --ref-seq "{input.fasta}" \
            --kmer-length {params.kmer} \
            --minimizer-length {params.minimizer}

        # simpleaf writes index to out_dir/index/,
        # but quant expects index files at parent dir; create symlinks
        for f in "{params.out_dir}"/index/piscem_idx.*; do
            base=$(basename "$f")
            [ -e "{params.out_dir}/$base" ] || ln -sf "index/$base" "{params.out_dir}/$base"
        done
        [ -f "{params.out_dir}/index/simpleaf_index.json" ] && \
            ln -sf "index/simpleaf_index.json" "{params.out_dir}/simpleaf_index.json"

        echo "Index built: {params.out_dir}"
    """
