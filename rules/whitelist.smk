# ==============================================================================
# whitelist.smk — Extract Cell Barcode Whitelist from GEX H5 Matrix
# ==============================================================================
# Filters low-quality droplets using hard thresholds (UMI > 1000,
# genes detected > 500), producing a cell barcode list validated by
# transcriptome activity.
# ==============================================================================

rule extract_whitelist:
    input:
        h5 = lambda wildcards: GROUPS[wildcards.group]["gex_h5"],
    output:
        wl_csv     = os.path.join(config["out_dir"], "{group}", "barcode_whitelist.csv"),
        wl_noheader = os.path.join(config["out_dir"], "{group}", "barcode_whitelist_noheader.txt"),
    params:
        script    = os.path.join(config["proj_dir"], "scripts", "filter_barcodes.py"),
        min_umi   = config["whitelist"]["min_umi"],
        min_genes = config["whitelist"]["min_genes"],
        lock_yaml = os.path.join(config["proj_dir"], "envs", "scp_analysis.lock.yaml"),
    log:
        os.path.join(config["log_dir"], "whitelist", "{group}.log"),
    benchmark:
        os.path.join(config["log_dir"], "benchmark", "whitelist_{group}.tsv"),
    threads: 1
    conda:
        os.path.join(config["proj_dir"], "envs", "scp_analysis.lock.yaml"),
    shell:"""
        set -euo pipefail
        exec &>> {log}
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
        python {params.script} \
            --h5 "{input.h5}" \
            --out-wl "{output.wl_csv}" \
            --out-noheader "{output.wl_noheader}" \
            --min-umi {params.min_umi} \
            --min-genes {params.min_genes}
    """
