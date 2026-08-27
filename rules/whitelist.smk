# ==============================================================================
# whitelist.smk — Extract Cell Barcode Whitelist from GEX Matrix
# ==============================================================================
# Filters low-quality droplets using hard thresholds (UMI > 1000,
# genes detected > 500), producing a cell barcode list validated by
# transcriptome activity.
#
# Supports three input formats:
#   .h5     — scprocess raw H5 (CSC sparse, h5py)  → filter_barcodes.py
#   .h5ad   — AnnData                             → inline Python
#   .h5mu   — MuData (reads mod['rna'])          → inline Python
# ==============================================================================

rule extract_whitelist:
    input:
        h5 = lambda wildcards: GROUPS[wildcards.group]["gex_h5"],
        # v0.2.2: pull the auto-derived translation table into the DAG so it is
        # downloaded on demand. Empty (no dependency) for chemistries that don't
        # translate (5') or when no translation table is configured.
        trans = (config["translation_table"]
                 if ((config.get("_chemistry") or {}).get("translation") and config.get("translation_table"))
                 else []),
    output:
        wl_csv     = os.path.join(config["out_dir"], "lanes", "{group}", "barcode_whitelist.csv"),
        wl_noheader = os.path.join(config["out_dir"], "lanes", "{group}", "barcode_whitelist_noheader.txt"),
    params:
        script    = os.path.join(config["proj_dir"], "scripts", "filter_barcodes.py"),
        min_umi   = config["whitelist"]["min_umi"],
        min_genes = config["whitelist"]["min_genes"],
        lock_yaml = os.path.join(config["proj_dir"], "envs", "scp_analysis.lock.yaml"),
        # v0.1.2: translation controlled by _chemistry spec.
        # True for dual-oligo systems (3' v3, cs1/cs2 bead capture) where
        # GEX barcodes (TruSeq) differ from Feature barcodes (Nextera).
        # False for single-oligo systems (5' v1/v2, soluble RT primer).
        do_translate = "true" if (config.get("_chemistry") or {}).get("translation") else "false",
        trans_script = os.path.join(config["proj_dir"], "scripts", "translate_barcodes.py"),
        trans_table  = config.get("translation_table", ""),
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

        GEX_FILE="{input.h5}"
        EXT="${{GEX_FILE##*.}}"

        if [[ "$EXT" == "h5ad" ]]; then
            echo "Detected h5ad format — extracting whitelist directly."
            python3 -c "
import anndata, sys
ad = anndata.read_h5ad('$GEX_FILE')
from scipy.sparse import issparse
X = ad.X
total_umi = X.sum(axis=1).A1 if issparse(X) else X.sum(axis=1)
genes_detected = (X>0).sum(axis=1).A1 if issparse(X) else (X>0).sum(axis=1)
qc = (total_umi >= {params.min_umi}) & (genes_detected >= {params.min_genes})
barcodes = [b.split('_')[1] if '_' in (b.decode() if isinstance(b, bytes) else b) else (b.decode() if isinstance(b, bytes) else b) for b in ad.obs_names[qc]]
with open('{output.wl_csv}', 'w') as f:
    f.write('barcode\\n')
    for b in barcodes: f.write(b + '\\n')
with open('{output.wl_noheader}', 'w') as f:
    for b in barcodes: f.write(b + '\\n')
print(f'Whitelist: {{len(barcodes)}} cells (≥{params.min_umi} UMI, ≥{params.min_genes} genes)')
"
        elif [[ "$EXT" == "h5mu" ]]; then
            echo "Detected h5mu format — extracting RNA modality."
            python3 -c "
import mudata as md, numpy as np
mdata = md.read_h5mu('$GEX_FILE')
ad = mdata.mod['rna']
from scipy.sparse import issparse
X = ad.X
total_umi = X.sum(axis=1).A1 if issparse(X) else X.sum(axis=1)
genes_detected = (X>0).sum(axis=1).A1 if issparse(X) else (X>0).sum(axis=1)
qc = (total_umi >= {params.min_umi}) & (genes_detected >= {params.min_genes})
barcodes = [b.split('_')[1] if '_' in (b.decode() if isinstance(b, bytes) else b) else (b.decode() if isinstance(b, bytes) else b) for b in ad.obs_names[qc]]
with open('{output.wl_csv}', 'w') as f:
    f.write('barcode\\n')
    for b in barcodes: f.write(b + '\\n')
with open('{output.wl_noheader}', 'w') as f:
    for b in barcodes: f.write(b + '\\n')
print(f'Whitelist: {{len(barcodes)}} cells (≥{params.min_umi} UMI, ≥{params.min_genes} genes)')
"
        else
            python {params.script} \
                --h5 "{input.h5}" \
                --out-wl "{output.wl_csv}" \
                --out-noheader "{output.wl_noheader}" \
                --min-umi {params.min_umi} \
                --min-genes {params.min_genes}
        fi

        # v0.1.2: Barcode translation (TO→FROM: TruSeq→Nextera).
        # Required for dual-oligo chemistries (3' v3) where GEX and Feature
        # libraries use different barcode formats on the same bead.
        # Skipped for single-oligo chemistries (5' v1/v2, soluble RT primer).
        DO_TRANSLATE="{params.do_translate}"
        if [ "$DO_TRANSLATE" = "true" ]; then
            echo "Translating barcodes (TO→FROM: TruSeq→Nextera)..."
            python3 "{params.trans_script}" \
                "{output.wl_noheader}" \
                --trans-table "{params.trans_table}" \
                --direction to_from
        fi
    """
