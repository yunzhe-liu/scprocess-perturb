# ==============================================================================
# integration.smk — GEX + guide-assignment AnnData export
# ==============================================================================
# This is the terminal step of the complete workflow. Existing MEX and
# assignment outputs remain unchanged; integration produces a new composite
# AnnData artifact for the selected mode.
#
# Input cell universe: concatenated per-group GEX cells.
# Cell keys: normalized 16mer + the same group suffix used by merge.smk.
# ============================================================================

import json
import shlex

INTEGRATION_CONFIG = config.get("integration", {})
INTEGRATION_MODE = INTEGRATION_CONFIG.get("mode", "construct")

OUTPUT_NAMES = {
    "guide_top1": "perturbation_adata_guide_top1.h5ad",
    "guide_full": "perturbation_adata_guide_full.h5ad",
    "construct": "perturbation_adata_construct.h5ad",
}


def _integration_pairs(values):
    return " ".join(
        f"--{label} {shlex.quote(f'{name}={path}') }"
        for label, name, path in values
    )


INTEGRATION_ASSIGNMENT = os.path.join(
    config["out_dir"], "assignment", _assignment_methods[0], "assignments.csv"
)


rule integrate_multimodal:
    input:
        gex = [GROUPS[group]["gex_h5"] for group in GROUPS],
        barcodes = os.path.join(
            config["out_dir"], "guide_matrix", "merged_barcodes.tsv.gz"
        ),
        assignment = INTEGRATION_ASSIGNMENT,
    output:
        adata = os.path.join(
            config["out_dir"], "integration", OUTPUT_NAMES[INTEGRATION_MODE]
        ),
    params:
        mode = INTEGRATION_MODE,
        script = os.path.join(config["proj_dir"], "scripts", "integrate_multimodal.py"),
        gex_args = lambda wildcards: _integration_pairs(
            [("gex", group, GROUPS[group]["gex_h5"]) for group in GROUPS]
        ),
        assignment = INTEGRATION_ASSIGNMENT,
        guide_csv = config.get("assignment", {}).get("guide_csv", ""),
        config_json = lambda wildcards: json.dumps({
            "dataset": INTEGRATION_CONFIG.get("dataset", "workflow"),
            "cell_key_mode": "auto",
            "cell_universe": "gex",
            "mode": INTEGRATION_MODE,
        }),
    log:
        os.path.join(
            config["log_dir"], "integration", f"{OUTPUT_NAMES[INTEGRATION_MODE]}.log"
        ),
    threads: 1
    conda:
        os.path.join(config["proj_dir"], "envs", "scp_analysis.lock.yaml")
    shell:"""
        set -euo pipefail
        exec &> {log}

        CONDA_BASE="$HOME/software/miniconda3"
        source "$CONDA_BASE/etc/profile.d/conda.sh"
        conda activate scp_analysis

        python3 "{params.script}" \\
            {params.gex_args} \\
            --assign "{params.assignment}" \\
            --barcodes "{input.barcodes}" \\
            --mode "{params.mode}" \\
            --guide-csv "{params.guide_csv}" \\
            --config-json '{params.config_json}' \\
            --out "{output.adata}"
    """
