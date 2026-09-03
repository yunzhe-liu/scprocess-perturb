# ==============================================================================
# integration.smk — GEX + guide-assignment AnnData export
# ==============================================================================
# This is the terminal step of the complete workflow. Existing MEX and
# assignment outputs remain unchanged; integration produces one canonical
# AnnData artifact.
#
# Input cell universe: concatenated per-group GEX cells.
# Cell keys: normalized 16mer + the same group suffix used by merge.smk.
# ============================================================================

import json
import shlex

INTEGRATION_CONFIG = config.get("integration", {})
GUIDE_DESIGN = config["assignment"]["guide_design"]
OUTPUT_NAME = "perturbation_adata.h5ad"


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
            config["out_dir"], "integration", OUTPUT_NAME
        ),
    params:
        guide_design = GUIDE_DESIGN,
        script = os.path.join(config["proj_dir"], "scripts", "integrate_multimodal.py"),
        gex_args = lambda wildcards: _integration_pairs(
            [("gex", group, GROUPS[group]["gex_h5"]) for group in GROUPS]
        ),
        assignment = INTEGRATION_ASSIGNMENT,
        guide_csv = config.get("assignment", {}).get("guide_csv") or config.get(
            "guide_csv", ""
        ),
        counts_source = INTEGRATION_CONFIG.get("counts_source", ""),
        normalized_source = INTEGRATION_CONFIG.get("normalized_source", ""),
        input_kind = INTEGRATION_CONFIG.get("input_kind", "auto"),
        counts_layer = INTEGRATION_CONFIG.get("counts_layer", "counts"),
        target_sum = INTEGRATION_CONFIG.get("target_sum", 10000.0),
        max_materialized_nnz = INTEGRATION_CONFIG.get(
            "max_materialized_nnz", 100000000
        ),
        stream_chunk_nnz = INTEGRATION_CONFIG.get("stream_chunk_nnz", 50000000),
        max_input_nnz = INTEGRATION_CONFIG.get("max_input_nnz", 5000000000),
        max_output_gb = INTEGRATION_CONFIG.get("max_output_gb", 250.0),
        min_free_disk_gb = INTEGRATION_CONFIG.get("min_free_disk_gb", 200.0),
        max_process_memory_gb = INTEGRATION_CONFIG.get(
            "max_process_memory_gb", 192.0
        ),
        config_json = lambda wildcards: json.dumps({
            "dataset": INTEGRATION_CONFIG.get("dataset", "workflow"),
            "cell_key_mode": "auto",
            "cell_universe": "gex",
            "assignment": {"guide_design": GUIDE_DESIGN},
            "integration": {
                "input_kind": INTEGRATION_CONFIG.get("input_kind", "auto"),
                "counts_layer": INTEGRATION_CONFIG.get("counts_layer", "counts"),
                "target_sum": INTEGRATION_CONFIG.get("target_sum", 10000.0),
                "counts_source": INTEGRATION_CONFIG.get("counts_source", ""),
                "normalized_source": INTEGRATION_CONFIG.get("normalized_source", ""),
            },
        }),
    log:
        os.path.join(
            config["log_dir"], "integration", f"{OUTPUT_NAME}.log"
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
            --guide-design "{params.guide_design}" \\
            --guide-csv "{params.guide_csv}" \\
            --counts-source "{params.counts_source}" \\
            --normalized-source "{params.normalized_source}" \\
            --input-kind "{params.input_kind}" \\
            --counts-layer "{params.counts_layer}" \\
            --target-sum "{params.target_sum}" \\
            --max-materialized-nnz "{params.max_materialized_nnz}" \\
            --stream-chunk-nnz "{params.stream_chunk_nnz}" \\
            --max-input-nnz "{params.max_input_nnz}" \\
            --max-output-gb "{params.max_output_gb}" \\
            --min-free-disk-gb "{params.min_free_disk_gb}" \\
            --max-process-memory-gb "{params.max_process_memory_gb}" \\
            --config-json '{params.config_json}' \\
            --out "{output.adata}"
    """
