# QC and standardization. Expression matrices and cell/gene universes are
# preserved; this rule validates contracts and adds a fixed status interface.

FINAL_CONFIG = config.get("finalization", {})
FINAL_METHOD = _status_method
FINAL_INPUT = os.path.join(config["out_dir"], "integration", "perturbation_adata.h5ad")
FINAL_STATUS = (
    os.path.join(
        config["out_dir"], "perturbation_status", FINAL_METHOD,
        "perturbation_status.tsv.gz",
    )
    if FINAL_METHOD != "none" else None
)


def _finalization_inputs(wildcards):
    inputs = {"adata": FINAL_INPUT}
    if FINAL_STATUS is not None:
        inputs["status"] = FINAL_STATUS
    return inputs


rule finalize_perturbation_adata:
    input:
        unpack(_finalization_inputs),
    output:
        adata = os.path.join(config["out_dir"], "final", "perturbation_adata.h5ad"),
        report = os.path.join(config["out_dir"], "final", "qc_report.json"),
    params:
        method = FINAL_METHOD,
        status_arg = (
            f'--status-table "{FINAL_STATUS}"' if FINAL_STATUS is not None else ""
        ),
        finalize_script = os.path.join(
            config["proj_dir"], "scripts", "finalize_perturbation_adata.py"
        ),
        validate_script = os.path.join(
            config["proj_dir"], "scripts", "validate_final_adata.py"
        ),
        scan_chunk_values = int(FINAL_CONFIG.get("scan_chunk_values", 10000000)),
        max_input_gb = float(FINAL_CONFIG.get("max_input_gb", 300.0)),
        min_free_disk_gb = float(FINAL_CONFIG.get("min_free_disk_gb", 100.0)),
    log:
        os.path.join(config["log_dir"], "finalization", "finalization.log"),
    threads: 1
    resources:
        mem_mb = int(FINAL_CONFIG.get("memory_gb", 16)) * 1024,
    conda:
        os.path.join(config["proj_dir"], "envs", "finalization.yaml")
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname '{log}')" "$(dirname '{output.adata}')"
        python3 "{params.finalize_script}" \
            --input "{input.adata}" \
            --status-method "{params.method}" \
            {params.status_arg} \
            --output "{output.adata}" \
            --report "{output.report}" \
            --scan-chunk-values "{params.scan_chunk_values}" \
            --max-input-gb "{params.max_input_gb}" \
            --min-free-disk-gb "{params.min_free_disk_gb}" \
            > "{log}" 2>&1
        python3 "{params.validate_script}" \
            --source "{input.adata}" \
            --final "{output.adata}" \
            --status-method "{params.method}" \
            --report "{output.report}" \
            --chunk-values "{params.scan_chunk_values}" \
            --merge-existing \
            >> "{log}" 2>&1
        """
