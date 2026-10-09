# Data Standardization and Validation. Expression matrices and cell/gene universes are
# preserved; this rule validates contracts and adds a fixed status interface.

import base64
import json

try:
    import snakemake as _snakemake_package
    PROVENANCE_SNAKEMAKE_VERSION = getattr(
        _snakemake_package, "__version__", "unknown"
    )
except ImportError:
    PROVENANCE_SNAKEMAKE_VERSION = "unknown"

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
FINAL_OUTPUT = os.path.join(config["out_dir"], "final", "perturbation_adata.h5ad")
FINAL_PROVENANCE = os.path.join(config["out_dir"], "final", "provenance.yaml")
FINAL_STATUS_MANIFEST = (
    os.path.join(
        config["out_dir"], "perturbation_status", FINAL_METHOD,
        "run_manifest.json",
    )
    if FINAL_METHOD != "none" else None
)
FINAL_STATUS_VALIDATION = (
    os.path.join(
        config["out_dir"], "perturbation_status", FINAL_METHOD,
        "validation.json",
    )
    if FINAL_METHOD != "none" else None
)


def _finalization_inputs(wildcards):
    inputs = {"adata": FINAL_INPUT}
    if FINAL_STATUS is not None:
        inputs["status"] = FINAL_STATUS
    return inputs


def _provenance_inputs(wildcards):
    inputs = {"adata": FINAL_OUTPUT}
    if FINAL_STATUS_MANIFEST is not None:
        inputs["status_manifest"] = FINAL_STATUS_MANIFEST
        inputs["status_validation"] = FINAL_STATUS_VALIDATION
    return inputs


def _json_b64(value):
    serialized = json.dumps(value, default=str, sort_keys=True).encode("utf-8")
    return base64.b64encode(serialized).decode("ascii")


PROVENANCE_GROUPS = {
    group: {
        "gex_h5": GROUPS[group]["gex_h5"],
        "sgRNA_r1": GROUPS[group]["_sgRNA_r1"],
        "sgRNA_r2": GROUPS[group]["_sgRNA_r2"],
    }
    for group in GROUPS
}


rule finalize_perturbation_adata:
    input:
        unpack(_finalization_inputs),
    output:
        adata = FINAL_OUTPUT,
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
        TEMP_OUTPUT="{output.adata}.partial"
        trap 'rm -f "$TEMP_OUTPUT"' EXIT
        python3 "{params.finalize_script}" \
            --input "{input.adata}" \
            --status-method "{params.method}" \
            {params.status_arg} \
            --output "$TEMP_OUTPUT" \
            --scan-chunk-values "{params.scan_chunk_values}" \
            --max-input-gb "{params.max_input_gb}" \
            --min-free-disk-gb "{params.min_free_disk_gb}" \
            > "{log}" 2>&1
        python3 "{params.validate_script}" \
            --source "{input.adata}" \
            --final "$TEMP_OUTPUT" \
            --status-method "{params.method}" \
            --chunk-values "{params.scan_chunk_values}" \
            >> "{log}" 2>&1
        mv "$TEMP_OUTPUT" "{output.adata}"
        trap - EXIT
        """


rule write_provenance:
    input:
        unpack(_provenance_inputs),
    output:
        yaml = FINAL_PROVENANCE,
    params:
        script = os.path.join(config["proj_dir"], "scripts", "write_provenance.py"),
        project_dir = config["proj_dir"],
        status_method = FINAL_METHOD,
        snakemake_version = PROVENANCE_SNAKEMAKE_VERSION,
        config_b64 = _json_b64(config),
        groups_b64 = _json_b64(PROVENANCE_GROUPS),
        status_args = (
            f'--status-manifest "{FINAL_STATUS_MANIFEST}" '
            f'--status-validation "{FINAL_STATUS_VALIDATION}"'
            if FINAL_STATUS_MANIFEST is not None else ""
        ),
    log:
        os.path.join(config["log_dir"], "finalization", "provenance.log"),
    threads: 1
    conda:
        os.path.join(config["proj_dir"], "envs", "finalization.yaml")
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname '{log}')" "$(dirname '{output.yaml}')"
        python3 "{params.script}" \
            --final-h5ad "{input.adata}" \
            --project-dir "{params.project_dir}" \
            --config-b64 "{params.config_b64}" \
            --groups-b64 "{params.groups_b64}" \
            --status-method "{params.status_method}" \
            --snakemake-version "{params.snakemake_version}" \
            {params.status_args} \
            --output "{output.yaml}" \
            > "{log}" 2>&1
        """
