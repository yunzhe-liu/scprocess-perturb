# Perturbation-status estimation. This stage is conditional and does not alter
# the integration artifact or any upstream rule.

STATUS_CONFIG = config.get("perturbation_status", {})
STATUS_METHOD = str(STATUS_CONFIG.get("method", "none")).lower()
STATUS_DIR = os.path.join(config["out_dir"], "perturbation_status", STATUS_METHOD)
STATUS_INPUT = os.path.join(config["out_dir"], "integration", "perturbation_adata.h5ad")


rule perturbation_status_estimation:
    input:
        adata = STATUS_INPUT,
    output:
        status = os.path.join(STATUS_DIR, "perturbation_status.tsv.gz"),
        validation = os.path.join(STATUS_DIR, "validation.json"),
        manifest = os.path.join(STATUS_DIR, "run_manifest.json"),
    params:
        method = STATUS_METHOD,
        perturbation_type = STATUS_CONFIG.get("perturbation_type", ""),
        shards = int(STATUS_CONFIG.get("shards", 4)),
        max_workers = int(STATUS_CONFIG.get("max_workers", 1)),
        min_target_cells = int(STATUS_CONFIG.get("min_target_cells", 3)),
        seed = int(STATUS_CONFIG.get("seed", 20260902)),
        script = os.path.join(
            config["proj_dir"], "scripts", "run_perturbation_status.py"
        ),
    log:
        os.path.join(config["log_dir"], "perturbation_status", f"{STATUS_METHOD}.log"),
    threads:
        lambda wildcards: max(1, int(STATUS_CONFIG.get("max_workers", 1)))
    resources:
        mem_mb = int(STATUS_CONFIG.get("memory_gb", 64)) * 1024,
    conda:
        os.path.join(
            config["proj_dir"], "envs",
            "mixscape.yaml" if STATUS_METHOD == "mixscape" else "ps.yaml"
        )
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname '{log}')" "{STATUS_DIR}"
        python3 "{params.script}" \
            --method "{params.method}" \
            --input "{input.adata}" \
            --output-dir "{STATUS_DIR}" \
            --perturbation-type "{params.perturbation_type}" \
            --shards "{params.shards}" \
            --max-workers "{params.max_workers}" \
            --min-target-cells "{params.min_target_cells}" \
            --seed "{params.seed}" \
            > "{log}" 2>&1
        """
