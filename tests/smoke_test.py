#!/usr/bin/env python3
"""Run a tiny, real FASTQ-to-final-H5AD workflow smoke test."""

from __future__ import annotations

import argparse
import subprocess
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import yaml
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
BARCODES = (
    "AAAAAAAAAAAAAAAA",
    "CCCCCCCCCCCCCCCC",
    "GGGGGGGGGGGGGGGG",
    "TTTTTTTTTTTTTTTT",
)
GUIDES = {
    "guide_A": "ACGTACGTACGTACGTACGT",
    "guide_B": "TGCATGCATGCATGCATGCA",
}


def _umi(number: int) -> str:
    alphabet = "ACGT"
    bases = []
    for _ in range(10):
        bases.append(alphabet[number % 4])
        number //= 4
    return "".join(bases)


def _write_fastqs(directory: Path) -> tuple[Path, Path]:
    r1 = directory / "smoke_R1_001.fastq"
    r2 = directory / "smoke_R2_001.fastq"
    r1_records = []
    r2_records = []
    read_number = 0
    for cell_index, barcode in enumerate(BARCODES):
        guide = GUIDES["guide_A" if cell_index < 2 else "guide_B"]
        for umi_index in range(6 + (cell_index % 2)):
            name = f"@smoke_{read_number}"
            read1 = barcode + _umi(cell_index * 32 + umi_index)
            r1_records.extend((name, read1, "+", "I" * len(read1)))
            r2_records.extend((name, guide, "+", "I" * len(guide)))
            read_number += 1
    r1.write_text("\n".join(r1_records) + "\n", encoding="utf-8")
    r2.write_text("\n".join(r2_records) + "\n", encoding="utf-8")
    return r1, r2


def _write_inputs(workdir: Path) -> Path:
    inputs = workdir / "inputs"
    outputs = workdir / "results"
    logs = workdir / "logs"
    inputs.mkdir(parents=True, exist_ok=True)
    _write_fastqs(inputs)

    counts = sparse.csr_matrix(
        np.asarray(
            [[10, 0, 2], [8, 1, 0], [0, 9, 2], [1, 7, 1]],
            dtype=np.int32,
        )
    )
    expression = ad.AnnData(
        X=counts,
        obs=pd.DataFrame(index=pd.Index(BARCODES, name="cell_id")),
        var=pd.DataFrame(index=pd.Index(["gene_A", "gene_B", "gene_C"])),
    )
    expression.write_h5ad(inputs / "expression.h5ad")

    pd.DataFrame(
        {
            "construct_id": ["construct_A", "construct_B"],
            "guide_id": list(GUIDES),
            "guide_sequence": list(GUIDES.values()),
            "target_label": ["gene_A", "gene_B"],
            "guide_position": [1, 1],
        }
    ).to_csv(inputs / "guide_library.csv", index=False)

    groups = {
        "groups": {
            "lane_01": {
                "group_id": "lane_01",
                "gex_h5": str(inputs / "expression.h5ad"),
                "sgRNA_fastq_dir": str(inputs),
                "sgRNA_r1_pattern": "*_R1_001.fastq",
                "sgRNA_r2_pattern": "*_R2_001.fastq",
            }
        }
    }
    groups_path = inputs / "groups.yaml"
    groups_path.write_text(yaml.safe_dump(groups, sort_keys=False), encoding="utf-8")

    config = {
        "proj_dir": str(ROOT),
        "out_dir": str(outputs),
        "log_dir": str(logs),
        "groups_file": str(groups_path),
        "tenx_chemistry": "5v1",
        "guide_csv": str(inputs / "guide_library.csv"),
        "guide_extraction": {"method": "simpleaf"},
        "simpleaf": {"af_home": ""},
        "whitelist": {"min_umi": 0, "min_genes": 0},
        "resources": {"simpleaf_quant_threads": 2, "simpleaf_index_threads": 2},
        "assignment": {
            "guide_design": "dual",
            "methods": ["umi_threshold"],
            "umi_threshold": {"threshold": 2},
        },
        "integration": {
            "input_kind": "auto",
            "cell_id_matching": "auto",
            "counts_layer": "counts",
            "target_sum": 10000,
            "min_free_disk_gb": 0,
            "max_process_memory_gb": 4,
        },
        "perturbation_status": {"method": "none", "feature_mode": "auto"},
        "finalization": {
            "memory_gb": 2,
            "scan_chunk_values": 1000,
            "max_input_gb": 1,
            "min_free_disk_gb": 0,
        },
    }
    config_path = inputs / "config.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    return config_path


def _validate(workdir: Path) -> None:
    final = workdir / "results" / "final" / "perturbation_adata.h5ad"
    provenance = workdir / "results" / "final" / "provenance.yaml"
    result = ad.read_h5ad(final)
    if result.n_obs != 4 or result.n_vars != 3:
        raise AssertionError(f"unexpected final shape: {result.shape}")
    if "counts" not in result.layers:
        raise AssertionError("final H5AD has no counts layer")
    required = {"assignment_structure", "perturbation_status_method"}
    if not required.issubset(result.obs.columns):
        raise AssertionError("final H5AD is missing required observation fields")
    record = yaml.safe_load(provenance.read_text(encoding="utf-8"))
    if record["workflow"].get("version") != "1.0.0":
        raise AssertionError("provenance does not record workflow version 1.0.0")
    if record["validation"].get("final_h5ad") != "PASS":
        raise AssertionError("final H5AD validation did not pass")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", type=Path, help="retain outputs in this directory")
    parser.add_argument("--cores", type=int, default=2)
    args = parser.parse_args()

    temporary = None
    if args.workdir is None:
        temporary = tempfile.TemporaryDirectory(prefix="scprocess-perturb-smoke-")
        workdir = Path(temporary.name)
    else:
        workdir = args.workdir.resolve()
        workdir.mkdir(parents=True, exist_ok=True)

    config = _write_inputs(workdir)
    subprocess.run(
        [
            "snakemake", "--use-conda", "--conda-frontend", "conda",
            "--cores", str(args.cores), "--configfile", str(config),
            "--rerun-incomplete",
        ],
        cwd=ROOT,
        check=True,
    )
    _validate(workdir)
    print(f"PASS: complete smoke workflow validated at {workdir}")
    if temporary is not None:
        temporary.cleanup()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
