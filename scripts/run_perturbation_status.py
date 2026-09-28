#!/usr/bin/env python3
"""Run one perturbation-status method against the integration H5AD contract."""

from __future__ import annotations

import argparse
import concurrent.futures
import gzip
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd
import yaml
import anndata as ad


ELIGIBLE_STRUCTURES = ["single_guide", "concordant_construct"]


def run(command: list[str]) -> None:
    subprocess.run(command, check=True)


def run_parallel(commands: list[list[str]], workers: int) -> None:
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(run, command) for command in commands]
        for future in futures:
            future.result()


def write_yaml(path: Path, value: dict) -> None:
    path.write_text(yaml.safe_dump(value, sort_keys=False), encoding="utf-8")


def boolean_series(values: pd.Series) -> pd.Series:
    if values.dtype == bool:
        return values
    result = values.astype(str).str.lower().map({"true": True, "false": False})
    if result.isna().any():
        raise ValueError("is_ntc contains invalid boolean values")
    return result.astype(bool)


def prepare_input(
    source: Path, work: Path, scripts: Path, dataset: str, feature_mode: str,
    method: str,
) -> tuple[Path, dict]:
    input_dir = work / "input"
    config = {
        "dataset": dataset,
        "method": method,
        "source_h5ad": str(source.resolve()),
        "input_dir": str(input_dir),
        "row_chunk_size": 512,
        "feature_mode": feature_mode,
        "selection": {
            "mode": "eligible",
            "eligible_assignment_structures": ELIGIBLE_STRUCTURES,
            "max_cells": 0,
        },
        "metadata_exclude_regex": r"^(guide_|guide_identity$|guide_ids$)",
    }
    config_path = work / "adapter.yaml"
    write_yaml(config_path, config)
    run([sys.executable, str(scripts / "adapt_input.py"), "--config", str(config_path)])
    adapter_manifest = json.loads((input_dir / "adapter_manifest.json").read_text())
    return input_dir, adapter_manifest


def run_mixscape(args, input_dir: Path, work: Path, scripts: Path) -> Path:
    metadata = pd.read_csv(input_dir / "metadata.tsv.gz", sep="\t", index_col=0)
    target_count = metadata.loc[~boolean_series(metadata["is_ntc"]), "target_label"].nunique()
    if target_count < 1:
        raise ValueError("Mixscape requires at least one eligible perturbed target")
    output = work / "method_output"
    config = {
        "module": "perturbation_status",
        "method": "Mixscape",
        "dataset": args.input.stem,
        "input_dir": str(input_dir),
        "output_dir": str(output),
        "seed": args.seed,
        "logging": {"directory": str(work / "resources"), "stage_log": str(work / "resources" / "stages.jsonl")},
        "checkpoints": {"enabled": True, "directory": str(work / "checkpoints")},
        "stop_after_checkpoint_b": True,
        "mixscape": {
            "input_assay": "RNA", "signature_assay": "PRTB", "signature_slot": "data",
            "classification_slot": "scale.data", "guide_class": "gene",
            "nt_cell_class": "non-targeting", "split_by": None, "num_neighbors": 20,
            "reduction": "pca", "labels": "gene", "nt_class_name": "non-targeting",
            "new_class_name": "mixscape_class", "min_de_genes": 10, "min_cells": 5,
            "de_assay": "RNA", "logfc_threshold": 0.25, "iter_num": 10,
            "prtb_type": args.perturbation_type, "nfeatures": 2000, "npcs": 30,
        },
    }
    config_path = work / "mixscape.yaml"
    write_yaml(config_path, config)
    run(["Rscript", str(scripts / "run_mixscape.R"), str(config_path)])
    shard_count = min(args.shards, int(target_count))
    shard_root = work / "mixscape_shards"
    shard_config = {
        "source_checkpoint_b": str(work / "checkpoints" / "checkpoint_B.rds"),
        "shard_root": str(shard_root), "dataset_prefix": args.input.stem,
        "source_h5ad": str(args.input.resolve()), "input_dir": str(input_dir),
        "n_shards": shard_count, "labels": "gene", "nt_class_name": "non-targeting",
        "selection": config["selection"] if "selection" in config else {
            "eligible_assignment_structures": ELIGIBLE_STRUCTURES
        },
        "seed": args.seed, "mixscape": config["mixscape"],
    }
    shard_config_path = work / "mixscape_shards.yaml"
    write_yaml(shard_config_path, shard_config)
    run(["Rscript", str(scripts / "build_mixscape_shards.R"), str(shard_config_path)])
    commands = [
        ["Rscript", str(scripts / "run_mixscape.R"), str(path)]
        for path in sorted((shard_root / "configs").glob("shard_*.yaml"))
    ]
    run_parallel(commands, min(args.max_workers, shard_count))
    run([
        sys.executable, str(scripts / "merge_mixscape_shards.py"),
        str(shard_root), str(input_dir / "metadata.tsv.gz"),
    ])
    return shard_root / "merged" / "mixscape_metadata.tsv.gz"


def run_ps(args, input_dir: Path, work: Path, scripts: Path) -> Path:
    metadata = pd.read_csv(input_dir / "metadata.tsv.gz", sep="\t", index_col=0)
    barcode = pd.DataFrame({
        "cell": metadata.index.astype(str),
        "barcode": metadata["guide_id"].astype(str),
        "gene": metadata["target_label"].astype(str),
    })
    barcode_path = input_dir / "barcode.tsv.gz"
    barcode.to_csv(barcode_path, sep="\t", index=False, compression="gzip")
    target_counts = barcode.loc[barcode["gene"] != "non-targeting", "gene"].value_counts()
    modeled_targets = int((target_counts >= args.min_target_cells).sum())
    shard_count = min(args.shards, modeled_targets)
    if shard_count < 1:
        raise ValueError("PS requires at least one target meeting min_target_cells")
    input_rds = work / "ps_input.rds"
    build_config = {
        "input_dir": str(input_dir), "seurat_rds": str(input_rds),
        "normalization": {"scale_factor": 10000},
    }
    build_path = work / "ps_build.yaml"
    write_yaml(build_path, build_config)
    run(["Rscript", str(scripts / "build_ps_seurat.R"), "--config", str(build_path)])
    ps_scripts = scripts / "ps"
    ps_work = work / "ps_shards"
    output = work / "method_output"
    ps_work.mkdir(parents=True, exist_ok=True)
    run([
        sys.executable, str(ps_scripts / "prepare_validation.py"), "--barcode",
        str(barcode_path), "--output-dir", str(ps_work), "--shards", str(shard_count),
        "--min-target-cells", str(args.min_target_cells),
    ])
    run([
        "Rscript", str(ps_scripts / "build_input_shards.R"), str(input_rds),
        str(barcode_path), str(ps_work), str(ps_work / "inputs"), str(shard_count),
    ])
    selection_commands = []
    for shard in range(1, shard_count + 1):
        name = f"shard_{shard:03d}"
        selection_commands.append([
            "Rscript", str(ps_scripts / "select_response_genes.R"),
            str(ps_work / "inputs" / f"{name}.rds"), str(barcode_path),
            str(ps_work / f"{name}_targets.txt"),
            str(ps_work / f"{name}_response_genes.tsv"), "non-targeting", "0.1",
        ])
    run_parallel(selection_commands, args.max_workers)
    manifest = json.loads((ps_work / "shard_manifest.json").read_text())
    run([
        sys.executable, str(ps_scripts / "merge_response_genes.py"),
        "--selection-dir", str(ps_work), "--output-dir", str(ps_work),
        "--expected-targets", str(manifest["targets"]),
    ])
    run([
        "Rscript", str(ps_scripts / "freeze_modeled_genes.R"), str(input_rds),
        str(ps_work / "global_response_genes.txt"), str(ps_work / "global_modeled_genes.txt"),
    ])
    score_commands = []
    for shard in range(1, shard_count + 1):
        name = f"shard_{shard:03d}"
        score_commands.append([
            "Rscript", str(ps_scripts / "run_shard.R"),
            str(ps_work / "inputs" / f"{name}.rds"), str(barcode_path),
            str(ps_work / f"{name}_targets.txt"), str(ps_work / "global_modeled_genes.txt"),
            str(ps_work / "scores" / name), str(args.seed + shard),
        ])
    run_parallel(score_commands, args.max_workers)
    run([
        sys.executable, str(ps_scripts / "merge_validate.py"), "--shard-dir",
        str(ps_work / "scores"), "--barcode", str(barcode_path),
        "--output-dir", str(output), "--excluded-targets", str(ps_work / "shard_manifest.json"),
    ])
    return output / "ps_metadata.tsv"


def write_skipped_status(method: str, metadata_path: Path, destination: Path, reason: str) -> dict:
    metadata = pd.read_csv(metadata_path, sep="\t", index_col=0)
    result = pd.DataFrame({
        "cell_id": metadata.index.astype(str),
        "target_label": metadata["target_label"].astype(str),
        "is_ntc": boolean_series(metadata["is_ntc"]),
        "assignment_structure": metadata["assignment_structure"].astype(str),
        "method": method,
        "score": pd.NA,
        "native_status": pd.NA,
        "scorable": False,
        "unscorable_reason": reason,
    })
    result.to_csv(destination, sep="\t", index=False, compression="gzip")
    return {"cells": len(result), "unique_cells": int(result.cell_id.nunique())}


def standardize(method: str, native_path: Path, metadata_path: Path, destination: Path) -> dict:
    frame = pd.read_csv(native_path, sep="\t", index_col=0 if method == "mixscape" else None)
    if method == "mixscape":
        if "cell_id" not in frame.columns:
            frame.index.name = "cell_id"
            frame = frame.reset_index()
        probability = next((c for c in frame.columns if c.startswith("mixscape_class_p_")), None)
        status = next((c for c in frame.columns if c.startswith("mixscape_class.global")), None)
        score = frame[probability] if probability else pd.Series(pd.NA, index=frame.index)
        label = frame[status] if status else pd.Series(pd.NA, index=frame.index)
        scorable = pd.Series(True, index=frame.index)
        reason = pd.Series(pd.NA, index=frame.index, dtype="object")
    else:
        score = frame["ps_score"]
        label = frame["ps_status"]
        scorable = frame["ps_status"] != "insufficient_target_cells"
        reason = frame["ps_status"].where(~scorable, pd.NA)
    canonical = pd.read_csv(metadata_path, sep="\t", index_col=0)
    structures = canonical["assignment_structure"].reindex(frame["cell_id"].astype(str))
    if structures.isna().any():
        raise ValueError("Method output contains cells absent from canonical metadata")
    ntc = frame["is_ntc"]
    if ntc.dtype != bool:
        ntc = ntc.astype(str).str.lower().map({"true": True, "false": False})
    if ntc.isna().any():
        raise ValueError("Method output contains invalid is_ntc values")
    result = pd.DataFrame({
        "cell_id": frame["cell_id"].astype(str),
        "target_label": frame["target_label"].astype(str),
        "is_ntc": ntc.astype(bool),
        "assignment_structure": structures.to_numpy(),
        "method": method,
        "score": score,
        "native_status": label,
        "scorable": scorable,
        "unscorable_reason": reason,
    })
    if result["cell_id"].duplicated().any():
        raise ValueError("Method output contains duplicate cell IDs")
    result.to_csv(destination, sep="\t", index=False, compression="gzip")
    return {"cells": len(result), "unique_cells": int(result.cell_id.nunique())}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--method", choices=("mixscape", "ps"), required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--perturbation-type", default="")
    parser.add_argument("--shards", type=int, default=4)
    parser.add_argument("--max-workers", type=int, default=1)
    parser.add_argument("--min-target-cells", type=int, default=3)
    parser.add_argument("--seed", type=int, default=20260902)
    parser.add_argument("--keep-work", action="store_true")
    parser.add_argument(
        "--feature-mode", choices=("auto", "gene", "usa_sa"), default="auto"
    )
    args = parser.parse_args()
    if args.method == "mixscape" and args.perturbation_type not in {"KO", "CRISPRa", "CRISPRi"}:
        parser.error("Mixscape requires --perturbation-type KO, CRISPRa, or CRISPRi")
    if not args.input.is_file():
        raise FileNotFoundError(args.input)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    work = args.output_dir / "work"
    work.mkdir(exist_ok=True)
    scripts = Path(__file__).resolve().parent / "perturbation_status"
    started = time.time()
    backed = ad.read_h5ad(args.input, backed="r")
    structures = backed.obs["assignment_structure"].astype(str)
    eligible_count = int(structures.isin(ELIGIBLE_STRUCTURES).sum())
    backed.file.close()
    destination = args.output_dir / "perturbation_status.tsv.gz"
    if eligible_count == 0:
        columns = [
            "cell_id", "target_label", "is_ntc", "assignment_structure",
            "method", "score", "native_status", "scorable", "unscorable_reason",
        ]
        pd.DataFrame(columns=columns).to_csv(
            destination, sep="\t", index=False, compression="gzip"
        )
        reason = "no_status_eligible_cells"
        print(
            f"[skip] method={args.method}; eligible_cells=0; runnable_targets=0; "
            f"reason={reason}", flush=True,
        )
        validation = {
            "status": "PASS", "method": args.method, "skipped": True,
            "cells": 0, "unique_cells": 0,
        }
        (args.output_dir / "validation.json").write_text(
            json.dumps(validation, indent=2) + "\n"
        )
        manifest = {
            "status": "PASS", "method": args.method,
            "input": str(args.input.resolve()),
            "eligible_assignment_structures": ELIGIBLE_STRUCTURES,
            "eligible_cells": 0, "perturbed_targets": 0, "runnable_targets": 0,
            "skipped": True, "skip_reason": reason,
            "feature_mode_requested": args.feature_mode,
            "feature_mode_resolved": "not_evaluated", "pse_feature_count": 0,
            "wall_seconds": time.time() - started,
        }
        (args.output_dir / "run_manifest.json").write_text(
            json.dumps(manifest, indent=2) + "\n"
        )
        if not args.keep_work:
            shutil.rmtree(work)
        return
    input_dir, adapter_manifest = prepare_input(
        args.input, work, scripts, args.input.stem, args.feature_mode, args.method
    )
    metadata = pd.read_csv(input_dir / "metadata.tsv.gz", sep="\t", index_col=0)
    non_ntc = metadata.loc[~boolean_series(metadata["is_ntc"]), "target_label"].value_counts()
    if args.method == "ps":
        runnable_targets = int((non_ntc >= args.min_target_cells).sum())
        skip_reason = "no_target_meets_min_target_cells"
    else:
        runnable_targets = int(non_ntc.size)
        skip_reason = "no_eligible_perturbed_target"
    skipped = runnable_targets == 0
    if skipped:
        print(
            f"[skip] method={args.method}; eligible_cells={len(metadata)}; "
            f"perturbed_targets={int(non_ntc.size)}; runnable_targets=0; "
            f"reason={skip_reason}", flush=True,
        )
        summary = write_skipped_status(
            args.method, input_dir / "metadata.tsv.gz", destination, skip_reason
        )
    else:
        native = run_mixscape(args, input_dir, work, scripts) if args.method == "mixscape" else run_ps(args, input_dir, work, scripts)
        summary = standardize(args.method, native, input_dir / "metadata.tsv.gz", destination)
    validation = {"status": "PASS", "method": args.method, "skipped": skipped, **summary}
    (args.output_dir / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    manifest = {
        "status": "PASS", "method": args.method, "input": str(args.input.resolve()),
        "eligible_assignment_structures": ELIGIBLE_STRUCTURES,
        "min_target_cells": args.min_target_cells, "shards_requested": args.shards,
        "max_workers_requested": args.max_workers, "wall_seconds": time.time() - started,
        "feature_mode_requested": args.feature_mode,
        "feature_mode_resolved": adapter_manifest["feature_mode_resolved"],
        "pse_feature_count": adapter_manifest["output_shape"][0],
        "eligible_cells": int(len(metadata)),
        "perturbed_targets": int(non_ntc.size),
        "runnable_targets": runnable_targets,
        "skipped": skipped,
        "skip_reason": skip_reason if skipped else None,
    }
    (args.output_dir / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    if not args.keep_work:
        shutil.rmtree(work)


if __name__ == "__main__":
    main()
