#!/usr/bin/env python3
"""Validate and finalize a perturbation-screen AnnData without changing expression."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


ELIGIBLE_STRUCTURES = {"single_guide", "concordant_construct"}
VALID_STRUCTURES = ELIGIBLE_STRUCTURES | {"mixed_construct"}
STATUS_COLUMNS = [
    "perturbation_status_method",
    "perturbation_status_score",
    "perturbation_status_label",
    "perturbation_status_scorable",
    "perturbation_status_reason",
]
REQUIRED_OBS = {
    "guide_id", "target_label", "perturbation_group", "is_ntc", "batch_id",
    "assignment_structure",
}


def decode(value):
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, np.generic):
        return value.item()
    return value


def read_array(dataset: h5py.Dataset) -> np.ndarray:
    if h5py.check_dtype(vlen=dataset.dtype) is not None or dataset.dtype.kind == "O":
        return np.asarray(dataset.asstr()[:], dtype=object)
    values = dataset[:]
    if values.dtype.kind == "S":
        return np.asarray([x.decode("utf-8") for x in values], dtype=object)
    return values


def read_column(item: h5py.Dataset | h5py.Group) -> np.ndarray:
    if isinstance(item, h5py.Group) and decode(item.attrs.get("encoding-type", "")) == "categorical":
        categories = np.asarray([decode(x) for x in read_array(item["categories"])], dtype=object)
        codes = np.asarray(item["codes"][:], dtype=np.int64)
        return np.asarray([None if code < 0 else categories[code] for code in codes], dtype=object)
    return np.asarray([decode(x) for x in read_array(item)], dtype=object)


def read_index(group: h5py.Group) -> np.ndarray:
    name = str(decode(group.attrs["_index"]))
    return np.asarray([str(x) for x in read_array(group[name])], dtype=object)


def matrix_shape(item: h5py.Dataset | h5py.Group) -> tuple[int, int]:
    if isinstance(item, h5py.Group):
        return tuple(map(int, item.attrs["shape"]))
    return tuple(map(int, item.shape))


def matrix_data(item: h5py.Dataset | h5py.Group) -> h5py.Dataset:
    return item["data"] if isinstance(item, h5py.Group) else item


def scan_values(dataset: h5py.Dataset, *, counts: bool, chunk: int) -> dict:
    total = int(dataset.size)
    nonfinite = negative = noninteger = explicit_zero = 0
    minimum = maximum = None
    row_width = max(1, int(np.prod(dataset.shape[1:]))) if dataset.ndim > 1 else 1
    step = max(1, chunk // row_width)
    first_dimension = dataset.shape[0] if dataset.ndim else 1
    for start in range(0, first_dimension, step):
        values = np.asarray(dataset[start : min(start + step, first_dimension)])
        flat = values.reshape(-1)
        finite = np.isfinite(flat)
        nonfinite += int((~finite).sum())
        valid = flat[finite]
        if valid.size:
            value_min, value_max = float(valid.min()), float(valid.max())
            minimum = value_min if minimum is None else min(minimum, value_min)
            maximum = value_max if maximum is None else max(maximum, value_max)
            negative += int((valid < 0).sum())
            explicit_zero += int((valid == 0).sum())
            if counts:
                noninteger += int((valid != np.floor(valid)).sum())
    return {
        "stored_values": total, "nonfinite": nonfinite, "negative": negative,
        "noninteger": noninteger, "explicit_zero": explicit_zero,
        "minimum": minimum, "maximum": maximum,
    }


def validate_input(path: Path, chunk: int) -> tuple[dict, dict[str, np.ndarray]]:
    with h5py.File(path, "r") as handle:
        required_root = {"X", "layers", "obs", "var", "uns"}
        missing_root = required_root.difference(handle.keys())
        if missing_root or "counts" not in handle["layers"]:
            raise ValueError(f"Missing H5AD components: {sorted(missing_root)}; counts layer required")
        obs = handle["obs"]
        missing_obs = REQUIRED_OBS.difference(obs.keys())
        if missing_obs:
            raise ValueError(f"Missing required obs fields: {sorted(missing_obs)}")
        x_shape = matrix_shape(handle["X"])
        counts_shape = matrix_shape(handle["layers/counts"])
        if x_shape != counts_shape:
            raise ValueError(f"X/counts shape mismatch: {x_shape} != {counts_shape}")
        cells, genes = read_index(obs), read_index(handle["var"])
        if len(cells) != x_shape[0] or len(genes) != x_shape[1]:
            raise ValueError("Matrix dimensions do not match obs/var indices")
        if len(set(cells)) != len(cells) or len(set(genes)) != len(genes):
            raise ValueError("Cell and gene identifiers must be unique")
        structures = read_column(obs["assignment_structure"])
        observed_structures = {str(x) for x in structures if x is not None}
        invalid = observed_structures.difference(VALID_STRUCTURES)
        if invalid:
            raise ValueError(f"Invalid assignment_structure values: {sorted(invalid)}")
        labels = read_column(obs["target_label"])
        ntc_raw = read_column(obs["is_ntc"])
        ntc = np.asarray([x is True or str(x).lower() == "true" for x in ntc_raw], dtype=bool)
        expected_ntc = np.asarray([x == "non-targeting" for x in labels], dtype=bool)
        if not np.array_equal(ntc, expected_ntc):
            raise ValueError("is_ntc is inconsistent with canonical target_label")
        x_scan = scan_values(matrix_data(handle["X"]), counts=False, chunk=chunk)
        counts_scan = scan_values(matrix_data(handle["layers/counts"]), counts=True, chunk=chunk)
        if x_scan["nonfinite"] or x_scan["negative"]:
            raise ValueError("X contains non-finite or negative values")
        if counts_scan["nonfinite"] or counts_scan["negative"] or counts_scan["noninteger"]:
            raise ValueError("layers['counts'] must contain finite, non-negative integer values")
        eligible = np.isin(structures, sorted(ELIGIBLE_STRUCTURES))
        report = {
            "shape": list(x_shape), "cells": len(cells), "genes": len(genes),
            "x": x_scan, "counts": counts_scan,
            "assignment_structure": {str(k): int(v) for k, v in zip(*np.unique(structures.astype(str), return_counts=True))},
            "eligible_cells": int(eligible.sum()), "noneligible_cells": int((~eligible).sum()),
            "ntc_cells": int(ntc.sum()),
        }
        arrays = {"cells": cells, "structures": structures, "eligible": eligible, "labels": labels, "is_ntc": ntc}
        return report, arrays


def build_status(method: str, status_path: Path | None, arrays: dict[str, np.ndarray]) -> tuple[pd.DataFrame, dict]:
    cells = pd.Index(arrays["cells"].astype(str), name="cell_id")
    eligible = pd.Series(arrays["eligible"], index=cells)
    final = pd.DataFrame(index=cells)
    final["perturbation_status_method"] = method
    final["perturbation_status_score"] = np.nan
    final["perturbation_status_label"] = pd.Series(pd.NA, index=cells, dtype="object")
    final["perturbation_status_scorable"] = False
    final["perturbation_status_reason"] = "not_run" if method == "none" else "not_eligible"
    if method == "none":
        return final, {"sidecar_rows": 0, "eligible_rows": int(eligible.sum()), "coverage": "not_run"}
    if status_path is None or not status_path.is_file():
        raise ValueError(f"Status table is required for method={method}")
    status = pd.read_csv(status_path, sep="\t")
    required = {
        "cell_id", "target_label", "is_ntc", "assignment_structure", "score",
        "native_status", "scorable", "unscorable_reason", "method",
    }
    missing = required.difference(status.columns)
    if missing:
        raise ValueError(f"Status table lacks columns: {sorted(missing)}")
    status["cell_id"] = status["cell_id"].astype(str)
    if status["cell_id"].duplicated().any():
        raise ValueError("Status table contains duplicate cell IDs")
    known = set(cells)
    unexpected = set(status["cell_id"]).difference(known)
    if unexpected:
        raise ValueError(f"Status table contains {len(unexpected)} cells absent from AnnData")
    eligible_ids = set(cells[eligible.to_numpy()])
    status_ids = set(status["cell_id"])
    missing_eligible = eligible_ids.difference(status_ids)
    noneligible_status = status_ids.difference(eligible_ids)
    if missing_eligible or noneligible_status:
        raise ValueError(
            f"Status coverage mismatch: missing eligible={len(missing_eligible)}, "
            f"included noneligible={len(noneligible_status)}"
        )
    methods = {str(x).lower() for x in status["method"].dropna().unique()}
    if methods != {method}:
        raise ValueError(f"Status method mismatch: expected {method}, observed {sorted(methods)}")
    status = status.set_index("cell_id").reindex(cells[eligible.to_numpy()])
    expected_positions = np.flatnonzero(arrays["eligible"])
    expected_labels = np.asarray(arrays["labels"], dtype=object)[expected_positions].astype(str)
    expected_structures = np.asarray(arrays["structures"], dtype=object)[expected_positions].astype(str)
    expected_ntc = np.asarray(arrays["is_ntc"], dtype=bool)[expected_positions]
    observed_ntc = status["is_ntc"]
    if observed_ntc.dtype != bool:
        observed_ntc = observed_ntc.astype(str).str.lower().map({"true": True, "false": False})
    if (
        not np.array_equal(status["target_label"].astype(str).to_numpy(), expected_labels)
        or not np.array_equal(status["assignment_structure"].astype(str).to_numpy(), expected_structures)
        or observed_ntc.isna().any()
        or not np.array_equal(observed_ntc.to_numpy(dtype=bool), expected_ntc)
    ):
        raise ValueError("Status table metadata is inconsistent with the input AnnData")
    scorable = status["scorable"]
    if scorable.dtype != bool:
        scorable = scorable.astype(str).str.lower().map({"true": True, "false": False})
    if scorable.isna().any():
        raise ValueError("Status scorable field contains invalid values")
    numeric_score = pd.to_numeric(status["score"], errors="coerce")
    if (~np.isfinite(numeric_score.loc[scorable])).any():
        raise ValueError("Scorable cells require finite status scores")
    if status.loc[scorable, "native_status"].isna().any():
        raise ValueError("Scorable cells require native_status")
    ids = status.index
    final.loc[ids, "perturbation_status_score"] = numeric_score.to_numpy()
    final.loc[ids, "perturbation_status_label"] = status["native_status"].to_numpy()
    final.loc[ids, "perturbation_status_scorable"] = scorable.to_numpy(dtype=bool)
    reasons = status["unscorable_reason"].astype("object")
    reasons.loc[scorable.to_numpy()] = pd.NA
    if reasons.loc[~scorable.to_numpy()].isna().any():
        raise ValueError("Unscorable eligible cells require an unscorable_reason")
    final.loc[ids, "perturbation_status_reason"] = reasons.to_numpy()
    return final, {
        "sidecar_rows": len(status), "eligible_rows": len(eligible_ids),
        "scorable_rows": int(scorable.sum()), "unscorable_rows": int((~scorable).sum()),
        "coverage": "exact eligible-cell coverage",
    }


def replace_item(group: h5py.Group, name: str) -> None:
    if name in group:
        del group[name]


def write_categorical(group: h5py.Group, name: str, values: pd.Series) -> None:
    replace_item(group, name)
    item = group.create_group(name)
    item.attrs["encoding-type"] = "categorical"
    item.attrs["encoding-version"] = "0.2.0"
    item.attrs["ordered"] = False
    nonmissing = [str(x) for x in values.dropna().unique()]
    categories = sorted(nonmissing)
    mapping = {value: index for index, value in enumerate(categories)}
    codes = np.asarray([-1 if pd.isna(x) else mapping[str(x)] for x in values], dtype=np.int32)
    string_dtype = h5py.string_dtype("utf-8")
    category_ds = item.create_dataset("categories", data=np.asarray(categories, dtype=object), dtype=string_dtype)
    category_ds.attrs["encoding-type"] = "string-array"
    category_ds.attrs["encoding-version"] = "0.2.0"
    code_ds = item.create_dataset("codes", data=codes, compression="gzip", compression_opts=4)
    code_ds.attrs["encoding-type"] = "array"
    code_ds.attrs["encoding-version"] = "0.2.0"


def write_array(group: h5py.Group, name: str, values: np.ndarray) -> None:
    replace_item(group, name)
    item = group.create_dataset(name, data=values, compression="gzip", compression_opts=4)
    item.attrs["encoding-type"] = "array"
    item.attrs["encoding-version"] = "0.2.0"


def write_scalar(group: h5py.Group, name: str, value) -> None:
    replace_item(group, name)
    if isinstance(value, str):
        item = group.create_dataset(name, data=value, dtype=h5py.string_dtype("utf-8"))
        item.attrs["encoding-type"] = "string"
    else:
        item = group.create_dataset(name, data=value)
        item.attrs["encoding-type"] = "numeric-scalar"
    item.attrs["encoding-version"] = "0.2.0"


def write_dict(group: h5py.Group, name: str, values: dict) -> None:
    replace_item(group, name)
    item = group.create_group(name)
    item.attrs["encoding-type"] = "dict"
    item.attrs["encoding-version"] = "0.1.0"
    for key, value in values.items():
        write_scalar(item, key, value)


def copy_and_annotate(source: Path, destination: Path, status: pd.DataFrame, method: str, input_report: dict) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + f".partial-{os.getpid()}")
    if temporary.exists():
        temporary.unlink()
    try:
        shutil.copyfile(source, temporary)
        with h5py.File(temporary, "r+") as handle:
            obs = handle["obs"]
            existing_order = [str(decode(x)) for x in obs.attrs["column-order"]]
            for column in STATUS_COLUMNS:
                if column in existing_order:
                    existing_order.remove(column)
            existing_order.extend(STATUS_COLUMNS)
            del obs.attrs["column-order"]
            obs.attrs.create(
                "column-order", np.asarray(existing_order, dtype=object),
                dtype=h5py.string_dtype("utf-8"),
            )
            write_categorical(obs, "perturbation_status_method", status["perturbation_status_method"])
            write_array(obs, "perturbation_status_score", status["perturbation_status_score"].to_numpy(dtype=np.float64))
            write_categorical(obs, "perturbation_status_label", status["perturbation_status_label"])
            write_array(obs, "perturbation_status_scorable", status["perturbation_status_scorable"].to_numpy(dtype=bool))
            write_categorical(obs, "perturbation_status_reason", status["perturbation_status_reason"])
            expression_contract = {
                "representation": "upstream normalized and log-transformed expression",
                "finalization_action": "validated and preserved; not recomputed or rescaled",
                "counts_layer": "counts",
                "normalization_provenance": "inherited from multimodal integration",
            }
            finalization = {
                "status_method": method,
                "eligible_assignment_structures": ",".join(sorted(ELIGIBLE_STRUCTURES)),
                "cell_filtering": "not_performed",
                "gene_filtering": "not_performed",
                "source_h5ad": str(source.resolve()),
            }
            write_dict(handle["uns"], "expression_contract", expression_contract)
            write_dict(handle["uns"], "finalization", finalization)
            handle.flush()
        os.replace(temporary, destination)
    finally:
        if temporary.exists():
            temporary.unlink()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--status-method", choices=("none", "mixscape", "ps"), default="none")
    parser.add_argument("--status-table", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--scan-chunk-values", type=int, default=10_000_000)
    parser.add_argument("--max-input-gb", type=float, default=300.0)
    parser.add_argument("--min-free-disk-gb", type=float, default=100.0)
    args = parser.parse_args()
    started = time.time()
    if args.output.resolve() == args.input.resolve():
        parser.error("Output must differ from input")
    input_gb = args.input.stat().st_size / (1024 ** 3)
    if input_gb > args.max_input_gb:
        raise RuntimeError(f"Input size {input_gb:.2f} GiB exceeds limit {args.max_input_gb:.2f} GiB")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    free_gb = shutil.disk_usage(args.output.parent).free / (1024 ** 3)
    required_gb = input_gb * 1.05 + args.min_free_disk_gb
    if free_gb < required_gb:
        raise RuntimeError(
            f"Insufficient disk: {free_gb:.2f} GiB free; {required_gb:.2f} GiB required"
        )
    input_report, arrays = validate_input(args.input, args.scan_chunk_values)
    status, status_report = build_status(args.status_method, args.status_table, arrays)
    copy_and_annotate(args.input, args.output, status, args.status_method, input_report)
    report = {
        "status": "PASS", "input": str(args.input.resolve()), "output": str(args.output.resolve()),
        "status_method": args.status_method, "input_validation": input_report,
        "status_validation": status_report, "expression_changed": False,
        "cell_filtering": "not_performed", "gene_filtering": "not_performed",
        "wall_seconds": time.time() - started,
        "resource_preflight": {
            "input_gb": input_gb, "free_disk_gb": free_gb,
            "required_disk_gb": required_gb, "max_input_gb": args.max_input_gb,
            "min_free_disk_gb": args.min_free_disk_gb,
        },
    }
    if args.report is not None:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        temporary_report = args.report.with_name(args.report.name + f".partial-{os.getpid()}")
        temporary_report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        os.replace(temporary_report, args.report)


if __name__ == "__main__":
    main()
