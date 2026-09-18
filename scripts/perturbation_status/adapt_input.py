#!/usr/bin/env python3
"""Convert a standardized H5AD file into a deterministic Seurat input directory."""

from __future__ import annotations

import argparse
import gzip
import json
import re
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import yaml


def decode_scalar(value: object) -> object:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, np.generic):
        return value.item()
    return value


def as_bool(value: object) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    return str(value).lower() == "true"


def read_h5_array(obj: h5py.Dataset) -> np.ndarray:
    if h5py.check_dtype(vlen=obj.dtype) is not None or obj.dtype.kind == "O":
        return np.asarray(obj.asstr()[:], dtype=object)
    values = obj[:]
    if values.dtype.kind == "S":
        return np.asarray([value.decode("utf-8") for value in values], dtype=object)
    return values


def read_obs_column(obj: h5py.Dataset | h5py.Group) -> list[object]:
    if isinstance(obj, h5py.Group) and obj.attrs.get("encoding-type", "") == "categorical":
        categories = [decode_scalar(value) for value in read_h5_array(obj["categories"])]
        codes = np.asarray(obj["codes"][:], dtype=np.int64)
        return [None if code < 0 else categories[code] for code in codes]
    return [decode_scalar(value) for value in read_h5_array(obj)]


def read_index(group: h5py.Group, attr_name: str) -> list[str]:
    index_name = decode_scalar(group.attrs[attr_name])
    return [str(value) for value in read_h5_array(group[index_name])]


def choose_cells(
    labels: list[object],
    is_ntc: list[object],
    assignment_structure: list[object],
    selection: dict[str, object],
) -> np.ndarray:
    label_array = np.asarray(labels, dtype=object)
    structure_array = np.asarray(assignment_structure, dtype=object)
    allowed_structures = selection.get(
        "eligible_assignment_structures",
        ["single_guide", "concordant_construct"],
    )
    eligible_mask = np.isin(structure_array, list(allowed_structures)) & pd.notna(label_array)
    eligible = np.flatnonzero(eligible_mask)
    mode = str(selection.get("mode", "eligible"))
    if mode in {"eligible", "all_labeled"}:
        selected = eligible.tolist()
    elif mode == "per_label":
        cells_per_label = int(selection["cells_per_label"])
        by_label: dict[str, list[int]] = {}
        for index in eligible:
            label = str(label_array[index])
            if len(by_label.setdefault(label, [])) < cells_per_label:
                by_label[label].append(int(index))
        labels = sorted(by_label)
        max_labels = int(selection.get("max_labels", 0))
        if max_labels > 0 and len(labels) > max_labels:
            control = "non-targeting"
            retained = [control] if control in by_label else []
            retained.extend(label for label in labels if label != control)
            labels = retained[:max_labels]
        selected = [index for label in labels for index in by_label[label]]
    else:
        raise ValueError(f"Unsupported cell selection mode: {mode}")

    max_cells = int(selection.get("max_cells", 0))
    if max_cells > 0 and len(selected) > max_cells:
        ntc = [index for index in selected if is_ntc[index] is True]
        non_ntc = [index for index in selected if is_ntc[index] is not True]
        selected = (ntc + non_ntc)[:max_cells]
    return np.asarray(sorted(selected), dtype=np.int64)


def write_matrix_market(
    h5_path: Path,
    output_path: Path,
    row_indices: np.ndarray,
    n_features: int,
    row_chunk_size: int,
) -> int:
    with h5py.File(h5_path, "r") as handle:
        matrix = handle["layers/counts"]
        indptr = np.asarray(matrix["indptr"][:], dtype=np.int64)
        starts = indptr[row_indices]
        ends = indptr[row_indices + 1]
        nnz = int(np.sum(ends - starts))
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(output_path, "wt", encoding="utf-8", compresslevel=1) as stream:
            stream.write("%%MatrixMarket matrix coordinate integer general\n")
            stream.write("% Generated from the standardized H5AD counts layer.\n")
            stream.write(f"{n_features} {len(row_indices)} {nnz}\n")
            indices = matrix["indices"]
            data = matrix["data"]
            for chunk_start in range(0, len(row_indices), row_chunk_size):
                chunk_rows = row_indices[chunk_start : chunk_start + row_chunk_size]
                run_start = 0
                while run_start < len(chunk_rows):
                    run_end = run_start + 1
                    while run_end < len(chunk_rows) and chunk_rows[run_end] == chunk_rows[run_end - 1] + 1:
                        run_end += 1
                    first_row = int(chunk_rows[run_start])
                    last_row = int(chunk_rows[run_end - 1])
                    data_start = int(indptr[first_row])
                    data_end = int(indptr[last_row + 1])
                    run_indices = np.asarray(indices[data_start:data_end], dtype=np.int64)
                    run_data = np.asarray(data[data_start:data_end])
                    for offset in range(run_end - run_start):
                        output_cell = chunk_start + run_start + offset
                        local_start = int(indptr[int(chunk_rows[run_start + offset])]) - data_start
                        local_end = int(indptr[int(chunk_rows[run_start + offset]) + 1]) - data_start
                        for feature_index, value in zip(
                            run_indices[local_start:local_end], run_data[local_start:local_end]
                        ):
                            stream.write(f"{int(feature_index) + 1} {output_cell + 1} {int(value)}\n")
                    run_start = run_end
    return nnz


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()
    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    source_h5ad = Path(config["source_h5ad"])
    output_dir = Path(config["input_dir"])
    selection = config.get("selection", {})
    row_chunk_size = int(config.get("row_chunk_size", 512))
    if not source_h5ad.is_file():
        raise FileNotFoundError(source_h5ad)

    with h5py.File(source_h5ad, "r") as handle:
        n_cells, n_features = map(int, handle["X"].attrs["shape"])
        obs_names = read_index(handle["obs"], "_index")
        feature_names = read_index(handle["var"], "_index")
        target_label = read_obs_column(handle["obs"]["target_label"])
        is_ntc = [as_bool(value) for value in read_obs_column(handle["obs"]["is_ntc"])]
        assignment_structure = read_obs_column(handle["obs"]["assignment_structure"])
        target_ntc = [value == "non-targeting" for value in target_label]
        if is_ntc != target_ntc:
            raise ValueError("is_ntc is inconsistent with target_label")
        row_indices = choose_cells(target_label, is_ntc, assignment_structure, selection)
        if row_indices.size == 0:
            raise ValueError("No eligible cells were selected")

        allowed_structures = set(
            selection.get(
                "eligible_assignment_structures",
                ["single_guide", "concordant_construct"],
            )
        )
        selected_structures = {assignment_structure[index] for index in row_indices}
        if not selected_structures.issubset(allowed_structures):
            raise ValueError("Selected cells contain a disallowed assignment structure")

        exclude_pattern = re.compile(str(config.get("metadata_exclude_regex", "^$")))
        required_metadata = {
            "guide_id",
            "target_label",
            "is_ntc",
            "batch_id",
        }
        metadata_columns = [
            name
            for name in handle["obs"].keys()
            if name in required_metadata or not exclude_pattern.search(str(name))
        ]
        obs = {}
        for name in metadata_columns:
            values = read_obs_column(handle["obs"][name])
            obs[name] = [values[index] for index in row_indices]

    output_dir.mkdir(parents=True, exist_ok=True)
    selected_cell_names = [obs_names[index] for index in row_indices]
    metadata = pd.DataFrame(obs, index=selected_cell_names)
    metadata.index.name = "cell_id"
    metadata["gene"] = metadata["target_label"]
    metadata.to_csv(output_dir / "metadata.tsv.gz", sep="\t", compression="gzip")
    features = pd.DataFrame(
        {"gene_id": feature_names, "gene_name": feature_names, "feature_type": "Gene Expression"}
    )
    features.to_csv(output_dir / "features.tsv.gz", sep="\t", header=False, index=False, compression="gzip")
    nnz = write_matrix_market(source_h5ad, output_dir / "counts.mtx.gz", row_indices, n_features, row_chunk_size)

    manifest = {
        "module": "M03",
        "method": "Mixscape",
        "dataset": config["dataset"],
        "source_h5ad": str(source_h5ad.resolve()),
        "source_shape": [n_cells, n_features],
        "output_shape": [int(n_features), int(row_indices.size)],
        "counts_nnz": nnz,
        "counts_source": "layers/counts",
        "selection": selection,
        "selected_cells": int(row_indices.size),
        "selected_ntc_cells": int(sum(is_ntc[index] is True for index in row_indices)),
        "selected_assignment_structures": sorted(
            {str(assignment_structure[index]) for index in row_indices}
        ),
        "metadata_columns": metadata_columns + ["gene"],
        "metadata_exclude_regex": config.get("metadata_exclude_regex", "^$"),
        "ordering": "sorted original H5AD row order",
    }
    (output_dir / "adapter_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
