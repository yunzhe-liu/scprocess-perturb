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
import scipy.sparse as sp
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


def resolve_feature_mode(feature_names: list[str], requested: str) -> tuple[str, list[str]]:
    """Resolve gene-level versus strict alevin-fry S/U/A block input."""
    if requested not in {"auto", "gene", "usa_sa"}:
        raise ValueError("feature_mode must be one of: auto, gene, usa_sa")
    exact_usa = False
    base_names: list[str] = []
    if len(feature_names) % 3 == 0:
        n_genes = len(feature_names) // 3
        suffixes = [name.rsplit("_", 1)[-1] for name in feature_names]
        bases = [name.rsplit("_", 1)[0] for name in feature_names]
        exact_usa = (
            suffixes == ["S"] * n_genes + ["U"] * n_genes + ["A"] * n_genes
            and bases[:n_genes] == bases[n_genes : 2 * n_genes] == bases[2 * n_genes :]
        )
        if exact_usa:
            base_names = bases[:n_genes]
    usa_suffix_count = sum(name.endswith(("_S", "_U", "_A")) for name in feature_names)
    suspected_usa = usa_suffix_count >= max(3, len(feature_names) // 10)
    if requested == "usa_sa" and not exact_usa:
        raise ValueError("feature_mode=usa_sa requires complete contiguous S/U/A feature blocks")
    if requested == "auto":
        if exact_usa:
            return "usa_sa", base_names
        if suspected_usa:
            raise ValueError("Input resembles USA features but fails the strict S/U/A block contract")
        return "gene", feature_names
    return requested, base_names if requested == "usa_sa" else feature_names


def aggregate_usa_block(
    matrix: h5py.Group,
    indptr: np.ndarray,
    source_start: int,
    source_end: int,
    selected_rows: np.ndarray,
    n_genes: int,
) -> sp.csr_matrix:
    data_start = int(indptr[source_start])
    data_end = int(indptr[source_end])
    raw = sp.csr_matrix(
        (
            np.asarray(matrix["data"][data_start:data_end]),
            np.asarray(matrix["indices"][data_start:data_end], dtype=np.int64),
            indptr[source_start : source_end + 1] - data_start,
        ),
        shape=(source_end - source_start, 3 * n_genes),
    )
    if not np.all(np.equal(raw.data, np.rint(raw.data))):
        raise ValueError("layers['counts'] contains non-integer values")
    selected = raw[selected_rows - source_start]
    aggregated = (selected[:, :n_genes] + selected[:, 2 * n_genes :]).tocsr()
    aggregated.sum_duplicates()
    aggregated.eliminate_zeros()
    return aggregated


def iter_usa_blocks(
    source_h5ad: Path,
    row_indices: np.ndarray,
    n_genes: int,
    row_chunk_size: int,
):
    with h5py.File(source_h5ad, "r") as handle:
        matrix = handle["layers/counts"]
        if decode_scalar(matrix.attrs.get("encoding-type", "")) != "csr_matrix":
            raise ValueError("USA S+A aggregation requires a CSR counts layer")
        indptr = np.asarray(matrix["indptr"][:], dtype=np.int64)
        n_source_cells = len(indptr) - 1
        for source_start in range(0, n_source_cells, row_chunk_size):
            source_end = min(source_start + row_chunk_size, n_source_cells)
            selected_rows = row_indices[
                (row_indices >= source_start) & (row_indices < source_end)
            ]
            if selected_rows.size:
                yield aggregate_usa_block(
                    matrix, indptr, source_start, source_end, selected_rows, n_genes
                )


def write_usa_sa_matrix(
    source_h5ad: Path,
    output_path: Path,
    row_indices: np.ndarray,
    n_genes: int,
    row_chunk_size: int,
) -> int:
    expected_nnz = sum(
        int(block.nnz)
        for block in iter_usa_blocks(source_h5ad, row_indices, n_genes, row_chunk_size)
    )
    written_nnz = 0
    written_cells = 0
    with gzip.open(output_path, "wt", encoding="utf-8", compresslevel=1) as stream:
        stream.write("%%MatrixMarket matrix coordinate integer general\n")
        stream.write("% Gene-level counts aggregated from USA features as S+A; U excluded.\n")
        stream.write(f"{n_genes} {len(row_indices)} {expected_nnz}\n")
        for block in iter_usa_blocks(source_h5ad, row_indices, n_genes, row_chunk_size):
            columns = np.repeat(
                np.arange(written_cells, written_cells + block.shape[0], dtype=np.int64) + 1,
                np.diff(block.indptr),
            )
            entries = np.column_stack(
                (block.indices.astype(np.int64) + 1, columns, block.data.astype(np.int64))
            )
            if entries.size:
                np.savetxt(stream, entries, fmt="%d %d %d")
            written_nnz += int(block.nnz)
            written_cells += block.shape[0]
    if written_cells != len(row_indices) or written_nnz != expected_nnz:
        raise RuntimeError("USA S+A write coverage mismatch")
    return written_nnz


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()
    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    source_h5ad = Path(config["source_h5ad"])
    output_dir = Path(config["input_dir"])
    selection = config.get("selection", {})
    row_chunk_size = int(config.get("row_chunk_size", 512))
    requested_feature_mode = str(config.get("feature_mode", "auto"))
    if not source_h5ad.is_file():
        raise FileNotFoundError(source_h5ad)

    with h5py.File(source_h5ad, "r") as handle:
        n_cells, n_features = map(int, handle["X"].attrs["shape"])
        obs_names = read_index(handle["obs"], "_index")
        feature_names = read_index(handle["var"], "_index")
        feature_mode, output_feature_names = resolve_feature_mode(
            feature_names, requested_feature_mode
        )
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
        {
            "gene_id": output_feature_names,
            "gene_name": output_feature_names,
            "feature_type": "Gene Expression",
        }
    )
    features.to_csv(output_dir / "features.tsv.gz", sep="\t", header=False, index=False, compression="gzip")
    if feature_mode == "usa_sa":
        nnz = write_usa_sa_matrix(
            source_h5ad, output_dir / "counts.mtx.gz", row_indices,
            len(output_feature_names), row_chunk_size,
        )
        counts_source = "layers/counts; gene-level S+A aggregation; U excluded"
    else:
        nnz = write_matrix_market(
            source_h5ad, output_dir / "counts.mtx.gz", row_indices,
            n_features, row_chunk_size,
        )
        counts_source = "layers/counts"

    manifest = {
        "module": "M03",
        "method": "Mixscape",
        "dataset": config["dataset"],
        "source_h5ad": str(source_h5ad.resolve()),
        "source_shape": [n_cells, n_features],
        "output_shape": [int(len(output_feature_names)), int(row_indices.size)],
        "counts_nnz": nnz,
        "counts_source": counts_source,
        "feature_mode_requested": requested_feature_mode,
        "feature_mode_resolved": feature_mode,
        "usa_aggregation": "S+A; U excluded" if feature_mode == "usa_sa" else None,
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
