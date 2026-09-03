#!/usr/bin/env python3
"""Fuse per-lane expression matrices with guide assignments.

The integration step is the terminal export for scprocess-perturb.  It reads:

  * one expression matrix per group (``group=path``),
  * optionally, the merged guide barcode list, and
  * one assignment CSV (standardized workflow schema or supported raw schema).

The output is a new AnnData file.  Existing GEX, MEX, and assignment files are
never modified.

Cell identity follows the same convention as ``rules/merge.smk``.  When the
expression ``obs`` contains ``lane`` and ``barcode_16mer``, the key is
``{16mer}-L{lane:02d}``; otherwise a normalized 16-base barcode receives the
group suffix (``-L01`` for a group ending in ``01``).  Already-merged
lane-prefixed IDs such as ``l1_AAAC...`` are reduced to the 16mer and repeated
keys are retained once in deterministic file order.  The GEX cells form the
output cell universe; assignment rows are left-joined onto that universe.  If
an already-merged input has IDs in the form ``{16mer}-{number}``, that suffix
is retained for batch-aware alignment.

Output layout:

  X                         log1p(CP10K) expression matrix, cells x genes
  layers["counts"]          original integer-valued counts
  obs.cell_key              canonical integration key
  obs.assignment_structure    single_guide/concordant_construct/mixed_construct
  obs.guide_count              number of unique valid guides per cell
  obs.construct_count          number of mapped constructs per cell
  obs.top_guide                top-ranked guide, for traceability only
  obs.resolved_construct       construct only when uniquely resolved
  obs.assigned_guide           compatibility alias of top_guide
  obs.assigned_construct       compatibility alias of resolved_construct
  obsm.guide_candidates        sparse native-score matrix for all guides
  obsm.construct_candidates    sparse construct-score matrix when mapped
  uns.assignment_score        score metadata for the assignment matrix
  uns.guide_candidates_guides / uns.construct_candidates_constructs
                              column order for candidate matrices
  uns.integration/manifest  provenance and alignment summary

Both the current standardized assignment schema and the older raw schema are
accepted.  The workflow itself supplies the standardized schema.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import gzip
import json
import re
import resource
import shutil
from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
from anndata import AnnData


_BARCODE_RE = re.compile(r"([ACGTN]{16})", re.IGNORECASE)
_LANE_PREFIXED_BARCODE_RE = re.compile(r"^l\d+_[ACGTN]{16}$", re.IGNORECASE)
_SUFFIXED_BARCODE_RE = re.compile(r"^[ACGTN]{16}-(\d+)$", re.IGNORECASE)
_CELL_SUFFIX_RE = re.compile(r"-L?0*(\d+)$", re.IGNORECASE)
_GUIDE_ID_SANITIZE_RE = re.compile(r"[^a-zA-Z0-9_-]")
_NONTARGETING_RE = re.compile(
    r"non[-_ ]?target(?:ing)?|negctrl|ntc|ntg\d+|control|intergenic",
    re.IGNORECASE,
)
GUIDE_DESIGN_CHOICES = ("single", "dual", "multi")
ASSIGNMENT_STRUCTURE_CHOICES = (
    "single_guide",
    "concordant_construct",
    "mixed_construct",
)


@dataclass
class ExpressionData:
    X: object
    counts: object | None
    obs: pd.DataFrame
    var: pd.DataFrame
    backing: object | None = None
    counts_origin: str = "unknown"
    normalized_origin: str = "unknown"


@dataclass
class AssignmentData:
    score_column: str
    score_type: str
    top1: pd.Series
    candidates: pd.DataFrame
    guides: list[str]


def _text(value) -> str:
    """Decode a scalar HDF5 value to text."""
    if isinstance(value, bytes):
        return value.decode()
    if isinstance(value, np.ndarray) and value.ndim == 0:
        return _text(value.item())
    return str(value)


def _text_array(values) -> np.ndarray:
    return np.asarray([_text(v) for v in values], dtype=str)


def _read_h5_dataframe(group: h5py.Group) -> pd.DataFrame:
    """Read the subset of AnnData dataframe encoding used by the inputs."""
    index_key = _text(group.attrs.get("_index", "_index"))
    order_attr = group.attrs.get("column-order", [])
    order = [_text(x) for x in np.asarray(order_attr).reshape(-1)]
    if not order:
        order = [k for k in group.keys() if k != index_key]

    def read_col(node):
        if isinstance(node, h5py.Group) and "categories" in node and "codes" in node:
            categories = _text_array(node["categories"][:])
            codes = node["codes"][:]
            return pd.Categorical.from_codes(codes, categories=categories)
        values = node[:]
        if values.dtype.kind in "SUO":
            return _text_array(values)
        return values

    index = read_col(group[index_key])
    data = {column: read_col(group[column]) for column in order if column in group}
    frame = pd.DataFrame(data)
    frame.index = pd.Index(np.asarray(index).astype(str))
    return frame


def _read_sparse_or_dense(node):
    if isinstance(node, h5py.Group):
        encoding = _text(node.attrs.get("encoding-type", ""))
        shape = tuple(int(x) for x in node.attrs["shape"])
        data = node["data"][:]
        indices = node["indices"][:]
        indptr = node["indptr"][:]
        if encoding == "csr_matrix":
            return sp.csr_matrix((data, indices, indptr), shape=shape)
        if encoding == "csc_matrix":
            return sp.csc_matrix((data, indices, indptr), shape=shape)
        raise ValueError(f"unsupported sparse encoding {encoding!r}")
    return sp.csr_matrix(node[:])


class _H5MatrixReader:
    """Small row-slice reader that keeps an h5ad matrix out of RAM."""

    def __init__(self, path: Path, node_path: str):
        self.handle = h5py.File(path, "r")
        self.node = self.handle[node_path]
        self.is_csr = isinstance(self.node, h5py.Group)
        self.shape = tuple(int(x) for x in self.node.attrs["shape"])
        self.dtype = (
            self.node["data"].dtype if self.is_csr else self.node.dtype
        )

    def __getitem__(self, item):
        if not self.is_csr:
            return self.node[item]
        if not isinstance(item, slice) or item.step not in (None, 1):
            raise TypeError("backed CSR reader supports contiguous row slices only")
        start, stop, _ = item.indices(self.shape[0])
        indptr = self.node["indptr"][start:stop + 1]
        if indptr.size == 0:
            return sp.csr_matrix((0, self.shape[1]), dtype=self.dtype)
        lo, hi = int(indptr[0]), int(indptr[-1])
        return sp.csr_matrix(
            (
                self.node["data"][lo:hi],
                self.node["indices"][lo:hi],
                indptr - lo,
            ),
            shape=(stop - start, self.shape[1]),
        )

    def close(self) -> None:
        self.handle.close()


def _h5ad_has_node(path: Path, node_path: str) -> bool:
    with h5py.File(path, "r") as handle:
        return node_path in handle


def _h5ad_sparse_stats(path: Path, node_path: str = "/X") -> tuple[int, int, int] | None:
    """Return rows, cols, nnz for a sparse h5ad node without reading data."""
    with h5py.File(path, "r") as handle:
        if node_path not in handle:
            return None
        node = handle[node_path]
        if not isinstance(node, h5py.Group) or "data" not in node:
            return None
        shape = tuple(int(x) for x in node.attrs["shape"])
        return shape[0], shape[1], int(node["data"].shape[0])


def _apply_memory_limit(max_memory_gb: float) -> None:
    if max_memory_gb <= 0:
        raise ValueError("--max-process-memory-gb must be positive")
    limit = int(max_memory_gb * (1024 ** 3))
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    if hard != resource.RLIM_INFINITY and limit > hard:
        limit = hard
    resource.setrlimit(resource.RLIMIT_AS, (limit, hard))


def preflight_resources(
    groups: list[tuple[str, Path]],
    output: Path,
    input_kind: str,
    counts_layer: str,
    max_input_nnz: int,
    max_output_gb: float,
    min_free_disk_gb: float,
    max_process_memory_gb: float,
) -> dict:
    """Fail before integration if a large run exceeds explicit safeguards."""
    input_nnz = 0
    estimated_output = 0
    stats = []
    for _, path in groups:
        result = _h5ad_sparse_stats(path)
        if result is None:
            continue
        n_rows, n_cols, x_nnz = result
        count_result = _h5ad_sparse_stats(path, f"/layers/{counts_layer}")
        count_nnz = count_result[2] if count_result is not None else x_nnz
        input_nnz += count_nnz
        if input_kind == "counts" and count_result is None:
            normalized_nnz = count_nnz
        else:
            normalized_nnz = x_nnz
        # The writer stores float32 data + int32 column indices for both CSR
        # matrices, plus int64 row pointers.  Add 25% for HDF5/metadata
        # overhead.
        matrix_bytes = (count_nnz + normalized_nnz) * (4 + 8)
        matrix_bytes += 2 * (n_rows + 1) * 8
        estimated_output += int(matrix_bytes * 1.25)
        stats.append({"path": str(path), "shape": [n_rows, n_cols], "nnz": count_nnz})

    if input_nnz > max_input_nnz:
        raise RuntimeError(
            f"resource preflight refused input nnz={input_nnz:,}; "
            f"limit={max_input_nnz:,}"
        )
    estimated_output_gb = estimated_output / 1024**3
    if estimated_output_gb > max_output_gb:
        raise RuntimeError(
            "resource preflight refused output estimate "
            f"{estimated_output_gb:.1f} GiB > "
            f"--max-output-gb={max_output_gb:.1f}"
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    free_bytes = shutil.disk_usage(output.parent).free
    required_bytes = estimated_output + int(min_free_disk_gb * (1024 ** 3))
    if estimated_output and free_bytes < required_bytes:
        raise RuntimeError(
            "resource preflight refused output: estimated "
            f"{estimated_output / 1024**3:.1f} GiB plus reserve "
            f"{min_free_disk_gb:.1f} GiB exceeds free disk "
            f"{free_bytes / 1024**3:.1f} GiB"
        )
    _apply_memory_limit(max_process_memory_gb)
    summary = {
        "input_nnz": input_nnz,
        "estimated_output_gb": estimated_output_gb,
        "free_disk_gb": free_bytes / 1024**3,
        "max_process_memory_gb": max_process_memory_gb,
        "inputs": stats,
    }
    print("[preflight] " + json.dumps(summary, sort_keys=True), flush=True)
    return summary


def _read_h5ad_group(
    group: h5py.Group,
    input_kind: str,
    counts_layer: str,
) -> ExpressionData:
    X = _read_sparse_or_dense(group["X"])
    layers = group.get("layers")
    counts = None
    counts_origin = "unavailable"
    normalized_origin = "input.X"
    if layers is not None and counts_layer in layers:
        counts = _read_sparse_or_dense(layers[counts_layer])
        counts_origin = f"input.layers[{counts_layer!r}]"
    elif input_kind == "counts":
        counts = X
        X = None
        counts_origin = "input.X"
        normalized_origin = "computed from input.X"
    elif input_kind == "auto":
        # An integer-valued X is an unambiguous raw-count input.  Do not
        # silently reinterpret a floating-point normalized h5ad as counts.
        if np.issubdtype(X.dtype, np.integer):
            counts = X
            X = None
            counts_origin = "input.X"
            normalized_origin = "computed from input.X"
        else:
            raise ValueError(
                "expression h5ad has no counts layer; provide --input-kind "
                "counts only when X contains raw counts"
            )
    elif input_kind != "standardized":
        raise ValueError(f"unsupported expression input kind {input_kind!r}")
    if counts is None:
        raise ValueError(
            "expression input must provide a counts layer or raw integer X"
        )
    obs = _read_h5_dataframe(group["obs"])
    var = _read_h5_dataframe(group["var"])
    matrix = counts if counts is not None else X
    if matrix.shape[0] != len(obs) or matrix.shape[1] != len(var):
        raise ValueError(
            f"expression dimensions {matrix.shape} do not match obs/var "
            f"({len(obs)}, {len(var)})"
        )
    if counts is not None and counts.shape != matrix.shape:
        raise ValueError("counts layer dimensions do not match expression X")
    return ExpressionData(
        X=None if X is None else X.tocsr(),
        counts=None if counts is None else counts.tocsr(),
        obs=obs,
        var=var,
        counts_origin=counts_origin,
        normalized_origin=normalized_origin,
    )


def _read_expression(
    path: Path,
    input_kind: str = "auto",
    counts_layer: str = "counts",
    counts_source: Path | None = None,
    normalized_source: Path | None = None,
) -> ExpressionData:
    """Read h5ad, h5mu RNA, or the scprocess raw H5 matrix."""
    suffix = path.suffix.lower()
    if suffix == ".h5ad":
        # Keep large sparse h5ad inputs backed.  This avoids materializing all
        # X data before the output writer can stream it to the new artifact.
        backed = ad.read_h5ad(path, backed="r")
        backings = [backed]

        def aligned(other, label: str) -> None:
            if other.shape != backed.shape:
                raise ValueError(f"{label}: shape differs from {path}")
            if not other.obs_names.equals(backed.obs_names):
                raise ValueError(f"{label}: obs_names are not identically ordered")
            if not other.var_names.equals(backed.var_names):
                raise ValueError(f"{label}: var_names are not identically ordered")

        if normalized_source is not None:
            normalized_backed = ad.read_h5ad(normalized_source, backed="r")
            aligned(normalized_backed, "normalized source")
            normalized = normalized_backed.X
            backings.append(normalized_backed)
            normalized_origin = f"{normalized_source.resolve()}:.X"
        elif input_kind == "counts":
            normalized = None
            normalized_origin = "computed from input.X"
        else:
            normalized = backed.X
            normalized_origin = "input.X"

        if counts_source is not None:
            counts_backed = ad.read_h5ad(counts_source, backed="r")
            aligned(counts_backed, "counts source")
            counts_node = f"/layers/{counts_layer}"
            if _h5ad_has_node(counts_source, counts_node):
                counts = _H5MatrixReader(counts_source, counts_node)
                counts_origin = f"{counts_source.resolve()}:layers[{counts_layer!r}]"
                backings.append(counts)
            elif input_kind == "counts" or np.issubdtype(counts_backed.X.dtype, np.integer):
                counts = counts_backed.X
                counts_origin = f"{counts_source.resolve()}:.X"
            else:
                raise ValueError(
                    f"{counts_source}: no layers[{counts_layer!r}] found and X "
                    "is not integer-valued"
                )
            backings.append(counts_backed)
        elif _h5ad_has_node(path, f"/layers/{counts_layer}"):
            counts = _H5MatrixReader(path, f"/layers/{counts_layer}")
            counts_origin = f"input.layers[{counts_layer!r}]"
            backings.append(counts)
        elif input_kind == "counts" or (
            input_kind == "auto" and np.issubdtype(backed.X.dtype, np.integer)
        ):
            counts = backed.X
            if normalized_source is None:
                normalized = None
                normalized_origin = "computed from input.X"
            counts_origin = "input.X"
        else:
            raise ValueError(
                f"{path}: no layers[{counts_layer!r}] found; provide a validated "
                "--counts-source or --input-kind counts for raw X"
            )
        return ExpressionData(
            X=normalized,
            counts=counts,
            obs=backed.obs.copy(),
            var=backed.var.copy(),
            backing=backings,
            counts_origin=counts_origin,
            normalized_origin=normalized_origin,
        )

    if suffix == ".h5mu":
        with h5py.File(path, "r") as handle:
            if "mod" not in handle or "rna" not in handle["mod"]:
                raise ValueError(f"{path}: h5mu does not contain mod/rna")
            return _read_h5ad_group(handle["mod"]["rna"], input_kind, counts_layer)

    if suffix == ".h5":
        with h5py.File(path, "r") as handle:
            if "matrix" not in handle:
                raise ValueError(f"{path}: expected matrix group")
            matrix = handle["matrix"]
            shape = tuple(int(x) for x in matrix["shape"][:])
            # scprocess raw H5 follows the 10x convention: genes x cells.
            X = sp.csc_matrix(
                (matrix["data"][:], matrix["indices"][:], matrix["indptr"][:]),
                shape=shape,
            ).T.tocsr()
            barcodes = _text_array(matrix["barcodes"][:])
            features = matrix.get("features")
            if features is None:
                raise ValueError(f"{path}: matrix/features is missing")
            feature_node = features.get("name", features.get("id"))
            if feature_node is None:
                raise ValueError(f"{path}: matrix/features needs name or id")
            feature_names = _text_array(feature_node[:])
        obs = pd.DataFrame(index=pd.Index(barcodes))
        var = pd.DataFrame(index=pd.Index(feature_names))
        return ExpressionData(
            X=None,
            counts=X,
            obs=obs,
            var=var,
            counts_origin="input.X",
            normalized_origin="computed from input.X",
        )

    raise ValueError(f"unsupported expression format {path}; expected .h5/.h5ad/.h5mu")


def _normalize_barcode(cell_id: str) -> str:
    match = _BARCODE_RE.search(str(cell_id).upper())
    if match is None:
        raise ValueError(f"cannot find a 16-base barcode in expression cell id {cell_id!r}")
    return match.group(1)


def _group_suffix(group: str) -> str:
    match = re.search(r"(\d+)$", group)
    return f"-L{match.group(1)}" if match else f"-{group}"


def _matrix_nnz(matrix) -> int:
    """Return nnz without materializing a backed sparse matrix."""
    if matrix is None:
        return 0
    if sp.issparse(matrix):
        return int(matrix.nnz)
    if hasattr(matrix, "nnz"):
        return int(matrix.nnz)
    return -1


def _materialize_csr(matrix) -> sp.csr_matrix:
    if sp.issparse(matrix):
        return matrix.tocsr()
    return sp.csr_matrix(matrix[:])


def _normalize_counts(matrix, target_sum: float) -> sp.csr_matrix:
    """Compute log1p(CP10K)-style normalization for an in-memory CSR matrix."""
    counts = matrix.tocsr(copy=True).astype(np.float32)
    totals = np.asarray(counts.sum(axis=1)).ravel().astype(np.float64)
    scale = np.divide(
        float(target_sum), totals,
        out=np.zeros_like(totals), where=totals > 0,
    )
    repeated = np.repeat(scale, np.diff(counts.indptr))
    counts.data = np.log1p(counts.data.astype(np.float64) * repeated).astype(
        np.float32
    )
    return counts


def _expression_cell_keys(obs: pd.DataFrame, group: str) -> list[str]:
    """Build keys from explicit lane metadata or the configured group suffix."""
    if {"lane", "barcode_16mer"}.issubset(obs.columns):
        lane_values = pd.to_numeric(obs["lane"], errors="raise").to_numpy(
            dtype=float
        )
        if (
            not np.isfinite(lane_values).all()
            or not np.equal(lane_values, np.floor(lane_values)).all()
            or (lane_values < 1).any()
        ):
            raise ValueError("expression obs.lane must contain positive integers")
        lanes = lane_values.astype(int)
        barcodes = obs["barcode_16mer"].astype(str)
        return [
            f"{_normalize_barcode(barcode)}-L{lane:02d}"
            for barcode, lane in zip(barcodes, lanes)
        ]

    original_ids = obs.index.astype(str).to_numpy()
    if original_ids.size and all(
        _LANE_PREFIXED_BARCODE_RE.fullmatch(cell_id) for cell_id in original_ids
    ):
        # Some already-merged h5ad files encode lane in the cell ID (for
        # example, Papalexi: l1_AAAC...).  The assignment table uses the
        # canonical 16mer-L## form, so retain the barcode and let the loader
        # deterministically collapse repeated 16mers across lanes below.
        return [_normalize_barcode(cell) for cell in original_ids]
    if original_ids.size and all(
        _SUFFIXED_BARCODE_RE.fullmatch(cell_id) for cell_id in original_ids
    ):
        # Already-merged 10x-style h5ad files commonly preserve a numeric
        # sample/lane suffix (for example, AAAC...-1).  Keep it so assignment
        # rows using the same IDs align to the correct batch.
        return [
            f"{_normalize_barcode(cell)}-"
            f"{_SUFFIXED_BARCODE_RE.fullmatch(cell).group(1)}"
            for cell in original_ids
        ]
    suffix = _group_suffix(group)
    return [f"{_normalize_barcode(cell)}{suffix}" for cell in original_ids]


def _parse_pairs(items: list[str], label: str) -> list[tuple[str, Path]]:
    pairs = []
    seen = set()
    for item in items:
        if "=" not in item:
            raise ValueError(f"--{label} expects name=path, got {item!r}")
        name, raw_path = item.split("=", 1)
        if not name or not raw_path:
            raise ValueError(f"--{label} contains an empty name or path: {item!r}")
        if name in seen:
            raise ValueError(f"duplicate {label} name {name!r}")
        seen.add(name)
        path = Path(raw_path)
        if not path.exists():
            raise FileNotFoundError(f"{label} input does not exist: {path}")
        pairs.append((name, path))
    return pairs


def _parse_path(item: str, label: str) -> Path:
    path = Path(item)
    if not path.exists():
        raise FileNotFoundError(f"{label} input does not exist: {path}")
    return path


def load_expression(
    groups: list[tuple[str, Path]],
    input_kind: str,
    counts_layer: str,
    target_sum: float,
    counts_source: Path | None = None,
    normalized_source: Path | None = None,
) -> tuple[ExpressionData, list[str]]:
    if not groups:
        raise ValueError("at least one --gex group=path is required")
    if (counts_source is not None or normalized_source is not None) and len(groups) != 1:
        raise ValueError("external counts/normalized sources require exactly one --gex input")

    matrices = []
    count_matrices = []
    counts_origins = []
    normalized_origins = []
    observations = []
    lane_prefixed_input = True
    first_var = None
    first_var_names = None
    backings = []
    for group, path in groups:
        data = _read_expression(
            path,
            input_kind=input_kind,
            counts_layer=counts_layer,
            counts_source=counts_source,
            normalized_source=normalized_source,
        )
        if data.backing is not None:
            if isinstance(data.backing, list):
                backings.extend(data.backing)
            else:
                backings.append(data.backing)
        var_names = data.var.index.astype(str).tolist()
        if first_var_names is None:
            first_var_names = var_names
            first_var = data.var.copy()
        elif var_names != first_var_names:
            raise ValueError(
                f"{path}: gene order differs from the first expression input; "
                "align genes before integration"
            )

        original_ids = data.obs.index.astype(str).to_numpy()
        lane_prefixed_input = lane_prefixed_input and bool(original_ids.size) and all(
            _LANE_PREFIXED_BARCODE_RE.fullmatch(cell_id) for cell_id in original_ids
        )
        cell_keys = _expression_cell_keys(data.obs, group)
        obs = data.obs.copy()
        obs["cell_key"] = cell_keys
        obs["source_group"] = group
        obs["source_cell_id"] = original_ids
        obs.index = pd.Index(cell_keys, name="cell_id")
        normalized = data.X
        counts = data.counts
        if normalized is None and counts is not None and data.backing is None:
            normalized = _normalize_counts(counts, target_sum)
            normalized_origins.append("computed from input.X")
        else:
            normalized_origins.append(data.normalized_origin)
        matrices.append(normalized)
        count_matrices.append(counts)
        counts_origins.append(data.counts_origin)
        observations.append(obs)

    obs = pd.concat(observations, axis=0, sort=False)
    if len(matrices) == 1:
        X = matrices[0]
        counts = count_matrices[0]
    else:
        # Per-lane inputs are normally the small/raw .h5 files.  Backed
        # multi-file concatenation cannot remain lazy, so materialize only
        # this non-Replogle path; the large merged h5ad path stays streamed.
        X = sp.vstack([
            matrix if matrix is not None else _normalize_counts(
                _materialize_csr(count), target_sum
            )
            for matrix, count in zip(matrices, count_matrices)
        ], format="csr")
        counts = sp.vstack(
            [_materialize_csr(count) for count in count_matrices], format="csr"
        )
    if not obs.index.is_unique:
        # Lane-prefixed merged h5ad inputs can contain the same 16mer in more
        # than one lane.  Keep the first row in file order, which is stable
        # and preserves the same cell universe used by the assignment table.
        if not lane_prefixed_input:
            raise ValueError("expression inputs produce duplicate integration cell keys")
        keep = ~obs.index.duplicated(keep="first")
        obs = obs.iloc[np.flatnonzero(keep)].copy()
        if X is not None:
            X = X[keep]
        counts = counts[keep]
    if not obs.index.is_unique:
        raise ValueError("expression inputs produce duplicate integration cell keys")
    if X is None:
        # This is only possible for a backed counts-only h5ad.  It is handled
        # by the streaming writer; keeping the source lazy is intentional.
        normalized_origin = "computed from input.X"
    else:
        normalized_origin = (
            normalized_origins[0] if len(set(normalized_origins)) == 1
            else "mixed: preserved input.X and computed from counts"
        )
    counts_origin = (
        counts_origins[0] if len(set(counts_origins)) == 1
        else "multiple input count sources"
    )
    return ExpressionData(
        X=X,
        counts=counts,
        obs=obs,
        var=first_var,
        backing=backings,
        counts_origin=counts_origin,
        normalized_origin=normalized_origin,
    ), list(obs.index)


def _pick_raw_score_column(columns) -> str:
    for column in ("prob_gaussian", "log_pval", "percent_counts", "UMI_counts"):
        if column in columns:
            return column
    raise ValueError(
        f"assignment CSV has no supported score column; found {list(columns)}"
    )


def _sanitize_guide_id(value: str) -> str:
    """Match guide ID sanitization used by feature_reference_adapter.py."""
    return _GUIDE_ID_SANITIZE_RE.sub("_", str(value).strip())


def _guide_aliases(value: str) -> list[str]:
    """Cover raw, comma-flattened, and feature-reference sanitized IDs."""
    raw = str(value).strip()
    return list(dict.fromkeys((raw, raw.replace(",", "_"), _sanitize_guide_id(raw))))


def load_guide_construct_map(path: Path) -> dict[str, str]:
    """Load a guide -> construct mapping from wide or normalized guide CSV."""
    frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    columns = set(frame.columns)
    recognized_columns = (
        {"sgID_A", "sgID_B"}.issubset(columns)
        or ("guide_id" in columns and "construct_id" in columns)
        or ("gRNA" in columns and any(
            column in columns for column in ("construct_id", "pair_id", "unique sgRNA pair ID")
        ))
    )
    if not recognized_columns:
        # The reference t2g files are headerless two-column tables.  Their
        # first column is the guide and their second column is the construct
        # identity; this is especially useful for single-guide systems where
        # guide-level construct IDs map each guide to itself.
        raw_frame = pd.read_csv(path, sep=None, engine="python", header=None, dtype=str,
                                keep_default_na=False)
        if raw_frame.shape[1] != 2:
            raise ValueError(
                f"{path}: unsupported guide mapping; expected a recognized CSV "
                "schema or a headerless two-column guide/construct table"
            )
        raw_frame.columns = ["guide_id", "construct_id"]
        frame = raw_frame
        columns = set(frame.columns)
    duplicate_column = "duplicated guide pair?"
    if duplicate_column in columns:
        duplicate = frame[duplicate_column].str.strip().str.lower().isin(
            {"true", "yes", "1"}
        )
        frame = frame.loc[~duplicate].copy()
    mapping = {}

    def add_mapping(guide: str, construct: str) -> None:
        previous = mapping.get(guide)
        if previous is not None and previous != construct:
            raise ValueError(
                f"{path}: guide {guide!r} maps to multiple constructs "
                f"({previous!r}, {construct!r}); construct mode needs a "
                "unique guide-to-construct mapping"
            )
        mapping[guide] = construct

    if {"sgID_A", "sgID_B"}.issubset(columns):
        construct_column = next(
            (
                column
                for column in ("construct_id", "pair_id", "unique sgRNA pair ID")
                if column in columns
            ),
            None,
        )
        if construct_column is None:
            raise ValueError(
                f"{path}: dual guide CSV needs construct_id, pair_id, or "
                "unique sgRNA pair ID"
            )
        for row in frame.to_dict(orient="records"):
            construct = str(row[construct_column]).strip()
            if not construct:
                continue
            for column in ("sgID_A", "sgID_B"):
                raw_guide = str(row[column]).strip()
                for guide in _guide_aliases(raw_guide):
                    if guide:
                        add_mapping(guide, construct)
        return mapping

    guide_column = (
        "guide_id" if "guide_id" in columns
        else ("gRNA" if "gRNA" in columns else None)
    )
    construct_column = next(
        (
            column
            for column in ("construct_id", "pair_id", "unique sgRNA pair ID")
            if column in columns
        ),
        None,
    )
    if guide_column is None or construct_column is None:
        raise ValueError(
            f"{path}: normalized guide CSV needs guide_id and construct_id "
            "(pair_id is also accepted)"
        )
    for row in frame.to_dict(orient="records"):
        raw_guide = str(row[guide_column]).strip()
        construct = str(row[construct_column]).strip()
        if raw_guide and construct:
            for guide in _guide_aliases(raw_guide):
                if guide:
                    add_mapping(guide, construct)
    return mapping


def _infer_target_label(value: str) -> str:
    """Conservative fallback for headerless guide→construct references."""
    value = str(value).strip()
    stripped = re.sub(r"g\d+$", "", value, flags=re.IGNORECASE)
    if stripped != value:
        return stripped or value
    # Headerless references used by dual-guide screens commonly retain the
    # target before the guide/strand suffix (for example AAAS_+_123-P1P2).
    # This fallback is only used when the library has no explicit target
    # column; explicit expression metadata always has precedence below.
    prefix = value.split("_", 1)[0]
    return prefix or value


def load_guide_target_map(path: Path) -> dict[str, str]:
    """Load optional guide→target metadata without requiring a new schema."""
    frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    columns = set(frame.columns)
    if {"sgID_A", "sgID_B"}.issubset(columns):
        guide_columns = ["sgID_A", "sgID_B"]
    elif "guide_id" in columns or "gRNA" in columns:
        guide_columns = ["guide_id" if "guide_id" in columns else "gRNA"]
    else:
        raw_frame = pd.read_csv(
            path, sep=None, engine="python", header=None, dtype=str,
            keep_default_na=False,
        )
        if raw_frame.shape[1] != 2:
            return {}
        raw_frame.columns = ["guide_id", "construct_id"]
        frame = raw_frame
        columns = set(frame.columns)
        guide_columns = ["guide_id"]

    target_column = next(
        (
            column for column in (
                "target_label", "target_gene_name", "target", "gene",
                "gene_symbol", "perturbation",
            ) if column in columns
        ),
        None,
    )
    target_map: dict[str, str] = {}

    def add(guide: str, target: str) -> None:
        guide = str(guide).strip()
        target = str(target).strip()
        if not guide:
            return
        target = target or _infer_target_label(guide)
        for alias in _guide_aliases(guide):
            previous = target_map.get(alias)
            if previous is not None and previous != target:
                raise ValueError(
                    f"{path}: guide {guide!r} maps to multiple targets"
                )
            target_map[alias] = target

    for row in frame.to_dict(orient="records"):
        if target_column:
            target = row.get(target_column, "")
        else:
            # A two-column guide/construct table does not provide target
            # metadata.  Derive the target from the guide identifier rather
            # than treating the construct ID itself as a target label.
            target = _infer_target_label(row.get("guide_id", ""))
        for column in guide_columns:
            add(row.get(column, ""), target)
    return target_map


def load_assignment(path: Path) -> AssignmentData:
    frame = pd.read_csv(path, keep_default_na=False)
    standardized = {"cell_barcode", "guide_id", "score"}.issubset(frame.columns)
    if standardized:
        cell_column, guide_column, score_column = "cell_barcode", "guide_id", "score"
        umi_column = "umi_count" if "umi_count" in frame.columns else None
        score_type = str(frame["score_type"].dropna().iloc[0]) if "score_type" in frame.columns and frame["score_type"].notna().any() else "native"
        rank_column = "rank" if "rank" in frame.columns else None
    else:
        if not {"cell", "gRNA"}.issubset(frame.columns):
            raise ValueError(
                f"{path}: expected standardized cell_barcode/guide_id or raw cell/gRNA columns"
            )
        cell_column, guide_column = "cell", "gRNA"
        score_column = _pick_raw_score_column(frame.columns)
        umi_column = "UMI_counts" if "UMI_counts" in frame.columns else None
        score_type = score_column
        rank_column = None

    work = pd.DataFrame({
        "cell_key": frame[cell_column].astype(str).str.strip(),
        "guide": frame[guide_column].astype(str).str.strip(),
        "weight": pd.to_numeric(frame[score_column], errors="raise"),
    })
    work["umi"] = (
        pd.to_numeric(frame[umi_column], errors="raise") if umi_column else 0
    )
    if rank_column:
        work["rank"] = pd.to_numeric(frame[rank_column], errors="raise")
    else:
        work["rank"] = np.nan
    work = work[(work["cell_key"] != "") & (work["guide"] != "")].reset_index(drop=True)

    if rank_column:
        top = work.sort_values(
            ["rank", "weight", "umi", "guide"],
            ascending=[True, False, False, True],
            kind="stable",
        )
    else:
        score_ascending = score_column == "log_pval"
        top = work.sort_values(
            ["weight", "umi", "guide"],
            ascending=[score_ascending, False, True],
            kind="stable",
        )
    top = top.drop_duplicates("cell_key", keep="first")
    top1 = pd.Series(top["guide"].to_numpy(), index=top["cell_key"].to_numpy())
    guides = sorted(work["guide"].unique().tolist())
    return AssignmentData(
        score_column=score_column,
        score_type=score_type,
        top1=top1,
        candidates=work[["cell_key", "guide", "weight"]],
        guides=guides,
    )


def _read_merged_barcodes(path: Path) -> set[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        return {line.strip() for line in handle if line.strip()}


def _candidate_matrix(assign: AssignmentData, cell_keys: list[str]) -> sp.csr_matrix:
    cell_pos = {key: i for i, key in enumerate(cell_keys)}
    guide_pos = {guide: i for i, guide in enumerate(assign.guides)}
    work = _unique_candidate_rows(assign)
    work = work[work["cell_key"].isin(cell_pos)].copy()
    if work.empty:
        return sp.csr_matrix((len(cell_keys), len(assign.guides)), dtype=np.float64)
    rows = work["cell_key"].map(cell_pos).to_numpy()
    cols = work["guide"].map(guide_pos).to_numpy()
    vals = work["weight"].to_numpy(dtype=np.float64)
    return sp.coo_matrix(
        (vals, (rows, cols)), shape=(len(cell_keys), len(assign.guides))
    ).tocsr()


def _align_assignment_cells(
    assignment: AssignmentData, cell_keys: list[str]
) -> AssignmentData:
    """Resolve assignment cell aliases against the GEX cell universe."""
    gex_keys = set(cell_keys)
    by_barcode: dict[str, list[str]] = {}
    for key in cell_keys:
        barcode = _normalize_barcode(key)
        by_barcode.setdefault(barcode, []).append(key)

    def resolve(raw_key: str) -> str:
        if raw_key in gex_keys:
            return raw_key
        barcode = _normalize_barcode(raw_key)
        candidates = by_barcode.get(barcode, [])
        raw_suffix = _CELL_SUFFIX_RE.search(raw_key)
        if _SUFFIXED_BARCODE_RE.fullmatch(raw_key) and raw_suffix:
            same_suffix = [
                key for key in candidates
                if (match := _CELL_SUFFIX_RE.search(key))
                and int(match.group(1)) == int(raw_suffix.group(1))
            ]
            return same_suffix[0] if len(same_suffix) == 1 else raw_key
        if len(candidates) == 1:
            return candidates[0]
        return raw_key

    candidates = assignment.candidates.copy()
    candidates["cell_key"] = candidates["cell_key"].map(resolve)
    top1 = assignment.top1.copy()
    top1.index = pd.Index([resolve(key) for key in top1.index], name=top1.index.name)
    return AssignmentData(
        score_column=assignment.score_column,
        score_type=assignment.score_type,
        top1=top1,
        candidates=candidates,
        guides=assignment.guides,
    )


def _unique_candidate_rows(assignment: AssignmentData) -> pd.DataFrame:
    """Keep one valid row per cell and guide for structure and score matrices."""
    return assignment.candidates.drop_duplicates(
        ["cell_key", "guide"], keep="first"
    )


def _classify_assignment_structure(
    assignment: AssignmentData,
    cell_keys: list[str],
    guide_design: str,
    guide_construct: dict[str, str] | None,
) -> tuple[pd.Series, pd.Series, pd.Series]:
    """Classify cells from all valid assignment candidates, not top-1 only."""
    rows = _unique_candidate_rows(assignment)
    rows = rows[rows["cell_key"].isin(set(cell_keys))]
    guide_count = rows.groupby("cell_key")["guide"].nunique()
    guide_count = guide_count.reindex(cell_keys, fill_value=0).astype("int64")

    construct_count = pd.Series(np.nan, index=cell_keys, dtype="float64")
    if guide_construct is not None:
        rows = rows.copy()
        rows["construct"] = rows["guide"].map(guide_construct)
        if rows["construct"].isna().any():
            missing = sorted(rows.loc[rows["construct"].isna(), "guide"].unique())
            raise ValueError(
                "construct library does not map all valid guides; missing: "
                + ", ".join(missing[:10])
                + (" ..." if len(missing) > 10 else "")
            )
        construct_count = (
            rows.groupby("cell_key")["construct"]
            .nunique()
            .reindex(cell_keys, fill_value=0)
            .astype("int64")
        )

    structure = pd.Series(pd.NA, index=cell_keys, dtype="string")
    if guide_design == "single":
        structure.loc[guide_count.eq(1)] = "single_guide"
        structure.loc[guide_count.gt(1)] = "mixed_construct"
    else:
        if guide_construct is None and guide_count.gt(1).any():
            raise ValueError(
                f"guide_design={guide_design!r} requires --guide-csv to classify "
                "multi-guide cells"
            )
        structure.loc[guide_count.eq(1)] = "single_guide"
        structure.loc[guide_count.gt(1) & construct_count.eq(1)] = (
            "concordant_construct"
        )
        structure.loc[guide_count.gt(1) & construct_count.gt(1)] = "mixed_construct"

    categories = list(ASSIGNMENT_STRUCTURE_CHOICES)
    structure = pd.Series(
        pd.Categorical(structure, categories=categories), index=cell_keys
    )
    return structure, guide_count, construct_count


def _construct_matrix(
    assign: AssignmentData,
    cell_keys: list[str],
    guide_construct: dict[str, str],
) -> tuple[sp.csr_matrix, list[str], int, str]:
    """Aggregate guide scores using the method's native score direction."""
    lower_is_better = (
        assign.score_column == "log_pval"
        or assign.score_type in {"log_pval", "neg_log_pval"}
    )
    aggregation = "min native log_pval" if lower_is_better else "max native score"
    cell_pos = {key: i for i, key in enumerate(cell_keys)}
    values = {}
    unmapped = 0
    for row in _unique_candidate_rows(assign).itertuples(index=False):
        cell_pos_i = cell_pos.get(row.cell_key)
        construct = guide_construct.get(row.guide)
        if cell_pos_i is None:
            continue
        if construct is None:
            unmapped += 1
            continue
        key = (cell_pos_i, construct)
        candidate = float(row.weight)
        if key not in values:
            values[key] = candidate
        elif lower_is_better:
            values[key] = min(values[key], candidate)
        else:
            values[key] = max(values[key], candidate)

    constructs = sorted({construct for _, construct in values})
    construct_pos = {construct: i for i, construct in enumerate(constructs)}
    if not values:
        return (
            sp.csr_matrix((len(cell_keys), 0), dtype=np.float64),
            constructs,
            unmapped,
            aggregation,
        )
    rows = [cell for cell, _ in values]
    cols = [construct_pos[construct] for _, construct in values]
    vals = list(values.values())
    matrix = sp.coo_matrix(
        (vals, (rows, cols)), shape=(len(cell_keys), len(constructs))
    ).tocsr()
    return matrix, constructs, unmapped, aggregation


def _string_series(values, index: pd.Index) -> pd.Series:
    return pd.Series(values, index=index, dtype="string")


def _existing_or_default(
    obs: pd.DataFrame,
    names: tuple[str, ...],
    index: pd.Index,
    default: str,
) -> pd.Series:
    for name in names:
        if name in obs.columns:
            values = obs[name].astype("string").reindex(index)
            values = values.mask(values.str.strip().eq(""))
            if values.notna().any():
                return values
    return _string_series([default] * len(index), index)


def _coalesce_metadata(
    obs: pd.DataFrame,
    names: tuple[str, ...],
    index: pd.Index,
) -> tuple[pd.Series, pd.Series]:
    """Coalesce metadata columns and report rows backed by source metadata."""
    values = pd.Series(pd.NA, index=index, dtype="string")
    supplied = pd.Series(False, index=index, dtype=bool)
    for name in names:
        if name not in obs.columns:
            continue
        candidate = obs[name].astype("string").reindex(index)
        candidate = candidate.mask(candidate.str.strip().eq(""))
        take = values.isna() & candidate.notna()
        values = values.mask(take, candidate)
        supplied = supplied | take
    return values, supplied


def _add_standardized_metadata(
    obs: pd.DataFrame,
    top_guide: pd.Series,
    resolved_construct: pd.Series,
    guide_count: pd.Series,
    construct_count: pd.Series,
    guide_target: dict[str, str],
    source_group: pd.Series | None = None,
) -> pd.DataFrame:
    """Add the small, stable M01 contract to the integration obs table."""
    result = obs.copy()
    index = result.index
    guide_id = top_guide.astype("string").reindex(index)
    target_from_guide = guide_id.map(guide_target).astype("string")
    target_label, target_supplied = _coalesce_metadata(
        result,
        (
            "target_label", "perturbation_name", "target_gene",
            "target_gene_name", "gene", "perturbation",
        ),
        index,
    )
    target_label = target_label.fillna(target_from_guide).fillna(guide_id)

    existing_ntc = result["is_ntc"].astype("boolean") if "is_ntc" in result else None
    # Source target metadata is authoritative (Norman uses
    # perturbation_name == control).  Only when it is absent do we inspect
    # the guide-derived target and guide ID.  In particular, never inspect a
    # resolved construct: dual-guide libraries often contain NegCtrl guides
    # inside otherwise non-control perturbations.
    ntc_text = target_label.copy()
    fallback_text = pd.concat([target_label, guide_id], axis=1).fillna("").astype(str).agg(" ".join, axis=1)
    ntc_text = ntc_text.where(target_supplied, fallback_text)
    inferred_ntc = ntc_text.str.contains(_NONTARGETING_RE, na=False)
    is_ntc = (
        existing_ntc.fillna(False).astype(bool)
        if existing_ntc is not None
        else inferred_ntc
    )
    target_label = target_label.mask(is_ntc, "non-targeting")

    result["guide_id"] = guide_id
    result["assigned_guide"] = guide_id
    result["assigned_construct_standard"] = resolved_construct.astype("string")
    result["assigned_construct"] = resolved_construct.astype("string")
    result["target_label"] = target_label.astype("string")
    result["is_ntc"] = is_ntc.astype(bool)

    # perturbation_group is the stable target-level grouping used by
    # downstream status estimation.  Construct identity remains available in
    # resolved_construct / assigned_construct.
    result["perturbation_group"] = target_label.astype("string")
    batch_id, _ = _coalesce_metadata(
        result,
        (
            "batch_id", "batch", "gemgroup", "lane", "replicate",
            "sample_id", "sample",
        ),
        index,
    )
    batch_id = batch_id.astype("string")
    # In a multi-group full workflow, the configured group is a meaningful
    # fallback batch when the source matrix has no batch column.  A single
    # standalone dataset remains explicitly unspecified rather than being
    # mislabeled with its dataset name.
    if source_group is not None:
        groups = source_group.astype("string").reindex(index)
        if groups.dropna().replace("", pd.NA).nunique() > 1:
            batch_id = batch_id.mask(batch_id.str.strip().eq(""), groups)
    result["batch_id"] = batch_id.mask(batch_id.str.strip().eq(""), "unspecified")
    result["guide_assignment_missing"] = guide_count.eq(0).astype(bool)
    return result


def integrate(
    groups: list[tuple[str, Path]],
    assignment: Path,
    merged_barcodes: Path | None,
    guide_design: str,
    config_snapshot: dict,
    guide_csv: Path | None = None,
    input_kind: str = "auto",
    counts_layer: str = "counts",
    target_sum: float = 10000.0,
    max_materialized_nnz: int = 100_000_000,
    counts_source: Path | None = None,
    normalized_source: Path | None = None,
) -> AnnData:
    if guide_design not in GUIDE_DESIGN_CHOICES:
        raise ValueError(
            f"unknown guide design {guide_design!r}; "
            f"choose from {GUIDE_DESIGN_CHOICES}"
        )
    expression, cell_keys = load_expression(
        groups, input_kind=input_kind, counts_layer=counts_layer,
        target_sum=target_sum,
        counts_source=counts_source,
        normalized_source=normalized_source,
    )
    merged_keys = (
        _read_merged_barcodes(merged_barcodes)
        if merged_barcodes is not None
        else set()
    )
    gex_keys = set(cell_keys)
    print(
        f"[gex] {expression.obs.shape[0]:,} cells x {expression.var.shape[0]:,} genes"
    )
    if merged_barcodes is not None:
        print(
            f"[merge] {len(merged_keys):,} guide-matrix cells; "
            f"{len(merged_keys & gex_keys):,} overlap GEX"
        )
    else:
        print("[merge] merged guide barcode list not supplied; overlap not computed")

    loaded = _align_assignment_cells(load_assignment(assignment), cell_keys)
    guide_construct = None
    if guide_csv is not None:
        guide_construct = load_guide_construct_map(guide_csv)
        if not guide_construct:
            raise ValueError(f"{guide_csv}: no guide-to-construct mappings found")
    elif guide_design in {"dual", "multi"}:
        raise ValueError(
            f"guide_design={guide_design!r} requires --guide-csv for construct mapping"
        )

    structure, guide_count, construct_count = _classify_assignment_structure(
        loaded, cell_keys, guide_design, guide_construct
    )
    obs = expression.obs.copy()
    obs["assignment_structure"] = structure.to_numpy()
    obs["guide_count"] = guide_count.to_numpy()
    obs["construct_count"] = construct_count.to_numpy()
    top_guide = pd.Series(
        loaded.top1.reindex(cell_keys).to_numpy(),
        index=obs.index,
        dtype="category",
    )
    obs["top_guide"] = top_guide
    # Keep the pre-canonical names for the downstream standardization adapter.
    obs["assigned_guide"] = top_guide.copy()

    adata = AnnData(X=expression.X, obs=obs, var=expression.var)
    # Keep backed source handles alive until the caller has written the output.
    # This attribute is private and is not serialized into the h5ad artifact.
    adata._scprocess_expression_backing = expression.backing
    adata.obs_names = pd.Index(cell_keys, name="cell_id")

    candidate_keys = ["guide_candidates"]
    adata.obsm["guide_candidates"] = _candidate_matrix(loaded, cell_keys)
    adata.uns["guide_candidates_guides"] = np.asarray(loaded.guides, dtype=object)
    print(
        f"[guides] guide_candidates: {adata.obsm['guide_candidates'].shape[0]:,}x"
        f"{adata.obsm['guide_candidates'].shape[1]:,} "
        f"nnz={adata.obsm['guide_candidates'].nnz:,}"
    )

    construct_aggregation = None
    resolved = pd.Series(pd.NA, index=cell_keys, dtype="string")
    if guide_construct is not None:
        adata.uns["guide_construct_map"] = guide_construct
        resolved = loaded.top1.reindex(cell_keys).map(guide_construct)
        resolved = resolved.where(construct_count.eq(1), pd.NA)
        matrix, constructs, unmapped, construct_aggregation = _construct_matrix(
            loaded, cell_keys, guide_construct
        )
        adata.obsm["construct_candidates"] = matrix
        adata.uns["construct_candidates_constructs"] = np.asarray(
            constructs, dtype=object
        )
        candidate_keys.append("construct_candidates")
        print(
            f"[construct] construct_candidates: {matrix.shape[0]:,}x"
            f"{matrix.shape[1]:,} "
            f"nnz={matrix.nnz:,}; unmapped={unmapped:,}"
        )
    resolved = pd.Series(
        pd.Categorical(resolved, categories=sorted(resolved.dropna().unique())),
        index=adata.obs.index,
    )
    adata.obs["resolved_construct"] = resolved
    adata.obs["assigned_construct"] = resolved.copy()
    guide_target = load_guide_target_map(guide_csv) if guide_csv is not None else {}
    adata.obs = _add_standardized_metadata(
        adata.obs,
        top_guide=adata.obs["top_guide"],
        resolved_construct=adata.obs["resolved_construct"],
        guide_count=guide_count,
        construct_count=construct_count,
        guide_target=guide_target,
        source_group=adata.obs.get("source_group"),
    )

    assignment_summary = {
        "score_column": loaded.score_column,
        "score_type": loaded.score_type,
        "n_rows": int(len(loaded.candidates)),
        "n_guides": int(len(loaded.guides)),
        "n_cells": int(loaded.candidates["cell_key"].nunique()),
    }
    adata.uns["assignment_score"] = {
        "column": loaded.score_column,
        "type": loaded.score_type,
    }
    adata.uns["config"] = json.loads(json.dumps(config_snapshot, default=str))
    adata.uns["integration"] = {
        "guide_design": guide_design,
        "cell_universe": "gex",
        "cell_key_rule": (
            "auto: explicit lane + barcode_16mer; lane-prefixed IDs; existing "
            "numeric suffix; otherwise normalized_16mer + merge group suffix"
        ),
        "construct_aggregation": construct_aggregation,
        "assignment_structure": list(ASSIGNMENT_STRUCTURE_CHOICES),
        "assignment": assignment_summary,
        "config": json.loads(json.dumps(config_snapshot, default=str)),
    }
    adata.uns["manifest"] = {
        "generated": _dt.datetime.now().isoformat(timespec="seconds"),
        "script": "integrate_multimodal.py",
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "assignment_rows": int(len(loaded.candidates)),
        "assignment_cells": int(loaded.candidates["cell_key"].nunique()),
        "obsm_candidate_keys": candidate_keys,
        "assignment_structure_counts": {
            str(key): int(value)
            for key, value in structure.value_counts(dropna=False).items()
        },
        "guide_design": guide_design,
        "construct_library_supplied": guide_construct is not None,
        "expression_input_kind": input_kind,
        "counts_layer": counts_layer,
        "counts_origin": expression.counts_origin,
        "normalized_origin": expression.normalized_origin,
        "target_sum": float(target_sum),
        "max_materialized_nnz": int(max_materialized_nnz),
        "counts_source": str(counts_source.resolve()) if counts_source else None,
        "normalized_source": (
            str(normalized_source.resolve()) if normalized_source else None
        ),
        "merged_barcode_count": len(merged_keys) if merged_barcodes is not None else None,
        "merged_barcode_overlap": (
            len(merged_keys & gex_keys) if merged_barcodes is not None else None
        ),
    }
    adata.uns["standardization"] = {
        "adapter": "scprocess-perturb multimodal integration",
        "adapter_contract": "M01",
        "counts_origin": expression.counts_origin,
        "normalized_origin": expression.normalized_origin,
        "target_sum": float(target_sum),
        "source_shape": [int(expression.obs.shape[0]), int(expression.var.shape[0])],
        "output_contract": {
            "X": "log1p(CP10K)",
            "layers/counts": "original integer-valued counts",
        },
    }
    adata._scprocess_expression_data = expression
    return adata


def _metadata_only(adata: AnnData) -> AnnData:
    """Create a small AnnData shell used by the streaming H5AD writer."""
    shell = AnnData(
        X=sp.csr_matrix(adata.shape, dtype=np.float32),
        obs=adata.obs.copy(),
        var=adata.var.copy(),
    )
    for key in adata.obsm.keys():
        shell.obsm[key] = adata.obsm[key]
    shell.uns = adata.uns.copy()
    return shell


def _append_csr(group, block: sp.csr_matrix, offset: int) -> int:
    """Append a CSR block and return its new cumulative nnz."""
    block = block.tocsr()
    block.sum_duplicates()
    block.eliminate_zeros()
    data = block.data.astype(np.float32, copy=False)
    if block.shape[1] > np.iinfo(np.int32).max:
        raise ValueError("CSR column count exceeds int32 index capacity")
    indices = block.indices.astype(np.int32, copy=False)
    old = group["data"].shape[0]
    new = old + data.size
    group["data"].resize((new,))
    group["indices"].resize((new,))
    group["data"][old:new] = data
    group["indices"][old:new] = indices
    end = offset + block.shape[0]
    group["indptr"][offset + 1:end + 1] = (
        block.indptr[1:].astype(np.int64) + old
    )
    return new


def _new_csr_group(parent, name: str, shape: tuple[int, int]):
    if name in parent:
        del parent[name]
    group = parent.create_group(name)
    group.attrs["encoding-type"] = "csr_matrix"
    group.attrs["encoding-version"] = "0.1.0"
    group.attrs["shape"] = np.asarray(shape, dtype=np.int64)
    group.create_dataset(
        "data", shape=(0,), maxshape=(None,), dtype=np.float32,
        chunks=(1_048_576,), compression="lzf", shuffle=True,
    )
    group.create_dataset(
        "indices", shape=(0,), maxshape=(None,), dtype=np.int32,
        chunks=(1_048_576,), compression="lzf", shuffle=True,
    )
    group.create_dataset(
        "indptr", shape=(shape[0] + 1,), dtype=np.int64,
    )
    group["indptr"][0] = 0
    return group


def _as_csr_block(matrix, start: int, stop: int) -> sp.csr_matrix:
    block = matrix[start:stop]
    return block if sp.issparse(block) else sp.csr_matrix(block)


def _stream_standardized_h5ad(
    adata: AnnData,
    expression: ExpressionData,
    output: Path,
    target_sum: float,
    chunk_nnz: int,
) -> None:
    """Write X and counts in bounded-memory CSR chunks."""
    if expression.counts is None:
        raise ValueError("standardized output requires an expression counts source")
    output.parent.mkdir(parents=True, exist_ok=True)
    partial = output.with_name(output.name + ".partial")
    if partial.exists():
        partial.unlink()
    shell = _metadata_only(adata)
    shell.write_h5ad(partial)
    try:
        with h5py.File(partial, "r+") as handle:
            if "X" in handle:
                del handle["X"]
            x_group = _new_csr_group(handle, "X", adata.shape)
            layers = handle.require_group("layers")
            counts_group = _new_csr_group(layers, "counts", adata.shape)
            n_rows = adata.n_obs
            start = 0
            while start < n_rows:
                stop = min(n_rows, start + 8192)
                while True:
                    counts_block = _as_csr_block(expression.counts, start, stop)
                    if expression.X is None:
                        x_block = _normalize_counts(counts_block, target_sum)
                    else:
                        x_block = _as_csr_block(expression.X, start, stop)
                    block_nnz = max(counts_block.nnz, x_block.nnz)
                    if block_nnz <= chunk_nnz or stop - start <= 1:
                        break
                    stop = start + max(
                        1, int((stop - start) * chunk_nnz / block_nnz)
                    )
                if np.any(counts_block.data < 0):
                    raise ValueError("expression counts contain negative values")
                if not np.all(np.equal(counts_block.data, np.floor(counts_block.data))):
                    raise ValueError("expression counts must be integer-valued")
                if x_block.data.size and not np.all(np.isfinite(x_block.data)):
                    raise ValueError("normalized expression contains non-finite values")
                _append_csr(counts_group, counts_block, start)
                _append_csr(x_group, x_block, start)
                print(f"[standardize] rows {stop:,}/{n_rows:,}", flush=True)
                start = stop
            # AnnData readers expect all row pointers to be populated.
            counts_group["indptr"][n_rows] = counts_group["data"].shape[0]
            x_group["indptr"][n_rows] = x_group["data"].shape[0]
        partial.replace(output)
    except Exception:
        if partial.exists():
            partial.unlink()
        raise


def write_standardized_output(
    adata: AnnData,
    expression: ExpressionData,
    output: Path,
    target_sum: float,
    max_materialized_nnz: int,
    stream_chunk_nnz: int,
) -> str:
    """Write the M01 contract, selecting memory or bounded CSR output."""
    counts_nnz = _matrix_nnz(expression.counts)
    normalized_nnz = _matrix_nnz(expression.X)
    total_nnz = counts_nnz + normalized_nnz if counts_nnz >= 0 and normalized_nnz >= 0 else -1
    use_stream = (
        expression.X is None
        or total_nnz < 0
        or total_nnz > max_materialized_nnz
    )
    if use_stream:
        _stream_standardized_h5ad(
            adata, expression, output, target_sum, stream_chunk_nnz
        )
        return "streamed_csr"

    counts = _materialize_csr(expression.counts).astype(np.float32, copy=False)
    if np.any(counts.data < 0) or not np.all(
        np.equal(counts.data, np.floor(counts.data))
    ):
        raise ValueError("expression counts must be non-negative integers")
    normalized = (
        _materialize_csr(expression.X).astype(np.float32, copy=False)
        if expression.X is not None
        else _normalize_counts(counts, target_sum)
    )
    if normalized.shape != counts.shape:
        raise ValueError("normalized expression and counts shapes differ")
    adata.X = normalized
    adata.layers["counts"] = counts
    output.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output)
    return "materialized"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gex", action="append", required=True, metavar="GROUP=PATH")
    parser.add_argument("--assign", required=True, metavar="PATH")
    parser.add_argument(
        "--guide-design", choices=GUIDE_DESIGN_CHOICES,
        help="assignment guide design; inherited from config in the workflow",
    )
    parser.add_argument(
        "--barcodes", default="", help="optional merged_barcodes.tsv.gz"
    )
    parser.add_argument("--guide-csv", default="", help="guide-to-construct library")
    parser.add_argument(
        "--counts-source", default="",
        help="optional aligned h5ad counts source for normalized external input",
    )
    parser.add_argument(
        "--normalized-source", default="",
        help="optional aligned h5ad source whose X is preserved as normalized input",
    )
    parser.add_argument(
        "--input-kind", choices=("auto", "counts", "standardized"), default="auto",
        help="expression input; auto detects raw integer X or layers[counts]",
    )
    parser.add_argument(
        "--counts-layer", default="counts",
        help="counts layer name for standardized h5ad/h5mu input",
    )
    parser.add_argument(
        "--target-sum", type=float, default=10000.0,
        help="target sum for computed log1p(CP10K) normalization",
    )
    parser.add_argument(
        "--max-materialized-nnz", type=int, default=100_000_000,
        help="maximum combined X/counts nnz before bounded CSR streaming",
    )
    parser.add_argument(
        "--stream-chunk-nnz", type=int, default=50_000_000,
        help="approximate CSR nnz limit for one streaming block",
    )
    parser.add_argument(
        "--max-input-nnz", type=int, default=5_000_000_000,
        help="hard input nnz limit checked before processing",
    )
    parser.add_argument(
        "--max-output-gb", type=float, default=250.0,
        help="hard estimated output-size limit",
    )
    parser.add_argument(
        "--min-free-disk-gb", type=float, default=200.0,
        help="required free-disk reserve after estimated output",
    )
    parser.add_argument(
        "--max-process-memory-gb", type=float, default=192.0,
        help="virtual-memory limit applied before reading expression data",
    )
    parser.add_argument("--out", required=True, help="output .h5ad path")
    parser.add_argument("--config-json", default="{}", help="JSON provenance snapshot")
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    groups = _parse_pairs(args.gex, "gex")
    assignment_path = _parse_path(args.assign, "assign")
    barcode_path = Path(args.barcodes) if args.barcodes else None
    if barcode_path is not None and not barcode_path.exists():
        raise FileNotFoundError(f"merged barcode file does not exist: {barcode_path}")
    try:
        config_snapshot = json.loads(args.config_json)
    except json.JSONDecodeError as exc:
        raise ValueError("--config-json must be valid JSON") from exc

    guide_design = args.guide_design
    if guide_design is None:
        guide_design = config_snapshot.get("guide_design")
    if guide_design is None:
        assignment_config = config_snapshot.get("assignment", {})
        guide_design = assignment_config.get("guide_design")
    if guide_design is None:
        raise ValueError(
            "guide design is required; provide --guide-design or "
            "assignment.guide_design in --config-json"
        )
    guide_csv = Path(args.guide_csv) if args.guide_csv else None
    counts_source = _parse_path(args.counts_source, "counts-source") if args.counts_source else None
    normalized_source = (
        _parse_path(args.normalized_source, "normalized-source")
        if args.normalized_source else None
    )
    output = Path(args.out)
    preflight = preflight_resources(
        groups=groups,
        output=output,
        input_kind=args.input_kind,
        counts_layer=args.counts_layer,
        max_input_nnz=args.max_input_nnz,
        max_output_gb=args.max_output_gb,
        min_free_disk_gb=args.min_free_disk_gb,
        max_process_memory_gb=args.max_process_memory_gb,
    )
    adata = integrate(
        groups=groups,
        assignment=assignment_path,
        merged_barcodes=barcode_path,
        guide_design=guide_design,
        config_snapshot=config_snapshot,
        guide_csv=guide_csv,
        input_kind=args.input_kind,
        counts_layer=args.counts_layer,
        target_sum=args.target_sum,
        max_materialized_nnz=args.max_materialized_nnz,
        counts_source=counts_source,
        normalized_source=normalized_source,
    )
    backend = write_standardized_output(
        adata=adata,
        expression=adata._scprocess_expression_data,
        output=output,
        target_sum=args.target_sum,
        max_materialized_nnz=args.max_materialized_nnz,
        stream_chunk_nnz=args.stream_chunk_nnz,
    )
    for backing in getattr(adata, "_scprocess_expression_backing", []) or []:
        if hasattr(backing, "close"):
            backing.close()
        elif hasattr(backing, "file"):
            backing.file.close()
    print(
        f"[done] wrote {output} ({adata.n_obs:,} x {adata.n_vars:,}, "
        f"guide_design={guide_design}, backend={backend})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
