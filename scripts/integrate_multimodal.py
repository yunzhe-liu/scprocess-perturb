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

  X                         expression matrix, cells x genes
  obs.cell_key              canonical integration key
  obs.assigned_guide          top-ranked guide
  obsm.guide_candidates       sparse native-score matrix (guide_full mode)
  obs.assigned_construct      top-ranked construct (construct mode)
  obsm.construct_candidates   sparse construct-score matrix (construct mode)
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
MODE_CHOICES = ("guide_top1", "guide_full", "construct")


@dataclass
class ExpressionData:
    X: object
    obs: pd.DataFrame
    var: pd.DataFrame
    backing: object | None = None


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


def _read_h5ad_group(group: h5py.Group) -> ExpressionData:
    X = _read_sparse_or_dense(group["X"])
    obs = _read_h5_dataframe(group["obs"])
    var = _read_h5_dataframe(group["var"])
    if X.shape[0] != len(obs) or X.shape[1] != len(var):
        raise ValueError(
            f"expression dimensions {X.shape} do not match obs/var "
            f"({len(obs)}, {len(var)})"
        )
    return ExpressionData(X=X.tocsr(), obs=obs, var=var)


def _read_expression(path: Path) -> ExpressionData:
    """Read h5ad, h5mu RNA, or the scprocess raw H5 matrix."""
    suffix = path.suffix.lower()
    if suffix == ".h5ad":
        # Keep large sparse h5ad inputs backed.  This avoids materializing all
        # X data before the output writer can stream it to the new artifact.
        backed = ad.read_h5ad(path, backed="r")
        return ExpressionData(
            X=backed.X,
            obs=backed.obs.copy(),
            var=backed.var.copy(),
            backing=backed,
        )

    if suffix == ".h5mu":
        with h5py.File(path, "r") as handle:
            if "mod" not in handle or "rna" not in handle["mod"]:
                raise ValueError(f"{path}: h5mu does not contain mod/rna")
            return _read_h5ad_group(handle["mod"]["rna"])

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
        return ExpressionData(X=X, obs=obs, var=var)

    raise ValueError(f"unsupported expression format {path}; expected .h5/.h5ad/.h5mu")


def _normalize_barcode(cell_id: str) -> str:
    match = _BARCODE_RE.search(str(cell_id).upper())
    if match is None:
        raise ValueError(f"cannot find a 16-base barcode in expression cell id {cell_id!r}")
    return match.group(1)


def _group_suffix(group: str) -> str:
    match = re.search(r"(\d+)$", group)
    return f"-L{match.group(1)}" if match else f"-{group}"


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


def load_expression(groups: list[tuple[str, Path]]) -> tuple[ExpressionData, list[str]]:
    if not groups:
        raise ValueError("at least one --gex group=path is required")

    matrices = []
    observations = []
    lane_prefixed_input = True
    first_var = None
    first_var_names = None
    backings = []
    for group, path in groups:
        data = _read_expression(path)
        if data.backing is not None:
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
        matrices.append(data.X)
        observations.append(obs)

    obs = pd.concat(observations, axis=0, sort=False)
    if len(matrices) == 1:
        X = matrices[0]
    else:
        X = sp.vstack(matrices, format="csr")
    if not obs.index.is_unique:
        # Lane-prefixed merged h5ad inputs can contain the same 16mer in more
        # than one lane.  Keep the first row in file order, which is stable
        # and preserves the same cell universe used by the assignment table.
        if not lane_prefixed_input:
            raise ValueError("expression inputs produce duplicate integration cell keys")
        keep = ~obs.index.duplicated(keep="first")
        obs = obs.iloc[np.flatnonzero(keep)].copy()
        X = X[keep]
    if not obs.index.is_unique:
        raise ValueError("expression inputs produce duplicate integration cell keys")
    return ExpressionData(X=X, obs=obs, var=first_var, backing=backings), list(obs.index)


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
                sanitized = _sanitize_guide_id(raw_guide)
                aliases = [raw_guide] if sanitized == raw_guide else [raw_guide, sanitized]
                for guide in aliases:
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
            sanitized = _sanitize_guide_id(raw_guide)
            aliases = [raw_guide] if sanitized == raw_guide else [raw_guide, sanitized]
            for guide in aliases:
                if guide:
                    add_mapping(guide, construct)
    return mapping


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
    work = assign.candidates[assign.candidates["cell_key"].isin(cell_pos)].copy()
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
    for row in assign.candidates.itertuples(index=False):
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


def integrate(
    groups: list[tuple[str, Path]],
    assignment: Path,
    merged_barcodes: Path | None,
    mode: str,
    config_snapshot: dict,
    guide_csv: Path | None = None,
) -> AnnData:
    if mode not in MODE_CHOICES:
        raise ValueError(f"unknown integration mode {mode!r}; choose from {MODE_CHOICES}")
    expression, cell_keys = load_expression(groups)
    merged_keys = (
        _read_merged_barcodes(merged_barcodes)
        if merged_barcodes is not None
        else set()
    )
    gex_keys = set(cell_keys)
    print(
        f"[gex] {expression.X.shape[0]:,} cells x {expression.X.shape[1]:,} genes"
    )
    if merged_barcodes is not None:
        print(
            f"[merge] {len(merged_keys):,} guide-matrix cells; "
            f"{len(merged_keys & gex_keys):,} overlap GEX"
        )
    else:
        print("[merge] merged guide barcode list not supplied; overlap not computed")

    loaded = _align_assignment_cells(load_assignment(assignment), cell_keys)
    obs = expression.obs.copy()
    obs["assigned_guide"] = pd.Series(
        loaded.top1.reindex(cell_keys).to_numpy(),
        index=obs.index,
        dtype="category",
    )

    adata = AnnData(X=expression.X, obs=obs, var=expression.var)
    # Keep backed source handles alive until the caller has written the output.
    # This attribute is private and is not serialized into the h5ad artifact.
    adata._scprocess_expression_backing = expression.backing
    adata.obs_names = pd.Index(cell_keys, name="cell_id")

    candidate_keys = []
    if mode == "guide_full":
        key = "guide_candidates"
        adata.obsm[key] = _candidate_matrix(loaded, cell_keys)
        adata.uns[f"{key}_guides"] = np.asarray(loaded.guides, dtype=object)
        candidate_keys.append(key)
        print(
            f"[guide_full] {key}: {adata.obsm[key].shape[0]:,}x"
            f"{adata.obsm[key].shape[1]:,} nnz={adata.obsm[key].nnz:,}"
        )

    construct_aggregation = None
    if mode == "construct":
        if guide_csv is None:
            raise ValueError("construct mode requires --guide-csv")
        guide_construct = load_guide_construct_map(guide_csv)
        if not guide_construct:
            raise ValueError(f"{guide_csv}: no guide-to-construct mappings found")
        adata.uns["guide_construct_map"] = guide_construct
        adata.obs["assigned_construct"] = pd.Series(
            adata.obs["assigned_guide"].map(guide_construct).to_numpy(),
            index=adata.obs.index,
            dtype="category",
        )
        key = "construct_candidates"
        matrix, constructs, unmapped, construct_aggregation = _construct_matrix(
            loaded, cell_keys, guide_construct
        )
        adata.obsm[key] = matrix
        adata.uns[f"{key}_constructs"] = np.asarray(constructs, dtype=object)
        candidate_keys.append(key)
        print(
            f"[construct] {key}: {matrix.shape[0]:,}x{matrix.shape[1]:,} "
            f"nnz={matrix.nnz:,}; unmapped={unmapped:,}"
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
        "mode": mode,
        "cell_universe": "gex",
        "cell_key_rule": (
            "auto: explicit lane + barcode_16mer; lane-prefixed IDs; existing "
            "numeric suffix; otherwise normalized_16mer + merge group suffix"
        ),
        "construct_aggregation": construct_aggregation,
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
        "merged_barcode_count": len(merged_keys) if merged_barcodes is not None else None,
        "merged_barcode_overlap": (
            len(merged_keys & gex_keys) if merged_barcodes is not None else None
        ),
    }
    return adata


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gex", action="append", required=True, metavar="GROUP=PATH")
    parser.add_argument("--assign", required=True, metavar="PATH")
    parser.add_argument(
        "--barcodes", default="", help="optional merged_barcodes.tsv.gz"
    )
    parser.add_argument(
        "--mode", "--assignment-mode", dest="mode", choices=MODE_CHOICES,
        required=True,
    )
    parser.add_argument("--guide-csv", default="", help="required for construct mode")
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

    guide_csv = Path(args.guide_csv) if args.guide_csv else None
    if args.mode == "construct" and guide_csv is None:
        raise ValueError("construct mode requires --guide-csv")
    adata = integrate(
        groups=groups,
        assignment=assignment_path,
        merged_barcodes=barcode_path,
        mode=args.mode,
        config_snapshot=config_snapshot,
        guide_csv=guide_csv,
    )
    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output)
    print(
        f"[done] wrote {output} ({adata.n_obs:,} x {adata.n_vars:,}, "
        f"mode={args.mode})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
