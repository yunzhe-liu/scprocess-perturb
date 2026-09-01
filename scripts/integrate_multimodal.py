#!/usr/bin/env python3
"""Fuse per-lane expression matrices with guide assignments.

The integration step is the terminal export for scprocess-perturb.  It reads:

  * one expression matrix per group (``group=path``),
  * the merged guide barcode list, and
  * one standardized assignment CSV per method (``method=path``).

The output is a new AnnData file.  Existing GEX, MEX, and assignment files are
never modified.

Cell identity follows the same convention as ``rules/merge.smk``: normalize a
16-base barcode from the expression matrix and append the group suffix
(``-L01`` for a group ending in ``01``).  The GEX cells form the output cell
universe; assignment rows are left-joined onto that universe.

Output layout:

  X                         expression matrix, cells x genes
  obs.cell_key              canonical integration key
  obs.assigned_<method>      top-ranked guide for every method
  obsm.guide_candidates_*    sparse native-score matrices (guide_full mode)
  obs.assigned_construct_*   top-ranked construct (construct mode)
  obsm.construct_candidates_* sparse construct-score matrices (construct mode)
  uns.method_weights        score column used by each candidate matrix
  uns.<candidate>_guides    column order for each candidate matrix
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
from anndata import AnnData


_BARCODE_RE = re.compile(r"([ACGTN]{16})", re.IGNORECASE)
_GUIDE_ID_SANITIZE_RE = re.compile(r"[^a-zA-Z0-9_-]")
MODE_CHOICES = ("guide_top1", "guide_full", "construct")


@dataclass
class ExpressionData:
    X: sp.spmatrix
    obs: pd.DataFrame
    var: pd.DataFrame


@dataclass
class AssignmentData:
    method: str
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
        with h5py.File(path, "r") as handle:
            return _read_h5ad_group(handle)

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


def load_expression(groups: list[tuple[str, Path]]) -> tuple[ExpressionData, list[str]]:
    if not groups:
        raise ValueError("at least one --gex group=path is required")

    matrices = []
    observations = []
    first_var = None
    first_var_names = None
    for group, path in groups:
        data = _read_expression(path)
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
        suffix = _group_suffix(group)
        cell_keys = [f"{_normalize_barcode(cell)}{suffix}" for cell in original_ids]
        if len(set(cell_keys)) != len(cell_keys):
            raise ValueError(f"{path}: duplicate integration cell keys within group {group!r}")

        obs = data.obs.copy()
        obs["cell_key"] = cell_keys
        obs["source_group"] = group
        obs["source_cell_id"] = original_ids
        obs.index = pd.Index(cell_keys, name="cell_id")
        matrices.append(data.X)
        observations.append(obs)

    obs = pd.concat(observations, axis=0, sort=False)
    if not obs.index.is_unique:
        raise ValueError("expression inputs produce duplicate integration cell keys")
    X = sp.vstack(matrices, format="csr")
    return ExpressionData(X=X, obs=obs, var=first_var), list(obs.index)


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
                guide = _sanitize_guide_id(row[column])
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
        guide = _sanitize_guide_id(row[guide_column])
        construct = str(row[construct_column]).strip()
        if guide and construct:
            add_mapping(guide, construct)
    return mapping


def load_assignment(method: str, path: Path) -> AssignmentData:
    frame = pd.read_csv(path)
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
        "cell_key": frame[cell_column].astype(str),
        "guide": frame[guide_column].astype(str),
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
        method=method,
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
    assignments: list[tuple[str, Path]],
    merged_barcodes: Path,
    mode: str,
    config_snapshot: dict,
    guide_csv: Path | None = None,
) -> AnnData:
    if mode not in MODE_CHOICES:
        raise ValueError(f"unknown integration mode {mode!r}; choose from {MODE_CHOICES}")
    expression, cell_keys = load_expression(groups)
    merged_keys = _read_merged_barcodes(merged_barcodes)
    gex_keys = set(cell_keys)
    print(
        f"[gex] {expression.X.shape[0]:,} cells x {expression.X.shape[1]:,} genes"
    )
    print(
        f"[merge] {len(merged_keys):,} guide-matrix cells; "
        f"{len(merged_keys & gex_keys):,} overlap GEX"
    )

    loaded = [load_assignment(method, path) for method, path in assignments]
    obs = expression.obs.copy()
    for assignment in loaded:
        obs[f"assigned_{assignment.method}"] = (
            assignment.top1.reindex(cell_keys).to_numpy()
        )

    adata = AnnData(X=expression.X, obs=obs, var=expression.var)
    adata.obs_names = pd.Index(cell_keys, name="cell_id")

    method_weights = {}
    candidate_keys = []
    if mode == "guide_full":
        for assignment in loaded:
            key = f"guide_candidates_{assignment.method}"
            adata.obsm[key] = _candidate_matrix(assignment, cell_keys)
            adata.uns[f"{key}_guides"] = np.asarray(assignment.guides, dtype=object)
            method_weights[assignment.method] = assignment.score_column
            candidate_keys.append(key)
            print(
                f"[guide_full] {key}: {adata.obsm[key].shape[0]:,}x"
                f"{adata.obsm[key].shape[1]:,} nnz={adata.obsm[key].nnz:,}"
            )

    construct_aggregation = {}
    if mode == "construct":
        if guide_csv is None:
            raise ValueError("construct mode requires --guide-csv")
        guide_construct = load_guide_construct_map(guide_csv)
        if not guide_construct:
            raise ValueError(f"{guide_csv}: no guide-to-construct mappings found")
        adata.uns["guide_construct_map"] = guide_construct
        for assignment in loaded:
            guide_top1 = assignment.top1.reindex(cell_keys)
            adata.obs[f"assigned_construct_{assignment.method}"] = guide_top1.map(
                guide_construct
            ).to_numpy()
            key = f"construct_candidates_{assignment.method}"
            matrix, constructs, unmapped, aggregation = _construct_matrix(
                assignment, cell_keys, guide_construct
            )
            adata.obsm[key] = matrix
            adata.uns[f"{key}_constructs"] = np.asarray(constructs, dtype=object)
            construct_aggregation[assignment.method] = aggregation
            method_weights[assignment.method] = assignment.score_column
            candidate_keys.append(key)
            print(
                f"[construct] {key}: {matrix.shape[0]:,}x{matrix.shape[1]:,} "
                f"nnz={matrix.nnz:,}; unmapped={unmapped:,}"
            )

    assignment_summary = {
        assignment.method: {
            "score_column": assignment.score_column,
            "score_type": assignment.score_type,
            "n_rows": int(len(assignment.candidates)),
            "n_guides": int(len(assignment.guides)),
            "n_cells": int(assignment.candidates["cell_key"].nunique()),
        }
        for assignment in loaded
    }
    adata.uns["method_weights"] = method_weights
    adata.uns["config"] = json.loads(json.dumps(config_snapshot, default=str))
    adata.uns["integration"] = {
        "mode": mode,
        "cell_universe": "gex",
        "cell_key_rule": "normalized_16mer + merge group suffix",
        "construct_aggregation": construct_aggregation if mode == "construct" else None,
        "assignment_methods": assignment_summary,
        "config": json.loads(json.dumps(config_snapshot, default=str)),
    }
    adata.uns["manifest"] = {
        "generated": _dt.datetime.now().isoformat(timespec="seconds"),
        "script": "integrate_multimodal.py",
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "assignment_methods": [method for method, _ in assignments],
        "obsm_candidate_keys": candidate_keys,
        "merged_barcode_count": len(merged_keys),
        "merged_barcode_overlap": len(merged_keys & gex_keys),
    }
    return adata


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gex", action="append", required=True, metavar="GROUP=PATH")
    parser.add_argument("--assign", action="append", required=True, metavar="METHOD=PATH")
    parser.add_argument("--barcodes", required=True, help="merged_barcodes.tsv.gz")
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
    assignments = _parse_pairs(args.assign, "assign")
    barcode_path = Path(args.barcodes)
    if not barcode_path.exists():
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
        assignments=assignments,
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
