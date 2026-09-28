#!/usr/bin/env python3
"""Read-only compatibility check for expression and assignment cell identifiers."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from collections import defaultdict
from pathlib import Path

import h5py
import numpy as np


BARCODE_RE = re.compile(r"([ACGTN]{16})", re.IGNORECASE)
SUFFIX_GROUP_RE = re.compile(r"-L?0*(\d+)$", re.IGNORECASE)
PREFIX_GROUP_RE = re.compile(r"^(?:lane|l)0*(\d+)[_-]", re.IGNORECASE)
CELL_COLUMNS = ("cell_barcode", "cell_id", "cell", "barcode")
BARCODE_COLUMNS = ("barcode_16mer", "cell_barcode_16mer")
GROUP_COLUMNS = ("lane", "lane_id", "gemgroup", "source_group")
FALLBACK_GROUP_COLUMNS = ("batch_id", "batch", "sample_id", "sample")


def text(value) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, np.generic):
        return str(value.item())
    return str(value)


def read_array(node: h5py.Dataset) -> list[str]:
    if h5py.check_dtype(vlen=node.dtype) is not None or node.dtype.kind == "O":
        return [str(value) for value in node.asstr()[:]]
    values = node[:]
    if values.dtype.kind == "S":
        return [value.decode("utf-8") for value in values]
    return [text(value) for value in values]


def read_obs_column(node: h5py.Dataset | h5py.Group) -> list[str]:
    if isinstance(node, h5py.Group) and "categories" in node and "codes" in node:
        categories = read_array(node["categories"])
        codes = np.asarray(node["codes"][:], dtype=np.int64)
        return ["" if code < 0 else categories[code] for code in codes]
    return read_array(node)


def expression_ids(path: Path, barcode_column: str | None, group_column: str | None):
    suffix = path.suffix.lower()
    with h5py.File(path, "r") as handle:
        if suffix == ".h5mu":
            obs = handle["mod/rna/obs"]
        elif suffix == ".h5ad":
            obs = handle["obs"]
        elif suffix == ".h5":
            ids = [text(value) for value in handle["matrix/barcodes"][:]]
            return ids, None, None, [], "matrix/barcodes"
        else:
            raise ValueError(f"unsupported expression format: {path.suffix}")

        index_name = text(obs.attrs.get("_index", "_index"))
        ids = read_obs_column(obs[index_name])
        available = sorted(str(name) for name in obs.keys() if name != index_name)
        selected_barcode = barcode_column or next(
            (name for name in BARCODE_COLUMNS if name in obs), None
        )
        selected_group = group_column or next(
            (name for name in GROUP_COLUMNS if name in obs), None
        )
        if selected_group is None:
            selected_group = next(
                (name for name in FALLBACK_GROUP_COLUMNS if name in obs), None
            )
        if selected_barcode and selected_barcode not in obs:
            raise ValueError(f"expression obs lacks barcode column {selected_barcode!r}")
        if selected_group and selected_group not in obs:
            raise ValueError(f"expression obs lacks group column {selected_group!r}")
        barcodes = read_obs_column(obs[selected_barcode]) if selected_barcode else None
        groups = read_obs_column(obs[selected_group]) if selected_group else None
        return ids, barcodes, groups, available, index_name


def open_text(path: Path):
    return gzip.open(path, "rt", newline="") if path.suffix == ".gz" else path.open(newline="")


def assignment_ids(path: Path, cell_column: str | None, group_column: str | None):
    with open_text(path) as stream:
        reader = csv.DictReader(stream)
        fields = reader.fieldnames or []
        selected_cell = cell_column or next((name for name in CELL_COLUMNS if name in fields), None)
        if selected_cell is None:
            raise ValueError(f"assignment table has no recognized cell column: {fields}")
        selected_group = group_column if group_column in fields else None
        if group_column and selected_group is None:
            raise ValueError(f"assignment table lacks group column {group_column!r}")
        ids: dict[str, str] = {}
        for row_number, row in enumerate(reader, start=2):
            cell = str(row.get(selected_cell, "")).strip()
            if not cell:
                raise ValueError(f"assignment row {row_number} has an empty cell ID")
            group = str(row.get(selected_group, "")).strip() if selected_group else ""
            previous = ids.get(cell)
            if previous is not None and previous != group:
                raise ValueError(f"assignment cell {cell!r} has conflicting group labels")
            ids[cell] = group
    return list(ids), [ids[cell] for cell in ids], fields, selected_cell, selected_group


def barcode(value: str) -> str | None:
    match = BARCODE_RE.search(str(value).upper())
    return match.group(1) if match else None


def normalized_group(value: str) -> str:
    raw = str(value).strip()
    numeric = re.fullmatch(r"(?:lane|l)?0*(\d+)", raw, flags=re.IGNORECASE)
    return f"lane:{int(numeric.group(1))}" if numeric else raw.lower()


def group_from_id(value: str) -> str | None:
    prefix = PREFIX_GROUP_RE.search(str(value))
    if prefix:
        return f"lane:{int(prefix.group(1))}"
    suffix = SUFFIX_GROUP_RE.search(str(value))
    if suffix:
        return f"lane:{int(suffix.group(1))}"
    return None


def build_keys(ids: list[str], barcodes: list[str] | None, groups: list[str] | None):
    barcode_values = [barcode(value) for value in (barcodes or ids)]
    group_values = []
    group_source = "none"
    if groups is not None:
        group_values = [normalized_group(value) if str(value).strip() else None for value in groups]
        group_source = "metadata"
    else:
        group_values = [group_from_id(value) for value in ids]
        if any(value is not None for value in group_values):
            group_source = "structured_cell_id"
    return barcode_values, group_values, group_source


def indexed_unique(keys: list[tuple[str, str] | str | None], ids: list[str]):
    values: dict[object, list[str]] = defaultdict(list)
    for key, cell in zip(keys, ids):
        if key is not None:
            values[key].append(cell)
    return {key: cells[0] for key, cells in values.items() if len(cells) == 1}, {
        key: cells for key, cells in values.items() if len(cells) > 1
    }


def check(args) -> dict:
    expression, expression_barcodes, expression_groups, obs_columns, index_name = expression_ids(
        args.expression, args.expression_barcode_column, args.expression_group_column
    )
    assignment, assignment_groups, assignment_columns, assignment_cell_column, assignment_group_column = assignment_ids(
        args.assignment, args.assignment_cell_column, args.assignment_group_column
    )
    if len(set(expression)) != len(expression):
        raise ValueError("expression cell IDs are not unique")

    expression_set = set(expression)
    exact = {cell for cell in assignment if cell in expression_set}
    remaining = [cell for cell in assignment if cell not in exact]
    fallback_enabled = expression_barcodes is not None or not exact

    expr_bc, expr_group, expr_group_source = build_keys(
        expression, expression_barcodes, expression_groups
    )
    asn_bc, asn_group, asn_group_source = build_keys(
        assignment, None, assignment_groups if assignment_group_column else None
    )
    asn_position = {cell: index for index, cell in enumerate(assignment)}

    expression_group_keys = [
        (bc, group) if bc and group else None for bc, group in zip(expr_bc, expr_group)
    ]
    assignment_group_keys = [
        (bc, group) if bc and group else None for bc, group in zip(asn_bc, asn_group)
    ]
    unique_expression_group, ambiguous_expression_group = indexed_unique(
        expression_group_keys, expression
    )
    unique_assignment_group, ambiguous_assignment_group = indexed_unique(
        assignment_group_keys, assignment
    )

    group_matches = set()
    used_expression = set(exact)
    if fallback_enabled:
        for cell in remaining:
            key = assignment_group_keys[asn_position[cell]]
            expression_cell = unique_expression_group.get(key)
            if (
                key is not None
                and key in unique_assignment_group
                and expression_cell is not None
                and expression_cell not in used_expression
            ):
                group_matches.add(cell)
                used_expression.add(expression_cell)

    remaining = [cell for cell in remaining if cell not in group_matches]
    unique_expression_bc, ambiguous_expression_bc = indexed_unique(expr_bc, expression)
    unique_assignment_bc, ambiguous_assignment_bc = indexed_unique(asn_bc, assignment)
    unique_matches = set()
    if fallback_enabled:
        for cell in remaining:
            key = asn_bc[asn_position[cell]]
            expression_cell = unique_expression_bc.get(key)
            if (
                key is not None
                and key in unique_assignment_bc
                and expression_cell is not None
                and expression_cell not in used_expression
            ):
                unique_matches.add(cell)
                used_expression.add(expression_cell)

    unresolved = [cell for cell in remaining if cell not in unique_matches]
    expression_by_barcode: dict[str, list[dict[str, str | None]]] = defaultdict(list)
    for cell, bc, group in zip(expression, expr_bc, expr_group):
        if bc:
            expression_by_barcode[bc].append({"cell_id": cell, "group": group})
    unresolved_details = []
    unresolved_without_expression_barcode = 0
    unresolved_with_expression_candidates = 0
    for cell in unresolved[:10]:
        position = asn_position[cell]
        unresolved_details.append(
            {
                "assignment_cell_id": cell,
                "barcode_16mer": asn_bc[position],
                "assignment_group": asn_group[position],
                "expression_candidates": expression_by_barcode.get(asn_bc[position], []),
            }
        )
    for cell in unresolved:
        candidates = (
            expression_by_barcode.get(asn_bc[asn_position[cell]], [])
            if fallback_enabled else []
        )
        if candidates:
            unresolved_with_expression_candidates += 1
        else:
            unresolved_without_expression_barcode += 1
    invalid_expression_barcodes = sum(value is None for value in expr_bc)
    invalid_assignment_barcodes = sum(value is None for value in asn_bc)
    lane_prefixed_expression = bool(expression) and all(
        PREFIX_GROUP_RE.search(cell) and barcode(cell) for cell in expression
    )
    pure_16mer_matches = sum(
        asn_bc[asn_position[cell]] in expression_by_barcode for cell in assignment
    )
    mixed_strategy_matches = len(exact) + len(group_matches) + len(unique_matches)
    if (
        lane_prefixed_expression
        and pure_16mer_matches > mixed_strategy_matches
        and not invalid_assignment_barcodes
    ):
        selected_strategy = "pure_16mer"
        verdict = "PASS"
    elif expression_barcodes is not None and expression_groups is not None:
        selected_strategy = "barcode_lane"
        verdict = "PASS" if not unresolved_with_expression_candidates else "NEEDS_GROUP_INPUT"
    elif exact:
        selected_strategy = "exact"
        verdict = "PASS"
    elif not unresolved_with_expression_candidates and not invalid_assignment_barcodes:
        selected_strategy = "unique_16mer"
        verdict = "PASS"
    elif unresolved_with_expression_candidates:
        selected_strategy = "unresolved"
        verdict = "NEEDS_GROUP_INPUT"
    else:
        selected_strategy = "unresolved"
        verdict = "FAIL"

    return {
        "verdict": verdict,
        "selected_strategy": selected_strategy,
        "expression": str(args.expression),
        "assignment": str(args.assignment),
        "expression_cells": len(expression),
        "assignment_cells": len(assignment),
        "exact_id_matches": len(exact),
        "fallback_matching_enabled": fallback_enabled,
        "pure_16mer_matches": pure_16mer_matches,
        "pure_16mer_duplicate_rows_removed": len(expression) - len(set(expr_bc)),
        "explicit_or_structured_group_matches": len(group_matches),
        "unique_16mer_matches": len(unique_matches),
        "unresolved_assignment_cells": len(unresolved),
        "unresolved_with_expression_candidates": unresolved_with_expression_candidates,
        "assignment_cells_absent_from_expression": unresolved_without_expression_barcode,
        "expression_cells_without_assignment": len(expression) - len(used_expression),
        "invalid_expression_barcodes": invalid_expression_barcodes,
        "invalid_assignment_barcodes": invalid_assignment_barcodes,
        "ambiguous_expression_16mers": len(ambiguous_expression_bc),
        "ambiguous_assignment_16mers": len(ambiguous_assignment_bc),
        "ambiguous_expression_group_keys": len(ambiguous_expression_group),
        "ambiguous_assignment_group_keys": len(ambiguous_assignment_group),
        "expression_index": index_name,
        "expression_barcode_column": args.expression_barcode_column or next(
            (name for name in BARCODE_COLUMNS if name in obs_columns), None
        ),
        "expression_group_column": args.expression_group_column or next(
            (name for name in GROUP_COLUMNS + FALLBACK_GROUP_COLUMNS if name in obs_columns), None
        ),
        "expression_group_source": expr_group_source,
        "assignment_cell_column": assignment_cell_column,
        "assignment_group_column": assignment_group_column,
        "assignment_group_source": asn_group_source,
        "unresolved_examples": unresolved[:10],
        "unresolved_details": unresolved_details,
        "ambiguous_expression_16mer_examples": list(ambiguous_expression_bc)[:10],
        "available_expression_obs_columns": obs_columns,
        "available_assignment_columns": assignment_columns,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expression", type=Path, required=True)
    parser.add_argument("--assignment", type=Path, required=True)
    parser.add_argument("--expression-barcode-column")
    parser.add_argument("--expression-group-column")
    parser.add_argument("--assignment-cell-column")
    parser.add_argument("--assignment-group-column")
    args = parser.parse_args()
    report = check(args)
    print(json.dumps(report, indent=2, ensure_ascii=False))
    raise SystemExit(0 if report["verdict"] == "PASS" else 2)


if __name__ == "__main__":
    main()
