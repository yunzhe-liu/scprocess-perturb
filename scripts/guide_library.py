#!/usr/bin/env python3
"""Normalize supported guide-library tables to one long-form contract."""

from __future__ import annotations

from pathlib import Path
import re

import pandas as pd


CANONICAL_COLUMNS = [
    "construct_id", "guide_id", "guide_sequence", "target_label", "guide_position"
]
GUIDE_ID_SANITIZE_RE = re.compile(r"[^a-zA-Z0-9_-]")


def sanitize_guide_id(value: str) -> str:
    """Apply the guide-ID normalization used by generated references."""
    return GUIDE_ID_SANITIZE_RE.sub("_", str(value).strip())


def guide_id_aliases(value: str) -> list[str]:
    """Return raw and generated-reference forms of one guide identifier."""
    raw = str(value).strip()
    return list(dict.fromkeys((raw, raw.replace(",", "_"), sanitize_guide_id(raw))))


def _first(columns: set[str], choices: tuple[str, ...]) -> str | None:
    return next((name for name in choices if name in columns), None)


def normalize_guide_library(path: str | Path, *, require_sequence: bool) -> pd.DataFrame:
    """Read a long or legacy-wide guide library into canonical long form.

    Shared guide IDs are retained as separate construct memberships. This
    adapter only normalizes representation; it does not change construct
    resolution semantics downstream.
    """
    path = Path(path)
    frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    columns = set(frame.columns)
    duplicate_col = "duplicated guide pair?"
    if duplicate_col in columns:
        duplicate = frame[duplicate_col].str.strip().str.lower().isin(
            {"true", "yes", "1"}
        )
        frame = frame.loc[~duplicate].copy()
        columns = set(frame.columns)

    target_col = _first(columns, (
        "target_label", "target_gene_name", "target", "gene", "gene_symbol",
        "perturbation",
    ))
    construct_col = _first(columns, (
        "construct_id", "pair_id", "unique sgRNA pair ID",
    ))

    rows: list[dict[str, str]] = []
    if {"sgID_A", "sgID_B"}.issubset(columns):
        sequence_a = _first(columns, ("targeting sequence A", "guide_sequence_A", "sequence_A"))
        sequence_b = _first(columns, ("targeting sequence B", "guide_sequence_B", "sequence_B"))
        if require_sequence and (sequence_a is None or sequence_b is None):
            raise ValueError(f"{path}: dual-guide reference generation requires both guide sequences")
        for record in frame.to_dict(orient="records"):
            for position, guide_col, sequence_col in (
                ("1", "sgID_A", sequence_a), ("2", "sgID_B", sequence_b)
            ):
                guide = str(record.get(guide_col, "")).strip()
                if guide:
                    rows.append({
                        "construct_id": str(record.get(construct_col, "")).strip() if construct_col else "",
                        "guide_id": guide,
                        "guide_sequence": str(record.get(sequence_col, "")).strip() if sequence_col else "",
                        "target_label": str(record.get(target_col, "")).strip() if target_col else "",
                        "guide_position": position,
                    })
    else:
        guide_col = _first(columns, ("guide_id", "gRNA", "sgRNA_id", "sgRNA"))
        if guide_col is None:
            raise ValueError(
                f"{path}: guide library requires guide_id or legacy sgID_A/sgID_B columns"
            )
        sequence_col = _first(columns, (
            "guide_sequence", "sequence", "targeting_sequence", "targeting sequence",
        ))
        position_col = _first(columns, ("guide_position", "position"))
        if require_sequence and sequence_col is None:
            raise ValueError(f"{path}: reference generation requires guide_sequence")
        for record in frame.to_dict(orient="records"):
            guide = str(record.get(guide_col, "")).strip()
            if guide:
                rows.append({
                    "construct_id": str(record.get(construct_col, "")).strip() if construct_col else "",
                    "guide_id": guide,
                    "guide_sequence": str(record.get(sequence_col, "")).strip() if sequence_col else "",
                    "target_label": str(record.get(target_col, "")).strip() if target_col else "",
                    "guide_position": str(record.get(position_col, "")).strip() if position_col else "",
                })

    result = pd.DataFrame(rows, columns=CANONICAL_COLUMNS)
    if result.empty:
        raise ValueError(f"{path}: guide library contains no usable guide rows")
    if require_sequence and result["guide_sequence"].str.strip().eq("").any():
        raise ValueError(f"{path}: every guide requires a non-empty guide_sequence")

    # A guide may belong to multiple constructs, but its physical sequence
    # must remain identical wherever it appears.
    sequence_counts = (
        result.loc[result["guide_sequence"].str.strip().ne("")]
        .groupby("guide_id")["guide_sequence"].nunique()
    )
    conflicting = sequence_counts[sequence_counts > 1]
    if not conflicting.empty:
        raise ValueError(
            f"{path}: guide IDs have conflicting sequences: "
            + ", ".join(conflicting.index.astype(str)[:10])
        )
    return result.drop_duplicates(CANONICAL_COLUMNS).reset_index(drop=True)
