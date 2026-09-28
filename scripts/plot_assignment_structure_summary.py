#!/usr/bin/env python3
"""Plot assignment-structure composition for a small set of datasets.

The plot uses the original assignment CSV for the total number of unique
assignment cells and the integration H5AD ``obs`` table for the three
assignment-structure counts.  The expression matrix is never loaded.

Example
-------
python scripts/plot_assignment_structure_summary.py \
  --dataset Papalexi_2021=/path/assignment.csv=/path/perturbation_adata.h5ad \
  --dataset Norman_2019=/path/assignment.csv=/path/perturbation_adata.h5ad \
  --dataset Replogle_2022=/path/assignment.csv=/path/perturbation_adata.h5ad \
  --output-dir /path/to/output
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


STRUCTURES = (
    "single_guide",
    "concordant_construct",
    "mixed_construct",
)
STRUCTURE_LABELS = {
    "single_guide": "single guide",
    "concordant_construct": "concordant construct",
    "mixed_construct": "mixed construct",
}
# Muted, colour-blind-friendly colours suitable for a publication figure.
COLORS = {
    "single_guide": "#4E79A7",
    "concordant_construct": "#59A14F",
    "mixed_construct": "#E15759",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dataset",
        action="append",
        required=True,
        metavar="NAME=ASSIGNMENT_CSV=INTEGRATION_H5AD",
        help="Dataset name and the two input paths; repeat exactly three times.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory for the PNG figure and summary CSV.",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=250_000,
        help="Rows per assignment CSV read chunk (default: 250000).",
    )
    return parser.parse_args()


def parse_dataset_spec(spec: str) -> tuple[str, Path, Path]:
    parts = spec.split("=", 2)
    if len(parts) != 3 or not all(parts):
        raise ValueError(
            "--dataset must have the form NAME=ASSIGNMENT_CSV=INTEGRATION_H5AD"
        )
    name, assignment, h5ad = parts
    return name, Path(assignment), Path(h5ad)


def count_assignment_cells(path: Path, chunk_size: int) -> int:
    """Count unique assignment cells without loading the whole CSV."""
    cells: set[str] = set()
    try:
        chunks = pd.read_csv(
            path,
            usecols=["cell"],
            dtype={"cell": "string"},
            keep_default_na=False,
            chunksize=chunk_size,
        )
        for chunk in chunks:
            cells.update(chunk["cell"].astype("string").str.strip().tolist())
    except ValueError as exc:
        raise ValueError(f"{path}: assignment CSV must contain a 'cell' column") from exc
    cells.discard("")
    return len(cells)


def decode_scalar(value: object) -> object:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, np.generic):
        return value.item()
    return value


def read_categorical_values(node: h5py.Dataset | h5py.Group) -> list[object]:
    """Read one AnnData obs column from either categorical or plain storage."""
    if isinstance(node, h5py.Group) and {"categories", "codes"}.issubset(node):
        categories = [decode_scalar(value) for value in node["categories"][()]]
        codes = node["codes"][()]
        return [
            categories[int(code)] if int(code) >= 0 else None
            for code in codes
        ]
    if isinstance(node, h5py.Dataset):
        return [decode_scalar(value) for value in node[()]]
    raise ValueError("unsupported AnnData obs column encoding")


def count_structures(path: Path) -> Counter[str]:
    """Count the three structure labels from H5AD obs without reading X."""
    with h5py.File(path, "r") as handle:
        if "obs" not in handle or "assignment_structure" not in handle["obs"]:
            raise ValueError(f"{path}: missing obs/assignment_structure")
        values = read_categorical_values(handle["obs"]["assignment_structure"])
    counts = Counter(str(value) for value in values if value is not None)
    return Counter({structure: counts.get(structure, 0) for structure in STRUCTURES})


def collect_records(specs: list[str], chunk_size: int) -> list[dict[str, object]]:
    records = []
    seen_names: set[str] = set()
    for spec in specs:
        name, assignment_path, h5ad_path = parse_dataset_spec(spec)
        if name in seen_names:
            raise ValueError(f"duplicate dataset name: {name}")
        seen_names.add(name)
        assignment_total = count_assignment_cells(assignment_path, chunk_size)
        structure_counts = count_structures(h5ad_path)
        classified_total = sum(structure_counts.values())
        unmatched = assignment_total - classified_total
        if unmatched < 0:
            raise ValueError(
                f"{name}: classified structure cells ({classified_total:,}) exceed "
                f"assignment cells ({assignment_total:,})"
            )
        denominator = classified_total
        if denominator == 0:
            raise ValueError(f"{name}: no classified assignment cells found")
        record = {
            "dataset": name,
            "assignment_cell_total": assignment_total,
            "single_guide": structure_counts["single_guide"],
            "concordant_construct": structure_counts["concordant_construct"],
            "mixed_construct": structure_counts["mixed_construct"],
            "classified_cell_total": classified_total,
            "unmatched_assignment_cells": unmatched,
        }
        for structure in STRUCTURES:
            record[f"{structure}_percent"] = (
                100.0 * structure_counts[structure] / denominator
            )
        records.append(record)
    return records


def plot_records(records: list[dict[str, object]], output_path: Path) -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 11,
            "axes.labelsize": 9,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    fig, ax = plt.subplots(figsize=(7.2, 4.8), dpi=300)
    x = np.arange(len(records), dtype=float)
    bottom = np.zeros(len(records), dtype=float)
    width = 0.40

    for structure in STRUCTURES:
        values = np.asarray(
            [record[f"{structure}_percent"] for record in records], dtype=float
        )
        bars = ax.bar(
            x,
            values,
            width=width,
            bottom=bottom,
            color=COLORS[structure],
            edgecolor="none",
            linewidth=0,
            label=STRUCTURE_LABELS[structure],
        )
        for bar, value, record in zip(bars, values, records):
            count = int(record[structure])
            if value >= 7.0:
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_y() + bar.get_height() / 2,
                    f"{count:,}",
                    ha="center",
                    va="center",
                    color="white",
                    fontsize=8.5,
                    fontweight="500",
                )
        bottom += values

    for xpos, record in zip(x, records):
        ax.text(
            xpos,
            104.0,
            f"N = {int(record['classified_cell_total']):,}",
            ha="center",
            va="bottom",
            color="#222222",
            fontsize=9,
            fontweight="500",
        )
    ax.set_title("Assignment structure composition across datasets", pad=16, loc="left")
    ax.set_ylabel("Assignment structure composition (%)")
    ax.set_ylim(0, 111)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_xticks(x)
    ax.set_xticklabels([record["dataset"] for record in records])
    ax.grid(axis="y", color="#B8B8B8", linewidth=0.55, alpha=0.45)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#4A4A4A")
    ax.spines["bottom"].set_color("#4A4A4A")
    ax.tick_params(axis="both", length=3, width=0.6, color="#4A4A4A")
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.14),
        ncol=3,
        frameon=False,
        handlelength=1.2,
        columnspacing=1.8,
        handletextpad=0.5,
    )
    fig.text(
        0.01,
        0.035,
        "Bar height is normalized among retained classified cells; labels show absolute cell counts. "
        "N is the final cell count retained after expression-assignment alignment.",
        ha="left",
        va="bottom",
        color="#5A5A5A",
        fontsize=7.5,
    )
    fig.text(
        0.01,
        0.012,
        "Norman: assignment features are construct-level GBC labels; per-sgRNA labels are unavailable, "
        "so its bar is not a guide-level concordance measurement.",
        ha="left",
        va="bottom",
        color="#5A5A5A",
        fontsize=7.5,
    )
    fig.subplots_adjust(left=0.12, right=0.98, top=0.86, bottom=0.27)
    fig.savefig(output_path, dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def write_summary(records: list[dict[str, object]], output_path: Path) -> None:
    fieldnames = list(records[0].keys())
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)


def main() -> None:
    args = parse_args()
    if len(args.dataset) != 3:
        raise ValueError("provide exactly three --dataset specifications")
    records = collect_records(args.dataset, args.chunk_size)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    plot_records(records, args.output_dir / "assignment_structure_summary.png")
    write_summary(records, args.output_dir / "assignment_structure_summary.csv")
    for record in records:
        print(
            f"{record['dataset']}: assignment={record['assignment_cell_total']:,}; "
            f"single={record['single_guide']:,}; "
            f"concordant={record['concordant_construct']:,}; "
            f"mixed={record['mixed_construct']:,}; "
            f"unmatched={record['unmatched_assignment_cells']:,}"
        )


if __name__ == "__main__":
    main()
