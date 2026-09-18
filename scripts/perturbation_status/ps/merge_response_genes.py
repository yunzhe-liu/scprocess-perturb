#!/usr/bin/env python3
"""Merge per-shard PS response-gene selections into one pooled union."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--selection-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--expected-targets", type=int, required=True)
    args = parser.parse_args()
    files = sorted(args.selection_dir.glob("shard_*_response_genes.tsv"))
    if not files:
        raise ValueError("No response-gene selections found")
    table = pd.concat([pd.read_csv(path, sep="\t") for path in files], ignore_index=True)
    if table["target_label"].nunique() != args.expected_targets:
        raise ValueError("Response-gene target coverage is incomplete")
    if table.duplicated(["target_label", "response_gene"]).any():
        raise ValueError("Duplicated target-response gene records")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    table.sort_values(["target_label", "response_gene"]).to_csv(
        args.output_dir / "response_gene_manifest.tsv", sep="\t", index=False
    )
    union = sorted(table["response_gene"].unique())
    (args.output_dir / "global_response_genes.txt").write_text("\n".join(union) + "\n", encoding="utf-8")
    summary = {"targets": int(table["target_label"].nunique()), "target_gene_pairs": len(table), "global_response_genes": len(union)}
    (args.output_dir / "response_gene_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
