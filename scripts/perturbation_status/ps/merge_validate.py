#!/usr/bin/env python3
"""Merge PS shards and compare them with the canonical pooled result."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--shard-dir", type=Path, required=True)
    parser.add_argument("--barcode", type=Path, required=True)
    parser.add_argument("--canonical", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--tolerance", type=float, default=1e-8)
    parser.add_argument("--excluded-targets", type=Path)
    args = parser.parse_args()
    files = sorted(args.shard_dir.glob("shard_*/ps_target_scores.tsv"))
    if not files:
        raise ValueError("No completed shard scores")
    merged = pd.concat([pd.read_csv(path, sep="\t") for path in files], ignore_index=True)
    if merged["cell_id"].duplicated().any():
        raise ValueError("Duplicated target cells")
    barcode = pd.read_csv(args.barcode, sep="\t")
    canonical_output = barcode[["cell", "gene"]].rename(columns={"cell": "cell_id", "gene": "target_label"})
    canonical_output["is_ntc"] = canonical_output["target_label"] == "non-targeting"
    canonical_output = canonical_output.merge(merged, on=["cell_id", "target_label"], how="left")
    canonical_output.loc[canonical_output["is_ntc"], "ps_score"] = 0.0
    excluded = set()
    if args.excluded_targets:
        excluded = set(json.loads(args.excluded_targets.read_text(encoding="utf-8")).get("excluded_targets", {}))
    canonical_output["ps_status"] = "modeled"
    canonical_output.loc[canonical_output["is_ntc"], "ps_status"] = "control"
    canonical_output.loc[canonical_output["target_label"].isin(excluded), "ps_status"] = "insufficient_target_cells"
    modeled_target = (~canonical_output["is_ntc"]) & (~canonical_output["target_label"].isin(excluded))
    if canonical_output.loc[modeled_target, "ps_score"].isna().any():
        raise ValueError("Target-cell shard coverage is incomplete")
    result = {"cells": len(canonical_output), "target_cells": int((~canonical_output["is_ntc"]).sum()),
              "modeled_target_cells": int(modeled_target.sum()), "unmodeled_target_cells": int((canonical_output["ps_status"] == "insufficient_target_cells").sum()),
              "excluded_targets": sorted(excluded), "ntc_cells": int(canonical_output["is_ntc"].sum()),
              "shards": len(files), "coverage_complete": True, "status": "PASS"}
    if args.canonical:
        canonical = pd.read_csv(args.canonical, sep="\t")
        canonical = canonical.loc[~canonical["is_ntc"], ["cell_id", "target_label", "ps_score"]]
        comparison = canonical.merge(merged, on=["cell_id", "target_label"], how="outer", suffixes=("_canonical", "_sharded"), indicator=True)
        if not (comparison["_merge"] == "both").all():
            raise ValueError("Shard coverage differs from reference result")
        difference = np.abs(comparison["ps_score_canonical"] - comparison["ps_score_sharded"])
        result["reference_comparison"] = {
            "max_abs_difference": float(difference.max()), "mean_abs_difference": float(difference.mean()),
            "pearson": float(np.corrcoef(comparison["ps_score_canonical"], comparison["ps_score_sharded"])[0, 1]),
            "exact_equal": bool(np.array_equal(comparison["ps_score_canonical"].to_numpy(), comparison["ps_score_sharded"].to_numpy()))
        }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    canonical_output.to_csv(args.output_dir / "ps_metadata.tsv", sep="\t", index=False)
    (args.output_dir / "merge_validation.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
