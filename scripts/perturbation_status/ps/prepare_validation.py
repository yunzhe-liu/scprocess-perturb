#!/usr/bin/env python3
"""Create balanced Papalexi target shards for pooled-equivalence validation."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--barcode", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--shards", type=int, default=4)
    parser.add_argument("--min-target-cells", type=int, default=3)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    barcode = pd.read_csv(args.barcode, sep="\t")
    counts = barcode.loc[barcode["gene"] != "non-targeting", "gene"].value_counts()
    excluded = {str(target): int(value) for target, value in counts[counts < args.min_target_cells].items()}
    counts = counts[counts >= args.min_target_cells]
    targets = [(target, int(counts[target])) for target in sorted(counts.index)]
    targets.sort(key=lambda item: (-item[1], item[0]))
    loads = [0] * args.shards
    assignments: list[list[tuple[str, int]]] = [[] for _ in range(args.shards)]
    for target, cells in targets:
        shard = min(range(args.shards), key=lambda index: (loads[index], index))
        assignments[shard].append((target, cells))
        loads[shard] += cells

    rows = []
    for shard, values in enumerate(assignments, start=1):
        path = args.output_dir / f"shard_{shard:03d}_targets.txt"
        path.write_text("\n".join(target for target, _ in values) + "\n", encoding="utf-8")
        rows.extend({"target_label": target, "cells": cells, "shard_id": shard} for target, cells in values)
    pd.DataFrame(rows).sort_values(["shard_id", "target_label"]).to_csv(
        args.output_dir / "target_manifest.tsv", sep="\t", index=False
    )
    manifest = {"shards": args.shards, "targets": len(targets), "target_cells": int(sum(loads)),
                "min_target_cells": args.min_target_cells, "excluded_targets": excluded,
                "excluded_target_cells": int(sum(excluded.values())), "shard_target_cell_loads": loads}
    (args.output_dir / "shard_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
