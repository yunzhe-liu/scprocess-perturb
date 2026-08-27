#!/usr/bin/env python3
"""
Standardize assignment CSV to unified schema.

Each assignment method produces a CSV with different columns and sorting
conventions.  This script normalises them to a single schema:

    cell_barcode, guide_id, umi_count, rank, score, score_type, method

Usage:
    python standardize_assignment.py \
        --input raw_assignments.csv \
        --output standardized.csv \
        --method pgmm_em
"""

import argparse, csv, time
from collections import defaultdict


# ── Method → sort + score mapping ───────────────────────────────────────

METHOD_SORT = {
    "pgmm_em": {
        "sort_keys": [("prob_gaussian", "desc"), ("UMI_counts", "desc")],
        "score_col": "prob_gaussian",
        "score_type": "prob_gaussian",
        "input_cols": ["cell", "gRNA", "UMI_counts", "prob_gaussian"],
    },
    "umi_threshold": {
        "sort_keys": [("UMI_counts", "desc")],
        "score_col": "UMI_counts",
        "score_type": "umi_count",
        "input_cols": ["cell", "gRNA", "UMI_counts"],
    },
    "fishash": {
        # Fed the raw fishash CSV directly (no top-K truncation in the pipeline).
        # All FDR-passing candidates are kept and ranked here by log_pval ASC
        # (more negative = more significant); per-cell selection happens later
        # in make_perturbation_obs.py.
        "sort_keys": [("log_pval", "asc")],
        "score_col": "log_pval",
        "score_type": "neg_log_pval",
        "input_cols": ["cell", "gRNA", "UMI_counts", "log_pval"],
    },
}


def main():
    parser = argparse.ArgumentParser(
        description="Standardize assignment CSV to unified schema")
    parser.add_argument("--input", required=True, help="Raw assignment CSV")
    parser.add_argument("--output", required=True, help="Output standardized CSV")
    parser.add_argument("--method", required=True,
                        choices=sorted(METHOD_SORT.keys()),
                        help="Assignment method name")
    args = parser.parse_args()

    mc = METHOD_SORT[args.method]
    t0 = time.time()

    # ── Load raw CSV ──
    print(f"Loading {args.input} …", end=' ', flush=True)
    pgmm = defaultdict(list)
    with open(args.input) as f:
        for row in csv.DictReader(f):
            cell = row.get("cell", "").strip()
            guide = row.get("gRNA", "").strip()
            if not cell or not guide:
                continue
            umi = int(float(row.get("UMI_counts", 0) or 0))
            # Build sort + score tuple
            sort_values = []
            for sk, _ in mc["sort_keys"]:
                val = row.get(sk, None)
                if val is not None:
                    sort_values.append(float(val))
                else:
                    sort_values.append(0.0)
            score_raw = row.get(mc["score_col"], None)
            score = float(score_raw) if score_raw else 0.0
            pgmm[cell].append((guide, umi, score, sort_values))

    # ── Sort per cell ──
    # Build comparison key respecting (sort_key, direction) pairs
    def make_sort_key(sort_vals):
        key_parts = []
        for (sk, direction), sv in zip(mc["sort_keys"], sort_vals):
            if direction == "desc":
                key_parts.append(-sv)
            else:
                key_parts.append(sv)
        return tuple(key_parts)

    for cell in pgmm:
        pgmm[cell].sort(key=lambda x: make_sort_key(x[3]))

    # ── Write unified schema ──
    n_cells = len(pgmm)
    n_rows = sum(len(v) for v in pgmm.values())
    print(f"{n_rows:,} rows, {n_cells:,} cells  [{time.time()-t0:.1f}s]")

    with open(args.output, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(["cell_barcode", "guide_id", "umi_count", "rank",
                     "score", "score_type", "method"])
        for cell, entries in pgmm.items():
            for rank_i, (guide, umi, score, _) in enumerate(entries):
                w.writerow([
                    cell, guide, umi, rank_i + 1,
                    f"{score:.6f}", mc["score_type"], args.method,
                ])

    print(f"  → {args.output}  [{time.time()-t0:.1f}s]")


if __name__ == "__main__":
    main()
