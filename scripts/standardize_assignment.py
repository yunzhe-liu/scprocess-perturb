#!/usr/bin/env python3
"""
Standardize assignment CSV to unified schema.

Each assignment method produces a CSV with different columns and sorting
conventions.  This script normalises them to a single schema:

    cell_barcode, guide_id, umi_count, rank, score, score_type, method,
    native_score, native_score_type

Usage:
    python standardize_assignment.py \
        --input raw_assignments.csv \
        --output standardized.csv \
        --method pgmm_em
"""

import argparse, csv, math, os, tempfile, time
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
        "score_type": "negative_log_pval",
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
    seen_pairs = set()
    with open(args.input) as f:
        reader = csv.DictReader(f)
        missing = set(mc["input_cols"]).difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing required columns: {sorted(missing)}")
        for line_number, row in enumerate(reader, start=2):
            cell = row.get("cell", "").strip()
            guide = row.get("gRNA", "").strip()
            if not cell or not guide:
                raise ValueError(f"Row {line_number}: cell and gRNA must be non-empty")
            pair = (cell, guide)
            if pair in seen_pairs:
                raise ValueError(f"Row {line_number}: duplicate cell-guide pair {pair}")
            seen_pairs.add(pair)
            umi_value = float(row["UMI_counts"])
            if not math.isfinite(umi_value) or umi_value < 0 or not umi_value.is_integer():
                raise ValueError(f"Row {line_number}: UMI_counts must be a non-negative integer")
            umi = int(umi_value)
            # Build sort + score tuple
            sort_values = []
            for sk, _ in mc["sort_keys"]:
                val = row.get(sk, None)
                value = float(val)
                if not math.isfinite(value):
                    raise ValueError(f"Row {line_number}: {sk} must be finite")
                sort_values.append(value)
            score_raw = row.get(mc["score_col"], None)
            if score_raw in (None, ""):
                raise ValueError(f"Row {line_number}: missing {mc['score_col']}")
            native_score = float(score_raw)
            if not math.isfinite(native_score):
                raise ValueError(f"Row {line_number}: score must be finite")
            score = -native_score if args.method == "fishash" else native_score
            if score < 0:
                raise ValueError(f"Row {line_number}: standardized score must be non-negative")
            pgmm[cell].append((guide, umi, score, native_score, sort_values))

    if not pgmm:
        raise ValueError("Assignment input contains no rows")

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
        pgmm[cell].sort(key=lambda x: make_sort_key(x[4]))

    # ── Write unified schema ──
    n_cells = len(pgmm)
    n_rows = sum(len(v) for v in pgmm.values())
    print(f"{n_rows:,} rows, {n_cells:,} cells  [{time.time()-t0:.1f}s]")

    output_dir = os.path.dirname(os.path.abspath(args.output))
    os.makedirs(output_dir, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=".assignments-", suffix=".csv", dir=output_dir)
    os.close(fd)
    with open(temporary, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(["cell_barcode", "guide_id", "umi_count", "rank",
                     "score", "score_type", "method", "native_score",
                     "native_score_type"])
        for cell, entries in pgmm.items():
            for rank_i, (guide, umi, score, native_score, _) in enumerate(entries):
                w.writerow([
                    cell, guide, umi, rank_i + 1,
                    f"{score:.6f}", mc["score_type"], args.method,
                    f"{native_score:.6f}", mc["score_col"],
                ])
    os.replace(temporary, args.output)

    print(f"  → {args.output}  [{time.time()-t0:.1f}s]")


if __name__ == "__main__":
    main()
