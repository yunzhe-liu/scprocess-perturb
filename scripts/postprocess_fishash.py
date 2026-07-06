#!/usr/bin/env python3
"""
Fishash post-processor — extracts top-K assignments per cell
from fishash's full (all-FDR-significant-pairs) CSV output.

Fishash assigns many guides per cell (gpC median 12–21 on typical data).
This script ranks them by log_pval (ASC, more negative = more significant)
and outputs a top-K CSV with standard columns.

Usage:
  python postprocess_fishash.py \
      --input fishash_raw.csv \
      --output assignments_topk.csv \
      --top-k 2
"""

import argparse, csv, os, time
from collections import defaultdict

def main():
    parser = argparse.ArgumentParser(description='Fishash top-K post-processor')
    parser.add_argument('--input', required=True, help='Fishash full assignments CSV')
    parser.add_argument('--output', required=True, help='Output top-K CSV path')
    parser.add_argument('--top-k', type=int, default=1, help='Guides per cell to keep (default: 1)')
    args = parser.parse_args()

    out_dir = os.path.dirname(args.output) or "."
    os.makedirs(out_dir, exist_ok=True)
    t0 = time.time()

    # Load fishash CSV — keep log_pval for ranking
    print(f"Loading {args.input} …", end=' ', flush=True)
    pgmm = defaultdict(list)
    with open(args.input) as f:
        for row in csv.DictReader(f):
            cell  = row.get('cell', '').strip()
            guide = row.get('gRNA', '').strip()
            if not cell or not guide:
                continue
            umi = int(float(row.get('UMI_counts', 0) or 0))
            lp  = float(row.get('log_pval', 0) or 0)
            pgmm[cell].append((guide, lp, umi))

    # Rank by log_pval ASC (more negative = more significant)
    for cell in pgmm:
        pgmm[cell].sort(key=lambda x: x[1])

    n_cells = len(pgmm)
    n_total = sum(len(v) for v in pgmm.values())
    print(f"{n_total:,} rows, {n_cells:,} cells  [{time.time()-t0:.1f}s]")

    # Output top-K — keep log_pval for downstream ranking
    rows = []
    for cell, entries in pgmm.items():
        for i in range(min(args.top_k, len(entries))):
            g, lp, umi = entries[i]
            rows.append([cell, g, str(umi), f"{lp:.6f}"])

    with open(args.output, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['cell', 'gRNA', 'UMI_counts', 'log_pval'])
        w.writerows(rows)

    n = len(rows)
    print(f"  top_{args.top_k}: {n:,} rows ({n/max(n_cells,1):.1f} guides/cell) → {args.output}")
    print(f"Done [{time.time()-t0:.1f}s]")


if __name__ == '__main__':
    main()
