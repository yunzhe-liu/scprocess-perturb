#!/usr/bin/env python3
"""
Standalone UMI Threshold Assignment.

Reads a merged MEX trio, assigns every (cell, guide) pair with UMI ≥ threshold.
Output matches PGMM EM CSV format (cell, gRNA, UMI_counts).

Usage:
  python run_umi_threshold.py \
      --input  /path/to/guide_matrix/ \
      --output /path/to/assignment/umi_threshold/_raw_assignments.csv \
      --threshold 3
"""

import argparse, os, sys, gzip, time, json, platform
import numpy as np
import scipy.io, scipy.sparse

try:
    import psutil
    _HAS_PSUTIL = True
except ImportError:
    _HAS_PSUTIL = False

try:
    import resource
    _HAS_RESOURCE = True
except ImportError:
    _HAS_RESOURCE = False

def get_memory_rss_mb():
    if _HAS_PSUTIL:
        return psutil.Process().memory_info().rss / (1024 * 1024)
    elif _HAS_RESOURCE:
        return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    return -1.0

def get_peak_memory_mb():
    if _HAS_PSUTIL:
        return psutil.Process().memory_info().rss / (1024 * 1024)
    elif _HAS_RESOURCE:
        return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    return -1.0

def get_system_info():
    info = {
        "hostname": platform.node(),
        "platform": platform.platform(),
        "python_version": sys.version,
        "numpy_version": np.__version__,
        "scipy_version": scipy.__version__,
        "cpu_count_logical": os.cpu_count(),
    }
    try:
        import subprocess
        out = subprocess.check_output(
            ["lscpu", "-p=cpu,core"], text=True, timeout=5)
        cores = set()
        for line in out.strip().split("\n"):
            if line.startswith("#"): continue
            parts = line.split(",")
            if len(parts) == 2: cores.add(parts[1])
        info["cpu_count_physical"] = len(cores)
    except Exception:
        info["cpu_count_physical"] = None
    if _HAS_PSUTIL:
        info["total_ram_gb"] = round(psutil.virtual_memory().total / (1024**3), 1)
    return info


def main():
    parser = argparse.ArgumentParser(
        description="Standalone UMI Threshold Assignment"
    )
    parser.add_argument("--input", required=True, help="Path to merged MEX directory")
    parser.add_argument("--output", required=True, help="Output CSV path")
    parser.add_argument("--threshold", type=int, required=True, help="UMI threshold (e.g. 3, 5, 10)")
    args = parser.parse_args()

    T0 = time.time()
    inp = args.input.rstrip("/")
    os.makedirs(os.path.dirname(args.output) if os.path.dirname(args.output) else ".", exist_ok=True)

    print(f"{'='*70}")
    print(f"UMI Threshold Assignment (standalone)")
    print(f"  Input:     {inp}")
    print(f"  Threshold: UMI ≥ {args.threshold}")
    print(f"  Output:    {args.output}")
    print(f"{'='*70}\n")

    # ── Stage 1: Load MEX ──
    print(f"[{time.time()-T0:.0f}s] Loading MEX …")
    t0 = time.time()
    mem_before = get_memory_rss_mb()

    with gzip.open(f"{inp}/merged_matrix.mtx.gz", "rt") as f:
        mtx = scipy.io.mmread(f).tocsc()
    with gzip.open(f"{inp}/merged_barcodes.tsv.gz", "rt") as f:
        barcodes = [l.strip() for l in f]
    with gzip.open(f"{inp}/merged_features.tsv.gz", "rt") as f:
        guide_names = [l.strip().split("\t")[0] for l in f]

    nc = mtx.shape[0]; ng = mtx.shape[1]
    mem_after = get_memory_rss_mb()
    load_t = time.time() - t0
    print(f"  Cells: {nc:,}  Guides: {ng:,}  NNZ: {mtx.nnz:,}  [{load_t:.1f}s]")

    # ── Stage 2: Apply threshold ──
    print(f"\n[{time.time()-T0:.0f}s] Applying UMI ≥ {args.threshold} filter …")
    t0 = time.time()
    mem_before2 = get_memory_rss_mb()

    # Convert to CSR for efficient row slicing (cell-level iteration)
    mtx_csr = mtx.tocsr()
    del mtx  # free CSC

    rows = []
    for cell_idx in range(nc):
        row = mtx_csr.getrow(cell_idx)
        if row.nnz == 0: continue
        cell_bc = barcodes[cell_idx]
        for j in range(row.nnz):
            umi = int(row.data[j])
            if umi >= args.threshold:
                g_idx = row.indices[j]
                rows.append((cell_bc, guide_names[g_idx], umi))

    filter_t = time.time() - t0
    mem_after2 = get_memory_rss_mb()
    print(f"  Filtered: {len(rows):,} assignments  [{filter_t:.1f}s]")

    # ── Stage 3: Write CSV ──
    print(f"\n[{time.time()-T0:.0f}s] Writing CSV …")
    t0 = time.time()

    # Sort by cell, then UMI descending
    import pandas as pd
    df = pd.DataFrame(rows, columns=["cell", "gRNA", "UMI_counts"])
    df = df.sort_values(["cell", "UMI_counts"], ascending=[True, False])
    df.to_csv(args.output, index=False)

    write_t = time.time() - t0
    out_size_mb = os.path.getsize(args.output) / (1024**2)
    print(f"  Wrote {len(df):,} rows  ({out_size_mb:.1f} MB)  [{write_t:.1f}s]")

    # ── Monitoring ──
    tt = time.time() - T0
    monitor = {
        "method": "umi_threshold",
        "parameters": {"threshold": args.threshold},
        "system": get_system_info(),
        "stages": {
            "load_mex": {
                "wall_s": round(load_t, 2),
                "mem_before_mb": round(mem_before, 1),
                "mem_after_mb": round(mem_after, 1),
                "ncells": nc, "nguides": ng, "nnz": int(mtx_csr.nnz),
            },
            "filter": {
                "wall_s": round(filter_t, 2),
                "mem_before_mb": round(mem_before2, 1),
                "mem_after_mb": round(mem_after2, 1),
                "assignments": len(rows),
            },
            "write": {
                "wall_s": round(write_t, 2),
                "output_size_mb": round(out_size_mb, 1),
            },
        },
        "summary": {
            "total_wall_s": round(tt, 2),
            "total_wall_min": round(tt / 60, 2),
            "peak_rss_mb": round(get_peak_memory_mb(), 1),
            "total_assignments": len(df),
            "cells_assigned": int(df["cell"].nunique()),
            "guides_detected": int(df["gRNA"].nunique()),
        },
    }

    # Match pgmm_em / fishash: monitoring.json in the output directory.
    mon_json = os.path.join(os.path.dirname(args.output) or ".", "monitoring.json")
    with open(mon_json, "w") as f:
        json.dump(monitor, f, indent=2, default=str)

    print(f"\n{'='*70}")
    print(f"UMI Threshold (standalone) DONE — t={args.threshold}")
    print(f"{'='*70}")
    print(f"  Assignments:        {len(df):>12,}")
    print(f"  Cells assigned:     {df['cell'].nunique():>12,}")
    print(f"  Guides detected:    {df['gRNA'].nunique():>12,}")
    print(f"  Total wall time:    {tt:.0f}s ({tt/60:.1f} min)")
    print(f"  Output:             {args.output}  ({out_size_mb:.1f} MB)")
    print(f"{'='*70}\n")


if __name__ == "__main__":
    main()
