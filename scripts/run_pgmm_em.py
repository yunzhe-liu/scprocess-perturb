#!/usr/bin/env python3
"""
PGMM assignment — Poisson-Gaussian Mixture Model via EM (MONITORED VERSION).

Identical EM algorithm to run_em_assignment.py, with additional monitoring:
  - Wall-clock time per stage
  - Peak RSS memory
  - System CPU count
  - Per-worker timing breakdown
  - Monitoring metadata saved as JSON alongside assignments.

Usage:
  python run_em_assignment_monitored.py \
      --input  /path/to/merged/ \
      --output /path/to/05_pgmm_em_assignment/{tool}/ \
      --tool-name cellranger \
      [--umi-threshold 3] [--prob-threshold 0.75] [--workers 16]
"""

import argparse, os, sys, gzip, time, json, platform, multiprocessing
from multiprocessing import Pool
import numpy as np
import pandas as pd
import scipy.io, scipy.sparse
from scipy.stats import poisson, norm
import warnings
warnings.filterwarnings("ignore")


# ══════════════════════════════════════════════════════════════════════════════
# Memory monitoring helper
# ══════════════════════════════════════════════════════════════════════════════
try:
    import resource
    _HAS_RESOURCE = True
except ImportError:
    _HAS_RESOURCE = False

try:
    import psutil
    _HAS_PSUTIL = True
except ImportError:
    _HAS_PSUTIL = False


def get_memory_rss_mb():
    """Return current RSS in MB using best available method."""
    if _HAS_PSUTIL:
        return psutil.Process().memory_info().rss / (1024 * 1024)
    elif _HAS_RESOURCE:
        return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    else:
        return -1.0


def get_peak_memory_mb():
    """Best-effort peak RSS so far, returns MB."""
    if _HAS_RESOURCE:
        # On Linux, ru_maxrss is in KB
        return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    elif _HAS_PSUTIL:
        return psutil.Process().memory_info().rss / (1024 * 1024)
    return -1.0


def get_system_info():
    """Collect system metadata for the monitoring report."""
    info = {
        "hostname": platform.node(),
        "platform": platform.platform(),
        "python_version": sys.version,
        "numpy_version": np.__version__,
        "scipy_version": scipy.__version__,
        "cpu_count_logical": os.cpu_count(),
    }
    try:
        # Physical cores via lscpu
        import subprocess
        out = subprocess.check_output(
            ["lscpu", "-p=cpu,core"], text=True, timeout=5
        )
        cores = set()
        for line in out.strip().split("\n"):
            if line.startswith("#"):
                continue
            parts = line.split(",")
            if len(parts) == 2:
                cores.add(parts[1])
        info["cpu_count_physical"] = len(cores)
    except Exception:
        info["cpu_count_physical"] = None
    if _HAS_PSUTIL:
        mem = psutil.virtual_memory()
        info["total_ram_gb"] = round(mem.total / (1024**3), 1)
        info["swap_total_gb"] = round(psutil.swap_memory().total / (1024**3), 1)
    return info


# ══════════════════════════════════════════════════════════════════════════════
# EM: 2-component Poisson-Gaussian mixture on log2(UMI)
# ══════════════════════════════════════════════════════════════════════════════
def em_poisson_gaussian(data, max_iter=200, tol=1e-4, n_init=5):
    """
    EM for 2-component Poisson-Gaussian mixture on log2(non-zero UMI).

    Constraint: Poisson mean < Gaussian mean (background < signal).

    Returns (w, lam, mu, sigma, converged) or None on failure.
    """
    N = len(data)
    if N < 2:
        return None

    best_ll = -np.inf
    best_params = None

    for init_i in range(n_init):
        rng = np.random.RandomState(42 + init_i)
        if init_i == 0:
            thresh = np.median(data)
        else:
            thresh = np.percentile(data, 40 + rng.uniform(-10, 10))
        low, high = data[data <= thresh], data[data > thresh]
        if len(low) < 2 or len(high) < 2:
            continue

        w = np.array([len(low) / N, len(high) / N])
        lam = low.mean()
        mu = high.mean()
        sigma = max(high.std(), 0.1)
        prev_ll = -np.inf
        converged = False

        for _ in range(max_iter):
            # E-step
            p_pois = poisson.pmf(np.round(data).astype(int), lam) + 1e-300
            p_gauss = norm.pdf(data, mu, sigma) + 1e-300
            resp = np.column_stack([w[0] * p_pois, w[1] * p_gauss])
            resp /= resp.sum(axis=1, keepdims=True)
            r0, r1 = resp[:, 0], resp[:, 1]
            s0, s1 = r0.sum(), r1.sum()
            if s0 < 1e-6 or s1 < 1e-6:
                break

            # M-step
            w = np.array([s0 / N, s1 / N])
            lam = np.average(data, weights=r0)
            mu = np.average(data, weights=r1)
            sigma = max(np.sqrt(np.average((data - mu) ** 2, weights=r1)), 0.05)

            # log-likelihood
            ll = np.sum(np.log(w[0] * p_pois + w[1] * p_gauss + 1e-300))
            if abs(ll - prev_ll) < tol:
                converged = True
                break
            prev_ll = ll

        if not converged or lam >= mu:
            continue
        if ll > best_ll:
            best_ll = ll
            best_params = (w, lam, mu, sigma)

    if best_params is None:
        return None
    return (*best_params, True)


# ══════════════════════════════════════════════════════════════════════════════
# Per-chunk worker
# ══════════════════════════════════════════════════════════════════════════════
def _worker(args):
    chunk_id, start_g, step, mtx_csc, barcodes, guide_names, \
        umi_th, prob_th, max_iter = args

    t0 = time.time()
    records = []
    stats = {"ok": 0, "fail": 0, "total_guides": 0}
    end_g = min(start_g + step, mtx_csc.shape[1])

    for g_idx in range(start_g, end_g):
        stats["total_guides"] += 1
        guide_name = guide_names[g_idx]
        col = mtx_csc[:, g_idx]
        if col.nnz < 2:
            stats["fail"] += 1
            continue

        data = col.data.astype(np.float64)
        log_data = np.log2(data)

        result = em_poisson_gaussian(log_data, max_iter=max_iter)
        if result is None:
            stats["fail"] += 1
            continue
        w, lam, mu, sigma, _ = result
        stats["ok"] += 1

        # P(Gaussian | UMI)
        p_pois = poisson.pmf(np.round(log_data).astype(int), lam) + 1e-300
        p_gauss = norm.pdf(log_data, mu, sigma) + 1e-300
        prob_gauss = w[1] * p_gauss / (w[0] * p_pois + w[1] * p_gauss)

        for cell_i, umi, prob in zip(col.indices, col.data, prob_gauss):
            records.append((barcodes[int(cell_i)], guide_name, int(umi), float(prob)))

    # ── Hard filters + flags ──
    df = pd.DataFrame(records, columns=["cell", "gRNA", "UMI_counts", "prob_gaussian"])
    df["pass_umi_filter"] = df["UMI_counts"] >= umi_th
    df["pass_prob_filter"] = df["prob_gaussian"] >= prob_th
    before_filter = len(df)
    df = df[df["pass_umi_filter"] & df["pass_prob_filter"]]
    df = df.sort_values(["cell", "prob_gaussian"], ascending=[True, False])

    stats["worker_wall_s"] = round(time.time() - t0, 2)
    stats["records_before_filter"] = before_filter
    stats["records_after_filter"] = len(df)

    return df, stats


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════
def main():
    parser = argparse.ArgumentParser(
        description="PGMM assignment — Poisson-Gaussian Mixture Model via EM (monitored)"
    )
    parser.add_argument("--input", required=True, help="Path to merged MEX directory")
    parser.add_argument("--output", required=True, help="Output CSV path")
    parser.add_argument("--umi-threshold", type=int, default=3)
    parser.add_argument("--prob-threshold", type=float, default=0.75)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--max-em-iter", type=int, default=200)
    args = parser.parse_args()

    # ── Monitoring bookkeeping ──
    monitor = {
        "method": "pgmm_em",
        "parameters": {
            "umi_threshold": args.umi_threshold,
            "prob_threshold": args.prob_threshold,
            "workers": args.workers,
            "max_em_iter": args.max_em_iter,
        },
        "stages": {},
        "system": get_system_info(),
    }
    T0 = time.time()
    mem_start_mb = get_memory_rss_mb()
    inp = args.input.rstrip("/")

    # Ensure output directory exists
    out_dir = os.path.dirname(args.output) or "."
    os.makedirs(out_dir, exist_ok=True)

    print(f"{'='*70}")
    print(f"PGMM EM Assignment — MONITORED RUN")
    print(f"  Tool:      pgmm_em")
    print(f"  Input:     {inp}")
    print(f"  Output:    {args.output}")
    print(f"  Workers:   {args.workers}")
    print(f"  UMI>={args.umi_threshold}  P(Gaussian)>={args.prob_threshold}")
    print(f"{'='*70}\n")

    # ── Stage 1: Load MEX ──
    stage_label = "load_mex"
    print(f"[{time.time()-T0:.0f}s] Stage 1/3: Loading MEX trio ...")
    t0 = time.time()
    mem_before = get_memory_rss_mb()

    with gzip.open(f"{inp}/merged_matrix.mtx.gz", "rt") as f:
        mtx = scipy.io.mmread(f).tocsc()
    with gzip.open(f"{inp}/merged_barcodes.tsv.gz", "rt") as f:
        barcodes = [l.strip() for l in f]
    with gzip.open(f"{inp}/merged_features.tsv.gz", "rt") as f:
        guide_names = [l.strip().split("\t")[0] for l in f]
    nc = mtx.shape[0]
    ng = mtx.shape[1]

    stage_elapsed = time.time() - t0
    mem_after = get_memory_rss_mb()
    monitor["stages"][stage_label] = {
        "wall_s": round(stage_elapsed, 2),
        "mem_before_mb": round(mem_before, 1),
        "mem_after_mb": round(mem_after, 1),
        "mem_delta_mb": round(mem_after - mem_before, 1),
        "input": inp,
        "ncells": nc,
        "nguides": ng,
        "nnz": int(mtx.nnz),
        "sparsity_pct": round((1 - mtx.nnz / (nc * ng)) * 100, 2),
        "mtx_size_mb": round(os.path.getsize(f"{inp}/merged_matrix.mtx.gz") / (1024**2), 1),
    }
    print(f"  Cells: {nc:,}  Guides: {ng:,}  NNZ: {mtx.nnz:,}  "
          f"Sparsity: {monitor['stages'][stage_label]['sparsity_pct']}%  [{stage_elapsed:.1f}s]")

    # ── Stage 2: Parallel EM ──
    stage_label = "em_fitting"
    print(f"\n[{time.time()-T0:.0f}s] Stage 2/3: Parallel EM fitting "
          f"({args.workers} workers, ~{ng} guides) ...")
    t0 = time.time()
    mem_before = get_memory_rss_mb()

    # Chunk guides across workers
    cs = ng // args.workers
    rm = ng % args.workers
    tasks = []
    s = 0
    for i in range(args.workers):
        st = cs + (1 if i < rm else 0)
        if st == 0:
            break
        tasks.append((i, s, st, mtx, barcodes, guide_names,
                      args.umi_threshold, args.prob_threshold, args.max_em_iter))
        s += st

    with Pool(processes=len(tasks)) as pool:
        results = pool.map(_worker, tasks)

    em_wall = time.time() - t0
    mem_after = get_memory_rss_mb()

    # Aggregate worker stats
    worker_wall_times = []
    worker_records_before = []
    worker_records_after = []
    tf = 0  # total fit OK
    tfl = 0  # total fit failed
    for _, stats in results:
        tf += stats["ok"]
        tfl += stats["fail"]
        worker_wall_times.append(stats.get("worker_wall_s", 0))
        worker_records_before.append(stats.get("records_before_filter", 0))
        worker_records_after.append(stats.get("records_after_filter", 0))

    monitor["stages"][stage_label] = {
        "wall_s": round(em_wall, 2),
        "wall_min": round(em_wall / 60, 2),
        "mem_before_mb": round(mem_before, 1),
        "mem_after_mb": round(mem_after, 1),
        "mem_delta_mb": round(mem_after - mem_before, 1),
        "guides_fit_ok": tf,
        "guides_failed": tfl,
        "guides_total": ng,
        "fit_success_rate_pct": round(tf / max(ng, 1) * 100, 2),
        "worker_wall_min": round(min(worker_wall_times), 1),
        "worker_wall_max": round(max(worker_wall_times), 1),
        "worker_wall_mean": round(np.mean(worker_wall_times), 1),
        "worker_wall_median": round(np.median(worker_wall_times), 1),
        "records_before_filter": sum(worker_records_before),
        "records_after_filter": sum(worker_records_after),
    }
    print(f"  Fit OK: {tf}  Failed: {tfl}  "
          f"Rate: {monitor['stages'][stage_label]['fit_success_rate_pct']}%")
    print(f"  EM wall time: {em_wall:.0f}s ({em_wall/60:.1f} min)")
    print(f"  Worker wall: min={min(worker_wall_times):.0f}s  "
          f"max={max(worker_wall_times):.0f}s  mean={np.mean(worker_wall_times):.0f}s")

    # ── Stage 3: Merge + Save ──
    stage_label = "merge_save"
    print(f"\n[{time.time()-T0:.0f}s] Stage 3/3: Merge & Save ...")
    t0 = time.time()
    mem_before = get_memory_rss_mb()

    all_dfs = [df for df, _ in results if len(df) > 0]
    merged = pd.concat(all_dfs, ignore_index=True) if all_dfs else pd.DataFrame(
        columns=["cell", "gRNA", "UMI_counts", "prob_gaussian",
                 "pass_umi_filter", "pass_prob_filter"])
    out_csv = args.output
    merged.to_csv(out_csv, index=False)

    stage_elapsed = time.time() - t0
    mem_after = get_memory_rss_mb()
    out_size_mb = os.path.getsize(out_csv) / (1024**2)
    monitor["stages"][stage_label] = {
        "wall_s": round(stage_elapsed, 2),
        "mem_before_mb": round(mem_before, 1),
        "mem_after_mb": round(mem_after, 1),
        "mem_delta_mb": round(mem_after - mem_before, 1),
        "output_csv": out_csv,
        "output_size_mb": round(out_size_mb, 1),
    }

    # ── Summary ──
    tt = time.time() - T0
    na = len(merged)
    nca = merged["cell"].nunique()
    ng2 = merged["gRNA"].nunique()
    gpc = merged.groupby("cell").size() if na > 0 else pd.Series(dtype=int)
    um = merged["UMI_counts"].median() if na > 0 else 0
    pm = merged["prob_gaussian"].median() if na > 0 else 0

    monitor["summary"] = {
        "total_wall_s": round(tt, 2),
        "total_wall_min": round(tt / 60, 2),
        "peak_rss_mb": round(get_peak_memory_mb(), 1),
        "mem_start_mb": round(mem_start_mb, 1),
        "mem_end_mb": round(get_memory_rss_mb(), 1),
        "total_assignments": int(na),
        "cells_assigned": int(nca),
        "cells_total": nc,
        "cell_recovery_pct": round(nca / max(nc, 1) * 100, 2),
        "guides_detected": int(ng2),
        "guides_per_cell_median": float(gpc.median()) if na > 0 else 0,
        "guides_per_cell_mean": round(float(gpc.mean()), 2) if na > 0 else 0,
        "guides_per_cell_max": int(gpc.max()) if na > 0 else 0,
        "cells_1_guide": int((gpc == 1).sum()) if na > 0 else 0,
        "cells_2_guides": int((gpc == 2).sum()) if na > 0 else 0,
        "cells_ge3_guides": int((gpc >= 3).sum()) if na > 0 else 0,
        "umi_median": float(um),
        "prob_gaussian_median": float(pm),
    }

    # ── Save monitoring JSON ──
    mon_json = os.path.join(out_dir, "monitoring.json")
    with open(mon_json, "w") as f:
        json.dump(monitor, f, indent=2, default=str)
    print(f"  Monitoring JSON saved: {mon_json}")

    print(f"\n{'='*70}")
    print(f"PGMM EM Assignment Complete")
    print(f"{'='*70}")
    print(f"  Input:                 {inp}")
    print(f"  Params:                UMI>={args.umi_threshold}  "
          f"P(Gaussian)>={args.prob_threshold}")
    print(f"  Guides fit OK:         {tf:>12,}")
    print(f"  Guides failed:         {tfl:>12,}")
    print(f"  Assignments:           {na:>12,}")
    print(f"  Cells assigned:        {nca:>12,}  ({nca/max(nc,1)*100:.1f}%)")
    print(f"  Guides detected:       {ng2:>12,}")
    if na > 0:
        print(f"  Guides/cell:           median={gpc.median():.0f}  "
              f"mean={gpc.mean():.2f}  max={gpc.max()}")
        print(f"  1 guide:               {(gpc==1).sum():>12,}")
        print(f"  2 guides:              {(gpc==2).sum():>12,}")
        print(f"  >=3 guides:            {(gpc>=3).sum():>12,}  "
              f"({(gpc>=3).sum()/max(nca,1)*100:.1f}%)")
    print(f"  UMI median:            {um:>12.0f}")
    print(f"  P(Gaussian) median:    {pm:>12.3f}")
    print(f"\n  Total wall time: {tt:.0f}s ({tt/60:.1f} min)")
    print(f"  Peak RSS:        {monitor['summary']['peak_rss_mb']:.0f} MB")
    print(f"  Assignment CSV:  {out_csv}  ({out_size_mb:.1f} MB)")
    print(f"  Monitoring JSON: {mon_json}")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
