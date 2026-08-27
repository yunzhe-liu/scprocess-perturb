#!/usr/bin/env python3
"""
Task 2: Global Barcode-Rank QC & Whitelist Generation

Performs barcode-rank (knee-plot) diagnostics on the merged 643K-cell
mRNA count matrix, determines thresholds via two strategies:
  1. Manual static:   UMI > 1000  &  Genes > 500
  2. Algorithm-based:  Knee point (max curvature) & Inflection point (max |y''|)
     via smoothing spline + 1st/2nd derivative analysis on log-log rank curve.

Outputs:
  - High-resolution QC plots with dual-threshold annotations
  - Three whitelist versions (manual / knee / inflection)
  - Summary statistics for comparative evaluation

Usage:
    conda activate scp_analysis
    python barcode_rank_qc.py
"""

import gc
import os
import sys
import time
import logging
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.interpolate import CubicSpline
from statsmodels.nonparametric.smoothers_lowess import lowess
import anndata as ad

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import FancyBboxPatch

# ─── Logging ────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)

# ─── Constants ──────────────────────────────────────────────────────────────
INPUT_H5AD = Path("/data/yunzliu/results/guide_extraction/merged/merged_48lanes_raw.h5ad")
OUTPUT_BASE = Path("/data/yunzliu/results/guide_extraction/whitelist")
MANUAL_UMI_THRESHOLD = 1000
MANUAL_GENE_THRESHOLD = 500
CHUNK_SIZE = 50000  # rows per chunk for metric computation
DEFAULT_EXCLUDED_FRACTION = 0.001  # exclude top/bottom 0.1% for spline fitting


def parse_args():
    p = argparse.ArgumentParser(description="Barcode-Rank QC and whitelist generation")
    p.add_argument("--input", type=Path, default=INPUT_H5AD)
    p.add_argument("--output-base", type=Path, default=OUTPUT_BASE)
    p.add_argument("--manual-umi", type=int, default=MANUAL_UMI_THRESHOLD)
    p.add_argument("--manual-gene", type=int, default=MANUAL_GENE_THRESHOLD)
    p.add_argument("--chunk-size", type=int, default=CHUNK_SIZE)
    p.add_argument("--no-plots", action="store_true", help="Skip plot generation")
    return p.parse_args()


# ═══════════════════════════════════════════════════════════════════════════════
# 1. Per-cell metric computation (chunked, memory-efficient)
# ═══════════════════════════════════════════════════════════════════════════════

def compute_cell_metrics_chunked(adata: ad.AnnData, chunk_size: int) -> pd.DataFrame:
    """
    Compute total_UMI and n_genes per barcode via chunked row-wise access.
    Works with both in-memory and backed AnnData.
    """
    n_cells = adata.n_obs
    total_umi = np.zeros(n_cells, dtype=np.float64)
    n_genes = np.zeros(n_cells, dtype=np.int32)

    log.info(f"Computing per-cell metrics for {n_cells:,} cells (chunk={chunk_size:,}) …")
    t0 = time.perf_counter()

    for start in range(0, n_cells, chunk_size):
        end = min(start + chunk_size, n_cells)
        # Extract chunk — backed CSRDataset supports slicing
        X_chunk = adata.X[start:end]
        if sp.issparse(X_chunk):
            X_chunk = X_chunk.tocsr()  # ensure CSR for fast row ops
            total_umi[start:end] = np.asarray(X_chunk.sum(axis=1)).ravel()
            n_genes[start:end] = np.diff(X_chunk.indptr)  # nnz per row for CSR
        else:
            total_umi[start:end] = X_chunk.sum(axis=1)
            n_genes[start:end] = (X_chunk > 0).sum(axis=1)

        if (start // chunk_size) % 10 == 0:
            progress = end / n_cells * 100
            log.info(f"  … {progress:.0f}% ({end:,}/{n_cells:,})")

    elapsed = time.perf_counter() - t0
    log.info(f"Metrics computed in {elapsed:.1f}s")

    df = pd.DataFrame(
        {"total_umi": total_umi, "n_genes": n_genes},
        index=adata.obs_names,
    )
    df["lane"] = adata.obs["lane"].values
    return df


# ═══════════════════════════════════════════════════════════════════════════════
# 2. BarcodeRanks algorithm
# ═══════════════════════════════════════════════════════════════════════════════

def barcode_ranks(
    total_umi: np.ndarray,
    exclude_frac: float = DEFAULT_EXCLUDED_FRACTION,
) -> dict:
    """
    Two-stage barcodeRanks: LOWESS denoising → Natural cubic spline.

    Stage 1 — LOWESS: non-parametric local regression (frac~2%), robust to
    outliers.  Irons out micro-oscillations while preserving the cliff edge.
    Stage 2 — Natural cubic spline (C², y''=0 at boundaries).  Fits the
    already-denoised trend for analytic derivative extraction.  Runge
    oscillations impossible since input data is free of high-freq noise.

    Parameters
    ----------
    total_umi : 1-D array of per-barcode total UMI counts
    exclude_frac : fraction of extreme points to exclude from fitting

    Returns
    -------
    dict with keys:
      ranks, sorted_umi, log_rank, log_umi,
      knee_rank, knee_umi, inflection_rank, inflection_umi,
      spline_y, spline_y1, spline_y2, curvature
    """
    n = len(total_umi)

    # Sort descending
    sorted_order = np.argsort(-total_umi)
    sorted_umi = total_umi[sorted_order]
    ranks = np.arange(1, n + 1)

    # Filter zero-UMI barcodes (they cause log10 issues)
    nonzero_mask = sorted_umi > 0
    ranks_nz = ranks[nonzero_mask]
    sorted_umi_nz = sorted_umi[nonzero_mask]

    log_rank = np.log10(ranks_nz)
    log_umi = np.log10(sorted_umi_nz)

    log.info(f"Nonzero UMI barcodes: {np.sum(nonzero_mask):,} / {n:,}")

    # Exclude extremes for robust spline fitting
    n_exclude = max(1, int(len(log_rank) * exclude_frac))
    fit_start = n_exclude
    fit_end = len(log_rank) - n_exclude

    x_fit = log_rank[fit_start:fit_end]
    y_fit = log_umi[fit_start:fit_end]

    # ══════════════════════════════════════════════════════════════════════
    # Stage 1: LOWESS — non-parametric local denoising (no global df assumption)
    # ══════════════════════════════════════════════════════════════════════
    N_LOWESS = 5000
    if len(x_fit) > N_LOWESS:
        idx_sub = np.linspace(0, len(x_fit) - 1, N_LOWESS, dtype=int)
        x_lowess = x_fit[idx_sub]
        y_lowess = y_fit[idx_sub]
    else:
        x_lowess = x_fit
        y_lowess = y_fit

    t0 = time.perf_counter()
    lowess_frac = 0.02
    lowess_result = lowess(y_lowess, x_lowess, frac=lowess_frac, return_sorted=True)
    x_smooth = lowess_result[:, 0]
    y_smooth = lowess_result[:, 1]
    log.info(f"LOWESS complete in {time.perf_counter()-t0:.2f}s "
             f"({len(x_smooth)} pts, frac={lowess_frac})")

    # ══════════════════════════════════════════════════════════════════════
    # Stage 2: Natural cubic spline on the denoised trend
    # ══════════════════════════════════════════════════════════════════════
    t0 = time.perf_counter()
    spline = CubicSpline(x_smooth, y_smooth, bc_type='natural', extrapolate=True)
    log.info(f"Natural cubic spline fitted in {time.perf_counter()-t0:.2f}s")

    x_eval = x_fit
    spline_y = spline(x_eval)
    spline_y1 = spline(x_eval, 1)
    spline_y2 = spline(x_eval, 2)

    curvature = np.abs(spline_y2) / (1 + spline_y1 ** 2) ** 1.5

    # ── Knee: max perpendicular distance from endpoint-connecting line ────
    x_first, y_first = x_eval[0], spline_y[0]
    x_last, y_last = x_eval[-1], spline_y[-1]
    a = y_last - y_first
    b = -(x_last - x_first)
    c = x_last * y_first - y_last * x_first
    denom = np.sqrt(a**2 + b**2)
    distances = np.abs(a * x_eval + b * spline_y + c) / denom

    n_eval = len(x_eval)
    search_lo = n_eval // 10
    search_hi = n_eval * 9 // 10
    knee_idx_local = search_lo + np.argmax(distances[search_lo:search_hi])
    knee_idx_global = fit_start + knee_idx_local
    knee_rank = ranks_nz[knee_idx_global]
    knee_umi = sorted_umi_nz[knee_idx_global]

    # ── Inflection: minimum 1st derivative (steepest slope) ───────────────
    inflection_idx_local = search_lo + np.argmin(spline_y1[search_lo:search_hi])
    inflection_idx_global = fit_start + inflection_idx_local
    inflection_rank = ranks_nz[inflection_idx_global]
    inflection_umi = sorted_umi_nz[inflection_idx_global]

    log.info(f"Knee point:       rank={knee_rank:,}  UMI={knee_umi:.0f}")
    log.info(f"Inflection point: rank={inflection_rank:,}  UMI={inflection_umi:.0f}")

    return {
        "ranks": ranks_nz,
        "sorted_umi": sorted_umi_nz,
        "log_rank": log_rank,
        "log_umi": log_umi,
        "knee_rank": knee_rank,
        "knee_umi": knee_umi,
        "inflection_rank": inflection_rank,
        "inflection_umi": inflection_umi,
        "fit_start": fit_start,
        "fit_end": fit_end,
        "x_eval": x_eval,
        "spline_y": spline_y,
        "spline_y1": spline_y1,
        "spline_y2": spline_y2,
        "curvature": curvature,
        "nonzero_mask": nonzero_mask,
    }


# ═══════════════════════════════════════════════════════════════════════════════
# 3. Whitelist generation
# ═══════════════════════════════════════════════════════════════════════════════

def generate_whitelists(
    metrics_df: pd.DataFrame,
    br_result: dict,
    manual_umi: int,
    manual_gene: int,
    output_base: Path,
) -> dict:
    """Generate three whitelist versions and save to disk."""

    total_cells = len(metrics_df)

    # ── v1: Manual threshold ──────────────────────────────────────────────
    mask_manual = (metrics_df["total_umi"] > manual_umi) & (metrics_df["n_genes"] > manual_gene)
    wl_manual = metrics_df.index[mask_manual]

    # ── v2: Knee threshold ───────────────────────────────────────────────
    knee_umi = br_result["knee_umi"]
    # Derive gene threshold from median n_genes near knee UMI
    knee_neighborhood = metrics_df[
        (metrics_df["total_umi"] >= knee_umi * 0.9) & (metrics_df["total_umi"] <= knee_umi * 1.1)
    ]
    knee_genes = int(knee_neighborhood["n_genes"].median()) if len(knee_neighborhood) > 0 else 0
    mask_knee = (metrics_df["total_umi"] > knee_umi) & (metrics_df["n_genes"] > knee_genes)
    wl_knee = metrics_df.index[mask_knee]

    # ── v3: Inflection threshold ──────────────────────────────────────────
    infl_umi = br_result["inflection_umi"]
    infl_neighborhood = metrics_df[
        (metrics_df["total_umi"] >= infl_umi * 0.9) & (metrics_df["total_umi"] <= infl_umi * 1.1)
    ]
    infl_genes = int(infl_neighborhood["n_genes"].median()) if len(infl_neighborhood) > 0 else 0
    mask_infl = (metrics_df["total_umi"] > infl_umi) & (metrics_df["n_genes"] > infl_genes)
    wl_infl = metrics_df.index[mask_infl]

    # ── Save ──────────────────────────────────────────────────────────────
    versions = {
        "v1_manual": (wl_manual, mask_manual, manual_umi, manual_gene,
                       output_base / "v1_manual"),
        "v2_knee": (wl_knee, mask_knee, knee_umi, knee_genes,
                     output_base / "v2_barcodeRanks_knee"),
        "v2_inflection": (wl_infl, mask_infl, infl_umi, infl_genes,
                           output_base / "v2_barcodeRanks_inflection"),
    }

    stats_list = []
    for version_name, (wl, mask, umi_th, gene_th, out_dir) in versions.items():
        out_dir.mkdir(parents=True, exist_ok=True)

        # Save whitelist (barcodes only, one per line, gzipped)
        wl_path = out_dir / f"whitelist_{version_name}.tsv.gz"
        pd.Series(wl).to_csv(wl_path, index=False, header=False, compression="gzip")
        fsize = wl_path.stat().st_size / 1024**2

        n_pass = mask.sum()
        stats_list.append({
            "version": version_name,
            "umi_threshold": umi_th,
            "gene_threshold": gene_th,
            "n_pass": n_pass,
            "pct_pass": round(n_pass / total_cells * 100, 2),
            "n_reject": total_cells - n_pass,
            "pct_reject": round((total_cells - n_pass) / total_cells * 100, 2),
            "whitelist_file": str(wl_path),
            "file_size_mb": round(fsize, 2),
        })

        log.info(f"  {version_name}: {n_pass:,} / {total_cells:,} "
                 f"({n_pass/total_cells*100:.1f}%) pass — "
                 f"UMI>{umi_th:.0f}, Genes>{gene_th} — saved {fsize:.1f} MB")

    stats_df = pd.DataFrame(stats_list)
    stats_path = output_base / "comparative" / "whitelist_comparison.csv"
    stats_path.parent.mkdir(parents=True, exist_ok=True)
    stats_df.to_csv(stats_path, index=False)

    return {vn: (wl, mask, umi_th, gene_th) for vn, (wl, mask, umi_th, gene_th, _) in versions.items()}


# ═══════════════════════════════════════════════════════════════════════════════
# 4. Visualization
# ═══════════════════════════════════════════════════════════════════════════════

def plot_barcode_rank(
    metrics_df: pd.DataFrame,
    br_result: dict,
    whitelists: dict,
    output_dir: Path,
):
    """Generate publication-quality barcode-rank plot with dual thresholds."""
    output_dir.mkdir(parents=True, exist_ok=True)

    log_rank = br_result["log_rank"]
    log_umi = br_result["log_umi"]
    ranks = br_result["ranks"]
    sorted_umi = br_result["sorted_umi"]
    x_eval = br_result["x_eval"]
    spline_y = br_result["spline_y"]
    spline_y1 = br_result["spline_y1"]
    spline_y2 = br_result["spline_y2"]
    curvature = br_result["curvature"]
    fit_start = br_result["fit_start"]
    knee_umi = br_result["knee_umi"]
    knee_rank = br_result["knee_rank"]
    infl_umi = br_result["inflection_umi"]
    infl_rank = br_result["inflection_rank"]
    manual_umi = whitelists["v1_manual"][2]

    # Determine which barcodes pass each threshold
    total_umi_all = metrics_df["total_umi"].values
    # Re-sort to match rank order
    sorted_order = np.argsort(-total_umi_all)
    pass_manual = np.isin(
        np.arange(len(metrics_df))[sorted_order],
        np.where(whitelists["v1_manual"][1])[0][np.argsort(-total_umi_all[np.where(whitelists["v1_manual"][1])[0]])]
    )
    # Simpler: compute mask in sorted order
    sorted_metrics = metrics_df.iloc[sorted_order]
    sorted_pass_manual = (sorted_metrics["total_umi"] > manual_umi) & (sorted_metrics["n_genes"] > 500)
    sorted_pass_knee = (sorted_metrics["total_umi"] > knee_umi) & (sorted_metrics["n_genes"] > whitelists["v2_knee"][3])
    sorted_pass_infl = (sorted_metrics["total_umi"] > infl_umi) & (sorted_metrics["n_genes"] > whitelists["v2_inflection"][3])

    # ── Figure 1: Main Knee Plot ─────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(14, 10))

    # Subsample for scatter (643K points is too dense)
    n_scatter = min(100000, len(log_rank))
    idx_scatter = np.linspace(0, len(log_rank) - 1, n_scatter, dtype=int)

    # Color by threshold category
    pass_manual_sc = sorted_pass_manual.values[idx_scatter] if hasattr(sorted_pass_manual, 'values') else sorted_pass_manual[idx_scatter]
    pass_knee_sc = sorted_pass_knee.values[idx_scatter] if hasattr(sorted_pass_knee, 'values') else sorted_pass_knee[idx_scatter]

    cat = np.full(len(idx_scatter), 0, dtype=int)
    # Category: 3=both, 2=manual only, 1=knee only, 0=neither
    cat[pass_manual_sc & pass_knee_sc] = 3
    cat[pass_manual_sc & ~pass_knee_sc] = 2
    cat[~pass_manual_sc & pass_knee_sc] = 1

    colors = {0: "#d0d0d0", 1: "#66c2a5", 2: "#fc8d62", 3: "#2b83ba"}
    labels = {0: "Both reject", 1: "Knee only", 2: "Manual only", 3: "Both pass"}
    for c in [0, 1, 2, 3]:
        mask_c = cat == c
        if mask_c.sum() > 0:
            ax.scatter(
                log_rank[idx_scatter][mask_c],
                log_umi[idx_scatter][mask_c],
                c=colors[c], s=0.5, alpha=0.4, rasterized=True, label=labels[c],
            )

    # Spline curve
    ax.plot(x_eval, spline_y, color="black", linewidth=2.0, label="Smoothing spline", zorder=5)

    # Knee point
    knee_log_rank = np.log10(knee_rank)
    knee_log_umi = np.log10(knee_umi)
    ax.axvline(x=knee_log_rank, color="#d7191c", linestyle="--", linewidth=1.5, alpha=0.8)
    ax.axhline(y=knee_log_umi, color="#d7191c", linestyle="--", linewidth=1.5, alpha=0.8)
    ax.scatter([knee_log_rank], [knee_log_umi], color="#d7191c", s=120, zorder=10,
               edgecolors="white", linewidths=1.5)
    ax.annotate(
        f"Knee (max distance)\nUMI={knee_umi:,.0f}\nRank={knee_rank:,}",
        xy=(knee_log_rank, knee_log_umi),
        xytext=(knee_log_rank + 0.3, knee_log_umi + 0.4),
        fontsize=9, color="#d7191c", fontweight="bold",
        arrowprops=dict(arrowstyle="->", color="#d7191c", lw=1.2),
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.85, edgecolor="#d7191c"),
    )

    # Inflection point
    infl_log_rank = np.log10(infl_rank)
    infl_log_umi = np.log10(infl_umi)
    ax.axvline(x=infl_log_rank, color="#fdae61", linestyle="-.", linewidth=1.5, alpha=0.8)
    ax.axhline(y=infl_log_umi, color="#fdae61", linestyle="-.", linewidth=1.5, alpha=0.8)
    ax.scatter([infl_log_rank], [infl_log_umi], color="#fdae61", s=120, zorder=10,
               edgecolors="white", linewidths=1.5, marker="s")
    ax.annotate(
        f"Inflection (min y')\nUMI={infl_umi:,.0f}\nRank={infl_rank:,}",
        xy=(infl_log_rank, infl_log_umi),
        xytext=(infl_log_rank - 1.8, infl_log_umi - 1.0),
        fontsize=9, color="#fdae61", fontweight="bold",
        arrowprops=dict(arrowstyle="->", color="#fdae61", lw=1.2),
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.85, edgecolor="#fdae61"),
    )

    # Manual threshold line
    manual_log_umi = np.log10(manual_umi)
    ax.axhline(y=manual_log_umi, color="#2c7bb6", linestyle=":", linewidth=2.0, alpha=0.9)
    ax.annotate(
        f"Manual UMI>{manual_umi:,}",
        xy=(log_rank[-1] * 0.98, manual_log_umi),
        xytext=(log_rank[-1] * 0.7, manual_log_umi + 0.15),
        fontsize=10, color="#2c7bb6", fontweight="bold", ha="right",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.85, edgecolor="#2c7bb6"),
    )

    ax.set_xlabel("log₁₀(Rank)", fontsize=14)
    ax.set_ylabel("log₁₀(Total UMI)", fontsize=14)
    ax.set_title("Global Barcode-Rank Plot (Knee Plot)\n"
                 f"643,038 barcodes × 115,983 genes | Knee UMI={knee_umi:,.0f} | Inflection UMI={infl_umi:,.0f} | Manual UMI>{manual_umi:,}",
                 fontsize=13, fontweight="bold")
    ax.legend(loc="lower left", fontsize=9, markerscale=8, framealpha=0.9)
    ax.grid(True, alpha=0.2)
    plt.tight_layout()

    plot_path = output_dir / "barcode_rank_kneeplot.png"
    fig.savefig(plot_path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    log.info(f"Main knee plot saved: {plot_path}")

    # ── Figure 2: Derivatives + Curvature ─────────────────────────────────
    fig, axes = plt.subplots(3, 1, figsize=(16, 14), sharex=True)

    # Panel A: Smoothed curve
    ax0 = axes[0]
    ax0.scatter(log_rank[::200], log_umi[::200], c="#cccccc", s=0.3, rasterized=True)
    ax0.plot(x_eval, spline_y, color="black", linewidth=2)
    ax0.axvline(x=knee_log_rank, color="#d7191c", linestyle="--", linewidth=1.5)
    ax0.axvline(x=infl_log_rank, color="#fdae61", linestyle="-.", linewidth=1.5)
    ax0.set_ylabel("log₁₀(Total UMI)", fontsize=12)
    ax0.set_title("A. Smoothing Spline Fit", fontsize=13, fontweight="bold", loc="left")
    ax0.legend(["Spline", "Knee", "Inflection"], fontsize=9, loc="lower left")
    ax0.grid(True, alpha=0.2)

    # Panel B: 1st derivative (slope)
    ax1 = axes[1]
    ax1.plot(x_eval, spline_y1, color="#2c7bb6", linewidth=1.8)
    ax1.axhline(y=-1.0, color="gray", linestyle=":", linewidth=1, alpha=0.5)
    ax1.axvline(x=knee_log_rank, color="#d7191c", linestyle="--", linewidth=1.5, alpha=0.7)
    ax1.axvline(x=infl_log_rank, color="#fdae61", linestyle="-.", linewidth=1.5, alpha=0.7)
    # Mark where y' crosses -1 (alternative knee definition)
    try:
        cross_idx = np.argmin(np.abs(spline_y1 + 1.0))
        ax1.scatter(x_eval[cross_idx], spline_y1[cross_idx], color="purple", s=60, zorder=10,
                    marker="D", edgecolors="white", linewidths=1)
        ax1.annotate(f"y'=-1 at rank≈{10**x_eval[cross_idx]:,.0f}",
                     xy=(x_eval[cross_idx], spline_y1[cross_idx]),
                     xytext=(x_eval[cross_idx]+0.2, spline_y1[cross_idx]+0.3),
                     fontsize=8, color="purple")
    except Exception:
        pass
    ax1.set_ylabel("1st Derivative y'(x)", fontsize=12)
    ax1.set_title("B. First Derivative (Slope)", fontsize=13, fontweight="bold", loc="left")
    ax1.grid(True, alpha=0.2)

    # Panel C: 2nd derivative + Curvature
    ax2 = axes[2]
    ax2.plot(x_eval, spline_y2, color="#d7191c", linewidth=1.5, alpha=0.7, label="y'' (2nd deriv)")
    ax2.plot(x_eval, curvature * 10, color="#fdae61", linewidth=1.8, label="Curvature κ (×10)")
    ax2.axhline(y=0, color="gray", linestyle=":", linewidth=0.8)
    ax2.axvline(x=knee_log_rank, color="#d7191c", linestyle="--", linewidth=1.5, alpha=0.7)
    ax2.axvline(x=infl_log_rank, color="#fdae61", linestyle="-.", linewidth=1.5, alpha=0.7)
    # Mark max curvature
    curv_max_idx = np.argmax(curvature[len(curvature)//20:]) + len(curvature)//20
    ax2.scatter(x_eval[curv_max_idx], curvature[curv_max_idx] * 10, color="#fdae61",
                s=100, zorder=10, marker="*", edgecolors="black", linewidths=0.8)
    ax2.set_xlabel("log₁₀(Rank)", fontsize=13)
    ax2.set_ylabel("2nd Derivative / Curvature", fontsize=12)
    ax2.set_title("C. Second Derivative (y'') & Curvature Reference",
                  fontsize=13, fontweight="bold", loc="left")
    ax2.legend(fontsize=9, loc="upper right")
    ax2.grid(True, alpha=0.2)

    plt.tight_layout()
    deriv_path = output_dir / "barcode_rank_derivatives.png"
    fig.savefig(deriv_path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    log.info(f"Derivative plot saved: {deriv_path}")


# ═══════════════════════════════════════════════════════════════════════════════
# 5. Summary report
# ═══════════════════════════════════════════════════════════════════════════════

def print_summary_report(
    metrics_df: pd.DataFrame,
    br_result: dict,
    whitelists: dict,
):
    """Print a formatted summary report to stdout."""
    total = len(metrics_df)
    knee_umi = br_result["knee_umi"]
    infl_umi = br_result["inflection_umi"]

    wl_manual, mask_m, umi_m, gene_m = whitelists["v1_manual"]
    wl_knee, mask_k, umi_k, gene_k = whitelists["v2_knee"]
    wl_infl, mask_i, umi_i, gene_i = whitelists["v2_inflection"]

    # Overlap analysis
    set_m = set(wl_manual)
    set_k = set(wl_knee)
    set_i = set(wl_infl)

    all_three = len(set_m & set_k & set_i)
    m_only = len(set_m - set_k - set_i)
    k_only = len(set_k - set_m - set_i)
    i_only = len(set_i - set_m - set_k)

    print("\n" + "=" * 72)
    print("  GLOBAL BARCODE-RANK QC — WHITELIST SUMMARY REPORT")
    print("=" * 72)
    print(f"  Total barcodes analyzed : {total:,}")
    print(f"  Total genes (features)  : 115,983")
    print(f"  Nonzero entries         : {metrics_df['total_umi'].sum():,.0f}")
    print("-" * 72)
    print(f"  {'Version':<25s} {'UMI>':>10s} {'Genes>':>8s} {'Pass':>10s} {'Pass%':>8s} {'Reject':>10s}")
    print("-" * 72)
    for name, (wl, mask, umi_th, gene_th) in whitelists.items():
        n_pass = len(wl)
        label = {"v1_manual": "Manual (static)", "v2_knee": "Knee (max distance from diagonal)",
                 "v2_inflection": "Inflection (min 1st deriv)"}[name]
        print(f"  {label:<25s} {umi_th:>10,.0f} {gene_th:>8,} {n_pass:>10,} {n_pass/total*100:>7.2f}% {total-n_pass:>10,}")
    print("-" * 72)
    print(f"  {'Knee UMI':<25s} {knee_umi:>10,.0f}")
    print(f"  {'Inflection UMI':<25s} {infl_umi:>10,.0f}")
    print(f"  {'Manual UMI threshold':<25s} {umi_m:>10,}")
    print("-" * 72)
    print(f"  Overlap (all 3 agree)   : {all_three:>10,}")
    print(f"  Manual-only             : {m_only:>10,}")
    print(f"  Knee-only               : {k_only:>10,}")
    print(f"  Inflection-only         : {i_only:>10,}")
    print("=" * 72)


# ═══════════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    args = parse_args()

    log.info("=" * 70)
    log.info("Task 2: Global Barcode-Rank QC & Whitelist Generation")
    log.info(f"  Input       : {args.input}")
    log.info(f"  Output base : {args.output_base}")
    log.info(f"  Manual thresholds: UMI>{args.manual_umi}, Genes>{args.manual_gene}")
    log.info("=" * 70)

    # ── Load data ─────────────────────────────────────────────────────────
    t_start = time.perf_counter()
    log.info("Loading merged AnnData (backed mode) …")
    adata = ad.read_h5ad(args.input, backed="r")
    log.info(f"  Shape: {adata.shape}")

    # ── 1. Compute per-cell metrics ───────────────────────────────────────
    metrics_df = compute_cell_metrics_chunked(adata, args.chunk_size)
    adata.file.close()  # close backed file
    del adata
    gc.collect()

    # ── 2. BarcodeRanks analysis ──────────────────────────────────────────
    log.info("Running barcodeRanks analysis …")
    br_result = barcode_ranks(metrics_df["total_umi"].values)

    # ── 3. Generate whitelists ────────────────────────────────────────────
    log.info("Generating whitelists …")
    whitelists = generate_whitelists(
        metrics_df, br_result,
        args.manual_umi, args.manual_gene,
        args.output_base,
    )

    # ── 4. Plots ──────────────────────────────────────────────────────────
    if not args.no_plots:
        log.info("Generating QC plots …")
        plot_barcode_rank(
            metrics_df, br_result, whitelists,
            args.output_base / "diagnostics",
        )

    # ── 5. Summary ────────────────────────────────────────────────────────
    print_summary_report(metrics_df, br_result, whitelists)

    t_total = time.perf_counter() - t_start
    log.info(f"\nTotal elapsed: {t_total/60:.1f} min")
    log.info("Done.")


if __name__ == "__main__":
    main()
