#!/usr/bin/env python3
"""
Generate per-lane whitelists (static, knee, inflection) by running
barcodeRanks (LOWESS + natural cubic spline) independently on each of
the 48 lanes from the merged mRNA count matrix.

Then consolidate with existing global whitelists into a single
comparison directory for downstream evaluation against published data.

Output directory:
  /data/yunzliu/results/guide_extraction/whitelist/all_versions/
    ├── global_static.tsv.gz        (symbolic link or copy)
    ├── global_knee.tsv.gz
    ├── global_inflection.tsv.gz
    ├── per_lane_static.tsv.gz
    ├── per_lane_knee.tsv.gz
    └── per_lane_inflection.tsv.gz
    └── thresholds_summary.csv      (all 6 versions' thresholds + stats)
"""

import gc
import os
import sys
import time
import logging
import argparse
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

# ─── Reuse barcodeRanks from existing module ────────────────────────────
sys.path.insert(0, str(Path(__file__).resolve().parent))
from barcode_rank_qc import barcode_ranks, compute_cell_metrics_chunked

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)

# ─── Paths ──────────────────────────────────────────────────────────────
INPUT_H5AD = Path("/data/yunzliu/results/guide_extraction/merged/merged_48lanes_raw.h5ad")
WHITELIST_BASE = Path("/data/yunzliu/results/guide_extraction/whitelist")
OUTPUT_DIR = WHITELIST_BASE / "all_versions"

MANUAL_UMI = 1000
MANUAL_GENE = 500
CHUNK_SIZE = 50000


def generate_per_lane_whitelists(metrics_df: pd.DataFrame) -> dict:
    """
    For each lane, run barcodeRanks independently to determine
    lane-specific knee and inflection UMI/gene thresholds.
    Returns dict with three merged whitelists (static, knee, inflection).
    """
    lanes = sorted(metrics_df["lane"].unique())
    log.info(f"Processing {len(lanes)} lanes independently …")

    results = {"static": [], "knee": [], "inflection": []}

    for lane in lanes:
        mask_lane = metrics_df["lane"] == lane
        df_lane = metrics_df[mask_lane]
        n_lane = len(df_lane)

        # ── Static threshold (per-lane, same rule) ───────────────────
        pass_static = df_lane[
            (df_lane["total_umi"] > MANUAL_UMI) & (df_lane["n_genes"] > MANUAL_GENE)
        ].index.tolist()
        results["static"].extend(pass_static)

        # ── BarcodeRanks on this lane ─────────────────────────────────
        br = barcode_ranks(df_lane["total_umi"].values)

        knee_umi = br["knee_umi"]
        infl_umi = br["inflection_umi"]

        # Derive gene thresholds from neighborhood
        knee_genes = _derive_gene_threshold(df_lane, knee_umi)
        infl_genes = _derive_gene_threshold(df_lane, infl_umi)

        # ── Knee whitelist ────────────────────────────────────────────
        pass_knee = df_lane[
            (df_lane["total_umi"] > knee_umi) & (df_lane["n_genes"] > knee_genes)
        ].index.tolist()
        results["knee"].extend(pass_knee)

        # ── Inflection whitelist ──────────────────────────────────────
        pass_infl = df_lane[
            (df_lane["total_umi"] > infl_umi) & (df_lane["n_genes"] > infl_genes)
        ].index.tolist()
        results["inflection"].extend(pass_infl)

        log.info(f"  {lane}: n={n_lane:,}  "
                 f"knee_UMI={knee_umi:.0f}  infl_UMI={infl_umi:.0f}  "
                 f"static={len(pass_static)}  knee={len(pass_knee)}  infl={len(pass_infl)}")

    return results


def _derive_gene_threshold(df_lane: pd.DataFrame, umi_threshold: float) -> int:
    """Derive gene threshold from median n_genes near the UMI threshold."""
    neighborhood = df_lane[
        (df_lane["total_umi"] >= umi_threshold * 0.9)
        & (df_lane["total_umi"] <= umi_threshold * 1.1)
    ]
    if len(neighborhood) > 0:
        return int(neighborhood["n_genes"].median())
    return 0


def save_whitelists(results: dict, output_dir: Path, total_cells: int) -> pd.DataFrame:
    """Save per-lane whitelists and return summary stats DataFrame."""
    output_dir.mkdir(parents=True, exist_ok=True)
    stats_rows = []

    for version_name, barcodes in results.items():
        fname = f"per_lane_{version_name}.tsv.gz"
        fpath = output_dir / fname

        # Validate uniqueness (barcodes already have lane suffixes from Task 1.2)
        barcodes_sorted = sorted(set(barcodes))
        if len(barcodes_sorted) != len(barcodes):
            log.warning(f"  {version_name}: duplicate barcodes detected ({len(barcodes)} → {len(barcodes_sorted)} unique)")

        pd.Series(barcodes_sorted, name="barcode").to_csv(
            fpath, index=False, header=False, compression="gzip"
        )
        fsize = fpath.stat().st_size / 1024**2

        n_pass = len(barcodes_sorted)
        stats_rows.append({
            "version": f"per_lane_{version_name}",
            "scope": "per_lane",
            "method": version_name,
            "n_pass": n_pass,
            "pct_pass": round(n_pass / total_cells * 100, 2),
            "file": str(fpath),
            "size_mb": round(fsize, 2),
        })
        log.info(f"  Saved per_lane_{version_name}: {n_pass:,} cells ({n_pass/total_cells*100:.1f}%) — {fsize:.1f} MB")

    return pd.DataFrame(stats_rows)


def collect_global_whitelists(input_base: Path, output_dir: Path, total_cells: int) -> pd.DataFrame:
    """Copy or symlink existing global whitelists into output_dir; return stats."""
    mapping = {
        "global_static": input_base / "v1_manual" / "whitelist_v1_manual.tsv.gz",
        "global_knee": input_base / "v2_barcodeRanks_knee" / "whitelist_v2_knee.tsv.gz",
        "global_inflection": input_base / "v2_barcodeRanks_inflection" / "whitelist_v2_inflection.tsv.gz",
    }

    stats_rows = []
    for version_name, src in mapping.items():
        dst = output_dir / f"{version_name}.tsv.gz"
        shutil.copy2(src, dst)
        fsize = dst.stat().st_size / 1024**2

        # Count entries
        n_pass = len(pd.read_csv(dst, header=None, compression="gzip"))
        method = version_name.replace("global_", "")
        stats_rows.append({
            "version": version_name,
            "scope": "global",
            "method": method,
            "n_pass": n_pass,
            "pct_pass": round(n_pass / total_cells * 100, 2),
            "file": str(dst),
            "size_mb": round(fsize, 2),
        })
        # Remove individual files (they're now in all_versions/)
        # Keep originals intact; just copy

    return pd.DataFrame(stats_rows)


def main():
    parser = argparse.ArgumentParser(description="Generate per-lane whitelists and consolidate all versions")
    parser.add_argument("--input", type=Path, default=INPUT_H5AD)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--whitelist-base", type=Path, default=WHITELIST_BASE)
    parser.add_argument("--chunk-size", type=int, default=CHUNK_SIZE)
    parser.add_argument("--skip-global-copy", action="store_true", help="Skip copying existing global whitelists")
    args = parser.parse_args()

    log.info("=" * 70)
    log.info("Per-Lane Whitelist Generation + All-Version Consolidation")
    log.info(f"  Input       : {args.input}")
    log.info(f"  Output dir  : {args.output_dir}")
    log.info("=" * 70)

    # ── Load data ───────────────────────────────────────────────────────
    t_start = time.perf_counter()
    log.info("Loading merged AnnData (backed mode) …")
    adata = ad.read_h5ad(args.input, backed="r")
    total_cells = adata.n_obs
    log.info(f"  Shape: {adata.shape}")

    # ── Compute per-cell metrics once ────────────────────────────────────
    metrics_df = compute_cell_metrics_chunked(adata, args.chunk_size)
    adata.file.close()
    del adata
    gc.collect()

    # ── Generate per-lane whitelists ─────────────────────────────────────
    log.info("Generating per-lane whitelists …")
    per_lane_results = generate_per_lane_whitelists(metrics_df)

    # ── Save per-lane whitelists ─────────────────────────────────────────
    stats_per_lane = save_whitelists(per_lane_results, args.output_dir, total_cells)

    # ── Collect existing global whitelists ───────────────────────────────
    if not args.skip_global_copy:
        log.info("Collecting existing global whitelists …")
        stats_global = collect_global_whitelists(args.whitelist_base, args.output_dir, total_cells)
    else:
        stats_global = pd.DataFrame()

    # ── Merge stats ──────────────────────────────────────────────────────
    all_stats = pd.concat([stats_global, stats_per_lane], ignore_index=True)
    stats_path = args.output_dir / "thresholds_summary.csv"
    all_stats.to_csv(stats_path, index=False)
    log.info(f"Summary saved: {stats_path}")

    # ── Final summary ────────────────────────────────────────────────────
    total_elapsed = time.perf_counter() - t_start
    log.info("\n" + "=" * 70)
    log.info("ALL WHITELISTS CONSOLIDATED")
    log.info(f"  Output directory: {args.output_dir}")
    log.info(f"  Total cells: {total_cells:,}")
    log.info("=" * 70)
    for _, row in all_stats.iterrows():
        log.info(f"  {row['version']:<28s} {row['n_pass']:>10,} ({row['pct_pass']:>6.2f}%)")
    log.info("=" * 70)
    log.info(f"Total elapsed: {total_elapsed/60:.1f} min")


if __name__ == "__main__":
    main()
