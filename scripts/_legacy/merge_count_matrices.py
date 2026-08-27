#!/usr/bin/env python3
"""
Task 1: Merge mRNA count matrices from 48 lanes of scprocess output.

Sub-tasks:
  1.1 - Stream-read 48 decontx HDF5 files, convert float64->float32
  1.2 - Add lane-specific suffix (-L01..-L48) to cell barcodes
  1.3 - vstack CSR matrices along cell axis, validate, serialize to .h5ad

Usage:
    conda activate scp_analysis
    python merge_count_matrices.py

Output:
    /data/yunzliu/results/guide_extraction/merged/merged_48lanes_raw.h5ad
"""

import gc
import os
import sys
import time
import logging
import argparse
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

# ─── Logging ────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)

# ─── Constants ──────────────────────────────────────────────────────────────
INPUT_ROOT = Path("/data/yunzliu/raw_fastq/full/scprocess_out")
OUTPUT_DIR = Path("/data/yunzliu/results/guide_extraction/merged")
OUTPUT_FILE = "merged_48lanes_raw.h5ad"
N_LANES = 48
DATE_STAMP = "2026-05-15"


def parse_args():
    parser = argparse.ArgumentParser(description="Merge 48-lane scprocess decontx HDF5 into single AnnData")
    parser.add_argument("--input-root", type=Path, default=INPUT_ROOT, help="Root dir of scprocess_out")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR, help="Output directory for merged .h5ad")
    parser.add_argument("--output-name", type=str, default=OUTPUT_FILE, help="Output filename")
    parser.add_argument("--n-lanes", type=int, default=N_LANES, help="Number of lanes")
    parser.add_argument("--date-stamp", type=str, default=DATE_STAMP, help="Date stamp in filenames")
    parser.add_argument("--dry-run", action="store_true", help="Only locate files, skip processing")
    return parser.parse_args()


def locate_files(input_root: Path, n_lanes: int, date_stamp: str) -> list[Path]:
    """Locate all 48 lane HDF5 files; raise if any missing."""
    files = []
    for i in range(1, n_lanes + 1):
        lane = f"lane{i:02d}"
        fpath = (
            input_root / lane / "output" / f"{lane}_ambient"
            / f"ambient_{lane}" / f"decontx_{lane}_{date_stamp}_filtered.h5"
        )
        if not fpath.is_file():
            raise FileNotFoundError(f"Missing: {fpath}")
        files.append(fpath)
    log.info(f"Located {len(files)}/{n_lanes} HDF5 files")
    return files


def read_one_lane(fpath: Path) -> tuple:
    """
    Read a single decontx HDF5 and return CSR components as float32.
    Returns: (data_f32, indices, indptr, shape, barcodes, gene_ids, gene_names)
    """
    with h5py.File(fpath, "r") as f:
        matrix = f["matrix"]
        data = matrix["data"][:].astype(np.float32, copy=False)
        indices = matrix["indices"][:]
        indptr = matrix["indptr"][:]
        shape = tuple(matrix["shape"][:])
        barcodes = matrix["barcodes"][:].astype(str)
        gene_ids = matrix["features"]["id"][:].astype(str)
        gene_names = matrix["features"]["name"][:].astype(str)
    return data, indices, indptr, shape, barcodes, gene_ids, gene_names


def merge_lanes(
    files: list[Path],
    n_lanes: int,
) -> tuple[sp.csr_matrix, np.ndarray, np.ndarray, np.ndarray]:
    """
    Stream-read 48 lane HDF5 files, add barcode suffixes,
    accumulate CSR blocks, then vstack into single sparse matrix.

    Returns:
        merged_csr  : scipy.sparse.csr_matrix (float32)
        all_barcodes: 1-D array of suffixed barcodes
        gene_ids    : 1-D array of ENSG IDs (from lane01 as reference)
        gene_names  : 1-D array of gene symbols
    """
    csr_blocks: list[sp.csr_matrix] = []
    barcode_blocks: list[np.ndarray] = []
    ref_gene_ids = None
    ref_gene_names = None
    total_nnz = 0

    for idx, fpath in enumerate(files, start=1):
        lane_tag = f"L{idx:02d}"
        t0 = time.perf_counter()

        data, indices, indptr, shape, barcodes, gene_ids, gene_names = read_one_lane(fpath)

        # Validate gene consistency against lane01
        if ref_gene_ids is None:
            ref_gene_ids = gene_ids
            ref_gene_names = gene_names
        else:
            if not np.array_equal(ref_gene_ids, gene_ids):
                log.error(f"Gene mismatch in lane {idx}! Cannot merge safely.")
                log.error(f"  Expected {len(ref_gene_ids)} genes, got {len(gene_ids)}")
                raise ValueError(f"Gene set mismatch at lane {idx:02d}")

        # Build CSC block (CellBender stores genes×cells in CSC format)
        csc = sp.csc_matrix((data, indices, indptr), shape=shape, dtype=np.float32)
        # Convert to CSR then transpose: genes×cells → cells×genes
        # (vstack along axis=0 requires uniform gene count across lanes)
        csr = csc.T.tocsr()
        csr_blocks.append(csr)
        del data, indices, indptr, csc  # free intermediates

        # Suffixed barcodes
        suffixed = np.array([f"{bc}-{lane_tag}" for bc in barcodes], dtype=str)
        barcode_blocks.append(suffixed)

        nnz_this = csr.nnz
        total_nnz += nnz_this
        elapsed = time.perf_counter() - t0
        log.info(
            f"  [{idx:02d}/{n_lanes}] {fpath.name} — "
            f"cells={csr.shape[0]:,}  genes={csr.shape[1]:,}  nnz={nnz_this:,}  "
            f"time={elapsed:.1f}s  mem_est={csr.data.nbytes/1024**2:.0f}MB"
        )
        gc.collect()

    # ── vstack all blocks ────────────────────────────────────────────────
    log.info(f"vstacking {len(csr_blocks)} CSR blocks ({total_nnz:,} total nnz) …")
    t0 = time.perf_counter()
    merged_csr = sp.vstack(csr_blocks, format="csr", dtype=np.float32)
    elapsed = time.perf_counter() - t0
    log.info(f"vstack complete in {elapsed:.1f}s — shape={merged_csr.shape}")

    # Free block list memory
    del csr_blocks
    gc.collect()

    all_barcodes = np.concatenate(barcode_blocks)
    del barcode_blocks
    gc.collect()

    return merged_csr, all_barcodes, ref_gene_ids, ref_gene_names


def build_and_save_anndata(
    merged_csr: sp.csr_matrix,
    barcodes: np.ndarray,
    gene_ids: np.ndarray,
    gene_names: np.ndarray,
    output_path: Path,
    files: list[Path],
) -> None:
    """Construct AnnData from merged CSR, add lane metadata, save to .h5ad."""

    log.info("Building AnnData object …")
    t0 = time.perf_counter()

    # ── var (genes) ──────────────────────────────────────────────────────
    var_df = pd.DataFrame(
        {"gene_id": gene_ids, "gene_name": gene_names},
        index=gene_ids,
    )
    var_df.index.name = "ensg_id"

    # ── obs (cells) ──────────────────────────────────────────────────────
    obs_df = pd.DataFrame(index=barcodes)

    # Parse lane from barcode suffix
    obs_df["lane"] = [bc.split("-")[-1] for bc in barcodes]
    obs_df["source_file"] = ""  # filled below if needed

    # ── AnnData (cells × genes) ──────────────────────────────────────────
    # AnnData expects obs × var; our CSR is (cells, genes)
    adata = ad.AnnData(
        X=merged_csr,
        obs=obs_df,
        var=var_df,
        dtype=np.float32,
    )

    # Attach provenance
    adata.uns["merge_info"] = {
        "n_input_files": len(files),
        "n_cells": merged_csr.shape[0],
        "n_genes": merged_csr.shape[1],
        "total_nnz": merged_csr.nnz,
        "sparsity_pct": round(merged_csr.nnz / (merged_csr.shape[0] * merged_csr.shape[1]) * 100, 3),
        "source": "scprocess decontx filtered (CellBender ambient RNA removal)",
        "date_merged": time.strftime("%Y-%m-%d %H:%M:%S"),
    }

    elapsed = time.perf_counter() - t0
    log.info(f"AnnData built in {elapsed:.1f}s — {adata}")

    # ── Save ─────────────────────────────────────────────────────────────
    output_path.parent.mkdir(parents=True, exist_ok=True)
    log.info(f"Saving to {output_path} …")
    t0 = time.perf_counter()
    adata.write_h5ad(output_path, compression="gzip")
    elapsed = time.perf_counter() - t0
    fsize_mb = output_path.stat().st_size / 1024**2
    log.info(f"Saved {fsize_mb:.0f} MB in {elapsed:.1f}s")


def main():
    args = parse_args()

    log.info("=" * 70)
    log.info("Task 1: Merge 48-lane scprocess mRNA count matrices")
    log.info(f"  Input root : {args.input_root}")
    log.info(f"  Output dir : {args.output_dir}")
    log.info(f"  Lanes      : {args.n_lanes}")
    log.info("=" * 70)

    # 1. Locate files
    files = locate_files(args.input_root, args.n_lanes, args.date_stamp)
    if args.dry_run:
        log.info("Dry-run complete. Exiting.")
        return

    # 2. Merge
    t_start = time.perf_counter()
    merged_csr, barcodes, gene_ids, gene_names = merge_lanes(files, args.n_lanes)

    # 3. Validate
    n_cells_expected = merged_csr.shape[0]
    n_genes_expected = len(gene_ids)
    assert len(barcodes) == n_cells_expected, "Barcode count mismatch!"
    log.info(f"Validation OK: {n_cells_expected:,} cells × {n_genes_expected:,} genes")
    log.info(f"  Nonzeros: {merged_csr.nnz:,}  Sparsity: {merged_csr.nnz/(n_cells_expected*n_genes_expected)*100:.3f}%")

    # 4. Build AnnData & save
    output_path = args.output_dir / args.output_name
    build_and_save_anndata(merged_csr, barcodes, gene_ids, gene_names, output_path, files)

    # 5. Summary
    t_total = time.perf_counter() - t_start
    log.info("=" * 70)
    log.info(f"Merge complete! Total elapsed: {t_total/60:.1f} min")
    log.info(f"Output: {output_path}")
    log.info(f"Cells: {n_cells_expected:,}  Genes: {n_genes_expected:,}")
    log.info("=" * 70)


if __name__ == "__main__":
    main()
