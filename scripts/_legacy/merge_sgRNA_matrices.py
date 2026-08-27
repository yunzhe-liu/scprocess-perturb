#!/usr/bin/env python3
"""
Merge per-lane sgRNA count matrices into a single unified count matrix.

Each lane's cell barcodes are suffixed with a lane identifier (e.g., "-L01")
to maintain consistency with mRNA count matrix conventions across lanes.

Input:
    Per-lane simpleaf quant output directories, each containing:
        quants_mat.mtx       — sparse count matrix (cells × guides)
        quants_mat_rows.txt  — cell barcodes (one per line)
        quants_mat_cols.txt  — guide feature IDs (one per line)

    Lanes are specified via a 3-column TSV (or the --list format):
        lane_id<tab>alevin_dir<tab>suffix
    e.g.:
        lane_01  /path/to/lane_01/simpleaf_quant/af_quant/alevin  -L01

Output:
    A standard 10x MEX trio in the output directory:
        matrix.mtx.gz     — vertically concatenated sparse matrix
        barcodes.tsv.gz   — cell barcodes with lane suffixes
        features.tsv.gz   — guide features (identical across lanes; taken from first)

Usage:
    # Via lane-list TSV
    python merge_sgRNA_matrices.py --lanes lanes.tsv --out merged/

    # Via explicit arguments (for quick ad-hoc merges)
    python merge_sgRNA_matrices.py \
        --id lane_01 --dir /path/to/lane_01/alevin --suffix -L01 \
        --id lane_02 --dir /path/to/lane_02/alevin --suffix -L02 \
        --out merged/
"""

import argparse
import gzip
import os
import sys
import time
from pathlib import Path

import numpy as np
from scipy.io import mmread, mmwrite
from scipy import sparse


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Merge per-lane sgRNA count matrices into a single unified matrix."
    )
    p.add_argument("--lanes", type=Path, default=None,
                   help="TSV file: lane_id<TAB>alevin_dir<TAB>suffix (no header)")
    p.add_argument("--id", action="append", dest="ids", default=[],
                   help="Lane identifier (repeatable, used with --dir and --suffix)")
    p.add_argument("--dir", action="append", dest="dirs", default=[],
                   help="Path to alevin output directory (repeatable)")
    p.add_argument("--suffix", action="append", dest="suffixes", default=[],
                   help="Barcode suffix for this lane, e.g. '-L01' (repeatable)")
    p.add_argument("--out", type=Path, required=True,
                   help="Output directory for merged MEX trio")
    p.add_argument("--prefix", type=str, default="merged",
                   help="Output file prefix (default: 'merged')")
    return p.parse_args()


def _find_file(directory: Path, candidates: list[str], desc: str) -> Path:
    """Try multiple candidate file names, supporting both compressed and uncompressed."""
    for name in candidates:
        p = directory / name
        if p.is_file():
            return p
        # Also try .gz variant
        p_gz = directory / (name + ".gz")
        if p_gz.is_file():
            return p_gz
    tried = ", ".join(candidates)
    raise FileNotFoundError(f"{desc} file not found in {directory}. Tried: {tried}")


def _read_lines(path: Path) -> list[str]:
    """Read lines from a file, handling gzip transparently."""
    if path.suffix == '.gz':
        import gzip
        with gzip.open(path, 'rt') as f:
            return [line.strip() for line in f]
    else:
        with open(path) as f:
            return [line.strip() for line in f]


def load_lane(alevin_dir: Path, suffix: str) -> tuple:
    """Load a single lane's quant output.  Supports both simpleaf (quants_mat.*)
    and hash_matcher (matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz) naming."""
    
    # Try simpleaf naming first, then hash_matcher naming
    mtx_path = _find_file(alevin_dir, 
        ["quants_mat.mtx", "matrix.mtx"], "MTX")
    rows_path = _find_file(alevin_dir, 
        ["quants_mat_rows.txt", "barcodes.tsv"], "barcodes")
    cols_path = _find_file(alevin_dir, 
        ["quants_mat_cols.txt", "features.tsv"], "features")

    mtx = mmread(str(mtx_path)).tocsc()
    barcodes = np.array(_read_lines(rows_path), dtype=str)
    features = np.array(_read_lines(cols_path), dtype=str)
    
    # Extract feature IDs (first column if TSV, whole line if single column)
    if features.size > 0 and '\t' in features[0]:
        features = np.array([f.split('\t')[0] for f in features], dtype=str)

    # Apply suffix
    if suffix:
        barcodes = np.array([f"{bc}{suffix}" for bc in barcodes], dtype=str)

    n_cells, n_guides = mtx.shape
    print(f"  Loaded: {alevin_dir} → {n_cells:,} cells × {n_guides} guides, "
          f"{mtx.nnz:,} non-zero, suffix='{suffix}'", file=sys.stderr)
    return mtx, barcodes, features


def merge_all(lanes: list[tuple[str, Path, str]], out_dir: Path, prefix: str) -> None:
    """Merge all lanes and write MEX trio.  Supports feature mismatch across lanes
    by taking the union of all features (missing features filled with 0)."""
    os.makedirs(out_dir, exist_ok=True)

    t0 = time.time()
    
    # ── Pass 1: collect all features across all lanes ──
    all_feature_ids = []
    feature_to_global = {}
    for lane_id, alevin_dir, suffix in lanes:
        # Just read features to build the global index
        cols_path = _find_file(alevin_dir, 
            ["quants_mat_cols.txt", "features.tsv"], "features")
        features = np.array(_read_lines(cols_path), dtype=str)
        if features.size > 0 and '\t' in features[0]:
            features = np.array([f.split('\t')[0] for f in features], dtype=str)
        for ft in features:
            if ft not in feature_to_global:
                feature_to_global[ft] = len(all_feature_ids)
                all_feature_ids.append(ft)
    n_global_features = len(all_feature_ids)
    print(f"Global feature set: {n_global_features} unique features", file=sys.stderr)

    # ── Pass 2: load matrices aligned to global feature space ──
    all_matrices = []
    all_barcodes = []
    total_cells = 0

    for lane_id, alevin_dir, suffix in lanes:
        mtx, barcodes, features = load_lane(alevin_dir, suffix)
        
        # Build local-to-global feature mapping
        local_to_global = np.array([feature_to_global.get(ft, 0) for ft in features], dtype=np.int32)
        
        # Remap matrix columns to global feature indices using COO format
        from scipy.sparse import csr_matrix
        mtx_coo = mtx.tocoo()
        global_cols = local_to_global[mtx_coo.col]
        remapped = csr_matrix(
            (mtx_coo.data, (mtx_coo.row, global_cols)),
            shape=(mtx.shape[0], n_global_features)
        ).tocsc()
        
        all_matrices.append(remapped)
        all_barcodes.append(barcodes)
        total_cells += mtx.shape[0]
        print(f"  {lane_id}: {mtx.shape[0]:,} cells, {len(features)}→{n_global_features} features, "
              f"{int(mtx.sum()):,} UMIs", file=sys.stderr)

    # Vertical stack
    print(f"\nMerging {len(lanes)} lanes...", file=sys.stderr)
    merged_mtx = sparse.vstack(all_matrices, format="csc")
    merged_barcodes = np.concatenate(all_barcodes)
    elapsed = time.time() - t0

    # Validate
    assert merged_mtx.shape[0] == total_cells, \
        f"Cell count mismatch: {merged_mtx.shape[0]} vs {total_cells}"
    assert merged_mtx.shape[1] == n_global_features
    assert len(merged_barcodes) == total_cells

    # Write MEX trio
    mtx_out = out_dir / f"{prefix}_matrix.mtx.gz"
    bc_out  = out_dir / f"{prefix}_barcodes.tsv.gz"
    ft_out  = out_dir / f"{prefix}_features.tsv.gz"

    print(f"Writing matrix ({merged_mtx.shape[0]:,} × {merged_mtx.shape[1]})...", file=sys.stderr)
    with gzip.open(mtx_out, "wb") as f:
        mmwrite(f, merged_mtx, field="integer")

    print(f"Writing barcodes ({len(merged_barcodes):,})...", file=sys.stderr)
    with gzip.open(bc_out, "wt") as f:
        f.write("\n".join(merged_barcodes) + "\n")

    print(f"Writing features ({len(all_feature_ids)})...", file=sys.stderr)
    with gzip.open(ft_out, "wt") as f:
        for feat in all_feature_ids:
            f.write(f"{feat}\t{feat}\tCRISPR Guide Capture\n")

    # Summary
    print(f"\n{'='*50}", file=sys.stderr)
    print(f"Merge complete ({elapsed:.1f}s)", file=sys.stderr)
    print(f"  Lanes:     {len(lanes)}", file=sys.stderr)
    print(f"  Cells:     {total_cells:,}", file=sys.stderr)
    print(f"  Features:  {n_global_features}", file=sys.stderr)
    print(f"  Non-zero:  {merged_mtx.nnz:,}", file=sys.stderr)
    print(f"  Output:    {out_dir}/", file=sys.stderr)
    for fname in [f"{prefix}_matrix.mtx.gz", f"{prefix}_barcodes.tsv.gz", f"{prefix}_features.tsv.gz"]:
        fpath = out_dir / fname
        if fpath.is_file():
            print(f"    {fname:40s} {fpath.stat().st_size/1024:.0f} KB", file=sys.stderr)


def main() -> None:
    args = parse_args()

    # Build lane list from either --lanes TSV or --id/--dir/--suffix triples
    lanes: list[tuple[str, Path, str]] = []

    if args.lanes:
        with open(args.lanes) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 2:
                    continue
                lane_id = parts[0].strip()
                alevin_dir = Path(parts[1].strip())
                suffix = parts[2].strip() if len(parts) >= 3 else f"-{lane_id}"
                lanes.append((lane_id, alevin_dir, suffix))

    if args.ids or args.dirs:
        ids = args.ids
        dirs = args.dirs
        suffixes = args.suffixes if args.suffixes else [f"-{i}" for i in ids]
        if len(ids) != len(dirs):
            print("ERROR: --id and --dir counts must match.", file=sys.stderr)
            sys.exit(1)
        if len(suffixes) < len(ids):
            suffixes += [f"-{ids[i]}" for i in range(len(suffixes), len(ids))]
        for i in range(len(ids)):
            lanes.append((ids[i], Path(dirs[i]), suffixes[i]))

    if not lanes:
        print("ERROR: No lanes specified. Use --lanes or --id/--dir/--suffix.", file=sys.stderr)
        sys.exit(1)

    merge_all(lanes, args.out, args.prefix)


if __name__ == "__main__":
    main()
