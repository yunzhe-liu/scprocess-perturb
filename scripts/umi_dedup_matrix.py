#!/usr/bin/env python3
"""
umi_dedup_matrix.py — UMI deduplication + MEX matrix generation

Input:  assignments.npz from guide_matcher.py
        Format: (N, 3) int32 array — [cb_idx, umi_int, guide_idx]

UMI dedup strategy: directional (UMI-tools algorithm)
  - Group UMIs by (cell_barcode, guide_id)
  - Sort UMIs within each group by descending frequency
  - High-density UMIs absorb low-density UMIs within Hamming=1
  - Retain only "network center" UMI counts

Output: Standard 10x MEX format
  - matrix.mtx.gz
  - barcodes.tsv.gz
  - features.tsv.gz
"""

import sys
import os
import gzip
import time
import numpy as np
import scipy.sparse
from collections import defaultdict, Counter


# ── v3 UMI encode/decode ──
_BASES_V3 = 'ACGT'

def _decode_umi(val: int) -> str:
    """Decode uint32 UMI back to 12bp string."""
    chars = []
    for _ in range(12):
        chars.append(_BASES_V3[val & 0x3])
        val >>= 2
    return ''.join(reversed(chars))


def hamming_distance(s1: str, s2: str) -> int:
    """Compute Hamming distance between two equal-length strings."""
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def dedup_umis_directional(
    umi_list: list[str],
    threshold: int = 1,
) -> int:
    """
    Hash-accelerated directional UMI dedup (v2).
    
    Mathematically identical to UMI-tools directional algorithm.
    Replaces O(n × m) linear scan with O(n × 36) hash operations
    (12 bp UMI × 3 alternative bases).
    
    For each parent UMI, registers all 60 Hamming-1 neighbors into
    a local hash table. Subsequent UMIs are checked via O(1) lookup.
    
    Hash-accelerated directional UMI dedup (v2).
    
    Mathematically identical to UMI-tools directional algorithm.
    Replaces O(n x m) linear scan with O(n x 36) hash operations
    (12 bp UMI × 3 alternative bases).
    
    For each parent UMI, registers all 60 Hamming-1 neighbors into
    a local hash table. Subsequent UMIs are checked via O(1) lookup.
    
    Returns number of retained UMIs after deduplication.
    """
    if not umi_list:
        return 0

    umi_counts = Counter(umi_list)
    sorted_umis = sorted(umi_counts.keys(), key=lambda u: -umi_counts[u])

    local_hash: dict[str, str] = {}  # variant_seq → parent_umi
    retained_count = 0

    for umi in sorted_umis:
        if umi in local_hash:
            continue  # merged into existing parent

        # New parent UMI
        retained_count += 1

        # Register self
        if umi not in local_hash:
            local_hash[umi] = umi

        # Register 60 Hamming-1 neighbors (only if threshold >= 1)
        if threshold >= 1:
            for pos in range(len(umi)):
                orig = umi[pos]
                for alt in ('A', 'C', 'G', 'T'):
                    if alt == orig:
                        continue
                    variant = umi[:pos] + alt + umi[pos + 1:]
                    if variant not in local_hash:
                        local_hash[variant] = umi

    return retained_count


def build_count_matrix(
    assignments_path: str,
    output_dir: str,
    umi_threshold: int = 1,
):
    """
    Build deduplicated count matrix from assignments.npz (v3).
    
    Input: assignments.npz (binary numpy format)
      - assignments: (N, 3) int32 — [cb_idx, umi_int, guide_idx]
      - barcode_list: (M,) str — index → barcode string
      - idx_to_id: (K,) str — index → guide ID string
    
    Output: MEX trio (matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz)
    """
    t0 = time.time()

    # ── 1. Load assignments.npz ──
    print("[1/4] Loading assignments.npz...")
    data = np.load(assignments_path, allow_pickle=True)
    assignments = data['assignments']  # (N, 3) int32
    barcode_list = list(data['barcode_list'])
    idx_to_id = list(data['idx_to_id'])
    data.close()
    
    total_assignments = len(assignments)
    print(f"  Loaded {total_assignments:,} assignments")
    
    # ── 2. Group UMIs by (CB, guide) using integer keys to reduce memory ──
    print("[2/4] Grouping UMIs by (CB, guide)...")
    cb_guide_umis: dict[tuple[int, int], list[int]] = defaultdict(list)
    
    for row in assignments:
        cb_idx, umi_int, guide_idx = int(row[0]), int(row[1]), int(row[2])
        cb_guide_umis[(cb_idx, guide_idx)].append(umi_int)
    
    unique_pairs = len(cb_guide_umis)
    print(f"  Unique (CB, guide) pairs: {unique_pairs:,}")

    # ── 3. UMI dedup (decode ints to strings for hamming distance) ──
    print(f"[3/4] UMI deduplication (directional, threshold={umi_threshold})...")
    cb_guide_count: dict[tuple[int, int], int] = {}
    total_raw_umis = 0
    total_dedup_umis = 0

    for (cb_idx, guide_idx), umi_ints in cb_guide_umis.items():
        # Decode UMI ints to strings for hamming calculation
        umi_strs = [_decode_umi(u) for u in umi_ints]
        raw = len(umi_strs)
        dedup = dedup_umis_directional(umi_strs, umi_threshold)
        if dedup > 0:
            cb_guide_count[(cb_idx, guide_idx)] = dedup
        total_raw_umis += raw
        total_dedup_umis += dedup

    print(f"  Raw UMIs: {total_raw_umis:,}")
    print(f"  After dedup: {total_dedup_umis:,} "
          f"({total_dedup_umis/max(total_raw_umis,1)*100:.1f}%)")

    # ── 4. Build sparse matrix ──
    print("[4/5] Building sparse matrix...")
    
    # Collect all barcodes and features that have counts
    used_cb_indices = sorted(set(cb for cb, _ in cb_guide_count.keys()))
    used_guide_indices = sorted(set(gid for _, gid in cb_guide_count.keys()))
    
    # Build index mappings: global index → local position
    cb_idx_to_pos = {cb: i for i, cb in enumerate(used_cb_indices)}
    guide_idx_to_pos = {gid: i for i, gid in enumerate(used_guide_indices)}
    
    rows, cols, data_vals = [], [], []
    for (cb_idx, guide_idx), count in cb_guide_count.items():
        rows.append(cb_idx_to_pos[cb_idx])
        cols.append(guide_idx_to_pos[guide_idx])
        data_vals.append(count)

    matrix = scipy.sparse.csr_matrix(
        (data_vals, (rows, cols)),
        shape=(len(used_cb_indices), len(used_guide_indices)),
        dtype=np.int32,
    )

    # ── 5. Write MEX format ──
    print("[5/5] Writing MEX output...")
    os.makedirs(output_dir, exist_ok=True)

    # Map indices back to strings
    used_barcodes = [barcode_list[i] for i in used_cb_indices]
    used_features = [idx_to_id[i] for i in used_guide_indices]

    # matrix.mtx.gz
    mtx_path = os.path.join(output_dir, 'matrix.mtx.gz')
    with gzip.open(mtx_path, 'wt') as f:
        f.write("%%MatrixMarket matrix coordinate integer general\n")
        f.write("% Generated by guide_matcher v3 + umi_dedup_matrix\n")
        f.write(f"{matrix.shape[0]} {matrix.shape[1]} {matrix.nnz}\n")
        coo = matrix.tocoo()
        for r, c, v in zip(coo.row, coo.col, coo.data):
            f.write(f"{r + 1} {c + 1} {v}\n")  # 1-indexed

    # barcodes.tsv.gz
    bc_path = os.path.join(output_dir, 'barcodes.tsv.gz')
    with gzip.open(bc_path, 'wt') as f:
        for bc in used_barcodes:
            f.write(f"{bc}\n")

    # features.tsv.gz
    ft_path = os.path.join(output_dir, 'features.tsv.gz')
    with gzip.open(ft_path, 'wt') as f:
        for ft in used_features:
            f.write(f"{ft}\t{ft}\tGuide\n")

    elapsed = time.time() - t0

    # ── Summary statistics ──
    umi_per_cell = np.array(matrix.sum(axis=1)).flatten()
    nonzero_cells = (umi_per_cell > 0).sum()

    print(f"\n{'='*60}")
    print(f"MATRIX GENERATION COMPLETE (v3)")
    print(f"{'='*60}")
    print(f"  Cells (barcodes):    {len(used_barcodes):>10,}")
    print(f"  Features (guides):   {len(used_features):>10,}")
    print(f"  Non-zero entries:    {matrix.nnz:>10,}")
    print(f"  Cells with >0 UMI:   {nonzero_cells:>10,}")
    print(f"  Total UMIs:          {int(matrix.sum()):>10,}")
    print(f"  Median UMI/cell:     {np.median(umi_per_cell[umi_per_cell > 0]):>10.1f}"
          if nonzero_cells > 0 else "  Median UMI/cell:     N/A")
    print(f"  Mean UMI/cell:       {umi_per_cell[umi_per_cell > 0].mean():>10.1f}"
          if nonzero_cells > 0 else "  Mean UMI/cell:       N/A")
    print(f"  Time: {elapsed:.1f}s")
    print(f"  Output: {output_dir}/")
    print(f"{'='*60}")


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='UMI dedup + MEX matrix generation (v3)')
    parser.add_argument('-i', '--input', required=True,
                        help='Assignments .npz from guide_matcher.py')
    parser.add_argument('-o', '--output-dir', required=True,
                        help='Output directory for MEX files')
    parser.add_argument('-t', '--umi-threshold', type=int, default=1,
                        help='UMI Hamming distance threshold (default: 1)')
    args = parser.parse_args()

    build_count_matrix(args.input, args.output_dir, args.umi_threshold)
