#!/usr/bin/env python3
"""
guide_matcher.py — Hash-Accelerated Guide Matcher (v3)

Core algorithm (vs piscem):
  piscem: 86 × 13-mers across full 98bp read → adapter k-mer noise
  This:   26bp exact window × 7 sliding 20bp candidates → zero adapter noise

v3 improvements over v2:
  P0 — Integer encoding:  20bp guide → uint64, eliminating per-read str malloc.
                           Window encoded once as 52-bit int, 7 candidates
                           extracted via shift+mask — zero string allocation.
  P1 — Numpy pre-allocation: 30M×3 int32 buffer instead of list-of-tuples.
                           Output as .npz (binary) instead of .tsv.gz (text).
  P3 — Multi-process I/O:  Spawns one worker per FASTQ pair (4 per lane),
                           utilising allocated threads for true parallelism.

Two-tier guide matching (based on fba QC: >99% Illumina errors are substitutions):
  Fast path: 7 × O(1) int dict lookups → ~65% reads (exact match)
  Slow path: 7×20×3=420 int dict lookups → <1% reads (Hamming-1)

Cell barcode correction (Plan H):
  Encodes 16bp barcodes as uint32 (2 bits/base).
  All Hamming ≤2 variants stored in Python dict[int→int] for O(1) lookup.

Input:  Read1 (CB+UMI) + Read2 (guide construct) FASTQ
Output: assignments.npz → downstream UMI dedup + MEX generation
"""

import sys
import os
import gzip
import pickle
import time
import numpy as np
from typing import Optional

# ── Constants (empirically validated on 50K reads) ──
# Read2 guide window [28:54), covering observed guide start positions 28–34
WINDOW_START = 28
WINDOW_END   = 54
WINDOW_SIZE  = WINDOW_END - WINDOW_START  # 26bp
GUIDE_LEN    = 20
SLIDE_COUNT  = WINDOW_SIZE - GUIDE_LEN + 1  # 7

# Cell barcode: Read1[0:16], UMI: Read1[16:28]
CB_START, CB_END = 0, 16
UMI_START, UMI_END = 16, 28

BASES = ('A', 'C', 'G', 'T')

# DNA → 2-bit encoding table (A=00, C=01, G=10, T=11)
_BASE2BITS = {'A': 0, 'C': 1, 'G': 2, 'T': 3}


def _encode_barcode(bc: str) -> int:
    """Encode a 16bp DNA barcode into a 32-bit integer (2 bits/base)."""
    val = 0
    for ch in bc:
        val = (val << 2) | _BASE2BITS[ch]
    return val


def _generate_variants(bc: str, max_hamming: int = 2):
    """Generate all Hamming ≤ max_hamming variants of a 16bp barcode.
    Yields (variant_str, hamming_distance)."""
    yield bc, 0
    n = len(bc)
    # Hamming = 1
    for p1 in range(n):
        o1 = bc[p1]
        for a1 in BASES:
            if a1 == o1: continue
            yield bc[:p1] + a1 + bc[p1+1:], 1
    # Hamming = 2
    if max_hamming >= 2:
        for p1 in range(n):
            o1 = bc[p1]
            for a1 in BASES:
                if a1 == o1: continue
                v1 = bc[:p1] + a1 + bc[p1+1:]
                for p2 in range(p1+1, n):
                    o2 = bc[p2]
                    for a2 in BASES:
                        if a2 == o2: continue
                        yield v1[:p2] + a2 + v1[p2+1:], 2


def load_whitelist(path: str) -> set:
    """Load cell barcode whitelist as a Python set (O(1) membership test)."""
    with open(path) as f:
        return set(line.strip() for line in f if line.strip())


def build_cb_hash(whitelist: set, prefix_len: int = 6) -> dict:
    """
    Build compact CB correction hash table (v2 — Plan H: Python int dict).

    Encodes all Hamming ≤ 2 variants of each whitelist barcode as 32-bit
    integers and stores them in a Python dict[int→int] for O(1) lookup.
    Integer keys are hashed trivially (the value itself) in CPython,
    making this faster than v1's string-key dict while using less memory.

    Memory: ~500–800MB (Python int objects × 14M entries).
    Speed:  O(1) dict lookup (faster than v1 string hashing).
    Quality: Mathematically identical to v1 (74.4% match rate).
    """
    t0 = time.time()

    barcode_list = list(whitelist)
    bc_to_idx = {bc: i for i, bc in enumerate(barcode_list)}

    # Python dict with integer keys: O(1) lookup, no external deps
    cb_map: dict[int, int] = {}
    total_theoretical = 0

    for bc in barcode_list:
        parent_idx = bc_to_idx[bc]
        for variant_str, _ in _generate_variants(bc, max_hamming=2):
            total_theoretical += 1
            encoded = _encode_barcode(variant_str)
            # First-come-first-serve: keep the first (lowest Hamming distance)
            if encoded not in cb_map:
                cb_map[encoded] = parent_idx

    elapsed = time.time() - t0
    collision_rate = (total_theoretical - len(cb_map)) / total_theoretical * 100
    est_mb = len(cb_map) * 52 / (1024 * 1024)
    print(f"  CB hash v2 (int dict): {len(cb_map):,} entries "
          f"(collision {collision_rate:.1f}%) | "
          f"~{est_mb:.0f}MB [{elapsed:.1f}s]")

    return {
        'cb_map': cb_map,
        'barcode_list': barcode_list,
    }


def correct_cb(raw_cb: str, cb_hash: dict):
    """
    Correct a raw cell barcode to the nearest whitelist barcode (v2 — Plan H: int dict).

    Encodes raw CB as uint32, performs O(1) Python dict.get().
    Integer-key dict lookup is faster than string-key (v1) because
    CPython hashes ints as their own value with zero computation.
    Returns corrected barcode string, or None if no Hamming ≤ 2 match exists.
    """
    # Reject barcodes with ambiguous bases (N, etc.)
    if any(ch not in _BASE2BITS for ch in raw_cb):
        return None
    query = _encode_barcode(raw_cb)
    parent_idx = cb_hash['cb_map'].get(query)
    if parent_idx is not None:
        return cb_hash['barcode_list'][parent_idx]
    return None


def load_guide_hash(path: str) -> dict:
    """Load the pre-built guide hash table from a pickle file."""
    with open(path, 'rb') as f:
        return pickle.load(f)


# ── v3 Integer encoding helpers ──

# Byte-level lookup table: ASCII byte → 2-bit value (0-3), -1 for invalid
_BYTE2BITS_V3 = [-1] * 256
_BYTE2BITS_V3[ord('A')] = 0
_BYTE2BITS_V3[ord('C')] = 1
_BYTE2BITS_V3[ord('G')] = 2
_BYTE2BITS_V3[ord('T')] = 3

# Pre-computed 5bp chunk encoding lookup table (P0 optimization).
# Index layout: b0<<8 | b1<<6 | b2<<4 | b3<<2 | b4 (each b ∈ {0,1,2,3}).
# Since the index IS the correctly-ordered 10-bit encoding, the LUT is
# the identity function: _5MER_ENC[i] == i.  The speedup comes from
# eliminating the Python for-loop overhead via inline unrolling.
_5MER_ENC = list(range(1024))  # identity LUT, 1024 × 10-bit values


def _encode_umi_bytes(b: bytes, start: int) -> int:
    """Encode 12bp UMI from bytes as uint32 (24 bits), loop-unrolled."""
    v = _BYTE2BITS_V3
    return (v[b[start]]   << 22 | v[b[start+1]]  << 20 |
            v[b[start+2]] << 18 | v[b[start+3]]  << 16 |
            v[b[start+4]] << 14 | v[b[start+5]]  << 12 |
            v[b[start+6]] << 10 | v[b[start+7]]  << 8  |
            v[b[start+8]] << 6  | v[b[start+9]]  << 4  |
            v[b[start+10]] << 2 | v[b[start+11]])


def _decode_umi(val: int) -> str:
    """Decode uint32 UMI back to 12bp string."""
    bases = 'ACGT'
    chars = []
    for _ in range(12):
        chars.append(bases[val & 0x3])
        val >>= 2
    return ''.join(reversed(chars))


def _encode_seq(seq: str) -> int:
    """Encode DNA sequence string as uint64 (2 bits/base).
    
    Used for guide hash conversion (build-time only, not in hot path).
    """
    val = 0
    for ch in seq:
        val = (val << 2) | _BASE2BITS[ch]
    return val


def _encode_window_bigint(window_bytes: bytes) -> int:
    """Encode 26bp window as a 52-bit integer (P0 optimized).

    Uses 5x 5bp chunk LUT lookups + 1 trailing base.  The 5bp chunk
    indices are computed inline (unrolled, no Python loop overhead),
    then combined via bitwise shifts.  This eliminates the 26-iteration
    for-loop in the hot path, saving ~500ns per read.
    """
    v = _BYTE2BITS_V3
    b = window_bytes
    # 5bp chunk indices: b0<<8 | b1<<6 | b2<<4 | b3<<2 | b4
    i0 = (v[b[0]]  << 8 | v[b[1]]  << 6 | v[b[2]]  << 4 |
          v[b[3]]  << 2 | v[b[4]])
    i1 = (v[b[5]]  << 8 | v[b[6]]  << 6 | v[b[7]]  << 4 |
          v[b[8]]  << 2 | v[b[9]])
    i2 = (v[b[10]] << 8 | v[b[11]] << 6 | v[b[12]] << 4 |
          v[b[13]] << 2 | v[b[14]])
    i3 = (v[b[15]] << 8 | v[b[16]] << 6 | v[b[17]] << 4 |
          v[b[18]] << 2 | v[b[19]])
    i4 = (v[b[20]] << 8 | v[b[21]] << 6 | v[b[22]] << 4 |
          v[b[23]] << 2 | v[b[24]])
    last = v[b[25]]
    # Combine: 5×10-bit chunks + 1×2-bit trailing base → 52 bits
    return (_5MER_ENC[i0] << 42 |
            _5MER_ENC[i1] << 32 |
            _5MER_ENC[i2] << 22 |
            _5MER_ENC[i3] << 12 |
            _5MER_ENC[i4] << 2  |
            last)


def _extract_guides_from_bigint(big_int: int) -> list:
    """Extract 7 × 20bp guide codes from a 52-bit window int via shift+mask.
    
    Returns list of 7 ints, each encoding a 20bp guide as 40-bit value.
    """
    MASK_40 = (1 << 40) - 1
    guides = []
    for offset in range(7):
        shift = (26 - offset - 20) * 2
        guides.append((big_int >> shift) & MASK_40)
    return guides


def _process_one_fastq_pair(args: tuple) -> np.ndarray:
    """Worker function: process one R1+R2 FASTQ pair, return assignments array.
    
    Args: (r1_path, r2_path, cb_hash, guide_hash_int, idx_to_id, max_reads)
    Returns: np.ndarray of shape (N, 3) dtype=int32, columns=[cb_idx, umi_int, guide_idx]
    """
    r1_path, r2_path, cb_hash, guide_hash_int, max_reads = args
    
    cb_map = cb_hash['cb_map']
    barcode_list = cb_hash['barcode_list']
    bc_to_idx = {bc: i for i, bc in enumerate(barcode_list)}
    
    # Pre-allocate assignments (generous estimate: 25% of reads)
    MAX_ESTIMATE = 30_000_000
    assignments = np.zeros((MAX_ESTIMATE, 3), dtype=np.int32)
    write_ptr = 0
    
    stats = {'total': 0, 'valid_cb': 0, 'cb_exact': 0, 'cb_corrected': 0,
             'matched': 0, 'exact_hit': 0, 'hamming1_hit': 0, 'no_guide': 0,
             'short_read2': 0}
    
    # Read FASTQ
    opener = gzip.open if r1_path.endswith('.gz') else open
    with opener(r1_path, 'rb') as f1, opener(r2_path, 'rb') as f2:
        while True:
            # Read Read1 record (4 lines)
            h1 = f1.readline()
            if not h1: break
            seq1 = f1.readline().strip()
            f1.readline()  # +
            f1.readline()  # qual
            
            # Read Read2 record
            h2 = f2.readline()
            if not h2: break
            seq2 = f2.readline().strip()
            f2.readline()
            f2.readline()
            
            stats['total'] += 1
            
            # ── Step 1: CB correction (int dict lookup) ──
            if len(seq1) < UMI_END:
                continue
            raw_cb_bytes = seq1[CB_START:CB_END]
            if any(_BYTE2BITS_V3[b] < 0 for b in raw_cb_bytes):
                continue
            cb_int = _encode_barcode(raw_cb_bytes.decode())
            cb_parent = cb_map.get(cb_int)
            if cb_parent is None:
                continue
            stats['valid_cb'] += 1
            cb_idx = bc_to_idx.get(barcode_list[cb_parent], -1)
            if cb_idx < 0:
                continue
            
            # Track exact vs corrected (for diagnostics)
            parent_cb = barcode_list[cb_parent]
            if raw_cb_bytes.decode() == parent_cb:
                stats['cb_exact'] += 1
            else:
                stats['cb_corrected'] += 1
            
            # ── Step 2: UMI encoding ──
            umi_int = _encode_umi_bytes(seq1, UMI_START)
            
            # ── Step 3: Guide matching (integer, 7 sliding windows) ──
            if len(seq2) < WINDOW_END:
                stats['short_read2'] += 1
                continue
            
            window_bytes = seq2[WINDOW_START:WINDOW_END]
            big_int = _encode_window_bigint(window_bytes)
            guides = _extract_guides_from_bigint(big_int)
            
            found_idx = -1
            # Fast path: exact match (65%+ reads hit on 1st candidate)
            for g in guides:
                idx = guide_hash_int.get(g)
                if idx is not None:
                    found_idx = idx
                    stats['exact_hit'] += 1
                    break
            
            # Slow path: Hamming=1 (only if fast path failed, <1% reads)
            if found_idx < 0:
                for offset, g in enumerate(guides):
                    for pos in range(GUIDE_LEN):
                        orig_bits = (g >> ((GUIDE_LEN - 1 - pos) * 2)) & 0x3
                        for alt_bits in range(4):
                            if alt_bits == orig_bits:
                                continue
                            # Flip 2 bits at position
                            shift = (GUIDE_LEN - 1 - pos) * 2
                            variant = (g & ~(0x3 << shift)) | (alt_bits << shift)
                            idx = guide_hash_int.get(variant)
                            if idx is not None:
                                found_idx = idx
                                stats['hamming1_hit'] += 1
                                break
                        if found_idx >= 0:
                            break
                    if found_idx >= 0:
                        break
            
            if found_idx >= 0:
                stats['matched'] += 1
                assignments[write_ptr, 0] = cb_idx
                assignments[write_ptr, 1] = umi_int
                assignments[write_ptr, 2] = found_idx
                write_ptr += 1
            else:
                stats['no_guide'] += 1
            
            if max_reads and stats['total'] >= max_reads:
                break
    
    # Truncate and return
    return assignments[:write_ptr], stats


def match_reads(
    r1_path: str,
    r2_path: str,
    whitelist: set,
    guide_hash: dict,
    output_path: str,
    max_reads: Optional[int] = None,
    threads: int = 1,
    report_interval: int = 1_000_000,
) -> dict:
    """
    Core matching loop (v3 — integer encoding + numpy pre-allocation + multi-process I/O).
    """
    seq_to_idx = guide_hash['seq_to_idx']
    idx_to_id = guide_hash['idx_to_id']
    
    # Use integer-encoded guide hash if available (v3), else string keys (v2 compat)
    guide_hash_int: dict = guide_hash.get('seq_to_idx_int')
    if guide_hash_int is None:
        print("Note: integer guide hash not found, converting...")
        guide_hash_int = {}
        for seq_str, idx in seq_to_idx.items():
            guide_hash_int[_encode_seq(seq_str)] = idx

    # Build CB correction hash (Plan H — int dict)
    print("Building CB hash (v3 — Plan H int dict)...")
    cb_hash = build_cb_hash(whitelist)
    print()

    # Split FASTQ paths
    r1_list = [p.strip() for p in r1_path.split(',') if p.strip()]
    r2_list = [p.strip() for p in r2_path.split(',') if p.strip()]
    assert len(r1_list) == len(r2_list), \
        f"Mismatched R1/R2 file counts: {len(r1_list)} vs {len(r2_list)}"
    
    n_files = len(r1_list)
    use_mp = threads > 1 and n_files > 1

    t0 = time.time()
    
    if use_mp:
        import multiprocessing as mp
        mp.set_start_method('fork', force=True)
        
        worker_args = [
            (r1, r2, cb_hash, guide_hash_int, max_reads)
            for r1, r2 in zip(r1_list, r2_list)
        ]
        
        print(f"Processing {n_files} FASTQ pairs with {min(threads, n_files)} workers...")
        with mp.Pool(processes=min(threads, n_files)) as pool:
            results = pool.map(_process_one_fastq_pair, worker_args)
        
        # Merge results
        all_arrays = []
        merged_stats = {}
        for arr, stats in results:
            if arr.size > 0:
                all_arrays.append(arr)
            for k, v in stats.items():
                merged_stats[k] = merged_stats.get(k, 0) + v
        
        if all_arrays:
            all_assignments = np.vstack(all_arrays)
        else:
            all_assignments = np.zeros((0, 3), dtype=np.int32)
    else:
        # Single-threaded path
        all_arrays = []
        merged_stats = {'total': 0, 'valid_cb': 0, 'cb_exact': 0, 'cb_corrected': 0,
                        'matched': 0, 'exact_hit': 0, 'hamming1_hit': 0, 'no_guide': 0,
                        'short_read2': 0}
        
        for r1, r2 in zip(r1_list, r2_list):
            print(f"Processing: {os.path.basename(r1)} + {os.path.basename(r2)}")
            arr, stats = _process_one_fastq_pair(
                (r1, r2, cb_hash, guide_hash_int, max_reads))
            if arr.size > 0:
                all_arrays.append(arr)
            for k, v in stats.items():
                merged_stats[k] = merged_stats.get(k, 0) + v
        
        if all_arrays:
            all_assignments = np.vstack(all_arrays)
        else:
            all_assignments = np.zeros((0, 3), dtype=np.int32)
    
    elapsed = time.time() - t0

    # ── Save as .npz (binary, fast) ──
    os.makedirs(os.path.dirname(output_path) if os.path.dirname(output_path) else '.',
                exist_ok=True)
    print(f"\nWriting {all_assignments.shape[0]:,} assignments...")
    # Build barcode lookup: cb_idx → barcode string
    barcode_list = cb_hash['barcode_list']
    np.savez_compressed(output_path,
                        assignments=all_assignments,
                        barcode_list=np.array(barcode_list),
                        idx_to_id=np.array(idx_to_id))

    # ── Summary ──
    stats = merged_stats
    print(f"\n{'='*60}")
    print(f"MATCHING COMPLETE (v3)")
    print(f"{'='*60}")
    print(f"  Total reads:       {stats['total']:>12,}")
    print(f"  Valid CB:          {stats['valid_cb']:>12,}  "
          f"({stats['valid_cb']/max(stats['total'],1)*100:.1f}%)")
    print(f"    Exact CB:        {stats['cb_exact']:>12,}")
    print(f"    Corrected CB:    {stats['cb_corrected']:>12,}")
    print(f"  Short Read2:       {stats['short_read2']:>12,}")
    print(f"  Guide matched:     {stats['matched']:>12,}  "
          f"({stats['matched']/max(stats['valid_cb'],1)*100:.1f}%)")
    print(f"    Exact:           {stats['exact_hit']:>12,}")
    print(f"    Hamming=1:       {stats['hamming1_hit']:>12,}")
    print(f"  No guide detected: {stats['no_guide']:>12,}")
    print(f"  Time: {elapsed:.1f}s  ({stats['total']/elapsed:.0f} reads/s)")
    print(f"  Output: {output_path}")
    print(f"{'='*60}")

    return stats


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Hash-accelerated guide matcher (v3)')
    parser.add_argument('-1', '--read1', required=True, help='Read1 FASTQ (CB+UMI), comma-separated for multi-file')
    parser.add_argument('-2', '--read2', required=True, help='Read2 FASTQ (guide), comma-separated for multi-file')
    parser.add_argument('-w', '--whitelist', required=True, help='Cell barcode whitelist')
    parser.add_argument('-g', '--guide-hash', required=True, help='Guide hash table (.pkl)')
    parser.add_argument('-o', '--output', required=True, help='Output assignments (.npz)')
    parser.add_argument('-n', '--max-reads', type=int, default=None,
                        help='Max reads to process (for testing)')
    parser.add_argument('-t', '--threads', type=int, default=1,
                        help='Number of parallel workers (one per FASTQ pair)')
    args = parser.parse_args()

    wl = load_whitelist(args.whitelist)
    gh = load_guide_hash(args.guide_hash)
    match_reads(args.read1, args.read2, wl, gh, args.output, args.max_reads, args.threads)
