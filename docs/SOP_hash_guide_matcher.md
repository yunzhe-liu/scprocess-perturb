# SOP: Hash-Accelerated Guide Matcher for Perturb-seq sgRNA Quantification

**Version**: 3.0  
**Date**: May 25, 2026  
**Status**: v3 — Integer encoding (P0), numpy pre-allocation (P1), multi-process I/O (P3).  
Validated on K562 Essential CRISPR Library (4,536 guides, 48 lanes, 625,911 cells)

---

## v3 Changelog (vs v2)

| Optimization | Description | Expected Gain |
|---|---|---|
| **P0 – Integer encoding** | 20bp guide → uint64, 26bp window → 52-bit int. Zero per-read str malloc. Guide hash keys are int, CPython hashes ints trivially. | 15% per-read |
| **P1 – Numpy pre-allocation** | 30M×3 int32 buffer replaces list-of-tuples. Output as .npz (binary) instead of .tsv.gz (text). Eliminates per-assignment malloc + TSV overhead. | 40% total |
| **P3 – Multi-process I/O** | One worker per FASTQ pair (4 per lane), fork-based COW memory. | Linear with threads |

**v3 per-lane wall time**: ~7-9 min (4 threads) vs v2's ~24 min → **~3× speedup**.  
**48-lane 48-core**: ~35 min vs v2's ~2.5h.  
**48-lane 16-core**: ~1.5h vs v2's ~6h.  

**Data format change**: Intermediate output changed from `assignments.tsv.gz` to `assignments.npz`. MEX output unchanged. Guide hash includes `seq_to_idx_int` alongside `seq_to_idx`.

---

## 1. Background and Diagnostic Foundations

### 1.1 The Problem: simpleaf Cell Recovery Failure

The `guide_extraction` workflow originally used **simpleaf** (piscem + alevin-fry v0.14.0)
for sgRNA guide quantification from Perturb-seq data. While simpleaf achieved a
reasonable mapping rate (66.9%–71.2%), the **cell recovery rate against the public
standard reference** (Cell Ranger-processed count matrix) stagnated at ~50%.
Diagnostic investigation revealed:

| Tool | Mapping Rate | Whitelist Coverage | Median UMI/Cell | Pearson r vs fba |
|------|:---:|:---:|:---:|:---:|
| simpleaf (k=13,m=9) | 71.2% | 45.0% | ~4 | 0.05 |
| **fba** (Levenshtein) | — | 68.5% (100K reads) | 4.0 | — |

The critical finding: **simpleaf's high mapping rate was driven by false-positive
k-mer matches from adapter sequences, not genuine guide detection.** The Pearson
correlation of r=0.05 between simpleaf and fba on the same cells confirmed that
the two tools were assigning UMIs to fundamentally different guides.

### 1.2 Root Cause: piscem k-mer Pseudoalignment Failure Mode

piscem (k=13, m=9) decomposes each 98bp Read2 into 86 overlapping 13-mers.
Only ~8 of these (positions 30–50) originate from the 20bp guide sequence.
The remaining ~78 13-mers come from the TSO adapter (positions 1–30) and
poly-A/scaffold tail (positions 51–98).

With 4,536 guide sequences and a 13-mer space of 4^13 ≈ 67 million, random
13-mer collisions between adapter-derived k-mers and guide index k-mers are
frequent. piscem's intersection-based scoring is thus dominated by noise,
causing **systematic guide misassignment**.

### 1.3 What fba Got Right

fba (Feature Barcoding Assays) demonstrated that **restricting the search to
the exact guide window** (Read2 positions 31–51, 20bp) and using Levenshtein
distance matching eliminated adapter noise entirely. fba achieved 68.5%
whitelist coverage with only 100K reads (1.9% of a full lane), compared to
simpleaf's 45.0% with 5.2M reads.

However, fba's BK-tree-based Levenshtein matching against 4,536 guides proved
impractically slow: >22 minutes per lane vs simpleaf's 8 seconds.

### 1.4 Design Principle

The Hash-Accelerated Guide Matcher combines:
- **fba's window-restricted search** (eliminating adapter noise)
- **Pre-computed hash tables** (O(1) lookup replacing BK-tree O(log N))
- **Hamming-distance error tolerance** (covering Illumina substitution-dominant errors)

---

## 2. Algorithm Description

### 2.1 Overview

```
┌──────────────────────────────────────────────────────────────────┐
│                    PREPROCESSING (once)                           │
│                                                                  │
│  build_guide_hash.py                                             │
│    Input:  guides.fasta (4,536 × 20bp)                           │
│    Output: guide_hash.pkl (~7 MB)                                │
│                                                                  │
│    For each guide sequence:                                      │
│      a. Store exact 20bp → guide_id mapping                      │
│      b. Generate 60 Hamming-1 neighbors → guide_id mapping       │
│    Total: ~276K entries                                          │
└──────────────────────────────────────────────────────────────────┘
                              ↓
┌──────────────────────────────────────────────────────────────────┐
│                    MATCHING (per lane)                            │
│                                                                  │
│  guide_matcher.py                                                │
│    Input:  Read1 FASTQ (CB+UMI), Read2 FASTQ (guide construct)   │
│            Cell barcode whitelist, guide_hash.pkl                 │
│    Output: assignments.tsv.gz (CB, UMI, guide_id)                │
│                                                                  │
│    For each read pair:                                           │
│      1. Extract raw CB = Read1[0:16]                             │
│      2. Correct CB to whitelist (Hamming ≤ 2)                    │
│      3. Extract UMI = Read1[16:28]                               │
│      4. Extract window = Read2[28:54] (26bp)                     │
│      5. Slide 20bp sub-window ×7 positions                       │
│      6. O(1) hash lookup → guide_id                              │
│      7. Fast path: exact match (7 lookups)                       │
│      8. Slow path: Hamming-1 variants (420 lookups)              │
└──────────────────────────────────────────────────────────────────┘
                              ↓
┌──────────────────────────────────────────────────────────────────┐
│                    UMI DEDUP + MATRIX (per lane)                  │
│                                                                  │
│  umi_dedup_matrix.py                                             │
│    Input:  assignments.tsv.gz                                    │
│    Output: MEX trio (matrix.mtx.gz, barcodes.tsv.gz,             │
│                      features.tsv.gz)                            │
│                                                                  │
│    1. Group assignments by (CB, guide_id)                        │
│    2. UMI-tools directional dedup (Hamming ≤ 1)                  │
│    3. Build sparse CSR matrix (cells × features)                 │
│    4. Write standard 10x MEX format                              │
└──────────────────────────────────────────────────────────────────┘
```

### 2.2 Key Parameters

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| **Read2 window** | `[28:54)` (26bp) | Covers guide start positions 28–34 (observed range from brute-force search on 50K reads) |
| **Guide length** | 20bp | All 4,536 guides in K562 essential library are exactly 20bp |
| **Guide match tolerance** | Hamming ≤ 1 | Illumina errors are >99% substitutions; Hamming=1 covers 60 variants per guide |
| **CB match tolerance** | Hamming ≤ 2 via prefix bucketing | L1: ~614K H1 hash (12MB) + L2: 4,096 prefix buckets (<200KB). Total <20MB (v2) |
| **CB prefix length** | 6 bp | Controls bucket granularity; 4^6=4096 buckets, avg 3 barcodes/bucket |
| **CB window** | `[0:16)` | 10x Chromium V3 chemistry: 16bp cell barcode |
| **UMI window** | `[16:28)` | 10x Chromium V3 chemistry: 12bp UMI |
| **UMI dedup** | Hash-accelerated directional (v2) | O(n×36) replacing (12bp UMI × 3 alt bases) O(n×m); mathematically identical to UMI-tools |
| **Guide hash size** | 276,450 entries (~7 MB) | 4,536 exact + 271,918 Hamming-1 variants |

### 2.3 Two-Tier Matching Strategy

**Fast path (exact match)** — covers ~99.6% of guide-matched reads:
- 7 hash lookups per read (sliding window positions)
- O(7) expected time per read
- Dominates execution: 2,838,726 / 2,849,915 = 99.6% of assignments

**Slow path (Hamming-1 fallback)** — covers ~0.4% of assignments:
- 7 × 20 × 3 = 420 hash lookups per read
- Only triggered when no exact match found
- Adds <1% overhead to total runtime

### 2.4 Cell Barcode Error Correction (v2: Prefix Bucketing)

Raw cell barcodes extracted from Read1[0:16] exhibit significant divergence
from the GEX-derived whitelist. Only 0.5% of reads have exact CB matches.

**v2 Design — Two-Level Prefix Bucketing (generalizable, no library-specific assumptions):**

```
Level 1 (Fast Cache): Exact + Hamming=1 hash
  → ~614K entries, ~12MB, O(1) lookup
  → Covers ~0.9% of reads

Level 2 (Prefix Buckets): First 6bp index → 4,096 buckets
  → Average 3 barcodes per bucket
  → Check self-bucket + 18 neighbor buckets (1bp prefix mutation)
  → Max ~60 candidates → real-time Hamming ≤ 2 verification
  → Covers ~73.5% of reads
```

This achieves 74.4% valid CB rate, while using **<20MB** memory (vs 640MB in v1).
The prefix-bucket design is fully general: it makes no assumptions about
specific mismatch positions or library preparation artifacts.

---

## 3. Performance Characteristics

### 3.1 Runtime (lane_01, 5.2M reads, 1 core, v2 estimated)

| Stage | v1 Time | v2 Time | Improvement |
|-------|:---:|:---:|:---:|
| Build guide hash | 0.2s | 0.2s | — |
| Build CB hash | 9.5s | **0.3s** | 32× |
| Match reads | 191.6s | **~90s** | 2× (faster CB lookup) |
| UMI dedup + matrix | 313.9s | **~30s** | 10× (hash-accelerated) |
| **Total (single lane)** | ~8.5 min | **~2 min** | 4× |
| **48 lanes (parallel)** | ~7h (limited by lane_45) | **~20 min** | 21× |

### 3.2 Memory

| Component | v1 | v2 |
|-----------|:---:|:---:|
| Guide hash | ~7 MB | ~7 MB |
| CB hash | 640 MB | **<20 MB** |
| Peak (matching) | ~800 MB | **~150 MB** |
| Peak (UMI dedup) | ~200 MB | ~200 MB |

### 3.3 Accuracy (vs simpleaf, same lane)

| Metric | simpleaf | Hash Matcher | Improvement |
|--------|:---:|:---:|:---:|
| Whitelist coverage | 45.0% | **98.4%** | 2.2× |
| Cells detected | 5,674 | **12,365** | 2.2× |
| Median UMI/cell | ~4 | **97.0** | 24× |
| Total UMIs | 94,062 | **2,142,753** | 22.8× |

*Note: simpleaf used 4 Illumina lanes; Hash Matcher used only 1 lane (L001).
With all 4 lanes, the Hash Matcher's advantage would be even larger.*

---

## 4. Generality and Applicability

### 4.1 When This Method is Appropriate

The Hash-Accelerated Guide Matcher is designed for **Feature Barcoding assays
where short, known sequences are embedded within longer sequencing reads
flanked by constant adapter sequences.** This includes:

- CRISPR sgRNA guide quantification (Perturb-seq, CROP-seq)
- Antibody-derived tag (ADT) quantification (CITE-seq, TotalSeq)
- Hashtag oligonucleotide (HTO) demultiplexing
- Any custom Feature Barcoding experiment with known barcode sequences

### 4.2 Adapting to New Libraries

To adapt the method for a new Feature Barcode library:

1. **Replace `guides.fasta`** with the new Feature Barcode sequences
2. **Adjust `WINDOW_START` and `WINDOW_END`** based on the barcode position in Read2
   (determined via brute-force search or `fba qc`)
3. **Adjust `GUIDE_LEN`** if barcodes have non-uniform lengths
4. **Update the whitelist** with GEX-validated cell barcodes

### 4.3 Limitations

- **Hamming-only error model**: Does not handle insertion/deletion errors.
  For platforms with significant indel rates (e.g., PacBio, Nanopore), the
  Levenshtein-based fba approach is more appropriate.
- **Fixed barcode length**: All Feature Barcodes must have the same length.
  Variable-length barcodes require extending the hash table or using the
  slower BK-tree fallback.
- **Memory for large whitelists**: The CB hash grows as O(N × K²) where N is
  the number of whitelist barcodes and K is the barcode length. For >100K
  barcodes, consider reducing max Hamming distance or using prefix-bucket
  indexing.

---

## 5. Relationship to Existing Tools

### 5.1 Inheritance from simpleaf

| Component | simpleaf | Hash Matcher | Inherited? |
|-----------|----------|-------------|:---:|
| Reference preparation | piscem index (k=13,m=9) | Hash table (exact + H1) | ✗ (simpler) |
| Read mapping | k-mer pseudoalignment | Windowed hash lookup | ✗ (different paradigm) |
| Whitelist handling | alevin-fry `--explicit-pl` | Hamming ≤ 2 hash correction | ✗ (tolerates errors) |
| UMI dedup | alevin-fry parsimony-gene | UMI-tools directional | ✗ (different algorithm) |
| Output format | MEX via alevin-fry | MEX via scipy | ✓ (same format) |
| Per-lane parallelism | Snakemake per-group | Snakemake per-group | ✓ (same architecture) |

### 5.2 Inheritance from fba

| Component | fba | Hash Matcher | Inherited? |
|-----------|-----|-------------|:---:|
| Window-restricted search | Read2[31:51] fixed | Read2[28:54] sliding | ✓ (extended) |
| Error model | Levenshtein (ins+del+sub) | Hamming (sub only) | ✗ (simplified) |
| CB matching | Levenshtein ≤ 2 | Hamming ≤ 2 hash | ✓ (same tolerance) |
| UMI dedup | UMI-tools directional | UMI-tools directional | ✓ (identical) |
| Speed optimization | BK-tree (O(log N)) | Hash table (O(1)) | ✗ (key innovation) |
| QC module | `fba qc` diagnostics | N/A (relies on fba qc for diagnostics) | ✗ (delegated) |

### 5.3 Algorithmic Lineage

```
Cell Ranger (10x Genomics)
  └─ Feature Barcode regex matching with (BC)/(5P(BC)) patterns
     └─ Public standard reference (our ground truth)

fba (Feature Barcoding Assays)
  └─ Window-restricted Levenshtein matching
     └─ Key insight: window restriction eliminates adapter noise
     └─ Inherited: UMI-tools directional dedup, CB error correction principle
     └─ Limitation: BK-tree O(log N) too slow for 4,536 guides

Hash-Accelerated Guide Matcher (this work)
  └─ Window-restricted exact + Hamming-1 hash matching
     └─ Extends fba's window approach with sliding window
     └─ Replaces BK-tree with O(1) pre-computed hash tables
     └─ Simplifies error model from Levenshtein to Hamming
        (justified by Illumina error profile: >99% substitutions)
     └─ Key innovation: pre-computation of all Hamming-1 guide variants
        and Hamming ≤ 2 CB variants enables O(1) lookup
```

---

## 6. Implementation Files

| File | Purpose |
|------|---------|
| `scripts/build_guide_hash.py` | One-time guide hash table construction |
| `scripts/guide_matcher.py` | Main read matching engine |
| `scripts/umi_dedup_matrix.py` | UMI dedup + MEX matrix generation |
| `rules/guide_quant.smk` | Snakemake rule for per-lane quantification |
| `config/config.yaml` | Configuration (method selection, parameters) |

---

## 7. Usage

### 7.1 Command Line (Single Lane)

```bash
# Step 1: Build guide hash (one-time, v3 includes int keys)
python scripts/build_guide_hash.py \
    /data/yunzliu/references/guides.fasta \
    /data/yunzliu/references/guide_hash_v3.pkl

# Step 2: Match reads (v3 — .npz binary output, multi-threaded)
python scripts/guide_matcher.py \
    -1 "lane_01_L001_R1.fastq.gz,lane_01_L002_R1.fastq.gz,lane_01_L003_R1.fastq.gz,lane_01_L004_R1.fastq.gz" \
    -2 "lane_01_L001_R2.fastq.gz,lane_01_L002_R2.fastq.gz,lane_01_L003_R2.fastq.gz,lane_01_L004_R2.fastq.gz" \
    -w barcode_whitelist_noheader.txt \
    -g /data/yunzliu/references/guide_hash_v3.pkl \
    -o assignments.npz \
    -t 4

# Step 3: UMI dedup + MEX (v3 — reads .npz)
python scripts/umi_dedup_matrix.py \
    -i assignments.npz \
    -o output_matrix/ \
    -t 1
```

### 7.2 Via Snakemake

```bash
# Set method in config/config.yaml:
#   guide_extraction:
#     method: hash_matcher  # or "simpleaf"

snakemake --cores 8
```

---

## 8. Validation Results

### 8.1 Single-Lane Test (lane_01, 1 Illumina sub-lane L001)

**Dataset**: K562 Essential CRISPR Library (Replogle et al., 2022)  
**Reference**: `K562_essential_raw_singlecell_01.h5ad` (310,385 cells, Cell Ranger-processed)

| Metric | Value |
|--------|-------|
| Whitelist coverage | 12,365 / 12,561 = **98.4%** |
| Guide-matched reads | 2,849,915 / 3,881,776 = **73.4%** |
| Total UMIs after dedup | 2,142,753 |
| Median UMI per cell | 97.0 |
| Unique guides detected | 4,480 / 4,536 |

### 8.2 48-Lane Full Run

| Metric | Value |
|--------|-------|
| Merged cells | 668,632 (with lane suffix) |
| Merged features | 4,532 / 4,536 |
| Total UMIs | 353,143,093 |
| Median UMI/cell | 253.0 |
| Mean UMI/cell | 528.2 |
| Per-lane cells (mean ± std) | 13,930 ± 4,226 |
| Per-lane UMIs (mean ± std) | 7.0M ± 2.2M |
