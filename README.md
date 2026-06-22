# sgprocess — Single-Cell Guide Extraction Workflow

**A Snakemake workflow for extracting sgRNA guide counts from Perturb-seq data.**

Integrates two quantification methods (simpleaf, HAM) under a common
configuration layer, and provides post-quantification merge and barcode
translation steps shared by both methods.

---

## Workflow Overview

```
                    ┌─────────────────────────────────────────┐
GEX matrix (H5) ───→ extract whitelist                        │
Guide FASTA ───────→ build index / hash table                │
sgRNA FASTQ ───────┤                                         │
                   │    ┌─ simpleaf: piscem + alevin-fry     │
                   ├──→ │─ HAM: hash match + dedup          ├──→ merge ──→ MEX
                   │    └───────────────────────────────────┘   + barcode
                   └─────────────────────────────────────────┘   translation
```

**Core design principle:** GEX mRNA signal (produced externally by alevin-fry-based
quantification) defines which droplets contain real cells. A whitelist of cell
barcodes passing GEX QC thresholds (minimum UMI, minimum genes) is extracted,
**translated from GEX format to Feature format** (see Barcode Translation below)
to match the sgRNA FASTQ barcodes, and fed to the guide quantification tools.
After per-lane quantification and cross-lane merge, barcodes are **translated
back from Feature format to GEX format** so the final output matches the mRNA
reference.

The workflow produces **standard MEX-formatted count matrices** (gzip-compressed
Market Exchange Format: `matrix.mtx.gz`, `barcodes.tsv.gz`, `features.tsv.gz`),
compatible with downstream tools (Seurat, Scanpy, scvi-tools).

---

## Installation & Prerequisites

### Conda environments

Two conda environments are required:

| Environment | Contains | Purpose |
|-------------|----------|---------|
| `scprocess` | Snakemake, HAM, numpy, scipy, h5py | Workflow orchestration, HAM quant, merge, analysis |
| `simpleaf` | simpleaf ≥ 0.24, piscem ≥ 0.19, alevin-fry ≥ 0.14 | simpleaf quant |

```bash
conda activate scprocess    # for Snakemake + HAM + merge
conda activate simpleaf     # for simpleaf quant
```

### Dependencies

- **Snakemake** ≥ 8.0 (workflow engine)
- **simpleaf** (piscem + alevin-fry) for k-mer pseudoalignment
- **HAM** (pip-installable from `/path/to/ham`) for hash-based matching
- **h5py, numpy, scipy** for whitelist extraction, matrix I/O, merge

### Workflow installation

```bash
git clone <sgprocess-url> /path/to/sgprocess
cd /path/to/sgprocess
conda env create -f envs/scprocess.lock.yaml   # HAM + Snakemake + analysis
conda env create -f envs/simpleaf.lock.yaml    # simpleaf only
```

---

## Repository Structure

```
sgprocess/
├── Snakefile                     ← Snakemake entry point
├── config/
│   ├── config.yaml               ← Default (ham)
│   ├── config_<variant>.yaml     ← Per-method/per-sweep configs
│   └── samples.yaml              ← Lane topology (48 groups)
├── rules/
│   ├── reference.smk             ← Guide FASTA → index / hash table
│   ├── whitelist.smk             ← GEX H5 → barcode whitelist
│   ├── quant.smk                 ← simpleaf quant
│   ├── guide_quant.smk           ← HAM quant (match + dedup)
│   └── merge.smk                 ← Per-lane merge + barcode translation
├── scripts/
│   ├── translate_barcodes.py     ← Feature ↔ GEX barcode translation
│   ├── build_guide_hash.py       ← HAM hash table builder (used by reference.smk)
│   ├── generate_per_lane_whitelists.py ← Whitelist extraction (used by whitelist.smk)
│   ├── filter_barcodes.py        ← Barcode QC thresholds
│   ├── feature_reference_adapter.py   ← Guide FASTA / t2g generation
│   ├── merge_count_matrices.py   ← Legacy merge (deprecated; ham merge supersedes)
│   ├── analyze_outputs.py        ← Post-run summary statistics
│   └── barcode_rank_qc.py        ← Knee-plot diagnostics
├── envs/
│   ├── simpleaf.lock.yaml        ← conda lock: simpleaf
│   └── scp_analysis.lock.yaml    ← conda lock: Python + Snakemake + HAM
└── profiles/
    └── local/                    ← Snakemake profile (local execution)
```

---

## Configuration

### `samples.yaml` — Lane Topology

Defines one group per physical 10x lane. Each group specifies:

```yaml
groups:
  lane_01:
    group_id: lane_01
    gex_h5: /path/to/scprocess/decontx_lane01.h5
    sgRNA_fastq_dir: /path/to/sgRNA_fastq/lane_01
    sgRNA_r1_pattern: "lane01_sgRNA_*_R1_001.fastq.gz"
    sgRNA_r2_pattern: "lane01_sgRNA_*_R2_001.fastq.gz"
```

### `config.yaml` — Global Settings

```yaml
guide_extraction:
  method: ham                     # "ham" | "simpleaf"

simpleaf:
  quant:
    chemistry: 10xv3              # 10x 3' v3 chemistry
    resolution: parsimony-gene    # UMI resolution mode
    use_knee: false               # true = --knee, false = --explicit-pl (whitelist)

ham:
  quant:
    cb_max_hamming: 1             # Hamming distance for CB correction

whitelist:
  min_umi: 1000                   # GEX UMI threshold
  min_genes: 500                  # GEX gene threshold

resources:
  threads: 12                     # Threads per lane
```

**Method switching:** Change `guide_extraction.method` to `simpleaf` or `ham`.

**Separate config files** are provided for each variant: k-mer sweeps
(`config_kmer_check_k15.yaml`, `config_kmer_check_k17.yaml`), knee mode
(`config_simpleaf_knee.yaml`), and minimizer sweeps.

---

## Snakemake Rule Sequence (Algorithm Order)

The workflow executes five Snakemake rules in dependency order:

### 1. `reference.smk` — Index / Hash Table Construction

**Input:** Guide FASTA file (one sequence per guide, 20 bp each).

**Output:** Method-specific lookup structure:
- **simpleaf:** piscem dense index (`simpleaf index --kmer-length K --minimizer-length M --ref-seq guides.fasta`)
- **HAM:** Python pickle hash table (`ham build-hash --fasta guides.fasta`)

This is a one-time step — the index is reused across all 48 lanes and all
parameter sweeps sharing the same k-mer and minimizer values.

### 2. `whitelist.smk` — GEX Whitelist Extraction

**Input:** Per-lane GEX count matrix (HDF5 from external quantification).

**Output:** One text file per lane containing cell barcodes that pass QC
thresholds (configured as `min_umi` and `min_genes`).

This ensures guide quantification is restricted to GEX-validated cells,
eliminating empty droplets and debris. The extracted whitelist is
automatically translated from GEX format to Feature format (TO→FROM)
so it matches the sgRNA FASTQ barcodes.

### 3a. `quant.smk` — simpleaf Quantification (per lane)

**Input:** sgRNA FASTQ, piscem index, whitelist.

**Process:**
1. `piscem map-sc` — k-mer pseudoalignment of R2 reads to guide index
   (`--geometry chromium_v3`, `--skipping-strategy permissive`)
2. `alevin-fry generate-permit-list` — restrict to whitelist barcodes
   (`--unfiltered-pl`, `--min-reads 1`)
3. `alevin-fry collate` — group mapped reads by cell barcode
4. `alevin-fry quant` — UMI deduplication and count matrix generation
   (`-r parsimony-gene`, `-t 12`)

Two cell-calling modes are supported: whitelist (`--explicit-pl`, uses
GEX-derived barcode list with no additional filtering) and knee
(`--knee`, automatic knee-detection on guide data alone).

### 3b. `guide_quant.smk` — HAM Quantification (per lane)

**Input:** sgRNA FASTQ, guide hash table, whitelist.

**Process:**
1. `ham match` — window-restricted hash lookup (positions 28–54 of R2)
   with integer-encoded guide sequences. Barcode correction at
   `--cb-max-hamming 1`. GEX whitelist restricts processing.
2. `ham dedup` — UMI-tools directional deduplication at Hamming
   distance 1, producing a per-lane MEX count matrix.

### 4. `merge.smk` — Cross-Lane Merge & Barcode Translation

**Input:** 48 per-lane MEX count matrices.

**Process:**
1. **Vertical merge:** All lane matrices are vertically stacked
   (`ham merge --lanes lane_list.tsv`). Lane suffixes (`-L{NN}`)
   are appended to barcodes for cross-lane uniqueness. Features
   must be identical across lanes.
2. **Barcode translation:** The 10x 3' v3 chemistry uses two
   oligo sets — TruSeq (GEX) and Nextera (Feature Barcoding).
   sgRNA FASTQ files carry Feature-format barcodes which differ
   from GEX-format barcodes at positions 7–8. The merge step
   translates Feature → GEX via the official 10x lookup table
   (`3M-february-2018.txt.gz`), ensuring output barcodes match
   the mRNA reference.

**Output:** A single merged MEX trio with GEX-format, lane-suffixed barcodes.
The merge step translates barcodes from Feature format back to GEX format
(FROM→TO), so the final output is directly comparable to the mRNA reference.

---

## Usage

### Basic Invocation

```bash
cd /path/to/sgprocess
conda activate scprocess

# HAM mode (default config)
snakemake --configfile config/config_final.yaml --cores 48

# simpleaf whitelist mode
snakemake --configfile config/config_simpleaf_final.yaml --cores 48

# simpleaf knee mode
snakemake --configfile config/config_simpleaf_knee.yaml --cores 48

# Dry-run (preview DAG)
snakemake --configfile config/config_final.yaml --cores 48 --dry-run
```

### Common Options

```
--dry-run           Preview execution plan without running
--cores N           Max parallel threads
--unlock            Remove stale locks after a crash
--rerun-incomplete  Re-run partially completed jobs
--latency-wait 60   Wait for filesystem latency (NFS-safe)
```

### Per-Lane Parallelism

The workflow is designed for 4 concurrent lanes × 12 threads = 48 cores total.
Concurrency is controlled by Snakemake's resource management (`resources:
threads=12` in config). No manual lane batching is needed — Snakemake
automatically schedules up to `--cores / threads_per_lane` lanes simultaneously.

### Output

```
{out_dir}/
├── lane_XX/
│   └── {quant_dir}/af_quant/alevin/
│       ├── quants_mat.mtx           ← per-lane sparse matrix
│       ├── quants_mat_rows.txt      ← cell barcodes (pre-merge)
│       └── quants_mat_cols.txt      ← guide feature IDs
└── merged/
    ├── merged_matrix.mtx.gz         ← cell × guide matrix
    ├── merged_barcodes.tsv.gz       ← barcodes (GEX format, -L{NN})
    └── merged_features.tsv.gz       ← guide feature IDs
```

---

## Barcode Translation

10x 3' v3 beads carry two oligo variants:
- **TruSeq Read 1** — captures poly(A) mRNA (Gene Expression)
- **Nextera Read 1** — captures Feature Barcoding libraries (CRISPR, cell surface proteins)

These differ by a complementary base swap at positions 7–8 of the 16 bp barcode.
The sgRNA library uses Nextera (Feature) chemistry, producing Feature-format
barcodes in the FASTQ files. The scprocess GEX pipeline uses TruSeq chemistry,
producing GEX-format barcodes. Translation between the two formats uses the
official 10x lookup table (`3M-february-2018.txt.gz`) distributed with 10x Genomics
software:

The workflow applies translation in two directions:
1. **FROM → TO:** Merge step — convert merged Feature-format barcodes to GEX format
   so output matches mRNA data.
2. **TO → FROM:** Whitelist step — convert GEX-format whitelist to Feature format
   so it matches sgRNA FASTQ barcodes.

Translation is performed by `scripts/translate_barcodes.py` using the
`--direction` flag (`from_to` or `to_from`). Translated barcodes that are not
found in the lookup table are retained unchanged, and original barcodes are
backed up in `*_from_backup.gz` files.

---

## Reference Files Required

| File | Description | Used by |
|------|-------------|---------|
| Guide FASTA | 20 bp guide spacer sequences, one per entry | All methods |
| t2g map | 2-column TSV: `guide_id → gene_id`, no header | simpleaf |
| GEX H5 (per lane) | scprocess decontX gene expression matrix | Whitelist extraction |
| 10x translation table | Feature ↔ GEX barcode mapping (3.7M entries) | Barcode translation |
| 10x barcode whitelist | 3M-february-2018 barcode list | simpleaf |
