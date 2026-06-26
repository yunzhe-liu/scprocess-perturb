# sgprocess — Single-Cell Guide Extraction Workflow

**A Snakemake workflow for extracting sgRNA guide counts from Perturb-seq data.**

Integrates two quantification methods (simpleaf — piscem + alevin-fry, and
[HAM](https://github.com/yunzhe-liu/ham) — Hash-Accelerated Matcher) under a
unified chemistry-configuration layer. A single `tenx_chemistry` setting
automatically resolves R1 geometry, barcode whitelist, translation behaviour,
and tool-specific parameters for all standard 10x direct-capture chemistries.

---

## Quick Start

```bash
git clone https://github.com/yunzhe-liu/sgprocess.git
cd sgprocess
conda env create -f envs/scp_analysis.lock.yaml
conda env create -f envs/simpleaf.lock.yaml

# Write your config (see Configuration section below), then:
conda activate scp_analysis
snakemake --configfile config/your_config.yaml --cores 48

# Preview execution plan
snakemake --configfile config/your_config.yaml --cores 1 --dry-run
```

Example configs are provided in `config/` — use them as starting points.

---

## Installation

### Conda environments

| Environment | Contains | Purpose |
|-------------|----------|---------|
| `scp_analysis` | Snakemake, HAM, numpy, scipy, h5py, anndata, mudata | Workflow orchestration, HAM quant, merge, analysis |
| `simpleaf` | simpleaf ≥ 0.24, piscem ≥ 0.19, alevin-fry ≥ 0.14 | simpleaf quant |

```bash
conda activate scp_analysis    # for Snakemake + HAM + merge
conda activate simpleaf        # for simpleaf quant
```

### Dependencies

- **Snakemake** ≥ 7.0
- **simpleaf** (piscem + alevin-fry) for k-mer pseudoalignment
- **[HAM](https://github.com/yunzhe-liu/ham)** — Hash-Accelerated Matcher for hash-based guide matching
- **Python packages:** numpy, scipy, h5py, anndata, mudata

---

## Repository Structure

```
sgprocess/
├── Snakefile                     ← Snakemake entry point
├── config/
│   ├── chemistry_spec.yaml       ← Chemistry parameters (single source of truth)
│   ├── config.yaml               ← Default config
│   ├── samples.yaml              ← Sample/lane topology
│   └── ...                       ← Additional config variants
├── rules/
│   ├── reference.smk             ← Guide FASTA → index / hash table + whitelist download
│   ├── whitelist.smk             ← GEX matrix → barcode whitelist (.h5 / .h5ad / .h5mu)
│   ├── quant.smk                 ← simpleaf quant
│   ├── guide_quant.smk           ← HAM quant (match + dedup)
│   └── merge.smk                 ← Per-lane merge + barcode translation
├── scripts/
│   ├── translate_barcodes.py     ← Feature ↔ GEX barcode translation
│   ├── build_guide_hash.py       ← HAM hash table builder
│   ├── filter_barcodes.py        ← Barcode QC thresholds (h5 format)
│   ├── feature_reference_adapter.py   ← Guide FASTA / t2g generation
│   └── ...
├── envs/
│   ├── simpleaf.lock.yaml
│   └── scp_analysis.lock.yaml
└── profiles/
    └── local/
```

---

## Configuration

You need **two files**: a global config and a samples file.

### 1. Global config (`config/your_config.yaml`)

All parameters with their meanings and defaults:

```yaml
# ── Project ──
proj_dir: /path/to/sgprocess                     # path to this repository
out_dir: /path/to/results                        # where output files are written
log_dir: /path/to/logs                           # where log files are written
samples_file: config/samples.yaml                # optional: path to samples file

# ── Chemistry (v0.1.2+) ──
# One field controls everything. See §Supported Chemistries below.
tenx_chemistry: "3v3"                            # required
#   Options: 3v3, 3v4, 3LT, 5v3, multiome, 5v1, 5v2, custom

# ── Optional: override individual chemistry fields ──
chemistry_overrides:                             # optional
  umi_len: 11                                    # example: override UMI length only

# ── Quantification method ──
guide_extraction:
  method: simpleaf                               # "simpleaf" | "hash_matcher"

# ── simpleaf settings (used when method=simpleaf) ──
simpleaf:
  af_home: /path/to/alevin-fry                   # required
  index:
    kmer_length: 15                              # default: 15
    minimizer_length: 11                         # default: 11
  quant:
    resolution: parsimony-gene                   # "parsimony-gene" | "cr-like"
    use_knee: false                              # true = knee calling, no whitelist needed

# ── HAM settings (used when method=hash_matcher) ──
hash_matcher:
  umi_threshold: 1                               # UMI Hamming threshold for dedup (default: 1)
  cb_max_hamming: 1                              # max Hamming for CB correction (default: 1)

# ── Whitelist QC thresholds ──
whitelist:
  min_umi: 1000                                  # minimum UMI count per cell
  min_genes: 500                                 # minimum genes detected per cell

# ── Reference files ──
references:
  guide_fasta: /path/to/guides.fasta             # required: guide protospacer sequences
  guide_t2g_2col: /path/to/t2g.tsv               # required for simpleaf: guide_id → gene_id, no header
  guide_hash: /path/to/guide_hash.pkl            # required for HAM: pre-built hash table
  sgRNA_index_dir: /path/to/piscem_index/        # required for simpleaf: piscem index directory
  whitelist_dir: /path/to/whitelist_cache/       # 10x whitelist cache directory
  guide_csv: /path/to/guide_library.csv          # optional: CSV for generating FASTA/t2g
  translation_table: /path/to/translation.txt    # required if translation=true

# ── Resource allocation ──
resources:
  simpleaf_quant_threads: 12                     # default: 12
  simpleaf_index_threads: 4                      # default: 4
  hash_quant_threads: 4                          # default: 4
```

**About `tenx_chemistry`:** This is the only chemistry field you need to set.
The workflow automatically derives all downstream values — simpleaf chemistry
string, barcode whitelist, translation behaviour, HAM chemistry, and UMI
length — from [config/chemistry_spec.yaml](config/chemistry_spec.yaml). See the
chemistry table below for what each value means.

If you omit `tenx_chemistry` entirely (pre-v0.1.2 configs), the workflow falls
back to reading manual chemistry fields from the config. This is supported but
deprecated.

### 2. Samples file

Defines one group per physical 10x lane. The `gex_h5` field auto-detects
format by extension:

```yaml
groups:
  lane_01:
    group_id: lane_01
    gex_h5: /path/to/expression_matrix.{h5,h5ad,h5mu}
    sgRNA_fastq_dir: /path/to/sgRNA_fastq/
    sgRNA_r1_pattern: "*.fastq.gz"               # glob pattern for R1 files
    sgRNA_r2_pattern: "*.fastq.gz"               # glob pattern for R2 files
```

Supported GEX formats:

| Extension | Format | How it's read |
|:---|:---|:---|
| `.h5` | scprocess raw H5 (CSC sparse) | `filter_barcodes.py` script |
| `.h5ad` | AnnData | Inline Python via `anndata.read_h5ad()` |
| `.h5mu` | MuData | Inline Python, reads `mod['rna']` |

### Supported Chemistries

| `tenx_chemistry` | 10x Kit | R1 layout | UMI | Whitelist | Translation | Guide capture |
|:---|:---|:---|:---:|:---|:---:|:---|
| `3v3` | 3' v3 / v3.1 | 28bp (16CB+12UMI) | 12bp | 3M-feb-2018 | Yes (TruSeq↔Nextera) | cs1/cs2 bead |
| `3v4` | 3' v4 (GEM-X) | 28bp (16CB+12UMI) | 12bp | 3M-feb-2018 | Yes | cs1/cs2 bead |
| `3LT` | 3' LT | 28bp (16CB+12UMI) | 12bp | 3M-feb-2018 | Yes | cs1/cs2 bead |
| `5v3` | 5' v3 | 28bp (16CB+12UMI) | 12bp | 3M-feb-2018 | Yes | cs1/cs2 bead |
| `multiome` | Multiome (GEX) | 28bp (16CB+12UMI) | 12bp | 3M-feb-2018 | Yes | cs1/cs2 bead |
| `5v1` | 5' v1.0 | 26bp (16CB+10UMI) | 10bp | 737K-aug-2016 | No (TruSeq only) | Soluble RT primer |
| `5v2` | 5' v2 | 26bp (16CB+10UMI) | 10bp | 737K-aug-2016 | No (TruSeq only) | Soluble RT primer |
| `custom` | User-defined | user-defined | user-defined | user-provided | user-defined | user-defined |

**What `tenx_chemistry` resolves automatically:**

| Parameter | Consumed by | Example for `5v1` |
|:---|:---|:---|
| simpleaf `--chemistry` | `quant.smk` | `1{b[16]u[10]x:}2{r:}` (explicit geometry) |
| Barcode whitelist | `whitelist.smk`, `quant.smk` | `737K-august-2016.txt` |
| Expected orientation | `quant.smk` | `rc` |
| Barcode translation | `whitelist.smk`, `merge.smk` | skipped (TruSeq only) |
| HAM `--chemistry` | `guide_quant.smk` | `10xv2-5p` |
| HAM `--umi-len` | `guide_quant.smk` | `10` |

**Unsupported:** Indirect capture methods (e.g. 3' v2 with GBC) where the
protospacer is absent from the guide FASTQ R2. This pipeline requires direct
capture.

### Chemistry Overrides

When you have standard hardware but a non-standard parameter (e.g. custom guide
length), override individual fields without leaving the standard chemistry path:

```yaml
tenx_chemistry: "5v1"
chemistry_overrides:
  umi_len: 11     # all other fields remain from the 5v1 spec entry
```

### Custom Chemistry

For completely non-standard hardware or custom barcode designs:

```yaml
tenx_chemistry: custom
custom_chemistry:
  af_chemistry: "1{b[14]u[8]x:}2{r:}"            # raw geometry for simpleaf
  whitelist: /path/to/custom_whitelist.txt
  expected_ori: fw
  translation: false
  ham_chemistry: custom                           # signals HAM to use manual params
  umi_len: 8

hash_matcher:
  custom_params:                                  # passed as individual CLI flags to HAM
    cb_start: 0
    cb_end: 14
    umi_start: 14
    umi_end: 22
    window_start: 8
    window_end: 28
    guide_len: 20
```

When `ham_chemistry: custom`, the workflow passes `--cb-start`, `--umi-start`,
`--window-start`, `--guide-len` etc. directly to `ham match` instead of using a
named chemistry.

---

## Usage

### Basic invocation

```bash
cd /path/to/sgprocess
conda activate scp_analysis

# Single config
snakemake --configfile config/your_config.yaml --cores 48

# Dry-run (check what will run without executing)
snakemake --configfile config/your_config.yaml --cores 1 --dry-run
```

### Common Snakemake options

```
--cores N           Max parallel threads
--dry-run (-n)      Preview execution plan without running
--unlock            Remove stale locks after a crash
--rerun-incomplete  Re-run partially completed jobs
--latency-wait 60   Wait for NFS filesystem latency (seconds)
```

### Switching quantification method

Change one field in your config:

```yaml
guide_extraction:
  method: simpleaf        # or: hash_matcher
```

The workflow selects the appropriate rules automatically — no other changes
needed.

---

## Inputs and Outputs

### What you need to provide

| File | Description | Required by |
|------|-------------|:---:|
| Guide FASTA | One protospacer sequence per guide entry | All methods |
| t2g map | 2-column TSV: `guide_id → gene_id`, no header | simpleaf |
| GEX matrix | Expression data per lane (`.h5`, `.h5ad`, or `.h5mu`) | All methods |
| Guide FASTQ | Paired sgRNA FASTQ files (R1 + R2) per lane | All methods |
| Translation table | Feature ↔ GEX barcode mapping (~3.7M entries) | 3' chemistries only |
| HAM guide hash | Pre-built pickle hash table | HAM only |
| piscem index | Pre-built piscem dense index | simpleaf only |

10x barcode whitelists (`3M-february-2018.txt`, `737K-august-2016.txt`) are
auto-downloaded by the workflow if not present in `whitelist_dir`.

### What the workflow produces

```
{out_dir}/
├── lane_XX/
│   └── {quant_subdir}/
│       └── (per-lane MEX matrix and barcode/feature lists)
└── merged/
    ├── merged_matrix.mtx.gz             ← cell × guide count matrix
    ├── merged_barcodes.tsv.gz           ← cell barcodes (GEX format, -L{NN} suffix)
    └── merged_features.tsv.gz           ← guide feature IDs
```

The merged output is in standard **MEX format** (Market Exchange Format,
gzip-compressed), directly loadable by Scanpy, Seurat, or any MEX-compatible
tool.

---

## Workflow Overview

```
                    ┌─────────────────────────────────────────┐
GEX matrix (H5/H5AD/H5MU) → extract whitelist                │
Guide FASTA ───────────────→ build index / hash table        │
sgRNA FASTQ ──────────────→┤                                 │
                            │   ┌─ simpleaf: piscem + alevin-fry
                            ├──→│─ HAM: hash match + dedup  ├──→ merge ──→ MEX
                            │   └───────────────────────────┘   + barcode
                            └───────────────────────────────────┘   translation
                                                                     (automatic)
```

**Core design:** GEX mRNA signal defines which droplets contain real cells.
A whitelist of cell barcodes passing QC thresholds (minimum UMI, minimum
genes) is extracted, automatically translated between barcode formats when
the chemistry requires it, and fed to the guide quantification tools.

---

## Snakemake Rule Sequence

1. **`reference.smk`** — Builds method-specific index/hash table from guide
   FASTA; downloads 10x whitelist files if not cached locally.
2. **`whitelist.smk`** — Extracts cell barcodes from GEX matrix (supports
   `.h5`, `.h5ad`, `.h5mu`). Applies TO→FROM translation when chemistry
   requires it.
3. **`quant.smk`** (simpleaf) or **`guide_quant.smk`** (HAM) — Quantifies
   guide counts per cell. Chemistry parameters are auto-resolved.
4. **`merge.smk`** — Stacks per-lane MEX matrices, appends lane suffixes.
   Applies FROM→TO translation when chemistry requires it.

---

## Barcode Translation

Two capture mechanisms exist in 10x chemistries. The workflow handles both
automatically through the `translation` field in `chemistry_spec.yaml`.

**Dual-oligo systems (3' v3/v4, 5' v3, multiome):** cs1/cs2 bead-borne
primers carry two oligo variants — TruSeq for mRNA (GEX) and Nextera for
Feature Barcoding (guide). These encode the same cell barcode in different
nucleotide formats. The workflow translates:
- **Whitelist:** GEX (TruSeq) → Feature (Nextera), so the whitelist
  matches guide FASTQ barcodes.
- **Merge:** Feature → GEX, so output barcodes match mRNA conventions.

**Single-oligo systems (5' v1/v2):** A soluble RT primer indexes guide cDNA
to the same TruSeq barcode as mRNA. Both libraries share identical barcode
formats. No translation is performed.

---

## Changelog

### v0.1.2 — 2026-06-26

- **`config/chemistry_spec.yaml`:** Central chemistry specification — one
  YAML file is the single source of truth for all chemistry parameters.
  Adding a new chemistry is ~7 lines of YAML.
- **`Snakefile`:** New `_resolve_chemistry()` function loads the spec and
  populates `config["_chemistry"]` before rules are included. The
  `samples_file` config key now allows overriding the samples file path.
- **`rules/whitelist.smk`:** TO→FROM barcode translation step added,
  automatically skipped for single-oligo chemistries.
- **`rules/merge.smk`:** Translation control reads from
  `_chemistry.translation` preferentially.
- **`rules/guide_quant.smk`:** Custom chemistry escape hatch fully
  implemented — `ham_chemistry: custom` unrolls individual position
  parameters as CLI flags.
- **`rules/quant.smk`:** Chemistry (including geometry overrides) resolved
  from `_chemistry` spec.
- **HAM:** `--chemistry custom` accepted; independent position flags
  (`--cb-start`, `--umi-start`, `--window-start`, `--guide-len`, etc.)
  added to CLI and Python API.
- **`5v2`** maps correctly to `10xv2-5p` (previously grouped with `3v2`).
- **`3v2`** removed from supported chemistries (indirect GBC capture, no
  protospacer in guide FASTQ R2).
- **`chemistry_overrides`** mechanism: override individual fields on top of
  a standard chemistry entry.
- **Backward compatibility:** Configs without `tenx_chemistry` continue to
  work unchanged.

### v0.1.1 — 2026-06-26

- **5' v1 chemistry support:** Added `5v1` chemistry mapping (10xv2-5p,
  737K whitelist, soluble RT primer, no translation).
- **`rules/reference.smk`:** `download_whitelists` rule auto-fetches
  `3M-february-2018.txt` and `737K-august-2016.txt` from 10x CDN.
- **`rules/whitelist.smk`:** Added `.h5ad` and `.h5mu` input support
  alongside the original `.h5` format.
- **`rules/merge.smk`:** `skip_translation` flag controls barcode
  translation at merge.
- **`rules/guide_quant.smk`:** HAM `--chemistry` and `--umi-len` flags
  passed from config.
- **Bug fixes:** `simpleaf set-paths` call added before quant; index path
  corrected to `.../index/piscem_idx`; hardcoded geometry in mapping.py
  replaced with `geometry_override` parameter.
