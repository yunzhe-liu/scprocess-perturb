# scprocess-perturb
## A Single-Cell Perturbation Screen Preprocessing Workflow

**From Perturb-seq FASTQ to per-cell perturbation calls.** Extracts sgRNA guide
counts via simpleaf (piscem + alevin-fry) or [HAM](https://github.com/yunzhe-liu/ham)
(Hash-Accelerated Matcher), then assigns each cell a perturbation identity via
PGMM EM, UMI threshold, or fishash. A single `tenx_chemistry` setting
auto-resolves R1 geometry, barcode whitelist, translation, and tool parameters
for all standard 10x direct-capture chemistries.

---

## Workflow Overview

```
                    sgRNA FASTQ                GEX matrix
                         │                         │
                         ▼                         ▼
                    ┌─────────────────────────────────────┐
                    │  guide quantification               │
                    │  simpleaf (piscem + alevin-fry)     │
                    │  or HAM (hash match + dedup)        │
                    │  per lane → per-lane MEX            │
                    └────────────────┬────────────────────┘
                                     │
                                     ▼
                    ┌─────────────────────────────────────┐
                    │  merge_matrices                     │
                    │  vertical stack + lane suffix        │
                    │  automatic barcode translation       │
                    └────────────────┬────────────────────┘
                                     │
                                     ▼
                              guide_matrix/
                       (MEX trio, cells × guides)
                                     │
                                     ▼
                    ┌─────────────────────────────────────┐
                    │  guide assignment                   │
                    │  pgmm_em / umi_threshold / fishash  │
                    │  MEX → unified CSV → perturbation   │
                    └────────────────┬────────────────────┘
                                     │
                                     ▼
                          perturbation_obs.csv
```

**Guide extraction** produces a per-lane (cells × guides) UMI count matrix in
standard MEX format under `guide_matrix/`. **Guide assignment** converts that
matrix into per-cell perturbation calls — `perturbation_obs.csv` — ready to
join into expression AnnData or MuData objects.

---

## Quick Start

```bash
git clone https://github.com/yunzhe-liu/scprocess-perturb.git
cd scprocess-perturb
conda env create -f envs/scp_analysis.lock.yaml
conda env create -f envs/simpleaf.lock.yaml

# Copy and edit the template config files:
#   config/config.yaml    — global settings (extraction + assignment + paths)
#   config/groups.yaml    — input file paths per lane
# Then:
conda activate scp_analysis
snakemake --configfile config/config.yaml --cores 48
```

By default the workflow runs `pgmm_em` assignment on `guide_design: dual`.
To skip assignment and produce only `guide_matrix/`, remove the
`assignment:` section from `config.yaml`.  To enable additional methods,
add them to `assignment.methods`.

---

## Guide Extraction

Guide extraction quantifies sgRNA UMIs per cell from raw FASTQ.  Two methods
are supported under a unified chemistry layer.

### Quantification methods

| Method | `config` key | Engine | Best for |
|--------|:---:|--------|----------|
| simpleaf | `simpleaf` | piscem (pseudoalignment) + alevin-fry (UMI dedup) | Standard Perturb-seq; fast; `parsimony-gene` mode handles multi-mapped reads |
| HAM | `hash_matcher` | Hash-accelerated exact matching + Hamming-distance error tolerance | Datasets where read-level guide assignment precision matters; lower false-positive rate |

Switch methods with one field:

```yaml
guide_extraction:
  method: simpleaf        # or: hash_matcher
```

The workflow selects the appropriate rules automatically.

### Per-lane quantification

Each group in `config/groups.yaml` corresponds to one physical 10x lane.
The workflow extracts a cell barcode whitelist from the GEX matrix, runs
quantification, and produces one MEX trio per lane.  Barcode translation
between TruSeq (GEX) and Nextera (feature) formats is applied automatically
for 3' dual-oligo chemistries; 5' single-oligo chemistries skip translation.

### Merge → Merged MEX

The `merge_matrices` rule vertically stacks per-lane MEX matrices, appends
`-L{NN}` lane suffixes to cell barcodes, and produces the final merged trio:

```
{out_dir}/guide_matrix/
├── merged_matrix.mtx.gz       ← cells × guides, gzip compressed
├── merged_barcodes.tsv.gz     ← cell barcodes (-L{NN} suffix)
└── merged_features.tsv.gz     ← guide feature IDs
```

This is standard MEX format, directly loadable by Scanpy, Seurat, or any
MEX-compatible tool.  It is also the input to guide assignment.

### Knee mode (no whitelist needed)

When a GEX matrix is not available, set `use_knee: true` to let simpleaf
call cells automatically from the UMI knee.  Only available with
method=simpleaf (HAM requires a whitelist).  Knee-called cells are not
filtered by GEX QC; counts will be noisier than whitelist mode.

```yaml
simpleaf:
  quant:
    use_knee: true
```

### Method parameters

```yaml
simpleaf:
  af_home: /path/to/alevin-fry
  index:
    kmer_length: 15
    minimizer_length: 11
  quant:
    resolution: parsimony-gene

hash_matcher:
  umi_threshold: 1
  cb_max_hamming: 1
```

---

## Guide Assignment

After guide extraction, the workflow can run per-cell guide assignment —
converting the merged MEX matrix into a `cell → perturbation` mapping.

**This is optional.**  If the `assignment:` section is absent from
`config.yaml` or `methods` is empty, assignment is skipped and only the
merged MEX trio is produced.

### Methods

| Method | Description | Strengths | Script |
|--------|-------------|-----------|--------|
| `pgmm_em` *(default)* | Poisson-Gaussian Mixture Model via **MAP-EM**. 2-component fit per guide on log₂(UMI) with weak priors (LogNormal on λ, Dirichlet on w) that stop the background rate collapsing to zero; assigns all (cell, guide) pairs with UMI ≥ 1 AND P(Gaussian) ≥ 0.75. | Calibrated confidence (`prob_gaussian`) plus an additive confidence layer (`lpo`, `lpo_pctl`, `guide_Sg`) and a per-guide `guide_qc.csv`; fast; dual-guide compatible. | `run_pgmm_em.py` |
| `umi_threshold` | Fixed integer UMI threshold (default t=3). No statistical model — all pairs with UMI ≥ threshold are assigned. | Fastest; no model assumptions; predictable behaviour. | `run_umi_threshold.py` |
| `fishash` | One-sided Fisher exact test on the (guides × cells) contingency table, with iterative Simpson's paradox correction and Guo & Sarkar FDR control. | No parametric assumptions; rigorous FDR; maximum cell recovery. | `run_fishash.R` |

All parameters are configurable with sensible defaults.

### Quick configuration

```yaml
assignment:
  guide_design: dual              # single | dual | multi
  methods:
    - pgmm_em                     # default
    # - umi_threshold              # baseline
    # - fishash                   # Fisher test, needs R + fishash
```

Multiple methods run independently; each produces its own output under
`{out_dir}/assignment/{method}/`.

### Guide design modes & guide_csv schemas

The `guide_design` field controls how top-N guides are mapped to a
perturbation identity.  The `guide_csv` format depends on the mode:

| `guide_design` | Top-N | Perturbation identity | `guide_csv` columns | Example |
|:---|:---:|---|---|---|
| `single` | top-1 | target gene symbol | `guide_id, gene` | `NFKBIA_g1, NFKBIA` |
| `dual` | top-2 | construct pair_id (if top-2 belong to same construct) | `sgID_A, sgID_B, gene, pair_id` | `AAAS_-, AAAS_+, AAAS, pair_2` |
| `multi` | all assigned guides | construct_id (if all assigned guides belong to same construct) | `guide_id, gene, construct_id` | `g1, BRCA1, construct_5` |

**Dual-guide edge case:** when a cell's top-2 guides belong to different
constructs, that cell is marked `ambiguous_pair` in `perturbation_obs.csv`.
The cell is not silently dropped.

**Single-guide without guide_csv:** if `guide_csv` is omitted, the guide_id
itself is used as the gene name.

### Method parameters

All parameters have defaults and can be overridden under the method name:

```yaml
assignment:
  pgmm_em:
    umi_threshold: 1       # minimum UMI
    prob_threshold: 0.75   # P(Gaussian) assignment cutoff
    workers: 16            # parallel EM fitting workers
    max_em_iter: 200       # EM convergence limit

  umi_threshold:
    threshold: 3           # UMI cutoff (3, 5, 10 supported)

  fishash:
    padj_cutoff: 0.05      # FDR level
    padj_method: GS        # Guo & Sarkar 2020 block-dependence method
    min_count: 2           # minimum UMI per (cell, guide)
    refit: 10              # Simpson's paradox correction iterations
```

### Unified output schema

All methods write through `standardize_assignment.py`, producing a single
consistent CSV format:

```
cell_barcode     guide_id    umi_count  rank  score   score_type       method
AAACCCAAG-L27    GENE_g1     45         1     0.98    prob_gaussian    pgmm_em
AAACCCAAG-L27    GENE_g2     38         2     0.91    prob_gaussian    pgmm_em
```

| Column | Description |
|--------|-------------|
| `cell_barcode` | Cell identifier (`16mer-L{NN}` format) |
| `guide_id` | sgRNA or guide feature ID |
| `umi_count` | Raw UMI count |
| `rank` | Per-cell rank (1 = best, per method's own scoring) |
| `score` | Raw score value used for ranking |
| `score_type` | `prob_gaussian` / `umi_count` / `neg_log_pval` |
| `method` | Assignment method name |

### Per-cell perturbation obs

`make_perturbation_obs.py` reduces the unified schema to a single-row-per-cell
`cell → perturbation` mapping — the delivery artefact for downstream integration:

```
cell_barcode     perturbation   n_guides_assigned  assignment_score  assignment_confidence  assignment_method
AAACCCAAG-L27    pair_123       2                   0.98              high                   pgmm_em
AAACCCAAG-L28    ambiguous_pair 2                   0.45              low                    pgmm_em
```

| Column | Description |
|--------|-------------|
| `perturbation` | Construct ID (dual/multi), gene symbol (single), or `ambiguous_pair`/`NA` |
| `n_guides_assigned` | Number of guides contributing to this call |
| `assignment_score` | Top-1 guide's score |
| `assignment_confidence` | `high` / `mid` / `low` |
| `assignment_method` | Method name for traceability |

Confidence thresholds (initial defaults, calibratable per dataset):
`pgmm_em`: high ≥ 0.90, mid ≥ 0.75; `umi_threshold`: high ≥ 10, mid ≥ 5;
`fishash`: high ≥ 15, mid ≥ 10.

Join `perturbation_obs.csv` into an expression AnnData `.obs` or MuData to
connect each cell's perturbation identity with its transcriptome.

### Assignment outputs

```
{out_dir}/assignment/{method}/
├── _raw_assignments.csv          ← raw per-(cell, guide) output (pre-standardize)
├── guide_qc.csv                  ← per-guide fit params + S_g (pgmm_em only)
├── assignments.csv               ← unified schema (all guides, ranked)
├── perturbation_obs.csv          ← per-cell perturbation call
└── monitoring.json               ← wall time, peak RSS, stats
```

Run logs are written under `{log_dir}/assignment/` (`{method}.log` for the
assignment step, `{method}_obs.log` for the per-cell call).

---

## Configuration

You need **two files** to run the workflow.

### 1. Global config (`config/config.yaml`)

```yaml
# ── Paths ──
proj_dir: /path/to/scprocess-perturb
out_dir: /path/to/results
log_dir: /path/to/logs

# ── Chemistry ──
# 3v3 | 3v4 | 3LT | multiome | 5v1 | 5v2 | 5v3 | custom
tenx_chemistry: "3v3"

# ── Guide library ──
# CSV used to generate FASTA + t2g for quantification and to map
# assigned guide IDs to genes/constructs.
guide_csv: /path/to/guide_library.csv

# ── Guide extraction ──
guide_extraction:
  method: simpleaf                    # simpleaf (default) | hash_matcher

simpleaf:
  af_home: /path/to/alevin-fry

# ── Guide assignment (v0.2.0+) ──
# Remove this entire section to produce guide_matrix/ only.
assignment:
  guide_design: dual                  # single | dual | multi
  methods:
    - pgmm_em

# ── Optional overrides ──
# All sections below are optional — defaults are auto-derived.

# Want to switch to HAM extraction?
# guide_extraction:
#   method: hash_matcher
# hash_matcher:
#   umi_threshold: 1
#   cb_max_hamming: 1

# Override auto-derived reference paths?
# references:
#   guide_fasta: /path/to/guides.fasta
#   guide_t2g_2col: /path/to/t2g.tsv
#   sgRNA_index_dir: /path/to/piscem_index/
#   guide_hash: /path/to/guide_hash.pkl
#   whitelist_dir: /path/to/whitelist_cache/

# Override whitelist QC thresholds?
# whitelist:
#   min_umi: 1000
#   min_genes: 500

# Override resource limits?
# resources:
#   simpleaf_quant_threads: 12
#   simpleaf_index_threads: 4
#   hash_quant_threads: 4
```

### 2. Input data topology (`config/groups.yaml`)

One group per physical 10x lane:

```yaml
groups:
  lane_01:
    group_id: lane_01
    gex_h5: /path/to/expression.h5ad       # .h5 / .h5ad / .h5mu
    sgRNA_fastq_dir: /path/to/sgRNA_fastq/
    sgRNA_r1_pattern: "*_R1_001.fastq.gz"
    sgRNA_r2_pattern: "*_R2_001.fastq.gz"
```

| Extension | Format | How it's read |
|:---|:---|:---|
| `.h5` | scprocess raw H5 (CSC sparse) | `filter_barcodes.py` |
| `.h5ad` | AnnData | `anndata.read_h5ad()` |
| `.h5mu` | MuData | reads `mod['rna']` |

The path to `groups.yaml` can be changed via `groups_file` in `config.yaml`.

**About `tenx_chemistry`:** this single field controls everything
chemistry-related.  The workflow derives downstream values (simpleaf chemistry,
whitelist, translation, HAM chemistry, UMI length) from
[config/chemistry_spec.yaml](config/chemistry_spec.yaml).  See the table below.

Configs without `tenx_chemistry` (pre-v0.1.2) fall back to reading manual
chemistry fields; this is supported but deprecated.

### Supported Chemistries

**Class A — 3' Direct Capture** · Bead-borne cs1/cs2 RT primer ·
Dual-oligo barcode (TruSeq + Nextera) · Translation required

| `tenx_chemistry` | 10x Kit | R1 layout | UMI | Whitelist | Translation | Validated |
|:---|:---|:---|:---:|:---|:---:|:---:|
| `3v3` | 3' v3 / v3.1 | 28bp (16CB + 12UMI) | 12bp | 3M-feb-2018 | Yes | ✅ |
| `3v4` | 3' v4 (GEM-X) | 28bp (16CB + 12UMI) | 12bp | 3M-3pgex-may-2023 | Yes | 💡 |
| `3LT` | 3' LT | 28bp (16CB + 12UMI) | 12bp | 3M-feb-2018 | Yes | 💡 |
| `multiome` | Multiome (GEX) | 28bp (16CB + 12UMI) | 12bp | 3M-feb-2018 | Yes | 💡 |

> `3v3` and `3LT` share identical library structure. `3v4` uses a distinct whitelist (`3M-3pgex-may-2023`) and translation file.

**Class B — 5' Direct Capture** · Soluble scaffold RT primer + barcoded TSO ·
Single-oligo barcode (TruSeq only) · No translation

| `tenx_chemistry` | 10x Kit | R1 layout | UMI | Whitelist | Translation | Validated |
|:---|:---|:---|:---:|:---|:---:|:---:|
| `5v1` | 5' v1.0 | 26bp (16CB + 10UMI) | 10bp | 737K-aug-2016 | No | ✅ |
| `5v2` | 5' v2 | 26bp (16CB + 10UMI) | 10bp | 737K-aug-2016 | No | 💡 |
| `5v3` | 5' v3 (GEM-X) | 28bp (16CB + 12UMI) | 12bp | 3M-5pgex-jan-2023 | No | 💡 |

> `5v3` shares the same single-oligo mechanism as `5v1`/`5v2` (no translation). Uses 12bp UMI and a distinct whitelist.

**Class C — Custom** · All parameters supplied by the user.

| `tenx_chemistry` | R1 layout | UMI | Whitelist | Translation | Validated |
|:---|:---|:---:|:---|:---:|:---:|
| `custom` | user-defined | user-defined | user-provided | user-defined | — |

> **`3v2` is explicitly absent.** It predates direct capture and uses GBC indirect capture — the guide sequence is not in the guide FASTQ R2.
> ✅ = validated on real data · 💡 = theoretically supported, pending validation

### Chemistry Overrides

Override individual fields without switching to `custom`:

```yaml
tenx_chemistry: "5v1"
chemistry_overrides:
  guide_len: 21
```

Available keys: `af_chemistry`, `whitelist`, `expected_ori`, `translation`,
`geometry_override`, `ham_chemistry`, `umi_len`.

### Custom Chemistry

```yaml
tenx_chemistry: custom
custom_chemistry:
  af_chemistry: "1{b[14]u[8]x:}2{r:}"
  whitelist: /path/to/whitelist.txt
  expected_ori: fw
  translation: false
  ham_chemistry: custom
  umi_len: 8

hash_matcher:
  custom_params:
    cb_start: 0; cb_end: 14
    umi_start: 14; umi_end: 22
    window_start: 8; window_end: 28
    guide_len: 20
```

---

## Installation

### Conda environments

| Environment | Contains | Purpose |
|-------------|----------|---------|
| `scp_analysis` | Snakemake, HAM, numpy, scipy, h5py, anndata, mudata | Orchestration, HAM quant, merge, assignment |
| `simpleaf` | simpleaf ≥ 0.24, piscem ≥ 0.19, alevin-fry ≥ 0.14 | simpleaf quant |

```bash
conda activate scp_analysis    # Snakemake + HAM + merge + assignment
conda activate simpleaf        # simpleaf quant
```

### Dependencies

- **Snakemake** ≥ 7.0
- **simpleaf** (piscem + alevin-fry)
- **[HAM](https://github.com/yunzhe-liu/ham)** — hash-accelerated guide matching
- **Python:** numpy, scipy, h5py, anndata, mudata, pandas
- **[fishash](https://github.com/jackkamm/fishash)** — R package (only needed if the `fishash` assignment method is enabled)

---

## Usage

```bash
cd /path/to/scprocess-perturb
conda activate scp_analysis
snakemake --configfile config/config.yaml --cores 48
```

Common Snakemake options:

```
--cores N           Max parallel threads
--dry-run (-n)      Preview execution plan
--unlock            Remove stale locks after a crash
--rerun-incomplete  Re-run partially completed jobs
--latency-wait 60   Wait for NFS latency (seconds)
```

---

## Inputs and Outputs

### What you need to provide

| Input | Format | Notes |
|-------|--------|-------|
| sgRNA FASTQ | Paired FASTQ (R1 + R2), per lane | R1 = cell barcode + UMI, R2 = protospacer sequence |
| GEX matrix | `.h5` / `.h5ad` / `.h5mu`, per lane | Used for cell barcode whitelist extraction |
| Guide library CSV | See [Guide Assignment](#guide-assignment) for schema | Auto-generates FASTA + t2g + piscem index. Provides guide→gene mapping for assignment. |

All other references (guide FASTA, t2g, piscem index, HAM hash table, 10x
whitelists, translation tables) are auto-generated or auto-downloaded by the
workflow under `{out_dir}/refs/`. The `3M-3pgex-may-2023` whitelist (3′ v4)
has no public mirror and must be copied manually from Cell Ranger ≥ 8.0.1.

### What the workflow produces

```
{out_dir}/
├── refs/                          ← auto-generated references (FASTA, index, whitelists)
├── lanes/
│   └── lane_XX/
│       └── {quant_subdir}/        ← per-lane MEX (intermediate)
├── guide_matrix/
│   ├── merged_matrix.mtx.gz       ← cells × guides UMI count matrix
│   ├── merged_barcodes.tsv.gz     ← cell barcodes (-L{NN} suffix)
│   └── merged_features.tsv.gz     ← guide feature IDs
└── assignment/
    └── {method}/
        ├── assignments.csv         ← unified schema (all guides, ranked)
        └── perturbation_obs.csv    ← per-cell perturbation call
```

The workflow also generates intermediate artefacts during reference setup:
t2g map, piscem index (simpleaf), guide hash table (HAM), 10x whitelists.

---

## Repository Structure

```
scprocess-perturb/
├── Snakefile
├── config/
│   ├── chemistry_spec.yaml       ← Chemistry parameters (all supported 10x kits)
│   ├── config.yaml               ← Global config template
│   └── groups.yaml               ← Input data topology
├── rules/
│   ├── reference.smk             ← Guide FASTA → index / hash + whitelist download
│   ├── whitelist.smk             ← GEX matrix → cell barcode whitelist
│   ├── quant.smk                 ← simpleaf quantification
│   ├── guide_quant.smk           ← HAM quantification
│   ├── merge.smk                 ← Per-lane merge + barcode translation
│   └── assignment.smk            ← Guide assignment
├── scripts/
│   ├── feature_reference_adapter.py ← Guide FASTA / t2g from CSV
│   ├── build_guide_hash.py       ← HAM hash table builder
│   ├── filter_barcodes.py        ← Cell barcode QC (h5)
│   ├── translate_barcodes.py     ← Feature ↔ GEX barcode translation
│   ├── preprocess_cutadapt.sh    ← Optional R2 trimming (preprocess.trimmed)
│   ├── parse_samples.py          ← groups.yaml helper for preprocessing
│   ├── run_pgmm_em.py            ← pgmm_em assignment
│   ├── run_umi_threshold.py      ← umi_threshold assignment
│   ├── run_fishash.R             ← fishash assignment (R)
│   ├── standardize_assignment.py ← Method CSV → unified schema
│   └── make_perturbation_obs.py  ← Unified schema → per-cell call
├── envs/
│   ├── simpleaf.lock.yaml
│   └── scp_analysis.lock.yaml
└── profiles/local/
```

---

## Snakemake Rule Sequence

1. **`reference.smk`** — Builds piscem index / HAM hash table from guide FASTA;
   downloads 10x whitelists if not cached.
2. **`whitelist.smk`** — Extracts cell barcodes from GEX matrix. Applies
   TO→FROM translation when chemistry requires it.
3. **`quant.smk`** (simpleaf) / **`guide_quant.smk`** (HAM) — Quantifies guide
   counts per cell per lane.
4. **`merge.smk`** — Vertically stacks per-lane MEX matrices, appends lane
   suffixes, applies FROM→TO translation.
5. **`assignment.smk`** — Runs guide assignment on merged MEX, standardises
   output, produces per-cell perturbation calls. Skipped when `assignment:`
   is absent from config.

---

## Barcode Translation

**Dual-oligo systems (3' v3/v4, 3LT, multiome):** Bead-borne cs1/cs2 primers
carry two oligo variants — TruSeq (mRNA) and Nextera (feature).  The workflow
translates: (1) GEX→Feature for whitelist matching, (2) Feature→GEX for merge
output. **Single-oligo systems (5' v1/v2/v3):** Soluble RT primer indexes
guide cDNA to the same TruSeq barcode as mRNA.  No translation.
