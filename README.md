# sgprocess — Single-Cell Guide Extraction Workflow

**A Snakemake workflow for extracting sgRNA guide counts from Perturb-seq data.**

Integrates two quantification methods (simpleaf piscem + alevin-fry, and HAM —
Hash-Accelerated Matcher) under a unified chemistry-configuration layer. A
single `tenx_chemistry` setting automatically resolves R1 geometry, barcode
whitelist, barcode translation behaviour, and tool-specific parameters for all
standard 10x chemistries.

---

## Quick Start

```bash
# Basic invocation
snakemake --configfile config/config.yaml --cores 48

# Dry-run (preview execution plan)
snakemake --configfile config/config.yaml --cores 12 --dry-run
```

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

**Core design principle:** GEX mRNA signal defines which droplets contain real
cells. A whitelist of cell barcodes passing GEX QC thresholds is extracted,
automatically translated between barcode formats when the chemistry requires it,
and fed to the guide quantification tools.

The workflow produces **standard MEX-formatted count matrices** (gzip-compressed
Market Exchange Format: `matrix.mtx.gz`, `barcodes.tsv.gz`, `features.tsv.gz`).

---

## Supported 10x Chemistries

Setting `tenx_chemistry` in your config is all that's needed. The workflow
automatically resolves R1 geometry, barcode whitelist, translation behaviour,
and tool-specific parameters from a central specification file
([config/chemistry_spec.yaml](config/chemistry_spec.yaml)).

| `tenx_chemistry` | 10x Kit | R1 | UMI | Whitelist | Translation | Guide capture |
|:---|:---|:---:|:---:|:---|:---:|:---|
| `3v3` / `3v4` / `3LT` | 3' v3 / v3.1 / v4 | 28bp | 12bp | 3M-feb-2018 | Yes (TruSeq↔Nextera) | cs1/cs2 bead primer |
| `5v3` / `multiome` | 5' v3 / Multiome | 28bp | 12bp | 3M-feb-2018 | Yes | cs1/cs2 bead primer |
| `5v1` | 5' v1.0 | 26bp | 10bp | 737K-aug-2016 | No (TruSeq only) | Soluble RT primer |
| `5v2` | 5' v2 | 26bp | 10bp | 737K-aug-2016 | No (TruSeq only) | Soluble RT primer |

**What the workflow resolves automatically from `tenx_chemistry`:**
- Piscem geometry string (including workarounds for known `__builtin` handler bugs)
- Barcode whitelist selection and auto-download from 10x CDN
- Whether Feature↔GEX barcode translation is needed at whitelist and merge steps
- HAM chemistry constants (UMI length, window position, guide length)

**Unsupported chemistries:** Indirect capture methods (e.g., 3' v2 with GBC)
where the protospacer is absent from the guide FASTQ R2 are not supported by
this workflow. The pipeline requires direct capture.

---

## Installation

### Conda environments

Two conda environments are required:

| Environment | Contains | Purpose |
|-------------|----------|---------|
| `scp_analysis` | Snakemake, HAM, numpy, scipy, h5py, anndata, mudata | Workflow orchestration, HAM quant, merge, analysis |
| `simpleaf` | simpleaf ≥ 0.24, piscem ≥ 0.19, alevin-fry ≥ 0.14 | simpleaf quant |

```bash
conda activate scp_analysis    # for Snakemake + HAM + merge
conda activate simpleaf        # for simpleaf quant
```

### Dependencies

- **Snakemake** ≥ 8.0
- **simpleaf** (piscem + alevin-fry) for k-mer pseudoalignment
- **HAM** — Hash-Accelerated Matcher for hash-based guide matching
- **Python packages:** numpy, scipy, h5py, anndata, mudata

### Workflow installation

```bash
git clone <url> /path/to/sgprocess
cd /path/to/sgprocess
conda env create -f envs/scp_analysis.lock.yaml   # HAM + Snakemake + analysis
conda env create -f envs/simpleaf.lock.yaml       # simpleaf only
```

---

## Repository Structure

```
sgprocess/
├── Snakefile                     ← Snakemake entry point
├── config/
│   ├── chemistry_spec.yaml       ← Single source of truth for chemistry parameters
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

### Chemistry — Single-Entry Resolution (v0.1.2)

The central configuration innovation in v0.1.2: **one field controls everything**.

```yaml
# All you need to set:
tenx_chemistry: "5v1"
```

From this single entry, the workflow automatically derives:

| Derived parameter | Value for `5v1` | Used by |
|:---|:---|:---|
| simpleaf chemistry | `1{b[16]u[10]x:}2{r:}` (explicit geometry) | `quant.smk` |
| Barcode whitelist | `737K-august-2016.txt` | `whitelist.smk`, `quant.smk` |
| Expected orientation | `rc` | `quant.smk` |
| Barcode translation | `false` (skip) | `whitelist.smk`, `merge.smk` |
| HAM chemistry | `10xv2-5p` | `guide_quant.smk` |
| UMI length | 10 bp | `guide_quant.smk` |

The mapping is defined in [config/chemistry_spec.yaml](config/chemistry_spec.yaml)
— a single YAML file that serves as the source of truth for all chemistry
parameters. Adding a new standard chemistry is a matter of adding a ~7-line
entry to this file.

### `samples.yaml` — Lane Topology

Defines one group per physical 10x lane:

```yaml
groups:
  lane_01:
    group_id: lane_01
    gex_h5: /path/to/expression_matrix.{h5,h5ad,h5mu}
    sgRNA_fastq_dir: /path/to/sgRNA_fastq/
    sgRNA_r1_pattern: "*.fastq.gz"
    sgRNA_r2_pattern: "*.fastq.gz"
```

The `gex_h5` field accepts three formats (auto-detected by extension):
- `.h5` — scprocess raw H5 matrix (CSC sparse, legacy format)
- `.h5ad` — AnnData (published data, pertpy output)
- `.h5mu` — MuData (reads `mod['rna']` modality)

### `config.yaml` — Global Settings

```yaml
tenx_chemistry: "3v3"             # ← Single chemistry entry (v0.1.2)

guide_extraction:
  method: simpleaf                # "simpleaf" | "hash_matcher"

simpleaf:
  af_home: /path/to/alevin-fry
  index:
    kmer_length: 15
    minimizer_length: 11
  quant:
    resolution: parsimony-gene
    use_knee: false
    # chemistry is auto-derived from tenx_chemistry

hash_matcher:
  umi_threshold: 1
  cb_max_hamming: 1
  # chemistry is auto-derived from tenx_chemistry

whitelist:
  min_umi: 1000
  min_genes: 500

references:
  guide_fasta: /path/to/guides.fasta
  guide_t2g_2col: /path/to/t2g.tsv
  guide_hash: /path/to/guide_hash.pkl        # HAM only
  sgRNA_index_dir: /path/to/simpleaf_index/   # simpleaf only
  whitelist_dir: /path/to/whitelist_cache/

resources:
  simpleaf_quant_threads: 12
  hash_quant_threads: 4
```

### Chemistry Overrides

For standard hardware with non-standard experimental parameters (e.g., custom
guide length), override individual fields without leaving the standard chemistry
path:

```yaml
tenx_chemistry: "5v1"
chemistry_overrides:
  umi_len: 11     # override only this field
```

All other parameters continue to be derived from the `5v1` entry in
`chemistry_spec.yaml`.

### Custom Chemistry

For completely non-standard hardware or custom barcode designs, use the
`custom` escape hatch:

```yaml
tenx_chemistry: custom
custom_chemistry:
  af_chemistry: "1{b[14]u[8]x:}2{r:}"
  whitelist: /path/to/custom_whitelist.txt
  expected_ori: fw
  translation: false
  ham_chemistry: custom
  umi_len: 8

hash_matcher:
  custom_params:
    cb_start: 0
    cb_end: 14
    umi_start: 14
    umi_end: 22
    window_start: 8
    window_end: 28
    guide_len: 20
```

When `ham_chemistry` is `custom`, `guide_quant.smk` passes individual position
flags (`--cb-start`, `--umi-start`, `--window-start`, `--guide-len`, etc.)
directly to HAM, bypassing HAM's built-in chemistry table.

### Legacy Config Compatibility

Configs that lack a `tenx_chemistry` field (pre-v0.1.2) continue to work
unchanged. The workflow detects the missing field and falls back to the
existing manual configuration paths (`simpleaf.quant.chemistry`,
`hash_matcher.chemistry`, `skip_translation`, etc.).

---

## Snakemake Rule Sequence

### 1. `reference.smk` — Index / Hash Table Construction + Whitelist Download

**Input:** Guide FASTA file (one sequence per guide).

**Output:** Method-specific index plus cached 10x whitelists:
- **simpleaf:** piscem dense index (`simpleaf index`)
- **HAM:** Python pickle hash table (`ham build-hash`)
- **Whitelists:** `3M-february-2018.txt` and `737K-august-2016.txt`, downloaded
  from 10x CDN if not already cached locally. Idempotent.

### 2. `whitelist.smk` — GEX Whitelist Extraction

**Input:** GEX expression matrix (`.h5`, `.h5ad`, or `.h5mu`).

**Output:** Cell barcodes passing QC (≥ UMI threshold, ≥ genes threshold).
Format is auto-detected by extension.

**v0.1.2:** After extraction, barcodes are automatically translated from GEX
format (TruSeq) to Feature format (Nextera) when the chemistry requires it
(dual-oligo systems: 3' v3/v4, 5' v3, multiome). For single-oligo chemistries
(5' v1/v2, soluble RT primer), translation is skipped.

### 3a. `quant.smk` — simpleaf Quantification

Piscem map-sc + alevin-fry permit-list → collate → quant. The chemistry
string (including any geometry override) and whitelist are resolved from
`tenx_chemistry`.

### 3b. `guide_quant.smk` — HAM Quantification

HAM match (with `--chemistry` + `--cb-max-hamming`) → UMI-tools directional
dedup → MEX matrix. Chemistry constants are auto-selected. For `custom`
chemistry, individual position parameters are passed as CLI flags.

### 4. `merge.smk` — Cross-Lane Merge & Barcode Translation

Vertically stacks per-lane MEX via `ham merge`, appends lane suffixes.
Feature↔GEX barcode translation is applied when the chemistry requires it
(3' chemistries: yes; 5' chemistries: no). For single-lane datasets, merge
degenerates to a no-op.

---

## Usage

```bash
cd /path/to/sgprocess
conda activate scp_analysis

# Standard invocation
snakemake --configfile config/config.yaml --cores 48

# Dry-run
snakemake --configfile config/config.yaml --cores 12 --dry-run
```

### Common Options

```
--dry-run           Preview execution plan without running
--cores N           Max parallel threads
--unlock            Remove stale locks after a crash
--rerun-incomplete  Re-run partially completed jobs
--latency-wait 60   Wait for filesystem latency (NFS-safe)
```

### Output

```
{out_dir}/
├── lane_XX/
│   └── {quant_dir}/
│       └── alevin/
│           ├── quants_mat.mtx           ← per-lane sparse matrix
│           ├── quants_mat_rows.txt      ← cell barcodes (pre-merge)
│           └── quants_mat_cols.txt      ← guide feature IDs
└── merged/
    ├── merged_matrix.mtx.gz             ← cell × guide matrix
    ├── merged_barcodes.tsv.gz           ← barcodes (GEX format, -L{NN})
    └── merged_features.tsv.gz           ← guide feature IDs
```

---

## Barcode Translation

Two capture mechanisms exist in 10x single-cell chemistries, and the workflow
handles both automatically:

**Dual-oligo systems (3' v3/v4, 5' v3, multiome):**
cs1/cs2 bead-borne primers carry two distinct oligo variants on the same bead:
- **TruSeq Read 1** — captures poly(A) mRNA → GEX library
- **Nextera Read 1** — captures Feature Barcoding → Guide library

The two variants encode different barcode formats for the same cell. The
workflow applies translation in two directions:
1. **Whitelist:** GEX-format (TruSeq) → Feature-format (Nextera), so the
   whitelist matches the guide FASTQ barcodes.
2. **Merge:** Feature-format → GEX-format, so output barcodes match mRNA
   reference conventions.

**Single-oligo systems (5' v1/v2):**
A soluble guide-specific RT primer indexes guide cDNA to the same TruSeq
barcode as mRNA. Both GEX and guide libraries share identical barcode formats.
No translation is performed at either point.

The `translation` field in `chemistry_spec.yaml` encodes this physical
difference; both translation points derive their behaviour from it
automatically.

---

## Reference Files Required

| File | Description | Used by |
|------|-------------|---------|
| Guide FASTA | Guide protospacer sequences, one per entry | All methods |
| t2g map | 2-column TSV: `guide_id → gene_id`, no header | simpleaf |
| GEX matrix (per lane) | Expression data (.h5 / .h5ad / .h5mu) | Whitelist extraction |
| 10x translation table | Feature ↔ GEX barcode mapping | Barcode translation (3' chemistries) |
| 10x barcode whitelists | 3M-feb-2018 (v3) + 737K-aug-2016 (5' v1) | Auto-downloaded by `reference.smk` |

---

## v0.1.2 Changelog — 2026-06-26

### Central Chemistry Specification (`config/chemistry_spec.yaml`)

All chemistry-dependent parameters are now defined in a single YAML file.
Previously, chemistry decisions were scattered across three files
(`scprocess_utils.py`, `mapping.py`, `quant.smk`) with if-elif chains and
hardcoded overrides. The new architecture:

- **Single source of truth:** One YAML entry per chemistry (~7 fields)
- **Single user-facing setting:** `tenx_chemistry` in config
- **Automatic derivation:** All downstream parameters (simpleaf chemistry,
  whitelist, translation, HAM chemistry, UMI length) resolved from the spec
- **Adding a chemistry:** ~7 lines in a YAML file; no code changes needed"
- **`custom` escape hatch:** For non-standard hardware, users provide their own
  parameter block; the workflow passes individual flags to HAM and simpleaf

### Snakefile Chemistry Resolution

New `_resolve_chemistry()` function in `Snakefile` loads the spec and populates
`config["_chemistry"]` before rules are included. It backfills legacy config
keys so pre-v0.1.2 configs (without `tenx_chemistry`) continue to work
unchanged.

### Automated Dual Translation Points (whitelist + merge)

**`whitelist.smk`:** After barcode extraction, a TO→FROM translation step
(TruSeq→Nextera) is conditionally executed based on the chemistry's
`translation` field. Previously, translation was only available at the merge
step and required manual configuration.

**`merge.smk`:** Translation control now preferentially reads from
`_chemistry.translation`, falling back to the legacy `skip_translation` config
key.

Both translation points are driven by the same `translation` field in
`chemistry_spec.yaml`, ensuring they can never be inconsistent.

### Custom Chemistry: Full Stack Implementation

The `custom` chemistry path described in the v0.1.1 README is now fully
implemented end-to-end:

- **`Snakefile`:** Detects `tenx_chemistry: custom`, reads `custom_chemistry`
  block, validates required fields
- **`guide_quant.smk`:** Detects `ham_chemistry: custom`, unrolls
  `hash_matcher.custom_params` into individual CLI flags
- **`ham/src/ham/cli.py`:** `--chemistry custom` now accepted; independent
  position flags added (`--cb-start`, `--cb-end`, `--umi-start`, `--umi-end`,
  `--window-start`, `--window-end`, `--guide-len`)
- **`ham/src/ham/matcher.py`:** `match_reads()` accepts optional `chem_cfg`
  dict; when `chemistry="custom"`, uses the provided dict directly instead of
  looking up the built-in `CHEMISTRY_CONFIGS` table

### Chemistry Coverage Refinements

- **`5v2`** now correctly maps to `10xv2-5p` (was incorrectly grouped with
  `3v2` under `10xv2` in some paths, which triggered piscem's known
  `chromium_v2` geometry crash). All 5' v2 parameters are identical to 5' v1
  (same R1 geometry, whitelist, capture mechanism).
- **`3v2`** (indirect GBC capture) explicitly removed from supported
  chemistries. The protospacer is absent from the guide FASTQ R2; this
  pipeline requires direct capture.
- **`chemistry_overrides`** mechanism added: users can override individual
  fields (e.g., `umi_len`, `guide_len`) on top of a standard chemistry entry
  without leaving the standard path.

### HAM Custom Chemistry: Two-Layer Architecture

The clean separation between workflow-layer and HAM-layer configuration is now
fully realised:

| Layer | Configuration | Responsibility |
|:---|:---|:---|
| Workflow (`chemistry_spec.yaml`) | chemistry → af_chemistry, whitelist, translation, ham_chemistry | What tools need to know |
| HAM (`CHEMISTRY_CONFIGS` in matcher.py) | ham_chemistry → cb/UMI/window positions, guide_len | Where to read on the physical read |

The only interaction point between the two layers is the `ham_chemistry` name
(`"10xv3"`, `"10xv2-5p"`, or `"custom"`), passed from workflow to HAM as a
single string. HAM does not need to know about whitelist choice or translation
requirements; the workflow does not need to know about CB/UMI window positions.

### Bug Fixes

- **`guide_quant.smk`:** UMI length is now read from `_chemistry.umi_len`
  rather than a hardcoded ternary on the chemistry name. This fixes incorrect
  UMI length for chemistries not named `10xv3` or `10xv2-5p` (e.g., future
  custom entries).
- **`scprocess_utils.py`:** Removed the erroneous `3v2` entry from
  `CHEMISTRY_SPEC` (indirect GBC capture is not supported).

### Backward Compatibility

All existing configs without `tenx_chemistry` continue to work. The
`_resolve_chemistry()` function detects the missing field, sets
`config["_chemistry"] = None`, and all rules fall back to their pre-v0.1.2
config paths. Configs can be migrated incrementally by adding a single
`tenx_chemistry` line.

---

## v0.1.1 Changelog — 2026-06-26

### 5' v1 Chemistry Support

Added support for 10x 5' v1.0 chemistry (soluble RT primer, TruSeq-only
barcodes, 737K whitelist) alongside the existing 3' v3.1 support.

### Whitelist Auto-Download (`rules/reference.smk`)

New `download_whitelists` rule fetches `3M-february-2018.txt` and
`737K-august-2016.txt` from 10x CDN if not cached locally. Idempotent.

### h5ad / h5mu Whitelist Extraction (`rules/whitelist.smk`)

The `extract_whitelist` rule now supports three GEX matrix formats: `.h5`
(scprocess CSC sparse, original), `.h5ad` (AnnData), `.h5mu` (MuData, reads
`mod['rna']`). Enables using published expression matrices directly as
whitelist sources.

### Automatic Barcode Translation Control (`rules/merge.smk`)

The merge rule's translation step is controlled by a `skip_translation` flag.
For 5' v1/v2 chemistries (TruSeq-only), translation is skipped.

### HAM Chemistry Support (`rules/guide_quant.smk`)

HAM guide quantification accepts a `--chemistry` flag. When set to
`10xv2-5p`, HAM switches to 5' v1 constants (10bp UMI, 19bp guide,
position 16–35 window).

### Custom Chemistry Escape Hatch (documented)

README documented the `custom` chemistry configuration pattern (full
implementation completed in v0.1.2).

### Bug Fixes

- **quant.smk:** Added `simpleaf set-paths` call before quant (fixes
  `simpleaf_info.json` missing in fresh environments).
- **quant.smk:** Index path corrected to `.../index/piscem_idx` (matches
  piscem's expected layout).
- **scprocess mapping.py:** Hardcoded geometry replaced with
  `geometry_override` parameter from `CHEMISTRY_SPEC`.
- **scprocess scprocess_utils.py:** Unified `CHEMISTRY_SPEC` table replaces
  scattered if-elif chains for chemistry resolution.
