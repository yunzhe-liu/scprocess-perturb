# Workflow details

This document collects method parameters, intermediate outputs, environments,
and implementation information that are not required for a first workflow run.

## Guide extraction details

### simpleaf

```yaml
simpleaf:
  af_home: /path/to/alevin-fry
  index:
    kmer_length: 15
    minimizer_length: 11
  quant:
    resolution: parsimony-gene
    use_knee: false
```

`use_knee: true` enables simpleaf UMI-knee cell calling. Knee-called cells have
not been restricted to the GEX QC cell set.

### HAM

```yaml
guide_extraction:
  method: hash_matcher

hash_matcher:
  umi_threshold: 1
  cb_max_hamming: 1
```

## Guide assignment details

### Methods and parameters

```yaml
assignment:
  guide_design: dual
  methods:
    - pgmm_em

  pgmm_em:
    umi_threshold: 1
    prob_threshold: 0.75
    workers: 16
    max_em_iter: 200

  umi_threshold:
    threshold: 3

  fishash:
    padj_cutoff: 0.05
    padj_method: GS
    min_count: 2
    refit: 10
```

Only the selected method is run. For a complete workflow,
`assignment.methods` must contain exactly one entry.

`pgmm_em` fits a two-component model independently for each guide and retains
cell-guide pairs meeting both the UMI and posterior-probability thresholds. It
also writes per-guide fit diagnostics.

`umi_threshold` applies a fixed UMI cutoff without fitting a statistical model.

`fishash` applies a one-sided Fisher test with iterative correction and FDR
control. Every candidate passing the configured filter is retained.

### Unified candidate schema

`standardize_assignment.py` converts native method output to:

```text
cell_barcode, guide_id, umi_count, rank, score, score_type, method
```

Ranks are calculated within each cell, but candidates are not truncated at a
fixed top-K during standardization.

### Per-cell summary

`make_perturbation_obs.py` additionally writes `perturbation_obs.csv`:

```text
cell_barcode, perturbation, n_guides_assigned,
assignment_score, assignment_confidence, assignment_method
```

This file is a reduced summary. Integration uses `assignments.csv` so that the
complete candidate set remains available for structure classification.

### Assignment outputs

```text
{out_dir}/assignment/{method}/
├── _raw_assignments.csv
├── assignments.csv
├── perturbation_obs.csv
├── monitoring.json
└── guide_qc.csv              pgmm_em only
```

## Multimodal integration details

Integration accepts `.h5`, `.h5ad`, and `.h5mu` expression inputs. For `.h5mu`,
the `rna` modality is used.

Important resource controls are configured under `integration`:

```yaml
integration:
  input_kind: auto
  counts_layer: counts
  # counts_source: /path/to/aligned_counts.h5ad
  # normalized_source: /path/to/aligned_normalized.h5ad
  target_sum: 10000
  max_materialized_nnz: 100000000
  stream_chunk_nnz: 50000000
  max_input_nnz: 5000000000
  max_output_gb: 250
  min_free_disk_gb: 200
  max_process_memory_gb: 192
```

For large sparse matrices, integration switches to bounded-memory streaming and
performs disk and memory preflight checks before writing.

Standalone integration can be invoked with:

```bash
python scripts/integrate_multimodal.py \
  --gex dataset=/path/to/expression.h5ad \
  --assign /path/to/assignments.csv \
  --guide-design dual \
  --guide-csv /path/to/guide_library.csv \
  --out /path/to/perturbation_adata.h5ad
```

## Perturbation-status details

Mixscape and PS are alternatives and are not run sequentially. Both consume a
temporary method input built from status-eligible cells. Outputs are standardized
to:

```text
cell_id, target_label, is_ntc, assignment_structure, method,
score, native_status, scorable, unscorable_reason
```

The stage writes:

```text
{out_dir}/perturbation_status/{method}/
├── perturbation_status.tsv.gz
├── validation.json
└── run_manifest.json
```

These files are intermediate inputs and execution records. Standardized status
fields are merged into the final H5AD by the last stage.

## Data validation and standardization details

The final stage scans `.X` and `layers["counts"]`, checks required metadata and
status coverage, writes the fixed status interface, and independently verifies
that expression and counts were preserved. It does not filter observations or
variables.

```yaml
finalization:
  memory_gb: 16
  scan_chunk_values: 10000000
  max_input_gb: 300
  min_free_disk_gb: 100
```

The internal configuration and rule retain the name `finalization`; the
user-facing workflow stage is Data validation and standardization.

## Conda environments

| Environment file | Purpose |
|---|---|
| `envs/scp_analysis.lock.yaml` | Snakemake orchestration, HAM, merge, assignment, and integration |
| `envs/simpleaf.lock.yaml` | simpleaf, piscem, and alevin-fry |
| `envs/mixscape.yaml` | Mixscape status estimation |
| `envs/ps.yaml` | PS status estimation |
| `envs/finalization.yaml` | Data validation and standardized output |

With `--use-conda`, Snakemake creates rule-specific environments as needed.

## Rule sequence

1. `reference.smk` builds the guide reference and chemistry resources.
2. `whitelist.smk` extracts or translates the cell whitelist.
3. `quant.smk` or `guide_quant.smk` quantifies guide UMIs per lane.
4. `merge.smk` combines lanes and harmonizes barcode representations.
5. `assignment.smk` produces canonical guide candidates and a per-cell summary.
6. `integration.smk` combines expression and assignment in AnnData.
7. `perturbation_status.smk` runs Mixscape or PS when selected.
8. `finalization.smk` validates and writes the final H5AD.

## Repository layout

```text
scprocess-perturb/
├── Snakefile
├── config/
│   ├── chemistry_spec.yaml
│   ├── config.yaml
│   └── groups.yaml
├── docs/
│   ├── assets/
│   ├── chemistry.md
│   └── workflow-details.md
├── rules/
│   ├── reference.smk
│   ├── whitelist.smk
│   ├── quant.smk
│   ├── guide_quant.smk
│   ├── merge.smk
│   ├── assignment.smk
│   ├── integration.smk
│   ├── perturbation_status.smk
│   └── finalization.smk
├── scripts/
│   ├── integrate_multimodal.py
│   ├── run_perturbation_status.py
│   ├── finalize_perturbation_adata.py
│   ├── validate_final_adata.py
│   └── perturbation_status/
├── envs/
└── tests/
```
