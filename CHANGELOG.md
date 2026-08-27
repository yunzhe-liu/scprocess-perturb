# Changelog

Archival version history for scprocess-perturb.

### v0.2.5 — 2026-08-12

- **Canonical `_raw_assignments.csv` order (pgmm_em):** rows are sorted
  `cell` asc, `UMI_counts` desc, `gRNA` asc before writing — a deterministic
  total order (the `gRNA` tie-break resolves the ~0.73% of cells with a tied top
  UMI). Pure reordering: no column values, retained set, or `guide_qc.csv`
  change. The redundant per-worker sort was removed.

### v0.2.4 — 2026-08-11

- **Docs/consistency cleanup:** removed the unused `--workers` flag and stale
  crispat references from `run_umi_threshold.py`; unified its monitoring output
  to `monitoring.json` (matching `pgmm_em` / `fishash`); dropped hardcoded
  fallback paths from `build_guide_hash.py`; corrected the assignment-outputs
  section in the README (log location, `_raw_assignments.csv`).
- **Removed unused scripts:** old analysis / param-sweep / merge helpers no
  longer used by the workflow.

### v0.2.3 — 2026-08-11

- **No top-K truncation in the pipeline:** `fishash` is no longer truncated by
  `postprocess_fishash.py`; its raw CSV is fed straight to
  `standardize_assignment.py`, which keeps and ranks every FDR-passing candidate.
  All three methods are now symmetric — `assignments.csv` holds all candidates,
  ranked; per-cell selection is deferred to `make_perturbation_obs.py`.
- **`rules/assignment.smk`:** removed the fishash post-processing branch and the
  `_TOP_K` map; the `run_assignment` shell is now run-method → standardize.
- **Multi guide_design fix:** `make_perturbation_obs.py` previously fell back to
  top-1 for `multi` (it silently considered only one guide). It now considers all
  ranked guides per cell (`TOP_GUIDES["multi"] = None`).
- Minor: simplified the dual-guide branch; refreshed stale docstrings.

### v0.2.2 — 2026-08-11

- **Translation table auto-derivation:** `tenx_chemistry` now fully resolves the
  RNA↔Feature translation table. New `translation_file` field in
  `chemistry_spec.yaml`; `Snakefile` derives `config["translation_table"]` into
  `{out_dir}/refs/whitelist_cache/`. Fixes Class A (3' dual-oligo) runs that
  previously failed on a fresh machine (empty/hardcoded translation path).
- **`rules/reference.smk`:** New on-demand `download_translation_table` rule
  (single-file, wildcard-based). Translation downloads removed from the orphan
  `download_whitelists` rule.
- **`rules/whitelist.smk` / `rules/merge.smk`:** Declare the translation table as
  an input so it is downloaded on demand; removed the hardcoded `/data/yunzliu`
  default from `merge.smk`.
- Class B (5') chemistries have `translation_file: null` and are unaffected.

### v0.2.1 — 2026-07-06

- **`rules/assignment.smk`:** New — runs guide assignment on merged MEX.
  Two rules: `run_assignment` (method → unified CSV) and
  `make_perturbation_obs` (CSV → per-cell perturbation call).
- **Assignment methods:** `pgmm_em` (PGMM EM, default), `umi_threshold`
  (simple UMI ≥ threshold), `fishash` (Fisher test, needs R).
- **Unified output schema:** `standardize_assignment.py` normalises each
  method to a single CSV format; `make_perturbation_obs.py` produces
  `perturbation_obs.csv`.
- **Three `guide_design` modes:** `single`, `dual`, `multi` with documented
  `guide_csv` schemas.
- **`config/config.yaml`:** New `assignment:` section.
- **`Snakefile`:** `rule all` dynamically extended with assignment targets.
- **Six assignment scripts** added to `scripts/`.

### v0.1.3 — 2026-06-30

- **5v3 chemistry corrected:** Moved from Class A to Class B (5' single-oligo).
- **3v4 whitelist corrected:** Changed to `3M-3pgex-may-2023`.
- **HAM `10xv2-5p-12umi` chemistry:** New entry for 5' v3.
- **`reference.smk` download overhaul:** teichlab mirror replaces broken CDN URLs.
- **README chemistry table** restructured to three-class format.

### v0.1.2 — 2026-06-26

- **`config/chemistry_spec.yaml`:** Central chemistry specification.
- **`Snakefile`:** `_resolve_chemistry()` populates config before rules.
- **`rules/whitelist.smk`:** TO→FROM barcode translation added.
- **HAM:** `--chemistry custom` accepted; independent position flags.
- **`chemistry_overrides`** mechanism.
- **Backward compatibility:** configs without `tenx_chemistry` unchanged.

### v0.1.1 — 2026-06-26

- **5' v1 chemistry support** (`5v1`).
- **`reference.smk`:** auto-downloads 10x whitelists.
- **`whitelist.smk`:** `.h5ad` and `.h5mu` input support.
- **Bug fixes:** `simpleaf set-paths` call; index path corrected.
