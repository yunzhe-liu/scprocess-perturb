# Design Decisions & Deferred Questions

A running log of non-obvious design decisions and questions that were
deliberately deferred. Each entry records the context, the options, and the
current status so the choice is not silently re-litigated or forgotten.

---

## D1 — pgmm_em top-1 ranking: `prob_gaussian` vs `UMI_counts` (DEFERRED)

**Status:** deferred — pipeline unchanged (top-1 stays `prob_gaussian`-based).
**Raised:** 2026-08-12.

### Context

`run_pgmm_em.py` writes `_raw_assignments.csv` in a canonical, deterministic
order: `cell` ascending, `UMI_counts` descending, `gRNA` ascending. The `gRNA`
tie-break makes `(cell, UMI_counts, gRNA)` a total order because `(cell, gRNA)`
is unique, so the row order is fully reproducible.

This canonical **file order is UMI-based** ("first row per cell" = highest-UMI
guide). But the workflow derives each cell's actual call through a different
path: `standardize_assignment.py` ranks pgmm_em by `prob_gaussian` descending
(then `UMI_counts`), and `make_perturbation_obs.py` takes that `rank`. So the
pipeline's top-1 is **`prob_gaussian`-based**, and it never depends on the
`_raw` file's row position.

### The finding (Replogle HAM, 579,295 cells)

- `(cell, gRNA)` is unique — 0 duplicate rows.
- Cells with a tied top `UMI_counts` across ≥2 guides: **0.73%** (4,233). These
  are the cells the `gRNA` tie-break makes deterministic.
- Cells where **UMI-top-1 ≠ `prob_gaussian`-top-1: 9.22%** (53,413). For these,
  the canonical file order and the pipeline's call disagree on which guide is
  "top-1".

### Options

1. **Keep `prob_gaussian` (current).** pgmm_em's per-cell call uses the model's
   native confidence score. The `_raw` file is UMI-sorted for external consumers,
   but the pipeline call comes from `standardize_assignment.py`'s rank, not from
   row position. The UMI-order file and the prob-order call differ for ~9% of
   cells, which must be documented (this entry).
2. **Switch to UMI.** Change `standardize_assignment.py`'s pgmm_em `sort_keys`
   to UMI-first so the pipeline top-1 matches the canonical file order. This
   changes ~9% of existing pgmm_em calls and moves away from `prob_gaussian` as
   the ranking quantity.

### Why deferred

The robustness goal ("don't rely on file row order for top-1") is **already
satisfied** — the pipeline ranks explicitly via `standardize_assignment.py`.
What remains is a genuine semantic choice ("is pgmm_em's top-1 the highest
`prob_gaussian` or the highest UMI?"), which is a deliberate design decision
rather than a side effect of the sort fix. Until then, **consumers must not
infer top-1 from `_raw_assignments.csv` row position.**
