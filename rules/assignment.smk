# ==============================================================================
# assignment.smk — Per-Cell Guide Assignment
# ==============================================================================
# Runs assignment methods on the merged (cells x guides) guide count matrix,
# standardises each method's output to a unified schema, and produces a
# per-cell perturbation call table for downstream multimodal integration.
#
# Methods supported:
#   pgmm_em        — PGMM EM (default, provides prob_gaussian confidence)
#   umi_threshold  — Simple UMI threshold (default t=3; fastest, no model)
#   fishash        — One-sided Fisher exact test with FDR control
#
# No top-K truncation happens in the pipeline: each method emits all candidates
# that pass its own filter, standardize_assignment.py ranks them into the unified
# schema, and make_perturbation_obs.py applies the guide_design reduction.
#
# Config keys used:
#   config["assignment"]["methods"]            — list of method names to run
#   config["assignment"]["guide_design"]       — "single" | "dual" | "multi"
#   config["assignment"]["guide_csv"]          — guide->gene/construct mapping CSV
#   config["assignment"][method]               — per-method parameter overrides
#
# Input:  {out_dir}/guide_matrix/merged_{matrix,barcodes,features}
# Output: {out_dir}/assignment/{method}/assignments.csv
#         {out_dir}/assignment/{method}/perturbation_obs.csv
# ==============================================================================

import os

SCRIPTS = os.path.join(config["proj_dir"], "scripts")

# ---- Method registry -------------------------------------------------------
# Each method defines:
#   runner              — shell command to execute
#   default_params      — fallback params if not specified in config
#
# All methods write a raw per-(cell, guide) CSV with every candidate that passes
# the method's own filter (no top-K truncation in the pipeline). Ranking and
# per-cell selection are deferred: standardize_assignment.py ranks all candidates
# into the unified schema, and make_perturbation_obs.py applies the guide_design
# reduction. This keeps the three methods symmetric and preserves the full
# fishash result (which is multi-guide-per-cell by construction) for downstream
# analysis, rather than baking an irreversible top-K choice into the output.

METHODS = {
    "pgmm_em": {
        "runner": f"python3 {SCRIPTS}/run_pgmm_em.py",
        "default_params": {
            "umi_threshold": 1,
            "prob_threshold": 0.75,
            "workers": 16,
            "max_em_iter": 200,
        },
    },
    "umi_threshold": {
        "runner": f"python3 {SCRIPTS}/run_umi_threshold.py",
        "default_params": {
            "threshold": 3,
        },
    },
    "fishash": {
        "runner": f"Rscript {SCRIPTS}/run_fishash.R",
        "default_params": {
            "padj_cutoff": 0.05,
            "padj_method": "GS",
            "min_count": 2,
            "refit": 10,
        },
    },
}

SELECTED = config.get("assignment", {}).get("methods", ["pgmm_em"])
GUIDE_DESIGN = config.get("assignment", {}).get("guide_design", "single")


def _method_flags(method_name):
    """Merge defaults + config overrides -> shell-safe CLI flags string."""
    mc = dict(METHODS[method_name].get("default_params", {}))
    mc.update(config.get("assignment", {}).get(method_name, {}))

    flags = []
    for key, val in mc.items():
        flag_key = key.replace("_", "-")
        if val is True:
            flags.append(f"--{flag_key}")
        elif val is False or val is None:
            pass
        elif isinstance(val, (int, float, str)):
            flags.append(f"--{flag_key} {val}")
    return " ".join(flags)


# ---- Rules -----------------------------------------------------------------

rule run_assignment:
    """Run one assignment method on merged MEX -> unified assignment CSV."""
    input:
        mtx  = os.path.join(config["out_dir"], "guide_matrix", "merged_matrix.mtx.gz"),
        bc   = os.path.join(config["out_dir"], "guide_matrix", "merged_barcodes.tsv.gz"),
        feat = os.path.join(config["out_dir"], "guide_matrix", "merged_features.tsv.gz"),
    output:
        csv = os.path.join(config["out_dir"], "assignment", "{method}",
                           "assignments.csv"),
    params:
        method       = "{method}",
        mex_dir      = os.path.join(config["out_dir"], "guide_matrix"),
        out_dir      = os.path.join(config["out_dir"], "assignment", "{method}"),
        raw_csv      = os.path.join(config["out_dir"], "assignment", "{method}",
                                    "_raw_assignments.csv"),
        runner       = lambda w: METHODS[w.method]["runner"],
        method_flags = lambda w: _method_flags(w.method),
    log:
        os.path.join(config["log_dir"], "assignment", "{method}.log"),
    threads: 1
    shell:"""
        set -euo pipefail
        mkdir -p "{params.out_dir}"
        exec &> {log}

        # Step 1: Run assignment method -> raw per-(cell, guide) CSV (all candidates)
        {params.runner} \
            --input  "{params.mex_dir}" \
            --output "{params.raw_csv}" \
            {params.method_flags}

        # Step 2: Standardize -> unified schema (ranks all candidates per cell)
        python3 {SCRIPTS}/standardize_assignment.py \
            --input  "{params.raw_csv}" \
            --output "{output.csv}" \
            --method "{params.method}"
    """


rule make_perturbation_obs:
    """Convert standardized assignment -> per-cell perturbation call."""
    input:
        csv = os.path.join(config["out_dir"], "assignment", "{method}",
                           "assignments.csv"),
    output:
        obs = os.path.join(config["out_dir"], "assignment", "{method}",
                           "perturbation_obs.csv"),
    params:
        guide_design = GUIDE_DESIGN,
        guide_csv_path = config.get("assignment", {}).get("guide_csv", ""),
    log:
        os.path.join(config["log_dir"], "assignment", "{method}_obs.log"),
    shell:"""
        python3 {SCRIPTS}/make_perturbation_obs.py \
            --input        "{input.csv}" \
            --output       "{output.obs}" \
            --guide-design "{params.guide_design}" \
            --guide-csv    "{params.guide_csv_path}" \
            --method       "{wildcards.method}"
    """
