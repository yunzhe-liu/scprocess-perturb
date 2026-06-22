# ==============================================================================
# merge.smk — Merge Per-Lane sgRNA Matrices into Unified Count Matrix
# ==============================================================================
# Vertically concatenates per-lane simpleaf quant output into a single
# count matrix. Cell barcodes are suffixed with lane identifiers (e.g. "-L01")
# to match mRNA count matrix conventions.
#
# Input:  Per-lane quants_mat.{mtx,rows,cols}
# Output: Single merged MEX trio — {prefix}_matrix.mtx.gz,
#         {prefix}_barcodes.tsv.gz, {prefix}_features.tsv.gz
# ==============================================================================

rule merge_matrices:
    input:
        mtx  = [os.path.join(config["out_dir"], g, QUANT_OUT_SUBDIR, QUANT_MTX_FILE) for g in GROUPS],
        rows = [os.path.join(config["out_dir"], g, QUANT_OUT_SUBDIR, QUANT_ROWS_FILE) for g in GROUPS],
        cols = [os.path.join(config["out_dir"], g, QUANT_OUT_SUBDIR, QUANT_COLS_FILE) for g in GROUPS],
    output:
        matrix   = os.path.join(config["out_dir"], "merged", "merged_matrix.mtx.gz"),
        barcodes = os.path.join(config["out_dir"], "merged", "merged_barcodes.tsv.gz"),
        features = os.path.join(config["out_dir"], "merged", "merged_features.tsv.gz"),
    params:
        out_dir     = os.path.join(config["out_dir"], "merged"),
        prefix      = "merged",
        lane_list   = os.path.join(config["out_dir"], "merged", ".lane_list.tsv"),
        groups_repr = "[{items}]".format(items=", ".join(repr(g) for g in GROUPS)),
        result_dir  = config["out_dir"],
        quant_subdir = QUANT_OUT_SUBDIR,
        quant_mtx    = QUANT_MTX_FILE,
    log:
        os.path.join(config["log_dir"], "merge", "merge.log"),
    threads: 1
    shell:"""
        set -euo pipefail
        exec &>> {log}

        # Activate analysis conda env (provides scipy, numpy for merge script)
        CONDA_BASE="$HOME/software/miniconda3"
        source "$CONDA_BASE/etc/profile.d/conda.sh"
        conda activate scp_analysis

        mkdir -p "{params.out_dir}"

        # Build lane-list TSV.
        # Suffix logic: extract trailing digits -> "-L{{NN}}"; fall back to group name.
        python3 -c "
import os, re
groups = {params.groups_repr}
result_dir = '{params.result_dir}'
quant_subdir = '{params.quant_subdir}'
with open('{params.lane_list}', 'w') as f:
    for g in groups:
        quant_dir = os.path.join(result_dir, g, quant_subdir)
        m = re.search(r'(\d+)$', g)
        suffix = f'-L{{m.group(1)}}' if m else f'-{{g}}'
        f.write(g + chr(9) + quant_dir + chr(9) + suffix + chr(10))
print('Lane list written: ' + str(len(groups)) + ' lanes')
"

        ham merge \
            --lanes "{params.lane_list}" \
            --out "{params.out_dir}" \
            --prefix "{params.prefix}"

        # Post-merge: translate sequencer format → 3v3 whitelist standard.
        # Cell Ranger internally applies this FROM→TO translation; since we
        # bypass Cell Ranger, we apply it here so output matches published data.
        echo "Translating barcodes (sequencer format → 3v3 whitelist)..."
        python3 "{config[proj_dir]}/scripts/translate_barcodes.py" \
            "{params.out_dir}/{params.prefix}_barcodes.tsv.gz" \
            --trans-table /data/yunzliu/scdata/cellranger_ref/cellranger_whitelist_translation_3v3.txt \
            --direction from_to
    """
