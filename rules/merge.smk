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
        mtx  = [os.path.join(config["out_dir"], "lanes", g, QUANT_OUT_SUBDIR, QUANT_MTX_FILE) for g in GROUPS],
        rows = [os.path.join(config["out_dir"], "lanes", g, QUANT_OUT_SUBDIR, QUANT_ROWS_FILE) for g in GROUPS],
        cols = [os.path.join(config["out_dir"], "lanes", g, QUANT_OUT_SUBDIR, QUANT_COLS_FILE) for g in GROUPS],
    output:
        matrix   = os.path.join(config["out_dir"], "guide_matrix", "merged_matrix.mtx.gz"),
        barcodes = os.path.join(config["out_dir"], "guide_matrix", "merged_barcodes.tsv.gz"),
        features = os.path.join(config["out_dir"], "guide_matrix", "merged_features.tsv.gz"),
    params:
        out_dir     = os.path.join(config["out_dir"], "guide_matrix"),
        prefix      = "merged",
        lane_list   = os.path.join(config["out_dir"], "guide_matrix", ".lane_list.tsv"),
        groups_repr = "[{items}]".format(items=", ".join(repr(g) for g in GROUPS)),
        result_dir  = config["out_dir"],
        quant_subdir = QUANT_OUT_SUBDIR,
        quant_mtx    = QUANT_MTX_FILE,
        # v0.1.2: prefer _chemistry.translation; fall back to legacy skip_translation
        skip_trans  = "true" if not (config.get("_chemistry") or {}).get("translation", not config.get("skip_translation", False)) else "false",
        trans_table = config.get("translation_table", "/data/yunzliu/scdata/cellranger_ref/cellranger_whitelist_translation_3v3.txt"),
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
        quant_dir = os.path.join(result_dir, 'lanes', g, quant_subdir)
        m = re.search(r'(\d+)$', g)
        suffix = f'-L{{m.group(1)}}' if m else f'-{{g}}'
        f.write(g + chr(9) + quant_dir + chr(9) + suffix + chr(10))
print('Lane list written: ' + str(len(groups)) + ' lanes')
"

        ham merge \
            --lanes "{params.lane_list}" \
            --out "{params.out_dir}" \
            --prefix "{params.prefix}"

        # Post-merge: barcode translation (Nextera→TruSeq for 3' v3 chemistry).
        # Skipped for 5' v1 (TruSeq-only) and single-lane datasets.
        if [ "{params.skip_trans}" = "true" ]; then
            echo "Barcode translation skipped (config: skip_translation=true)."
        else
            echo "Translating barcodes (sequencer format → whitelist)..."
            python3 "{config[proj_dir]}/scripts/translate_barcodes.py" \
                "{params.out_dir}/{params.prefix}_barcodes.tsv.gz" \
                --trans-table "{params.trans_table}" \
                --direction from_to
        fi
    """
