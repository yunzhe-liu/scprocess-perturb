# ==============================================================================
# reference.smk — Guide Reference Preparation & piscem Index Building
# ==============================================================================
# Rule 1: generate_guide_reference — one-time: CSV → FASTA + t2g
# Rule 2: build_sgRNA_index        — one-time: FASTA → piscem index
# ==============================================================================

# ---------------------------------------------------------------------------
# Rule: generate_guide_reference
# Decomposes wide-format dual-sgRNA library into single-guide FASTA + t2g map.
# Skips if outputs already exist (rule is idempotent).
# ---------------------------------------------------------------------------
rule generate_guide_reference:
    input:
        csv = config["references"]["guide_csv"],
    output:
        fasta = config["references"]["guide_fasta"],
        t2g   = config["references"]["guide_t2g_2col"],
    params:
        adapter  = os.path.join(config["proj_dir"], "scripts", "feature_reference_adapter.py"),
        lock_yaml = os.path.join(config["proj_dir"], "envs", "scp_analysis.lock.yaml"),
    conda:
        os.path.join(config["proj_dir"], "envs", "scp_analysis.lock.yaml"),
    shell:"""
        set -euo pipefail
        CONDA_BASE="$HOME/software/miniconda3"
        ENV_NAME="scp_analysis"
        source "$CONDA_BASE/etc/profile.d/conda.sh"
        if [ -d "$CONDA_BASE/envs/$ENV_NAME" ]; then
            conda activate "$ENV_NAME"
        else
            echo "Creating conda environment '$ENV_NAME' from lock file..."
            conda env create -f "{params.lock_yaml}"
            conda activate "$ENV_NAME"
        fi
        # Skip if outputs already exist (adapter is idempotent; avoid redundant runs)
        if [ -f "{output.fasta}" ] && [ -f "{output.t2g}" ]; then
            echo "Guide reference files exist, skipping generation."
            exit 0
        fi
        python {params.adapter} --csv {input.csv} \
            --out-fasta {output.fasta} \
            --out-t2g {output.t2g}
    """

# ---------------------------------------------------------------------------
# Rule: build_sgRNA_index
# Builds piscem index (k=15, m=11) from guides.fasta.
# The index is globally shared and reused by all sample groups.
# ---------------------------------------------------------------------------
rule build_sgRNA_index:
    input:
        fasta = config["references"]["guide_fasta"],
    output:
        # Sentinel file for the piscem index (pre-built; skipped if exists)
        ctab = os.path.join(config["references"]["sgRNA_index_dir"], "index", "piscem_idx.ctab"),
    params:
        out_dir   = config["references"]["sgRNA_index_dir"],
        kmer      = config["simpleaf"]["index"]["kmer_length"],
        minimizer = config["simpleaf"]["index"]["minimizer_length"],
        af_home   = config["simpleaf"]["af_home"],
        lock_yaml = os.path.join(config["proj_dir"], "envs", "simpleaf.lock.yaml"),
    threads: config["resources"]["simpleaf_index_threads"]
    conda:
        os.path.join(config["proj_dir"], "envs", "simpleaf.lock.yaml"),
    shell:"""
        set -euo pipefail
        CONDA_BASE="$HOME/software/miniconda3"
        ENV_NAME="simpleaf"
        source "$CONDA_BASE/etc/profile.d/conda.sh"
        if [ -d "$CONDA_BASE/envs/$ENV_NAME" ]; then
            conda activate "$ENV_NAME"
        else
            echo "Creating conda environment '$ENV_NAME' from lock file..."
            conda env create -f "{params.lock_yaml}"
            conda activate "$ENV_NAME"
        fi
        export ALEVIN_FRY_HOME="{params.af_home}"
        mkdir -p "{params.out_dir}"

        simpleaf index \
            --output "{params.out_dir}" \
            --threads {threads} \
            --ref-seq "{input.fasta}" \
            --kmer-length {params.kmer} \
            --minimizer-length {params.minimizer}

        # simpleaf writes index to out_dir/index/,
        # but quant expects index files at parent dir; create symlinks
        for f in "{params.out_dir}"/index/piscem_idx.*; do
            base=$(basename "$f")
            [ -e "{params.out_dir}/$base" ] || ln -sf "index/$base" "{params.out_dir}/$base"
        done
        [ -f "{params.out_dir}/index/simpleaf_index.json" ] && \
            ln -sf "index/simpleaf_index.json" "{params.out_dir}/simpleaf_index.json"

        echo "Index built: {params.out_dir}"
    """

# ---------------------------------------------------------------------------
# Rule: download_whitelists
# Downloads 10x barcode whitelists and translation tables if not cached.
#
# Whitelist files (barcode lists for CB correction):
#   737K-august-2016.txt      — 3' v2, 5' v1/v2           (teichlab mirror)
#   3M-february-2018.txt      — 3' v3/v3.1, 3LT, multiome  (teichlab mirror)
#   3M-5pgex-jan-2023.txt     — 5' v3 (GEM-X)              (teichlab mirror)
#   3M-3pgex-may-2023.txt     — 3' v4 (GEM-X)              (Cell Ranger only)
#
# Translation files (RNA ↔ Feature barcode mapping for dual-oligo chemistries):
#   translation_3M-february-2018.txt     — 3' v3/v3.1       (teichlab mirror)
#   translation_3M-3pgex-may-2023.txt    — 3' v4 (GEM-X)    (teichlab mirror)
#
# All files are cached under the whitelist_dir specified in config.
# Files that already exist locally are skipped.
# ---------------------------------------------------------------------------
rule download_whitelists:
    output:
        wl_737   = os.path.join(config["references"].get("whitelist_dir", config["references"]["sgRNA_index_dir"]), "737K-august-2016.txt"),
        wl_3M    = os.path.join(config["references"].get("whitelist_dir", config["references"]["sgRNA_index_dir"]), "3M-february-2018.txt"),
        wl_5pgex = os.path.join(config["references"].get("whitelist_dir", config["references"]["sgRNA_index_dir"]), "3M-5pgex-jan-2023.txt"),
        wl_3pgex = os.path.join(config["references"].get("whitelist_dir", config["references"]["sgRNA_index_dir"]), "3M-3pgex-may-2023.txt"),
        tr_3M    = os.path.join(config["references"].get("whitelist_dir", config["references"]["sgRNA_index_dir"]), "translation_3M-february-2018.txt"),
        tr_3pgex = os.path.join(config["references"].get("whitelist_dir", config["references"]["sgRNA_index_dir"]), "translation_3M-3pgex-may-2023.txt"),
    threads: 1
    shell:"""
        set -euo pipefail
        WL_DIR=$(dirname "{output.wl_3M}")
        mkdir -p "$WL_DIR"
        MIRROR="https://teichlab.github.io/scg_lib_structs/data/10X-Genomics"

        # 737K-august-2016 (5' v1/v2)
        if [ ! -f "{output.wl_737}" ]; then
            echo "Downloading 737K-august-2016 whitelist..."
            wget -q -O "{output.wl_737}.gz" "$MIRROR/737K-august-2016.txt.gz"
            gunzip "{output.wl_737}.gz"
            echo "  Saved: {output.wl_737}"
        else
            echo "737K-august-2016 whitelist already cached."
        fi

        # 3M-february-2018 (3' v3/v3.1, 3LT, multiome)
        if [ ! -f "{output.wl_3M}" ]; then
            echo "Downloading 3M-february-2018 whitelist..."
            wget -q -O "{output.wl_3M}.gz" "$MIRROR/3M-february-2018.txt.gz"
            gunzip "{output.wl_3M}.gz"
            echo "  Saved: {output.wl_3M}"
        else
            echo "3M-february-2018 whitelist already cached."
        fi

        # 3M-5pgex-jan-2023 (5' v3, GEM-X)
        if [ ! -f "{output.wl_5pgex}" ]; then
            echo "Downloading 3M-5pgex-jan-2023 whitelist..."
            wget -q -O "{output.wl_5pgex}.gz" "$MIRROR/3M-5pgex-jan-2023.txt.gz"
            gunzip "{output.wl_5pgex}.gz"
            echo "  Saved: {output.wl_5pgex}"
        else
            echo "3M-5pgex-jan-2023 whitelist already cached."
        fi

        # 3M-3pgex-may-2023 (3' v4, GEM-X) — no public mirror
        if [ ! -f "{output.wl_3pgex}" ]; then
            echo "NOTE: 3M-3pgex-may-2023 whitelist has no public download URL."
            echo "  This file is bundled with Cell Ranger ≥ 8.0.1."
            echo "  Copy it from:"
            echo "    cellranger-8.x.x/lib/python/cellranger/barcodes/3M-3pgex-may-2023.txt.gz"
            echo "  to: {output.wl_3pgex}"
            echo "  Then re-run the workflow. (3v4 chemistry only; other chemistries unaffected.)"
        else
            echo "3M-3pgex-may-2023 whitelist already cached."
        fi

        # Translation: 3M-february-2018 (3' v3/v3.1)
        if [ ! -f "{output.tr_3M}" ]; then
            echo "Downloading translation_3M-february-2018..."
            wget -q -O "{output.tr_3M}.gz" "$MIRROR/translation_3M-february-2018.txt.gz"
            gunzip "{output.tr_3M}.gz"
            echo "  Saved: {output.tr_3M}"
        else
            echo "translation_3M-february-2018 already cached."
        fi

        # Translation: 3M-3pgex-may-2023 (3' v4, GEM-X)
        if [ ! -f "{output.tr_3pgex}" ]; then
            echo "Downloading translation_3M-3pgex-may-2023..."
            wget -q -O "{output.tr_3pgex}.gz" "$MIRROR/translation_3M-3pgex-may-2023.txt.gz"
            gunzip "{output.tr_3pgex}.gz"
            echo "  Saved: {output.tr_3pgex}"
        else
            echo "translation_3M-3pgex-may-2023 already cached."
        fi

        echo "Whitelist download complete."
    """
