# ==============================================================================
# guide_quant.smk — Hash-Accelerated Guide Matcher Quantification
# ==============================================================================
# Alternative to simpleaf quant. Uses window-restricted hash lookup with
# Hamming-distance error tolerance, replacing piscem k-mer pseudoalignment.
#
# Pipeline per group:
#   1. ham match  →  hits.npz (binary int32 array)
#   2. ham dedup  → MEX trio (matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz)
#
# Requires:
#   - Pre-built guide hash table (build_guide_hash.py, one-time)
#   - GEX cell barcode whitelist (from whitelist.smk)
# ==============================================================================

# ---------------------------------------------------------------------------
# Rule: build_guide_hash — one-time guide hash table construction
# ---------------------------------------------------------------------------
rule build_guide_hash:
    input:
        fasta = config["references"]["guide_fasta"],
    output:
        hash_file = config["references"]["guide_hash"],
    params:
        script = os.path.join(config["proj_dir"], "scripts", "build_guide_hash.py"),
    log:
        os.path.join(config["log_dir"], "hash_build.log"),
    threads: 1
    shell:"""
        set -euo pipefail
        exec &>> {log}
        python3 {params.script} {input.fasta} {output.hash_file}
    """


# ---------------------------------------------------------------------------
# Rule: hash_guide_quant — per-group guide quantification
# ---------------------------------------------------------------------------
rule hash_guide_quant:
    input:
        r1_files = lambda wildcards: _resolve_reads(wildcards, "r1"),
        r2_files = lambda wildcards: _resolve_reads(wildcards, "r2"),
        wl       = os.path.join(config["out_dir"], "{group}", "barcode_whitelist_noheader.txt"),
        hash_file = config["references"]["guide_hash"],
    output:
        matrix   = os.path.join(config["out_dir"], "{group}", "guide_quant", "matrix", "matrix.mtx.gz"),
        barcodes = os.path.join(config["out_dir"], "{group}", "guide_quant", "matrix", "barcodes.tsv.gz"),
        features = os.path.join(config["out_dir"], "{group}", "guide_quant", "matrix", "features.tsv.gz"),
    params:
        out_dir        = os.path.join(config["out_dir"], "{group}", "guide_quant"),
        umi_threshold  = config.get("hash_matcher", {}).get("umi_threshold", 1),
        cb_max_hamming = config.get("hash_matcher", {}).get("cb_max_hamming", 1),
        # v0.1.2: chemistry resolved from _chemistry spec preferentially,
        # falling back to the raw hash_matcher.chemistry config value.
        ham_chemistry  = (config.get("_chemistry") or {}).get("ham_chemistry") or config.get("hash_matcher", {}).get("chemistry", "10xv3"),
        umi_len        = (config.get("_chemistry") or {}).get("umi_len") or (10 if config.get("hash_matcher", {}).get("chemistry", "10xv3") == "10xv2-5p" else 12),
        # Custom chemistry support: when ham_chemistry == "custom", unroll
        # hash_matcher.custom_params as individual CLI flags.
        is_custom      = ((config.get("_chemistry") or {}).get("ham_chemistry") == "custom") or (config.get("hash_matcher", {}).get("chemistry") == "custom"),
        custom_cb_start    = config.get("hash_matcher", {}).get("custom_params", {}).get("cb_start", 0),
        custom_cb_end      = config.get("hash_matcher", {}).get("custom_params", {}).get("cb_end", 16),
        custom_umi_start   = config.get("hash_matcher", {}).get("custom_params", {}).get("umi_start", 16),
        custom_umi_end     = config.get("hash_matcher", {}).get("custom_params", {}).get("umi_end", 28),
        custom_window_start = config.get("hash_matcher", {}).get("custom_params", {}).get("window_start", 28),
        custom_window_end  = config.get("hash_matcher", {}).get("custom_params", {}).get("window_end", 54),
        custom_guide_len   = config.get("hash_matcher", {}).get("custom_params", {}).get("guide_len", 20),
        reads1 = lambda wildcards: _join_fastq(_resolve_reads(wildcards, "r1")),
        reads2 = lambda wildcards: _join_fastq(_resolve_reads(wildcards, "r2")),
    log:
        os.path.join(config["log_dir"], "guide_quant", "{group}.log"),
    benchmark:
        os.path.join(config["log_dir"], "benchmark", "guide_quant_{group}.tsv"),
    threads: config.get("resources", {}).get("hash_quant_threads", 4)
    shell:"""
        set -euo pipefail
        exec &>> {log}

        echo "=== HAM Guide Quant: {wildcards.group} ==="
        mkdir -p "{params.out_dir}"

        # Step 1: Match reads (HAM — integer encoding + numpy + multi-process)
        echo "[1/2] Matching reads (HAM, {threads} threads)..."
        IS_CUSTOM="{params.is_custom}"
        if [ "$IS_CUSTOM" = "True" ]; then
            echo "  Mode: custom chemistry"
            echo "  cb_start={params.custom_cb_start} cb_end={params.custom_cb_end}"
            echo "  umi_start={params.custom_umi_start} umi_end={params.custom_umi_end}"
            echo "  window_start={params.custom_window_start} window_end={params.custom_window_end}"
            echo "  guide_len={params.custom_guide_len}"
            ham match \
                -1 "{params.reads1}" \
                -2 "{params.reads2}" \
                -w "{input.wl}" \
                -g "{input.hash_file}" \
                -o "{params.out_dir}/hits.npz" \
                -t {threads} \
                --cb-max-hamming {params.cb_max_hamming} \
                --chemistry custom \
                --cb-start {params.custom_cb_start} \
                --cb-end {params.custom_cb_end} \
                --umi-start {params.custom_umi_start} \
                --umi-end {params.custom_umi_end} \
                --window-start {params.custom_window_start} \
                --window-end {params.custom_window_end} \
                --guide-len {params.custom_guide_len}
        else
            echo "  Chemistry: {params.ham_chemistry}"
            ham match \
                -1 "{params.reads1}" \
                -2 "{params.reads2}" \
                -w "{input.wl}" \
                -g "{input.hash_file}" \
                -o "{params.out_dir}/hits.npz" \
                -t {threads} \
                --cb-max-hamming {params.cb_max_hamming} \
                --chemistry {params.ham_chemistry}
        fi

        # Step 2: UMI dedup + matrix (HAM)
        echo "[2/2] UMI dedup + matrix generation..."
        ham dedup \
            -i "{params.out_dir}/hits.npz" \
            -o "{params.out_dir}/matrix" \
            -t {params.umi_threshold} \
            --umi-len {params.umi_len}

        echo "  Done. Output: {params.out_dir}/matrix/"
    """


# ── Helpers (shared with quant.smk) ──
def _join_fastq(files: list) -> str:
    """Join FASTQ file paths with comma for tools accepting multi-file input."""
    return ",".join(files)


def _resolve_reads(wildcards, read: str) -> list:
    """Return trimmed FASTQ if available, otherwise raw FASTQ."""
    gcfg = GROUPS[wildcards.group]
    trimmed = gcfg.get(f"_sgRNA_{read}_trimmed", [])
    if trimmed:
        return trimmed
    return gcfg.get(f"_sgRNA_{read}", [])
