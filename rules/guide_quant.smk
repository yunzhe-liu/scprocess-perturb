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
        ham match \
            -1 "{params.reads1}" \
            -2 "{params.reads2}" \
            -w "{input.wl}" \
            -g "{input.hash_file}" \
            -o "{params.out_dir}/hits.npz" \
            -t {threads} \
            --cb-max-hamming {params.cb_max_hamming}

        # Step 2: UMI dedup + matrix (HAM)
        echo "[2/2] UMI dedup + matrix generation..."
        ham dedup \
            -i "{params.out_dir}/hits.npz" \
            -o "{params.out_dir}/matrix" \
            -t {params.umi_threshold}

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
