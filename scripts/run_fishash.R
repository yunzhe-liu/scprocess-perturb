#!/usr/bin/env Rscript
#
# run_fishash.R — fishash runner (standardized CLI)
#
# Reads a merged MEX trio, runs fishash(), unpacks assignments to
# standard CSV format (cell, gRNA, UMI_counts, log_pval, odds_ratio_regularized),
# and writes monitoring JSON alongside the CSV.
#
# Usage:
#   Rscript run_fishash.R \
#       --input  merged/ \
#       --output path/to/assignments.csv \
#       --padj-cutoff 0.05 \
#       --padj-method GS \
#       --min-count 2 \
#       --refit 10

# ── Parse args ────────────────────────────────────────────────────────
library(optparse)

option_list <- list(
    make_option("--input", type="character", help="MEX trio directory"),
    make_option("--output", type="character", help="Output CSV path"),
    make_option("--padj-cutoff", type="double", default=0.05),
    make_option("--padj-method", type="character", default="GS"),
    make_option("--min-count", type="integer", default=2),
    make_option("--min-frac", type="double", default=0),
    make_option("--refit", type="integer", default=10),
    make_option("--exclude-empty", type="logical", default=TRUE)
)
parser <- OptionParser(option_list=option_list)
args <- parse_args(parser)

if (is.null(args$input) || is.null(args$output)) {
    print_help(parser)
    quit(status=1)
}

mex_dir  <- sub("/$", "", args$input)
out_csv  <- args$output
out_dir  <- dirname(out_csv)

suppressPackageStartupMessages({
    library(Matrix)
    library(S4Vectors)
    library(SummarizedExperiment)
    library(fishash)
    library(jsonlite)
})

# ── Monitoring helpers ────────────────────────────────────────────────
T0 <- proc.time()["elapsed"]
monitor <- list(
    method    = "fishash",
    parameters = list(
        padj_cutoff    = args$`padj-cutoff`,
        padj_method    = args$`padj-method`,
        min_count      = args$`min-count`,
        min_frac       = args$`min-frac`,
        refit          = args$refit,
        exclude_empty  = args$`exclude-empty`
    ),
    stages = list(),
    system  = list(
        hostname        = Sys.info()["nodename"],
        r_version       = paste(R.version$major, R.version$minor, sep="."),
        fishash_version = as.character(packageVersion("fishash")),
        cpu_count       = as.integer(parallel::detectCores())
    )
)

dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

cat(sprintf("%s\n", paste(rep("=", 70), collapse="")))
cat(sprintf("fishash — One-sided Fisher exact test + Simpson paradox correction\n"))
cat(sprintf("  MEX dir:       %s\n", mex_dir))
cat(sprintf("  Output:        %s\n", out_csv))
cat(sprintf("  padj_cutoff:   %.2f\n", args$`padj-cutoff`))
cat(sprintf("  padj_method:   %s\n",  args$`padj-method`))
cat(sprintf("  min_count:     %d\n",  args$`min-count`))
cat(sprintf("  refit:         %d\n",  args$refit))
cat(sprintf("%s\n\n", paste(rep("=", 70), collapse="")))

# ── Stage 1: Load MEX ─────────────────────────────────────────────────
cat(sprintf("[%.0fs] Stage 1/3: Loading MEX trio …\n", proc.time()["elapsed"] - T0))
t1 <- proc.time()["elapsed"]

mtx_path     <- file.path(mex_dir, "merged_matrix.mtx.gz")
barcodes_path <- file.path(mex_dir, "merged_barcodes.tsv.gz")
features_path <- file.path(mex_dir, "merged_features.tsv.gz")

counts <- readMM(gzfile(mtx_path))
counts <- t(counts)
counts <- as(counts, "CsparseMatrix")

barcodes <- read.table(barcodes_path, header=FALSE, stringsAsFactors=FALSE)[,1]
features <- read.table(features_path, header=FALSE, stringsAsFactors=FALSE)[,1]

dimnames(counts) <- list(features, barcodes)

n_cells  <- ncol(counts)
n_guides <- nrow(counts)
nnz      <- length(counts@x)
sparsity <- round((1 - nnz / (n_cells * n_guides)) * 100, 2)

load_t <- proc.time()["elapsed"] - t1
monitor$stages$load_mex <- list(
    wall_s    = round(load_t, 2),
    ncells    = n_cells,
    nguides   = n_guides,
    nnz       = nnz,
    sparsity_pct = sparsity
)

cat(sprintf("  Cells: %d  Guides: %d  NNZ: %d  Sparsity: %.2f%%  [%.1fs]\n",
            n_cells, n_guides, nnz, sparsity, load_t))

# ── Stage 2: Run fishash ───────────────────────────────────────────────
cat(sprintf("\n[%.0fs] Stage 2/3: Running fishash …\n", proc.time()["elapsed"] - T0))
t2 <- proc.time()["elapsed"]

res <- fishash(
    counts,
    padj_cutoff   = args$`padj-cutoff`,
    padj_method   = args$`padj-method`,
    min_count     = args$`min-count`,
    min_frac      = args$`min-frac`,
    refit         = args$refit,
    exclude_empty = args$`exclude-empty`
)

fishash_t <- proc.time()["elapsed"] - t2
n_iter    <- metadata(res)$num_iter
logp_cut  <- metadata(res)$log_pval_cutoff

monitor$stages$fishash <- list(
    wall_s          = round(fishash_t, 2),
    wall_min        = round(fishash_t / 60, 2),
    num_iter        = as.integer(n_iter),
    log_pval_cutoff = logp_cut
)

cat(sprintf("  Wall: %.0fs (%.1f min)  Iterations: %d\n",
            fishash_t, fishash_t/60, n_iter))

# ── Stage 3: Unpack assignments → standard CSV ─────────────────────────
cat(sprintf("\n[%.0fs] Stage 3/3: Unpacking assignments (with scores) …\n",
            proc.time()["elapsed"] - T0))
t3 <- proc.time()["elapsed"]

demux_type <- colData(res)$demux_type
assignment <- colData(res)$assignment

n_singlet  <- sum(demux_type == "singlet")
n_doublet  <- sum(demux_type == "doublet")
n_unknown  <- sum(demux_type == "unknown")

logpval_full  <- assay(res, "log_pval")
or_full       <- assay(res, "odds_ratio_regularized")

mat_assigned <- assay(res, "assigned")
mat_assigned <- as(mat_assigned, "TsparseMatrix")

n_assigned <- length(mat_assigned@i)
if (n_assigned > 0) {
    guide_names <- rownames(counts)
    cell_names  <- colnames(counts)

    idx          <- cbind(mat_assigned@i + 1, mat_assigned@j + 1)
    cells_vec    <- cell_names[idx[, 2]]
    guides_vec   <- guide_names[idx[, 1]]
    umi_vec      <- as.integer(counts[idx])
    logpval_vec  <- logpval_full[idx]
    or_vec       <- or_full[idx]

    df <- data.frame(
        cell                    = cells_vec,
        gRNA                    = guides_vec,
        UMI_counts              = umi_vec,
        log_pval                = round(logpval_vec, 6),
        odds_ratio_regularized  = round(or_vec, 6),
        stringsAsFactors        = FALSE
    )
} else {
    df <- data.frame(
        cell                    = character(),
        gRNA                    = character(),
        UMI_counts              = integer(),
        log_pval                = numeric(),
        odds_ratio_regularized  = numeric(),
        stringsAsFactors        = FALSE
    )
}

df <- df[order(df$cell, -df$UMI_counts), ]

write.csv(df, out_csv, row.names=FALSE, quote=FALSE)

unpack_t <- proc.time()["elapsed"] - t3
out_size_mb <- file.info(out_csv)$size / (1024^2)
n_assign <- nrow(df)
n_cells_assigned <- length(unique(df$cell))

monitor$stages$unpack <- list(
    wall_s         = round(unpack_t, 2),
    output_csv     = out_csv,
    output_size_mb = round(out_size_mb, 1)
)

# ── Summary ────────────────────────────────────────────────────────────
tt <- proc.time()["elapsed"] - T0
gpc <- table(df$cell)

monitor$summary <- list(
    total_wall_s          = round(tt, 2),
    total_wall_min        = round(tt / 60, 2),
    total_assignments     = n_assign,
    cells_assigned        = n_cells_assigned,
    cells_total           = n_cells,
    cell_recovery_pct     = round(n_cells_assigned / n_cells * 100, 2),
    guides_detected       = length(unique(df$gRNA)),
    guides_per_cell_median = as.numeric(median(gpc)),
    guides_per_cell_mean   = round(mean(gpc), 2),
    guides_per_cell_max    = as.integer(max(gpc)),
    cells_1_guide         = as.integer(sum(gpc == 1)),
    cells_2_guides        = as.integer(sum(gpc == 2)),
    cells_ge3_guides      = as.integer(sum(gpc >= 3)),
    cells_singlet         = n_singlet,
    cells_doublet         = n_doublet,
    cells_unknown         = n_unknown,
    fishash_iterations    = as.integer(n_iter),
    umi_median            = as.numeric(median(df$UMI_counts))
)

# ── Save monitoring ────────────────────────────────────────────────────
mon_json <- file.path(out_dir, "monitoring.json")
write_json(monitor, mon_json, pretty=TRUE, auto_unbox=TRUE, digits=10)

cat(sprintf("\n%s\n", paste(rep("=", 70), collapse="")))
cat(sprintf("fishash DONE\n"))
cat(sprintf("%s\n", paste(rep("=", 70), collapse="")))
cat(sprintf("  Assignments:        %12d\n", n_assign))
cat(sprintf("  Cells assigned:     %12d  (%.1f%%)\n", n_cells_assigned,
            n_cells_assigned / n_cells * 100))
cat(sprintf("  Guides detected:    %12d\n", monitor$summary$guides_detected))
cat(sprintf("  Guides/cell:        median=%.0f  mean=%.2f  max=%d\n",
            median(gpc), mean(gpc), max(gpc)))
cat(sprintf("  Singlet: %10d  Doublet: %10d  Unknown: %10d\n",
            n_singlet, n_doublet, n_unknown))
cat(sprintf("  1g: %d  2g: %d  >=3g: %d\n",
            sum(gpc==1), sum(gpc==2), sum(gpc>=3)))
cat(sprintf("  Fishash iterations: %d / %d\n", n_iter, args$refit))
cat(sprintf("  Total wall time:    %.0fs (%.1f min)\n", tt, tt/60))
cat(sprintf("  Output:             %s  (%.1f MB)\n", out_csv, out_size_mb))
cat(sprintf("  Monitoring:         %s\n", mon_json))
cat(sprintf("%s\n\n", paste(rep("=", 70), collapse="")))
