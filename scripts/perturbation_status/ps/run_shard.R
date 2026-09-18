#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(scMAGeCK)
  library(Seurat)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) stop("Usage: run_shard.R <input.rds> <barcode.tsv.gz> <targets.txt> <genes.txt> <output_dir> <seed>")
input_rds <- args[[1]]
barcode_path <- args[[2]]
target_path <- args[[3]]
gene_path <- args[[4]]
output_dir <- args[[5]]
seed <- as.integer(args[[6]])
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
targets <- readLines(target_path)
response_genes <- readLines(gene_path)
barcode <- read.delim(gzfile(barcode_path), stringsAsFactors = FALSE)
barcode <- barcode[barcode$gene %in% c(targets, "non-targeting"), c("cell", "barcode", "gene")]
object <- readRDS(input_rds)
object <- subset(object, cells = barcode$cell)
object <- object[, barcode$cell]
DefaultAssay(object) <- "RNA"
set.seed(seed)

started <- Sys.time()
timing <- system.time({
  regression <- scmageck_lr(
    BARCODE = barcode,
    RDS = object,
    NEGCTRL = "non-targeting",
    SELECT_GENE = response_genes,
    LABEL = "targetgene_lr",
    PERMUTATION = 100,
    SAVEPATH = NULL,
    LAMBDA = 0,
    GENE_FRAC = 0,
    SLOT = "data"
  )
  tr_x <- as.matrix(regression$regression_matrix$Xmat)
  tr_y <- as.matrix(regression$regression_matrix$Ymat)
  tr_score <- as.matrix(regression[[1]][, -1, drop = FALSE])
  scores <- scMAGeCK:::scmageck_optim_core(tr_x, tr_y, tr_score, scale_factor = 3, lambda = 0)
  scores <- scores / 3
  for (target in targets) {
    maximum <- max(scores[, target])
    if (maximum >= 0.01) scores[, target] <- scores[, target] / maximum
  }
})
cells <- rownames(scores)
labels <- barcode$gene[match(cells, barcode$cell)]
target_rows <- which(labels != "non-targeting")
metadata <- data.frame(
  cell_id = cells[target_rows],
  target_label = labels[target_rows],
  ps_score = scores[cbind(target_rows, match(labels[target_rows], colnames(scores)))],
  stringsAsFactors = FALSE
)
if (anyNA(metadata$ps_score) || any(!is.finite(metadata$ps_score))) stop("Invalid shard scores")
write.table(metadata, file.path(output_dir, "ps_target_scores.tsv.partial"), sep = "\t", quote = FALSE, row.names = FALSE)
if (!file.rename(file.path(output_dir, "ps_target_scores.tsv.partial"), file.path(output_dir, "ps_target_scores.tsv"))) stop("Atomic rename failed")
summary <- list(targets = length(targets), target_cells = nrow(metadata), ntc_cells = sum(labels == "non-targeting"),
                response_genes = length(response_genes), elapsed_seconds = unname(timing[["elapsed"]]),
                started_at = format(started, "%Y-%m-%dT%H:%M:%S%z"), status = "PASS")
jsonlite::write_json(summary, file.path(output_dir, "summary.json"), pretty = TRUE, auto_unbox = TRUE)
