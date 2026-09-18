#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(scMAGeCK)
  library(Seurat)
  library(jsonlite)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2 || args[[1]] != "--config") stop("Usage: run_ps.R --config <config.yaml>")
config <- yaml::read_yaml(args[[2]])
params <- config$parameters
dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(as.integer(config$seed))

barcode <- read.delim(gzfile(config$barcode_table), check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("cell", "barcode", "gene") %in% colnames(barcode))) stop("Invalid barcode table")
if (anyDuplicated(barcode$cell)) stop("PS input requires one target-level assignment per cell")
targets <- sort(setdiff(unique(barcode$gene), config$non_target_control))
if (config$target_mode != "all") stop("Only target_mode=all is supported by the canonical pooled runner")

object <- readRDS(config$input_rds)
if (!identical(barcode$cell, colnames(object))) stop("Barcode and Seurat cell order differ")
DefaultAssay(object) <- "RNA"

old_directory <- getwd()
on.exit(setwd(old_directory), add = TRUE)
setwd(config$output_dir)
started <- Sys.time()
timing <- system.time({
  result <- scmageck_eff_estimate(
    rds_object = object,
    bc_frame = barcode,
    perturb_gene = targets,
    non_target_ctrl = config$non_target_control,
    perturb_target_gene = params$perturb_target_gene,
    scale_factor = as.numeric(params$scale_factor),
    target_gene_min = as.integer(params$target_gene_min),
    target_gene_max = as.integer(params$target_gene_max),
    assay_for_cor = params$assay_for_cor,
    subset_rds = as.logical(params$subset_rds),
    scale_score = as.logical(params$scale_score),
    lambda = as.numeric(params$lambda),
    background_correction = as.logical(params$background_correction),
    use_perturb_exp = as.logical(params$use_perturb_exp),
    logfc.threshold = as.numeric(params$logfc_threshold)
  )
})
finished <- Sys.time()

score_matrix <- result$eff_matrix
if (is.null(dim(score_matrix))) score_matrix <- matrix(score_matrix, ncol = 1, dimnames = list(names(score_matrix), targets))
if (!all(targets %in% colnames(score_matrix))) stop("PS matrix is missing targets")
cell_ids <- rownames(score_matrix)
cell_targets <- barcode$gene[match(cell_ids, barcode$cell)]
is_ntc <- cell_targets == config$non_target_control
assigned_score <- numeric(length(cell_ids))
target_rows <- which(!is_ntc)
assigned_score[target_rows] <- score_matrix[cbind(target_rows, match(cell_targets[target_rows], colnames(score_matrix)))]
metadata <- data.frame(cell_id = cell_ids, target_label = cell_targets, is_ntc = is_ntc,
                       ps_score = assigned_score, stringsAsFactors = FALSE)
if (anyNA(metadata$target_label) || any(!is.finite(metadata$ps_score))) stop("Invalid canonical PS scores")

temporary_scores <- "ps_metadata.tsv.partial"
write.table(metadata, temporary_scores, sep = "\t", quote = FALSE, row.names = FALSE)
if (!file.rename(temporary_scores, "ps_metadata.tsv")) stop("Score rename failed")
saveRDS(score_matrix, "ps_score_matrix.rds.partial", compress = FALSE)
if (!file.rename("ps_score_matrix.rds.partial", "ps_score_matrix.rds")) stop("Matrix rename failed")
saveRDS(result$target_gene_search_result, "target_gene_search_result.rds.partial", compress = FALSE)
if (!file.rename("target_gene_search_result.rds.partial", "target_gene_search_result.rds")) stop("Target-gene rename failed")

target_summary <- do.call(rbind, lapply(targets, function(target) {
  values <- metadata$ps_score[metadata$target_label == target]
  data.frame(target_label = target, cells = length(values), selected_response_genes = length(result$target_gene_search_result[[target]]$target_gene_list),
             min = min(values), median = median(values), mean = mean(values), max = max(values), sd = sd(values))
}))
write.table(target_summary, "target_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

run_summary <- list(
  module = config$module,
  method = config$method,
  dataset = config$dataset,
  execution_semantics = "official pooled multi-target estimation",
  started_at = format(started, "%Y-%m-%dT%H:%M:%S%z"),
  finished_at = format(finished, "%Y-%m-%dT%H:%M:%S%z"),
  elapsed_seconds = unname(timing[["elapsed"]]),
  user_seconds = unname(timing[["user.self"]]),
  system_seconds = unname(timing[["sys.self"]]),
  cells = nrow(metadata),
  target_cells = sum(!metadata$is_ntc),
  ntc_cells = sum(metadata$is_ntc),
  targets = length(targets),
  response_gene_union = ncol(result$optimization_matrix$tr_y),
  optimization_dimensions = lapply(result$optimization_matrix, dim),
  parameters = params,
  package_versions = list(scMAGeCK = as.character(packageVersion("scMAGeCK")), presto = as.character(packageVersion("presto")), Seurat = as.character(packageVersion("Seurat"))),
  status = "PASS"
)
jsonlite::write_json(run_summary, "run_summary.json", pretty = TRUE, auto_unbox = TRUE, null = "null")
cat(jsonlite::toJSON(run_summary, auto_unbox = TRUE, null = "null"), "\n")
