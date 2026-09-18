#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(yaml)
  library(jsonlite)
})

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 1) stop("Usage: build_target_shards.R <config.yaml>")
config <- read_yaml(arguments[[1]])
source_checkpoint <- normalizePath(config$source_checkpoint_b, mustWork = TRUE)
shard_root <- normalizePath(config$shard_root, mustWork = FALSE)
for (directory in c(shard_root, file.path(shard_root, "checkpoints"), file.path(shard_root, "configs"), file.path(shard_root, "outputs"), file.path(shard_root, "resources"))) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

atomic_write <- function(writer, path) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary), add = TRUE)
  writer(temporary)
  if (!file.rename(temporary, path)) stop("Could not atomically rename: ", path)
}

object <- readRDS(source_checkpoint)
if (!inherits(object, "Seurat")) stop("source_checkpoint_b is not a Seurat object")
if (!all(c("RNA", "PRTB") %in% names(object@assays))) stop("Checkpoint B must contain RNA and PRTB assays")
metadata <- object[[]]
required_metadata <- c(config$labels, "is_ntc", "batch_id", "assignment_structure")
if (!all(required_metadata %in% colnames(metadata))) stop("Checkpoint B is missing required metadata")
if (!all(metadata$assignment_structure %in% config$selection$eligible_assignment_structures)) stop("Checkpoint B contains ineligible cells")

nt_cells <- rownames(metadata)[as.logical(metadata$is_ntc)]
if (length(nt_cells) == 0) stop("Checkpoint B contains no eligible NTC cells")
labels <- as.character(metadata[[config$labels]])
names(labels) <- rownames(metadata)
if (!all(labels[nt_cells] == config$nt_class_name)) stop("NTC labels do not match nt_class_name")
target_cells <- rownames(metadata)[!as.logical(metadata$is_ntc)]
target_labels <- labels[target_cells]
if (anyNA(target_labels) || any(target_labels == config$nt_class_name)) stop("Invalid perturbed target labels")

target_table <- as.data.frame(table(target_labels), stringsAsFactors = FALSE)
colnames(target_table) <- c("target_label", "n_cells")
target_table$n_batch_pairs <- vapply(target_table$target_label, function(target) {
  length(unique(metadata[target_cells[target_labels == target], "batch_id"]))
}, integer(1))
cell_scale <- median(target_table$n_cells)
pair_scale <- median(target_table$n_batch_pairs)
target_table$weight <- target_table$n_cells / cell_scale + target_table$n_batch_pairs / pair_scale
target_table <- target_table[order(-target_table$weight, -target_table$n_cells, target_table$target_label), ]

n_shards <- as.integer(config$n_shards)
load_cells <- rep(0, n_shards)
load_pairs <- rep(0, n_shards)
assignment <- integer(nrow(target_table))
for (index in seq_len(nrow(target_table))) {
  score <- load_cells / sum(target_table$n_cells) + load_pairs / sum(target_table$n_batch_pairs)
  shard_id <- which.min(score)
  assignment[[index]] <- shard_id
  load_cells[[shard_id]] <- load_cells[[shard_id]] + target_table$n_cells[[index]]
  load_pairs[[shard_id]] <- load_pairs[[shard_id]] + target_table$n_batch_pairs[[index]]
}
target_table$shard_id <- assignment
target_table <- target_table[order(target_table$shard_id, target_table$target_label), ]

for (shard_id in seq_len(n_shards)) {
  shard_name <- sprintf("shard_%03d", shard_id)
  shard_targets <- target_table$target_label[target_table$shard_id == shard_id]
  keep_cells <- Cells(object)[Cells(object) %in% c(nt_cells, names(labels)[labels %in% shard_targets])]
  shard_object <- subset(object, cells = keep_cells)
  command_name <- grep("CalcPerturbSig", Command(shard_object), value = TRUE)
  if (length(command_name) != 1 || !(command_name %in% Tool(shard_object))) stop("CalcPerturbSig contract failed for ", shard_name)
  if (!all(c("RNA", "PRTB") %in% names(shard_object@assays))) stop("Assay contract failed for ", shard_name)
  if (!all(Cells(shard_object) == keep_cells)) stop("Cell order contract failed for ", shard_name)
  checkpoint_path <- file.path(shard_root, "checkpoints", paste0(shard_name, ".rds"))
  atomic_write(function(path) saveRDS(shard_object, path, compress = FALSE), checkpoint_path)

  logging_directory <- file.path(shard_root, "logs", shard_name)
  dir.create(logging_directory, recursive = TRUE, showWarnings = FALSE)
  shard_config <- list(
    module = "M03", method = "Mixscape", dataset = paste0(config$dataset_prefix, "_", shard_name),
    execution_mode = "checkpoint_only", source_h5ad = config$source_h5ad, input_dir = config$input_dir,
    output_dir = file.path(shard_root, "outputs", shard_name), resume_checkpoint = checkpoint_path,
    selection = config$selection, seed = as.integer(config$seed) + shard_id,
    logging = list(directory = logging_directory, stage_log = file.path(logging_directory, "stage_resources.jsonl")),
    checkpoints = list(enabled = FALSE, directory = file.path(shard_root, "checkpoints_unused", shard_name)),
    mixscape = config$mixscape
  )
  atomic_write(function(path) write_yaml(shard_config, path), file.path(shard_root, "configs", paste0(shard_name, ".yaml")))
  rm(shard_object)
  gc()
}

atomic_write(function(path) write.table(target_table, path, sep = "\t", quote = FALSE, row.names = FALSE), file.path(shard_root, "target_manifest.tsv"))
summary <- list(source_checkpoint_b = source_checkpoint, n_shards = n_shards, n_cells = length(Cells(object)), n_ntc = length(nt_cells), n_targets = nrow(target_table), shard_cells = as.list(load_cells + length(nt_cells)), shard_target_cells = as.list(load_cells), shard_batch_target_pairs = as.list(load_pairs), seed = as.integer(config$seed))
atomic_write(function(path) write(toJSON(summary, auto_unbox = TRUE, pretty = TRUE), path), file.path(shard_root, "shard_summary.json"))
