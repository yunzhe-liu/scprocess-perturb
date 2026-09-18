#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Seurat))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) stop("Usage: build_input_shards.R <input.rds> <barcode.tsv.gz> <target_dir> <output_dir> <shard_count>")
input_rds <- args[[1]]
barcode_path <- args[[2]]
target_dir <- args[[3]]
output_dir <- args[[4]]
shard_count <- as.integer(args[[5]])
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

barcode <- read.delim(gzfile(barcode_path), stringsAsFactors = FALSE)
object <- readRDS(input_rds)
for (shard_id in seq_len(shard_count)) {
  targets <- readLines(file.path(target_dir, sprintf("shard_%03d_targets.txt", shard_id)))
  cells <- barcode$cell[barcode$gene %in% c(targets, "non-targeting")]
  shard <- subset(object, cells = cells)
  shard <- shard[, cells]
  destination <- file.path(output_dir, sprintf("shard_%03d.rds", shard_id))
  temporary <- paste0(destination, ".partial")
  saveRDS(shard, temporary, compress = FALSE)
  if (!file.rename(temporary, destination)) stop("Atomic input-shard rename failed")
  rm(shard)
  gc()
}
