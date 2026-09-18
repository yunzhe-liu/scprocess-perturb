#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2 || args[[1]] != "--config") stop("Usage: build_ps_seurat.R --config <config.yaml>")
config <- yaml::read_yaml(args[[2]])
input_dir <- config$input_dir

counts <- Matrix::readMM(gzfile(file.path(input_dir, "counts.mtx.gz")))
counts <- as(counts, "dgCMatrix")
features <- read.delim(gzfile(file.path(input_dir, "features.tsv.gz")), header = FALSE, stringsAsFactors = FALSE)[[1]]
metadata <- read.delim(gzfile(file.path(input_dir, "metadata.tsv.gz")), row.names = 1, check.names = FALSE)

if (nrow(counts) != length(features) || ncol(counts) != nrow(metadata)) stop("Input dimensions do not agree")
if (anyDuplicated(features) || anyDuplicated(rownames(metadata))) stop("Feature or cell identifiers are duplicated")
rownames(counts) <- features
colnames(counts) <- rownames(metadata)

object <- CreateSeuratObject(counts = counts, assay = "RNA", meta.data = metadata)
object <- NormalizeData(
  object,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = as.numeric(config$normalization$scale_factor),
  verbose = FALSE
)
# scMAGeCK 0.99.1 passes a legacy multi-value layer default through
# GetAssayData(). A v3 Assay preserves official helper compatibility under
# Seurat v5 without changing the installed scMAGeCK source.
rna_counts <- LayerData(object, assay = "RNA", layer = "counts")
rna_data <- LayerData(object, assay = "RNA", layer = "data")
legacy_rna <- CreateAssayObject(counts = rna_counts)
legacy_rna <- SetAssayData(legacy_rna, layer = "data", new.data = rna_data)
object[["RNA"]] <- legacy_rna
rm(rna_counts, rna_data, legacy_rna)
gc()
DefaultAssay(object) <- "RNA"

temporary <- paste0(config$seurat_rds, ".partial")
saveRDS(object, temporary, compress = FALSE)
if (!file.rename(temporary, config$seurat_rds)) stop("Atomic RDS rename failed")
