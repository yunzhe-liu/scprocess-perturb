#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Usage: freeze_modeled_genes.R <input.rds> <raw_union.txt> <output.txt>")
object <- readRDS(args[[1]])
raw_union <- readLines(args[[2]])
counts <- GetAssayData(object, assay = "RNA", layer = "counts")
genes <- rownames(counts)[rownames(counts) %in% raw_union]
keep <- Matrix::rowSums(counts[genes, , drop = FALSE] != 0) >= ncol(counts) * 0.01
modeled <- genes[keep]
writeLines(modeled, args[[3]])
cat(length(modeled), "modeled genes\n")
