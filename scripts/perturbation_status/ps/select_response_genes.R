#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(scMAGeCK)
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) stop("Usage: select_response_genes.R <input.rds> <barcode.tsv.gz> <targets.txt> <output.tsv> <NTC label> <logFC threshold>")
input_rds <- args[[1]]
barcode_path <- args[[2]]
target_path <- args[[3]]
output_path <- args[[4]]
ntc_label <- args[[5]]
logfc_threshold <- as.numeric(args[[6]])

targets <- readLines(target_path)
barcode <- read.delim(gzfile(barcode_path), stringsAsFactors = FALSE)
barcode <- barcode[barcode$gene %in% c(targets, ntc_label), c("cell", "barcode", "gene")]
object <- readRDS(input_rds)
object <- subset(object, cells = barcode$cell)
object <- object[, barcode$cell]
DefaultAssay(object) <- "RNA"

records <- lapply(targets, function(target) {
  selected <- scMAGeCK:::select_target_gene(
    rds_object = object,
    bc_frame = barcode,
    perturb_gene = target,
    non_target_ctrl = ntc_label,
    assay_for_expcor = "RNA",
    min_gene_num = 10,
    max_gene_num = 500,
    logfc.threshold = logfc_threshold
  )$target_gene_list
  data.frame(target_label = target, response_gene = selected, stringsAsFactors = FALSE)
})
result <- do.call(rbind, records)
temporary <- paste0(output_path, ".partial")
write.table(result, temporary, sep = "\t", quote = FALSE, row.names = FALSE)
if (!file.rename(temporary, output_path)) stop("Atomic response-gene rename failed")
