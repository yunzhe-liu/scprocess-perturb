#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: run_mixscape.R <config.yaml>")
}

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(future)
  library(jsonlite)
  library(yaml)
})

config <- yaml::read_yaml(args[[1]])
input_dir <- normalizePath(config$input_dir, mustWork = TRUE)
output_dir <- normalizePath(config$output_dir, mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
params <- config$mixscape

plan("sequential")
seed <- if (is.null(config$seed)) 20260902 else as.integer(config$seed)
set.seed(seed)

resource_dir <- if (!is.null(config$logging$directory)) {
  normalizePath(config$logging$directory, mustWork = FALSE)
} else {
  file.path(output_dir, "resources")
}
dir.create(resource_dir, recursive = TRUE, showWarnings = FALSE)
stage_log <- if (!is.null(config$logging$stage_log)) {
  config$logging$stage_log
} else {
  file.path(resource_dir, "stage_resources.jsonl")
}
checkpoint_dir <- if (!is.null(config$checkpoints$directory)) {
  config$checkpoints$directory
} else {
  file.path(output_dir, "checkpoints")
}
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
resume_checkpoint <- config$resume_checkpoint

read_proc_status <- function() {
  lines <- tryCatch(readLines("/proc/self/status", warn = FALSE), error = function(e) character())
  out <- list()
  for (line in lines) {
    key <- sub(":.*", "", line)
    if (key %in% c("VmRSS", "VmSize", "VmPeak", "VmHWM")) {
      value <- as.numeric(sub(".*:[[:space:]]*([0-9]+).*", "\\1", line)) * 1024
      out[[key]] <- value
    }
  }
  out
}

read_meminfo <- function() {
  lines <- tryCatch(readLines("/proc/meminfo", warn = FALSE), error = function(e) character())
  wanted <- c("MemTotal", "MemAvailable", "MemFree", "Cached", "Slab", "AnonPages", "SwapTotal", "SwapFree")
  out <- list()
  for (line in lines) {
    key <- sub(":.*", "", line)
    if (key %in% wanted) {
      out[[key]] <- as.numeric(sub(".*:[[:space:]]*([0-9]+).*", "\\1", line)) * 1024
    }
  }
  out
}

record_stage <- function(stage, event, object = NULL, extra = list()) {
  gc_snapshot <- gc()
  object_sizes <- list()
  if (!is.null(object)) {
    object_sizes$seurat_object_bytes <- as.numeric(object.size(object))
    object_sizes$assays <- lapply(names(object@assays), function(name) {
      list(name = name, bytes = as.numeric(object.size(object[[name]])))
    })
  }
  record <- c(
    list(
      timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE),
      epoch = as.numeric(Sys.time()),
      stage = stage,
      event = event,
      process = read_proc_status(),
      meminfo = read_meminfo(),
      gc = unname(as.data.frame(gc_snapshot)),
      object_sizes = object_sizes
    ),
    extra
  )
  write(toJSON(record, auto_unbox = TRUE, null = "null"), file = stage_log, append = TRUE)
}

save_checkpoint <- function(label, object) {
  gc()
  record_stage(paste0("checkpoint_", label), "start", object)
  final_path <- file.path(checkpoint_dir, paste0("checkpoint_", label, ".rds"))
  temporary_path <- paste0(final_path, ".tmp-", Sys.getpid())
  saveRDS(object, temporary_path, compress = FALSE)
  if (!file.rename(temporary_path, final_path)) {
    stop("Could not atomically rename checkpoint: ", final_path)
  }
  record_stage(paste0("checkpoint_", label), "end", object, list(path = final_path))
}

if (is.null(resume_checkpoint)) {
  record_stage("input_read", "start")
  counts <- readMM(gzfile(file.path(input_dir, "counts.mtx.gz")))
  counts <- as(counts, "dgCMatrix")
  features <- read.delim(
    gzfile(file.path(input_dir, "features.tsv.gz")),
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE,
    col.names = c("gene_id", "gene_name", "feature_type")
  )
  metadata <- read.delim(
    gzfile(file.path(input_dir, "metadata.tsv.gz")),
    header = TRUE,
    sep = "\t",
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  rownames(counts) <- make.unique(features$gene_id)
  colnames(counts) <- rownames(metadata)
  if (!identical(colnames(counts), rownames(metadata))) {
    stop("Counts columns and metadata row names are not aligned")
  }
  required <- c("gene", "target_label", "is_ntc", "guide_id", "batch_id")
  if (!all(required %in% colnames(metadata))) {
    stop("Required Mixscape metadata columns are missing")
  }
  metadata$gene <- as.character(metadata$gene)
  metadata$target_label <- as.character(metadata$target_label)
  metadata$is_ntc <- as.logical(metadata$is_ntc)
  metadata$batch_id <- as.character(metadata$batch_id)

  object <- CreateSeuratObject(
    counts = counts,
    meta.data = metadata,
    assay = "RNA",
    project = paste0(config$dataset, "_mixscape"),
    min.cells = 0,
    min.features = 0
  )
  rm(counts, features, metadata)
  gc()
  record_stage("input_read", "end", object)

  record_stage("NormalizeData", "start", object)
  object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  record_stage("NormalizeData", "end", object)

  record_stage("FindVariableFeatures", "start", object)
  object <- FindVariableFeatures(
    object,
    selection.method = "vst",
    nfeatures = min(as.integer(params$nfeatures), nrow(object)),
    verbose = FALSE
  )
  record_stage("FindVariableFeatures", "end", object)

  record_stage("ScaleData", "start", object)
  object <- ScaleData(object, features = VariableFeatures(object), verbose = FALSE)
  record_stage("ScaleData", "end", object)

  record_stage("RunPCA", "start", object)
  npcs <- min(as.integer(params$npcs), ncol(object) - 1, length(VariableFeatures(object)))
  object <- RunPCA(object, features = VariableFeatures(object), npcs = npcs, verbose = FALSE)
  record_stage("RunPCA", "end", object)

  if ("scale.data" %in% Layers(object[["RNA"]])) {
    LayerData(object[["RNA"]], layer = "scale.data") <- NULL
  }
  gc()
  record_stage("RunPCA_cleanup", "end", object)
  if (isTRUE(config$checkpoints$enabled %||% TRUE)) {
    save_checkpoint("A", object)
  }

  record_stage("CalcPerturbSig", "start", object)
  object <- CalcPerturbSig(
    object = object,
    assay = params$input_assay,
    slot = params$signature_slot,
    features = VariableFeatures(object),
    gd.class = params$guide_class,
    nt.cell.class = params$nt_cell_class,
    split.by = params$split_by,
    num.neighbors = as.integer(params$num_neighbors),
    reduction = params$reduction,
    ndims = npcs,
    new.assay.name = params$signature_assay,
    verbose = FALSE
  )
  record_stage("CalcPerturbSig", "end", object)
  gc()
  record_stage("CalcPerturbSig_cleanup", "end", object)
  if (isTRUE(config$checkpoints$enabled %||% TRUE)) {
    save_checkpoint("B", object)
  }
} else {
  resume_checkpoint <- normalizePath(resume_checkpoint, mustWork = TRUE)
  record_stage("resume_load", "start")
  object <- readRDS(resume_checkpoint)
  if (!inherits(object, "Seurat")) {
    stop("Resume checkpoint is not a Seurat object: ", resume_checkpoint)
  }
  if (!(params$signature_assay %in% names(object@assays))) {
    stop("Resume checkpoint is missing assay: ", params$signature_assay)
  }
  record_stage("resume_load", "end", object, list(path = resume_checkpoint))
  gc()
}

if (isTRUE(config$stop_after_checkpoint_b)) {
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
  quit(save = "no", status = 0)
}

record_stage("RunMixscape", "start", object)
object <- RunMixscape(
  object = object,
  assay = params$signature_assay,
  slot = params$classification_slot,
  labels = params$labels,
  nt.class.name = params$nt_class_name,
  new.class.name = params$new_class_name,
  min.de.genes = as.integer(params$min_de_genes),
  min.cells = as.integer(params$min_cells),
  de.assay = params$de_assay,
  logfc.threshold = as.numeric(params$logfc_threshold),
  iter.num = as.integer(params$iter_num),
  split.by = params$split_by,
  prtb.type = params$prtb_type,
  verbose = FALSE
)
record_stage("RunMixscape", "end", object)

scores <- object[[]]
write.table(
  scores,
  file = gzfile(file.path(output_dir, "mixscape_metadata.tsv.gz")),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)
saveRDS(object, file.path(output_dir, "mixscape_result.rds"), compress = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
