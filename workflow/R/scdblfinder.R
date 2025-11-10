#!/usr/bin/env Rscript
# Doublet detection with scDblFinder (Germain et al. 2021, F1000Research).
#
# Input : HDF5 written by ipfatlas.io.write_sparse_h5 (genes x cells CSC counts,
#         plus obs/library_id).
# Output: CSV with barcode, scDblFinder.score, scDblFinder.class per cell.
#
# scDblFinder simulates artificial doublets from the data, then trains a
# gradient-boosted classifier (xgboost) on cell-level features of real vs
# artificial doublets. It must be run per capture (library), because doublets
# can only form between cells of the same 10x channel; `samples=` does that.
#
# Usage: Rscript scdblfinder.R <input.h5> <output.csv> [n_threads]
suppressPackageStartupMessages({
  library(rhdf5)
  library(Matrix)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(BiocParallel)
})

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 2)
h5 <- args[1]
out <- args[2]
threads <- if (length(args) >= 3) as.integer(args[3]) else 4L
set.seed(0)

shape <- h5read(h5, "shape")
counts <- new("dgCMatrix",
  i = as.integer(h5read(h5, "indices")),
  p = as.integer(h5read(h5, "indptr")),
  x = as.numeric(h5read(h5, "data")),
  Dim = as.integer(shape)
)
rownames(counts) <- make.unique(as.character(h5read(h5, "genes")))
colnames(counts) <- as.character(h5read(h5, "barcodes"))
library_id <- as.character(h5read(h5, "obs/library_id"))
H5close()
message(sprintf("loaded %d genes x %d cells from %d libraries",
                nrow(counts), ncol(counts), length(unique(library_id))))

sce <- SingleCellExperiment(assays = list(counts = counts))
sce$library_id <- library_id

# Windows has no fork(): use a socket cluster there, forked workers elsewhere
bp <- if (.Platform$OS.type == "windows") SnowParam(workers = threads, RNGseed = 0) else
  MulticoreParam(workers = threads, RNGseed = 0)

# dbr = NULL lets scDblFinder use the 10x expected rate (~1 % per 1,000 cells)
sce <- scDblFinder(sce, samples = "library_id", BPPARAM = bp)

res <- data.frame(
  barcode = colnames(sce),
  library_id = sce$library_id,
  scDblFinder_score = sce$scDblFinder.score,
  scDblFinder_class = as.character(sce$scDblFinder.class)
)
write.csv(res, out, row.names = FALSE)
tab <- table(res$library_id, res$scDblFinder_class)
print(round(100 * prop.table(tab, 1), 1))
message(sprintf("doublets called: %d / %d (%.1f %%)", sum(res$scDblFinder_class == "doublet"),
                nrow(res), 100 * mean(res$scDblFinder_class == "doublet")))
