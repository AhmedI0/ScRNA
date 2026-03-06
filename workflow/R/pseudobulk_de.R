#!/usr/bin/env Rscript
# Pseudobulk differential expression with edgeR (quasi-likelihood) and DESeq2.
#
# Why pseudobulk? Cells from the same donor are not independent replicates.
# Testing cells as replicates (e.g. Wilcoxon per cell) inflates false
# discoveries massively (Squair et al. 2021, Nat Commun 12:5692). Summing raw
# counts per donor x cell type and using bulk RNA-seq GLMs with the donor as the
# unit of replication controls the false-discovery rate.
#
# Input directory (written by ipfatlas.pseudobulk.write_inputs):
#   <key>__counts.csv.gz   genes x pseudobulk samples (integer counts)
#   <key>__meta.csv        one row per sample: sample, condition, covariates...
#   index.csv              key, celltype
# Usage:
#   Rscript pseudobulk_de.R <indir> <outdir> <covariates|none> <case> <control> [methods]
#   e.g. Rscript pseudobulk_de.R pb_in pb_out sex IPF Control edgeR,DESeq2
suppressPackageStartupMessages({
  library(edgeR)
  library(DESeq2)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
indir <- args[1]; outdir <- args[2]
covars <- if (args[3] == "none") character(0) else strsplit(args[3], ",")[[1]]
case <- args[4]; control <- args[5]
methods <- if (length(args) >= 6) strsplit(args[6], ",")[[1]] else c("edgeR", "DESeq2")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
set.seed(0)

index <- fread(file.path(indir, "index.csv"))
summary_rows <- list()

for (k in seq_len(nrow(index))) {
  key <- index$key[k]; ct <- index$celltype[k]
  counts <- as.matrix(read.csv(gzfile(file.path(indir, paste0(key, "__counts.csv.gz"))),
                               row.names = 1, check.names = FALSE))
  meta <- read.csv(file.path(indir, paste0(key, "__meta.csv")), row.names = 1, check.names = FALSE)
  meta <- meta[colnames(counts), , drop = FALSE]
  meta$condition <- factor(meta$condition, levels = c(control, case))

  # drop covariates that are constant (or confounded) within this cell type
  use <- covars[vapply(covars, function(v) length(unique(meta[[v]])) > 1, logical(1))]
  for (v in use) meta[[v]] <- factor(meta[[v]])
  form <- as.formula(paste("~", paste(c(use, "condition"), collapse = " + ")))
  design <- model.matrix(form, data = meta)
  if (qr(design)$rank < ncol(design)) {          # e.g. site fully confounded
    use <- character(0); form <- ~ condition; design <- model.matrix(form, data = meta)
  }
  coef <- paste0("condition", case)
  n_case <- sum(meta$condition == case); n_ctrl <- sum(meta$condition == control)

  # gene filter shared by both methods (edgeR::filterByExpr, Chen et al. 2016)
  y <- DGEList(counts = counts, samples = meta)
  keep <- filterByExpr(y, design = design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- normLibSizes(y)  # TMM

  row <- data.table(celltype = ct, key = key, n_case = n_case, n_control = n_ctrl,
                    n_genes_tested = sum(keep), design = deparse(form))

  if ("edgeR" %in% methods) {
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design, robust = TRUE)
    qlf <- glmQLFTest(fit, coef = coef)
    tt <- topTags(qlf, n = Inf, sort.by = "none")$table
    res <- data.table(gene = rownames(tt), logFC = tt$logFC, logCPM = tt$logCPM, F = tt$F,
                      PValue = tt$PValue, FDR = tt$FDR)
    res[, stat := sign(logFC) * sqrt(F)]  # signed statistic for enrichment analysis
    fwrite(res[order(PValue)], file.path(outdir, paste0(key, "__edgeR.csv")))
    row[, edgeR_n_fdr05 := sum(res$FDR < 0.05)]
    row[, edgeR_n_up := sum(res$FDR < 0.05 & res$logFC > 0)]
    row[, edgeR_n_down := sum(res$FDR < 0.05 & res$logFC < 0)]
  }

  if ("DESeq2" %in% methods) {
    dds <- DESeqDataSetFromMatrix(countData = round(counts[keep, ]), colData = meta, design = form)
    dds <- DESeq(dds, quiet = TRUE)
    rn <- resultsNames(dds)
    cname <- rn[grepl(paste0("^condition_", case), rn)]
    r <- results(dds, name = cname, alpha = 0.05)
    shr <- tryCatch(lfcShrink(dds, coef = cname, type = "apeglm", quiet = TRUE),
                    error = function(e) NULL)
    res <- data.table(gene = rownames(r), baseMean = r$baseMean, log2FC = r$log2FoldChange,
                      lfcSE = r$lfcSE, stat = r$stat, pvalue = r$pvalue, padj = r$padj,
                      log2FC_shrunk = if (is.null(shr)) NA_real_ else shr$log2FoldChange)
    fwrite(res[order(pvalue)], file.path(outdir, paste0(key, "__DESeq2.csv")))
    row[, DESeq2_n_fdr05 := sum(res$padj < 0.05, na.rm = TRUE)]
  }
  summary_rows[[k]] <- row
  message(sprintf("[%d/%d] %-32s case=%d control=%d genes=%d", k, nrow(index), ct, n_case, n_ctrl, sum(keep)))
}
fwrite(rbindlist(summary_rows, fill = TRUE), file.path(outdir, "summary.csv"))
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
message("pseudobulk DE finished")
