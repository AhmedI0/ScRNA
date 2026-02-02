#!/usr/bin/env Rscript
# Differential cell-type proportions with propeller (speckle; Phipson et al. 2022).
#
# propeller computes per-sample cell-type proportions, applies a logit (or
# arcsin-sqrt) variance-stabilising transform and tests group differences with
# empirical-Bayes moderated linear models (limma). The donor is the replicate.
#
# Input CSV: one row per cell with columns celltype, donor_id, condition, and
# optional covariates (constant within donor).
# Usage: Rscript propeller.R <cells.csv> <out.csv> <case> <control> [covariates|none]
suppressPackageStartupMessages({
  library(speckle)
  library(limma)
  library(data.table)
})
args <- commandArgs(trailingOnly = TRUE)
cells <- fread(args[1]); out <- args[2]; case <- args[3]; control <- args[4]
covars <- if (length(args) >= 5 && args[5] != "none") strsplit(args[5], ",")[[1]] else character(0)

props <- getTransformedProps(clusters = cells$celltype, sample = cells$donor_id, transform = "logit")
donor <- unique(cells[, c("donor_id", "condition", covars), with = FALSE])
donor <- donor[match(colnames(props$Proportions), donor$donor_id)]
donor$condition <- factor(donor$condition, levels = c(control, case))
for (v in covars) donor[[v]] <- factor(donor[[v]])

form <- as.formula(paste("~ 0 + condition", if (length(covars)) paste("+", paste(covars, collapse = " + ")) else ""))
design <- model.matrix(form, data = donor)
colnames(design) <- make.names(colnames(design))
contr <- makeContrasts(contrasts = paste0("condition", case, " - condition", control), levels = design)
res <- propeller.ttest(props, design = design, contrasts = contr, robust = TRUE, trend = FALSE, sort = TRUE)
res <- data.table(celltype = rownames(res), res)
fwrite(res, out)
prop_wide <- as.data.frame.matrix(props$Proportions)
fwrite(data.table(celltype = rownames(prop_wide), prop_wide), sub("\\.csv$", "_proportions.csv", out))
print(res[, .(celltype, PropRatio, P.Value, FDR)])
