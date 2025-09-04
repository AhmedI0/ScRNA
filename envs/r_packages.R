# Install the R/Bioconductor packages used by workflow/R/*.R
#   Rscript envs/r_packages.R
# Packages are installed into the first entry of .libPaths() (user library).
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 4L)

cran <- c("BiocManager", "data.table", "jsonlite", "optparse", "Matrix")
bioc <- c(
  "SingleCellExperiment", "rhdf5", "BiocParallel",
  "scDblFinder",   # doublet detection (Germain et al. 2021)
  "edgeR", "limma", "DESeq2", "apeglm",  # pseudobulk differential expression
  "speckle"        # propeller: differential cell-type proportions (Phipson et al. 2022)
)

missing_cran <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran)) install.packages(missing_cran)

missing_bioc <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc)) BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)

ok <- vapply(c(cran, bioc), requireNamespace, logical(1), quietly = TRUE)
print(ok)
cat("R", R.version.string, "| Bioconductor", as.character(BiocManager::version()), "\n")
if (!all(ok)) stop("missing packages: ", paste(names(ok)[!ok], collapse = ", "))
writeLines(capture.output(sessionInfo()), "envs/lock/R_sessionInfo.txt")
