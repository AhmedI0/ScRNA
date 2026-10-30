# Running the pipeline

The analysis is a chain of 15 notebooks backed by a small Python package
(`src/ipfatlas`) and three R scripts. Every notebook reads the files written by
the previous ones, so it can be re-run on its own once its inputs exist.

```
                    GEO (GSE135893, GSE136831)            ENA FASTQ (8 libraries)
                               │                                   │
                 download_data.py (SHA-256 manifest)     linux_upstream/ (STARsolo,
                               │                          CellBender, pySCENIC)
 ┌─────────────────────────────┴───────────────────────┐           │
 │ 01 acquisition → 02 QC + doublets → 03 normalisation │           │
 │ → 04 integration (Harmony vs scVI, scIB) → 05 annot. │           │
 └──────┬──────────────┬──────────────┬─────────────────┘           │
        │              │              │                             │
 06 composition   07 pseudobulk DE   09 trajectory ──► scenic_input.loom ──┘
 (propeller,      (edgeR, DESeq2)    (PAGA, DPT,              │
  scCODA, Milo)        │             CellRank 2)              ▼
        │        08 pathways / TFs        │            14 upstream re-processing
        │        10 communication         │            15 SCENIC regulons
        ▼              ▼                  ▼
 11 validation preprocessing → 12 replication → 13 main figures
```

## 1. Set up (Windows)

```powershell
# conda env on J:, CUDA PyTorch, scverse stack, kernel "scrna-ipf"
powershell -ExecutionPolicy Bypass -File envs\setup_windows.ps1
# R/Bioconductor packages into the R library at K:\R_setup
& "K:\R_setup\R\R-4.6.1\bin\x64\Rscript.exe" envs\r_packages.R
```

Every later PowerShell session needs the environment *activated*, because MKL
on Windows loads DLLs lazily through `PATH`:

```powershell
. .\envs\activate_env.ps1
```

On Linux/macOS: `bash run_pipeline.sh` creates the environment the same way.

## 2. Run everything

```powershell
powershell -ExecutionPolicy Bypass -File run_pipeline.ps1 -SkipSetup
```

This runs the unit tests, downloads the data (resumable, size and SHA-256
checked) and executes all notebooks with papermill. Executed notebooks (with
outputs) are written to `notebooks/*.ipynb`; logs go to `logs/notebooks/`.

Single notebooks or a resume point:

```powershell
python workflow\scripts\run_notebooks.py 05 06      # selected notebooks
python workflow\scripts\run_notebooks.py --from 07  # 07 and everything after it
```

Notebooks 14 and 15 need the results of the Linux package and are skipped
automatically until `data/linux_results/results_to_return/` exists.

If the Linux machine has no NVIDIA GPU, CellBender is skipped there and the raw
droplet matrices come back instead; run it on this PC's GPU before Notebook 14
(own environment, `envs/cellbender.yml`, because CellBender pins Jupyter < 7;
`cellbender_compat.py` lets CellBender 0.3.2 write its checkpoint under PyTorch 2):

```powershell
python workflow\scripts\run_cellbender.py
```

## 3. Hardware and run time

Measured on an i7-12700K (12 cores), 128 GB RAM and an RTX 4070 (12 GB):

| notebook | main cost | time |
|---|---|---|
| 01 | reading the 220k-barcode matrix | ~1 min |
| 02 | scDblFinder (R) on 139k cells, Scrublet per library | ~25 min |
| 03 | HVG selection, PCA, UMAP | ~5 min |
| 04 | scVI training (GPU), Harmony, scIB benchmark | ~20 min |
| 05 | Leiden, CellTypist, marker scoring | ~10 min |
| 06 | propeller, scCODA (NUTS), Milo | ~15 min |
| 07 | edgeR + DESeq2 for ~35 cell types | ~10 min |
| 08–10 | decoupler, CellRank, LIANA | ~30 min |
| 11–12 | validation cohort (240k cells), scVI, CellTypist, replication | ~45 min |

Peak memory is ~40 GB (notebooks 02 and 11). With less RAM, subsample the
validation cohort in Notebook 11.

## 4. Where things go

| path | content | in git? |
|---|---|---|
| `config/` | all parameters (`config.yaml`), marker genes, dataset URLs | yes |
| `src/ipfatlas/` | reusable functions (QC, doublets, pseudobulk, plotting …) | yes |
| `workflow/R/` | scDblFinder, edgeR/DESeq2 and propeller scripts | yes |
| `notebooks/py/` | notebook sources (Jupytext percent format) | yes |
| `notebooks/*.ipynb` | executed notebooks with outputs | yes |
| `results/figures/` | every figure as SVG + 300-dpi JPG + 300-dpi PNG | yes |
| `results/tables/` | every table as CSV + XLSX (tables over 50,000 rows as `.csv.gz` only) | yes |
| `data/` | raw downloads and intermediate `.h5ad` files | no (re-created) |
| `linux_upstream/` | self-contained Linux package (FASTQ → counts, CellBender, pySCENIC) | yes |

## 5. Changing parameters

All thresholds live in `config/config.yaml` (QC MADs and floors, number of
HVGs, scVI architecture, DE thresholds …). After changing a value, re-run the
notebook that uses it and everything downstream, e.g. for QC:

```powershell
python workflow\scripts\run_notebooks.py --from 02
```

## 6. Tests

```powershell
python -m pytest -q tests
```

The tests cover the pieces that are easiest to get silently wrong: MAD outlier
flags, gene-group annotation, pseudobulk sums and minimum-cell filters, the CLR
transform and the 300-dpi figure export.
