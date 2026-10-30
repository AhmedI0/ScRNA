# Methods

Written in the style of a manuscript methods section. Every parameter is read from
[`config/config.yaml`](../config/config.yaml) unless a notebook is cited; numbers of
cells and donors are those produced by the notebooks in `notebooks/`. References are
listed in [references.md](references.md).

## Data

**Discovery cohort.** Processed count matrices, barcodes, features and cell metadata
of GSE135893 (Habermann *et al.* 2020) were downloaded from GEO. The cohort contains
220,213 barcodes from 43 libraries of 30 donors (10 controls, 12 IPF, 3 NSIP,
2 chronic hypersensitivity pneumonitis, 2 sarcoidosis, 1 unclassified ILD) profiled
with 10x Genomics 5′ v1 chemistry at three sites. Site is partly confounded with
diagnosis: one site contributed controls only and one disease samples only
([data_provenance.md](data_provenance.md)). Non-IPF ILDs were kept for clustering and
annotation and excluded from IPF-versus-control statistics.

**Validation cohort.** GSE136831 (Adams *et al.* 2020; 10x 3′ v2) was downloaded as the
authors' filtered matrix with cell metadata. Control and IPF donors were retained and
cells annotated by the authors as multiplets were removed (COPD donors were not used).

Downloads are resumable, size-checked against the GEO listing and recorded with SHA-256
checksums (`workflow/scripts/download_data.py`).

## Quality control (Notebook 02)

QC metrics were computed with scanpy (total UMIs, detected genes, percentage of
mitochondrial and ribosomal counts, fraction of counts in the 20 most expressed genes).
Outliers were defined per library as cells more than 4 median absolute deviations (MADs)
from the library median of log(total UMIs) or log(detected genes). The mitochondrial
threshold was the library median + 3 MADs, bounded to 15–25 %. Hard floors of
500 genes and 1,000 UMIs were applied after a gene-bin analysis showed that barcodes
with 300–500 genes resembled ambient-dominated droplets. The top-20-gene fraction was
reported but not used for filtering, because it removed plasma, goblet, club and AT2
cells, whose transcriptomes are dominated by a few secretory transcripts. Every QC
criterion was checked against the authors' cell-type labels before it was adopted.
Genes detected in fewer than 10 cells and libraries with fewer than 200 cells were
removed. Donor sex was inferred from *XIST* and chromosome-Y gene expression and
checked against the metadata.

Doublets were called per library with scDblFinder (R/Bioconductor; cluster-based
artificial doublets, expected rate from the number of cells captured). Scrublet (expected
rate 0.06, automatic threshold) was run as a second method and its concordance with
scDblFinder is reported. scDblFinder calls were used to remove doublets.

## Normalisation, feature selection and dimensionality reduction (Notebook 03)

Counts were normalised by size factors (total counts scaled to the median library size)
and log1p-transformed (shifted logarithm; Ahlmann-Eltze & Huber 2023). Analytic Pearson
residuals were computed for comparison. 3,000 highly variable genes (HVGs) were selected
with the seurat_v3 method on raw counts in a batch-aware way: only libraries with
≥500 cells contributed to the per-batch ranking. Mitochondrial, ribosomal, haemoglobin,
sex-chromosome and *MALAT1*/*NEAT1* genes were excluded from HVGs. Cell-cycle phase was
scored with the Tirosh *et al.* gene sets and was not regressed out. Principal component
analysis used 50 components.

## Integration (Notebook 04)

Library was treated as the batch. Two methods were compared: Harmony (θ = 2 on the
50 PCs) and scVI (2 layers of 128 hidden units, 30 latent dimensions, dropout 0.1,
negative-binomial likelihood, up to 300 epochs with early stopping, batch size 1,024,
trained on raw HVG counts). Embeddings were benchmarked with scib-metrics on a
40,000-cell stratified subsample using the authors' cell-type labels for the
biological-conservation metrics and library for the batch-correction metrics
(aggregate = 0.6 × biological conservation + 0.4 × batch correction). Harmony had the
higher aggregate score (0.636 vs 0.608 for scVI and 0.602 for uncorrected PCA) and was
used for all downstream steps of the discovery cohort.

## Clustering and annotation (Notebook 05)

Annotation was hierarchical. **Level 1:** a 15-nearest-neighbour graph was built from
the integrated embedding with exact Euclidean neighbours (scikit-learn
`KNeighborsTransformer`), so that the graph and the Leiden clustering (resolution 1.0)
are reproducible. Each cluster was assigned to the epithelial, immune, endothelial or
stromal compartment by the mean fraction of its cells detecting pan-compartment
markers ([config/markers.yaml](../config/markers.yaml)). Detection fractions were used
instead of control-gene-subtracted scores, which are distorted by total RNA content
(RNA-rich alveolar macrophages received negative immune scores). **Level 2:** each
compartment was re-embedded (2,000 HVGs, 30 PCs, Harmony by library) and re-clustered
(Leiden, resolution 1.0). **Level 3:** every subcluster with ≥ 500 cells was
re-embedded (20 PCs) and re-clustered once more (resolution 0.5); rare populations
such as aberrant basaloid cells and plasmacytoid dendritic cells (which share the
secretory programme of plasma cells) separated only at this level. Harmony was run on
the CPU with a single BLAS/OpenMP thread: on CUDA, or with multithreaded k-means, two
identical runs differed in the third decimal, which was enough to move the Leiden
boundaries of rare states.

Each final subcluster was labelled by explicit rules recorded in evidence tables:
(i) the best-scoring curated marker set (scanpy `score_genes`); (ii) when the margin to
the second-best set was < 0.3, the CellTypist majority label (Human Lung Cell Atlas
model; purity ≥ 0.6), but only to break the tie between the two best marker sets or for
types defined only by the reference (monocyte-derived macrophages), because a healthy
reference pulls disease states towards the nearest healthy label; otherwise the best
marker set was kept, unless its score was < 0.15 (no programme detected; all
unambiguous calls scored ≥ 0.19), in which case the subcluster was left unresolved
and removed; (iii) literature
definitions of transitional and disease states: transitional AT2 cells co-express the
AT2 and AT1 programmes at similar levels (both scores ≥ 0.5) (Strunz *et al.* 2020; Kobayashi *et al.* 2020),
aberrant basaloid cells are KRT17⁺KRT5⁻ epithelial cells (KRT17 detected in ≥ 50 % and
KRT5 in < 20 % of cells) with the aberrant basaloid programme among the two best
marker sets (Habermann *et al.* 2020; Adams *et al.* 2020), and proliferating
macrophages combine a macrophage label with a high proliferation score. Subclusters
were removed when their pan-compartment markers pointed to another compartment
(except mesothelium, which co-expresses epithelial keratins and stromal genes), when
> 15 % of cells co-expressed markers of two lineages (*PTPRC*/*EPCAM*/*CLDN5*/*COL1A2*
pairs; heterotypic doublets missed by scDblFinder), or when they contained < 50 cells.
Concordance with the authors' annotation is reported but was never used for a decision.

## Differential abundance (Notebook 06)

Cell-type counts per donor were tested for IPF versus control (12 vs 10 donors) with
propeller (speckle; logit transform, robust empirical Bayes, sex as covariate), both
across all cells and within each compartment, which is less sensitive to
compartment-level dissociation bias. scCODA (pertpy; NUTS sampling, 20,000 samples after
5,000 warm-up, automatic reference cell type, expected FDR 0.1) was run within
compartments. A chain with a mean acceptance rate outside 0.6–0.95 was restarted with
the next random seed (at most three tries); one of four compartments needed a second
seed. Because scCODA replaces zero counts by 0.5, numpyro's argument validation, which
rejects non-integer totals, was disabled for the sampler. A change was called robust if propeller FDR < 0.05 within its compartment
and scCODA reported a credible effect. Milo (pertpy; 30-neighbour graph, 5 % of cells
sampled as neighbourhood indices, design `~ sex + diagnosis`, PyDESeq2 solver, spatial FDR)
tested differential abundance without cluster boundaries. A Vanderbilt-only analysis
(the only site that contributed both groups) was used as a site sensitivity analysis.

## Differential expression (Notebook 07)

Raw counts were summed per donor and cell type (pseudobulk). Pseudobulk samples with
fewer than 10 cells were removed and cell types were tested only when at least three
donors per group remained. Genes were filtered with `edgeR::filterByExpr` and libraries normalised by TMM. The primary
test was edgeR quasi-likelihood (`glmQLFit` with `robust = TRUE`) with the design
`~ sex + diagnosis`; DESeq2 (Wald test, apeglm-shrunken log2 fold changes) was run on
the same genes as a second method. Significance was FDR < 0.05 (Benjamini–Hochberg);
"strong" DEGs additionally had |log2 FC| > 0.585. Site sensitivity was assessed by
adding site as a covariate (the design falls back to `~ diagnosis` where it is not of full
rank; the design used is recorded per cell type) and by repeating the
test in Vanderbilt donors only.

## Pathway and transcription-factor activity (Notebook 08)

Activities were estimated with decoupler (univariate linear model, ULM) from the signed
square root of the edgeR QL F statistic of each cell type: PROGENy pathway footprints (top 500 target genes
per pathway), CollecTRI TF regulons and MSigDB Hallmark gene sets. Per-cell PROGENy
scores were computed on log-normalised expression for visualisation. P values were
adjusted by Benjamini–Hochberg across cell types and sources.

## Alveolar epithelial trajectories (Notebook 09)

AT2, transitional AT2, AT1, aberrant basaloid and basal cells were subsampled to at
most 1,800 cells per type (all cells of rarer types were kept). PAGA summarised
connectivity between cell types, and diffusion pseudotime (DPT) was rooted in the AT2
cell with the highest AT2 marker score. CellRank 2 combined a pseudotime kernel on DPT
(soft threshold scheme, weight 0.8) with a connectivity kernel (0.2); the direction
therefore rests on the established role of AT2 cells as alveolar progenitors. A CytoTRACE
kernel was evaluated and rejected: it ranked aberrant basaloid cells as the most potent
cells of the lineage, and every population, AT1 included, was absorbed into the aberrant
basaloid state. Six macrostates were computed with GPCCA (Schur decomposition by the
Brandts method, because PETSc/SLEPc are unavailable on Windows) as an unsupervised check
of the end points. The whole AT1, aberrant basaloid and basal populations were set as
absorbing states, and fate probabilities of AT2 and transitional AT2 cells were computed
with the direct solver. Fate probabilities of transitional
AT2 cells were compared between IPF and control donors with a two-sided Mann–Whitney
test on donor means. Lineage drivers were genes correlated with the aberrant basaloid
fate. A 50-gene aberrant basaloid signature was derived by paired pseudobulk analysis
(aberrant basaloid versus alveolar epithelial cells of the same donors; edgeR QL with
`~ donor_id + group`; FDR < 0.01, log2 FC > 2, logCPM > 3, top 50 by log2 FC).

## Cell–cell communication (Notebook 10)

Ligand–receptor interactions were scored per condition with LIANA+ `rank_aggregate`
(consensus resource; ligand and receptor expressed in ≥10 % of the cells of their type;
1,000 permutations; at most 1,500 cells per type per condition). Differential
interactions were derived with `liana.multi.df_to_lr` from the pseudobulk edgeR results,
combining the change of ligand and receptor expression between IPF and control. Axes
targeting fibroblast states and the summed changes per sender and receiver were
reported.

## Replication in the validation cohort (Notebooks 11–12)

The validation cohort was processed with the same QC rules, normalisation, HVG
selection and scVI integration (batch = library). A CellTypist logistic-regression model
was trained on the discovery labels (up to 2,000 cells per type, shared genes) and used
to annotate the validation cells with majority voting over Leiden over-clusters at
CellTypist's default resolution for the number of cells (30 for > 200,000 cells); coarser
clusters let abundant neighbours outvote rare populations. Fine labels were mapped to a
common level shared with the authors' annotation, and agreement was measured on the
labels the model can predict by the adjusted Rand index, accuracy and balanced accuracy
(mean recall over cell types). The aberrant basaloid signature was
scored with `scanpy.tl.score_genes` among the authors' epithelial cells, and its ability
to identify the authors' aberrant basaloid cells was measured by the area under the ROC
curve, with *KRT17* alone as a single-gene baseline. Composition changes (propeller) and
pseudobulk DE (edgeR QL, `~ sex + diagnosis`) were repeated in the validation cohort.
Replication of DE was quantified by the sign concordance of discovery DEGs (one-sided
binomial test against 50 %). PROGENy activities were computed from the validation edgeR
statistics as in Notebook 08, and replication was assessed per pathway by the Spearman
correlation of activities across cell types and the sign agreement where the discovery
activity had FDR < 0.05. Discordant pathways were examined with Notebook 07's
site-adjusted and Vanderbilt-only contrasts and by their leading genes (PROGENy weight ×
discovery statistic, summed over cell types).

## Upstream re-processing and gene-regulatory networks (Linux; Notebooks 14–15)

FASTQ files of eight discovery libraries (4 control, 4 IPF) were re-processed on a
Linux workstation ([linux_upstream/](../linux_upstream/)) with STARsolo (GRCh38-3.0.0
reference, 737K-august-2016 whitelist, Cell Ranger-compatible UMI and cell-calling
options, spliced/unspliced counts). CellBender 0.3.2 `remove-background` (FPR 0.01,
150 epochs, expected cells from STARsolo's EmptyDrops call, max(3 × cells, 25,000)
droplets) removed ambient RNA from the raw droplet matrices on a CUDA GPU (RTX 4070).
STARsolo and the GEO matrices were compared per library; ambient contamination was
estimated per annotated cell type, and marker counts in the cell types of interest
were compared before and after CellBender. RNA velocity (scVelo 0.3.4, stochastic
model) was computed for the alveolar lineage as an exploratory analysis on the atlas
UMAP (genes with ≥ 10 shared counts, 2,000 highly variable genes, moments over 30
neighbours in 30 PCs). pySCENIC 0.12.1 (GRNBoost2, cisTarget ranking databases hg38
v10, AUCell) was run three times with different seeds on the 5,892 lineage cells and
regulons recovered in at least two runs were kept. GRNBoost2 ran through a small
driver that shares the expression matrix between worker processes; it calls the same
arboreto regressions with the same seeds as pySCENIC's
`arboreto_with_multiprocessing.py` and gave identical adjacencies in a per-gene check.
Regulon specificity scores (Suo *et al.* 2018) were computed per cell type and regulon
activity was correlated with the aberrant basaloid fate probability.

## Software and reproducibility

Analyses used Python 3.11 (scanpy, anndata, scvi-tools, harmonypy, scib-metrics,
CellTypist, pertpy, decoupler, LIANA+, CellRank, scVelo) and R (scDblFinder, edgeR,
DESeq2, apeglm, speckle). Exact versions are in [`envs/lock/`](../envs/lock/). Random
seeds are set in the configuration; GPU training of scVI is not bit-for-bit
deterministic, and annotation is rule-based so that labels do not depend on cluster
numbering. CellTypist training (stochastic gradient descent on multithreaded BLAS) is
not bit-for-bit reproducible either: two complete runs of Notebook 12 differed by 0.002
in balanced accuracy, by ≤ 0.01 in DE sign concordance and by ≤ 0.04 in the pathway
correlations, and no replication call changed. Figures were exported at 300 dpi
as SVG, JPG and PNG. Unit tests (`pytest`) and linting (`ruff`) run on every push
(GitHub Actions).
