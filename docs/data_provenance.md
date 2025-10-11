# Data provenance

All data are public. Nothing is redistributed in this repository; every file is
downloaded by `workflow/scripts/download_data.py`, which checks the remote size
and records a SHA-256 checksum in `data/raw/download_manifest.tsv`.

## Discovery cohort — GSE135893

| item | value |
|---|---|
| Study | Habermann AC, Gutierrez AJ, Bui LT, *et al.* Single-cell RNA sequencing reveals profibrotic roles of distinct epithelial and mesenchymal lineages in pulmonary fibrosis. *Sci Adv* 2020;6:eaba1972. doi:10.1126/sciadv.aba1972 |
| GEO / SRA | GSE135893 · BioProject PRJNA560492 · SRP218543 |
| Tissue | Lung biopsies / explants: non-fibrotic donor lungs (controls) and pulmonary fibrosis lungs |
| Sites | Vanderbilt University (VU), Norton Thoracic Institute (NTI), Donor Network of Arizona (DNA, controls only) |
| Chemistry | 10x Genomics Chromium Single Cell 5′ (v1); Illumina HiSeq 4000 / NovaSeq 6000 / NextSeq 500 |
| Processing by authors | Cell Ranger 3.0.2, reference GRCh38-3.0.0; Seurat v3; cells with < 1,000 genes or > 25 % mitochondrial reads removed |
| Files used | `GSE135893_matrix.mtx.gz` (33,694 genes × 220,213 barcodes, **before** author QC), `_genes.tsv.gz`, `_barcodes.tsv.gz`, `_IPF_metadata.csv.gz` (114,396 annotated cells), `_family.soft.gz` |
| Design | 43 libraries from 30 donors: 10 control, 12 IPF, 3 NSIP, 2 cHP, 2 sarcoidosis, 1 unclassifiable ILD |

Notes
* The donor `VUHD071` (metadata) and `VUHD71` (GEO sample title) are the same person; harmonised to `VUHD71`.
* Several donors contributed two libraries (e.g. `ILD59-1`, `ILD59-2`); the **donor** is the statistical replicate.
* Site is partly confounded with disease (DNA: controls only; NTI: disease only). See `docs/methods.md`, *Sensitivity analyses*.

## Validation cohort — GSE136831

| item | value |
|---|---|
| Study | Adams TS, Schupp JC, Poli S, *et al.* Single-cell RNA-seq reveals ectopic and aberrant lung-resident cell populations in idiopathic pulmonary fibrosis. *Sci Adv* 2020;6:eaba1983. doi:10.1126/sciadv.aba1983 |
| Sites | Yale School of Medicine / Brigham and Women's Hospital |
| Chemistry | 10x Genomics Chromium Single Cell 3′ (v2) |
| Files used | `GSE136831_RawCounts_Sparse.mtx.gz` (45,947 genes × 312,928 cells, author-filtered ≥ 1,000 UMIs), `_AllCells.GeneIDs.txt.gz`, `_AllCells.cellBarcodes.txt.gz`, `_AllCells.Samples.CellType.MetadataTable.txt.gz` |
| Design | 78 donors: 28 control, 32 IPF, 18 COPD (COPD excluded here) |

## Raw reads (Linux upstream package)

FASTQ files for 8 single-library donors of GSE135893 (4 control, 4 IPF; ~96 GB) are
fetched from the European Nucleotide Archive with MD5 verification; accessions,
URLs and checksums are in `linux_upstream/config/samples.tsv`.

| resource | source |
|---|---|
| Reference genome + annotation | 10x Genomics `refdata-cellranger-GRCh38-3.0.0` (identical to the original study) |
| Barcode whitelist | `737K-august-2016.txt` (10x 3′ v2 / 5′ v1), MD5 `f62a276e262fdd85262a889d0f48556b` |
| cisTarget databases | Aerts lab, hg38 refseq_r80, motif collection v10 (clustered), 500 bp/100 bp and 10 kb/10 kb windows |

## Prior-knowledge resources

| resource | used for | access |
|---|---|---|
| Human Lung Cell Atlas CellTypist model (`Human_Lung_Atlas.pkl`, v2) | reference-based annotation | `celltypist.models.download_models` |
| PROGENy (top 500 genes / pathway) | pathway activity | OmniPath via `decoupler.op.progeny` |
| CollecTRI | TF activity | OmniPath via `decoupler.op.collectri` |
| MSigDB Hallmarks | programme enrichment | OmniPath via `decoupler.op.hallmark` |
| LIANA consensus ligand–receptor resource | cell–cell communication | bundled with `liana` |
