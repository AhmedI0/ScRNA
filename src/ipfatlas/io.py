"""Read the public count matrices into AnnData with harmonised metadata.

Raw UMI counts are returned in ``.X`` (float32, CSR). Notebooks copy them to
``layers["counts"]`` once barcodes have been filtered, to avoid holding two
copies of the unfiltered matrix.

Harmonised ``.obs`` columns (both cohorts)
------------------------------------------
library_id   sequencing library (10x channel) - technical batch
donor_id     patient / donor - the biological replicate
diagnosis    Control, IPF, NSIP, cHP, Sarcoidosis, ILD, COPD
status       Control / Disease
site         recruiting centre
cohort       GSE accession
author_celltype, author_compartment   labels published by the original authors
"""
from __future__ import annotations

import gzip
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

from .utils import log, timer


def _read_mtx(path: Path) -> sp.csr_matrix:
    """Read a (gzipped) MatrixMarket file and return a CSR matrix (genes x cells)."""
    with timer(f"reading {path.name}"):
        with gzip.open(path, "rb") if path.suffix == ".gz" else open(path, "rb") as fh:
            m = scipy.io.mmread(fh)
    return sp.csr_matrix(m)


def _cells_by_genes(m_genes_by_cells: sp.spmatrix) -> sp.csr_matrix:
    x = sp.csr_matrix(m_genes_by_cells.T)
    x.data = x.data.astype(np.float32)
    x.sort_indices()
    return x


# --------------------------------------------------------------------------- #
# Discovery cohort: Habermann et al. 2020 (GSE135893)
# --------------------------------------------------------------------------- #
def read_habermann(raw_dir: Path) -> ad.AnnData:
    """GSE135893: Cell Ranger 3.0.2 filtered barcodes of 43 libraries (pre author-QC).

    The GEO matrix holds 220,213 barcodes; the authors' metadata annotates the
    114,396 cells that passed their QC. We keep *all* barcodes and perform our
    own QC; the author labels are stored for later comparison only.
    """
    d = Path(raw_dir) / "GSE135893"
    genes = pd.read_csv(d / "GSE135893_genes.tsv.gz", header=None, sep="\t")[0].astype(str)
    barcodes = pd.read_csv(d / "GSE135893_barcodes.tsv.gz", header=None, sep="\t")[0].astype(str)
    x = _cells_by_genes(_read_mtx(d / "GSE135893_matrix.mtx.gz"))
    adata = ad.AnnData(X=x, obs=pd.DataFrame(index=barcodes.values), var=pd.DataFrame(index=genes.values))
    adata.var_names_make_unique()
    adata.var["gene_symbol"] = genes.values

    # library prefix: everything before the last "_" (e.g. "ILD59-1_AAACCTG...")
    adata.obs["library_id"] = adata.obs_names.str.rsplit("_", n=1).str[0]

    meta = pd.read_csv(d / "GSE135893_IPF_metadata.csv.gz", index_col=0)
    lib = (
        meta.groupby("orig.ident")[["Sample_Name", "Diagnosis", "Status", "Sample_Source"]]
        .first()
        .rename(columns={"Sample_Name": "donor_id", "Diagnosis": "diagnosis", "Status": "status", "Sample_Source": "site"})
    )
    # the same donor is written VUHD071 in one place and VUHD71 elsewhere
    lib["donor_id"] = lib["donor_id"].replace({"VUHD071": "VUHD71"})
    unknown = sorted(set(adata.obs["library_id"]) - set(lib.index))
    if unknown:
        log.warning(f"{len(unknown)} libraries without metadata are dropped: {unknown}")
        adata = adata[adata.obs["library_id"].isin(lib.index)].copy()
    adata.obs = adata.obs.join(lib, on="library_id")
    adata.obs["author_celltype"] = meta["celltype"].reindex(adata.obs_names).values
    adata.obs["author_compartment"] = meta["population"].reindex(adata.obs_names).values
    adata.obs["author_qc_pass"] = adata.obs_names.isin(meta.index)
    adata.obs["cohort"] = "GSE135893"
    adata.obs["chemistry"] = "10x 5' v1"
    for c in ["library_id", "donor_id", "diagnosis", "status", "site", "author_celltype", "author_compartment"]:
        adata.obs[c] = adata.obs[c].astype("category")
    log.info(f"Habermann: {adata.n_obs:,} barcodes x {adata.n_vars:,} genes, "
             f"{adata.obs.library_id.nunique()} libraries, {adata.obs.donor_id.nunique()} donors")
    return adata


# --------------------------------------------------------------------------- #
# Validation cohort: Adams et al. 2020 (GSE136831)
# --------------------------------------------------------------------------- #
ADAMS_TO_PROJECT = {
    # Adams "Manuscript_Identity" -> harmonised label used for cross-cohort comparison
    "ATII": "AT2", "ATI": "AT1", "Aberrant_Basaloid": "Aberrant basaloid", "Basal": "Basal",
    "Club": "Secretory (club/goblet)", "Goblet": "Secretory (club/goblet)", "Ciliated": "Ciliated",
    "PNEC": "PNEC", "Ionocyte": "Ionocyte", "Mesothelial": "Mesothelial",
    "Fibroblast": "Fibroblasts", "Myofibroblast": "Myofibroblasts", "SMC": "Smooth muscle", "Pericyte": "Pericytes",
    "VE_Capillary_A": "Capillary EC", "VE_Capillary_B": "Capillary EC", "VE_Arterial": "Arterial EC",
    "VE_Venous": "Venous EC", "VE_Peribronchial": "Venous EC", "Lymphatic": "Lymphatic EC",
    "Macrophage_Alveolar": "Alveolar macrophages", "Macrophage": "Monocyte-derived macrophages",
    "cMonocyte": "Classical monocytes", "ncMonocyte": "Non-classical monocytes",
    "cDC1": "cDC", "cDC2": "cDC", "DC_Mature": "cDC", "DC_Langerhans": "cDC", "pDC": "pDC",
    "Mast": "Mast cells", "B": "B cells", "B_Plasma": "Plasma cells",
    "T": "CD4 T cells", "T_Regulatory": "CD4 T cells", "T_Cytotoxic": "CD8 T cells",
    "NK": "NK cells", "ILC_A": "ILC", "ILC_B": "ILC",
    "Multiplet": "Multiplet",
}


def read_adams(raw_dir: Path) -> ad.AnnData:
    """GSE136831: 312,928 cells (author-filtered), 78 donors (IPF, COPD, control)."""
    d = Path(raw_dir) / "GSE136831"
    genes = pd.read_csv(d / "GSE136831_AllCells.GeneIDs.txt.gz", sep="\t")
    barcodes = pd.read_csv(d / "GSE136831_AllCells.cellBarcodes.txt.gz", header=None)[0].astype(str)
    meta = pd.read_csv(d / "GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz", sep="\t")
    meta = meta.set_index("CellBarcode_Identity")
    if not (meta.index == barcodes.values).all():
        meta = meta.reindex(barcodes.values)
    x = _cells_by_genes(_read_mtx(d / "GSE136831_RawCounts_Sparse.mtx.gz"))
    var = pd.DataFrame(
        {"gene_ids": genes["Ensembl_GeneID"].values, "gene_symbol": genes["HGNC_EnsemblAlt_GeneID"].values},
        index=genes["HGNC_EnsemblAlt_GeneID"].astype(str).values,
    )
    adata = ad.AnnData(X=x, obs=pd.DataFrame(index=barcodes.values), var=var)
    adata.var_names_make_unique()
    obs = pd.DataFrame(index=adata.obs_names)
    obs["library_id"] = meta["Library_Identity"].astype(str).values
    obs["donor_id"] = meta["Subject_Identity"].astype(str).values
    obs["diagnosis"] = meta["Disease_Identity"].astype(str).values
    obs["status"] = np.where(obs["diagnosis"] == "Control", "Control", "Disease")
    obs["site"] = "Yale/BWH"
    obs["author_celltype"] = meta["Manuscript_Identity"].astype(str).values
    obs["author_subclass"] = meta["Subclass_Cell_Identity"].astype(str).values
    obs["author_compartment"] = meta["CellType_Category"].astype(str).values
    obs["harmonised_celltype"] = obs["author_celltype"].map(ADAMS_TO_PROJECT).fillna("Other")
    obs["cohort"] = "GSE136831"
    obs["chemistry"] = "10x 3' v2"
    for c in obs.columns:
        obs[c] = obs[c].astype("category")
    adata.obs = obs
    log.info(f"Adams: {adata.n_obs:,} cells x {adata.n_vars:,} genes, {obs.donor_id.nunique()} donors")
    return adata


def write_sparse_h5(adata: ad.AnnData, path: Path, layer: str = "counts", obs_cols: list[str] | None = None) -> Path:
    """Write counts as a CSC (genes x cells) HDF5 file that R can read without anndata.

    A cells x genes CSR matrix has exactly the same ``data/indices/indptr`` as a
    genes x cells CSC matrix, so R reconstructs a ``dgCMatrix`` directly.
    """
    import h5py

    x = adata.layers[layer] if layer else adata.X
    x = sp.csr_matrix(x)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as f:
        f.create_dataset("data", data=x.data.astype(np.float64), compression="gzip", compression_opts=1)
        f.create_dataset("indices", data=x.indices.astype(np.int32), compression="gzip", compression_opts=1)
        f.create_dataset("indptr", data=x.indptr.astype(np.int64))
        f.create_dataset("shape", data=np.array([adata.n_vars, adata.n_obs], dtype=np.int64))
        f.create_dataset("genes", data=np.array(adata.var_names, dtype="S"))
        f.create_dataset("barcodes", data=np.array(adata.obs_names, dtype="S"))
        for c in obs_cols or []:
            f.create_dataset(f"obs/{c}", data=np.array(adata.obs[c].astype(str), dtype="S"))
    return path
