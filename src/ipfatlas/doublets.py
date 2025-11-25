"""Doublet detection: scDblFinder (R, primary) and Scrublet (Python, secondary).

Heterotypic doublets (two different cell types captured in one droplet) form
spurious "intermediate" clusters that can be mistaken for transitional states -
a serious risk in this project, whose focus is a transitional epithelial state.
Doublets are therefore removed *before* clustering.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

from .io import write_sparse_h5
from .utils import log, run_r


def run_scrublet(adata: AnnData, batch_key: str, expected_rate: float = 0.06, seed: int = 0) -> pd.DataFrame:
    """Scrublet (Wolock et al. 2019) per library on raw counts."""
    tmp = AnnData(X=adata.layers["counts"].copy(), obs=adata.obs[[batch_key]].copy())
    tmp.var_names = adata.var_names
    sc.pp.scrublet(tmp, batch_key=batch_key, expected_doublet_rate=expected_rate, random_state=seed)
    return tmp.obs[["doublet_score", "predicted_doublet"]].rename(
        columns={"doublet_score": "scrublet_score", "predicted_doublet": "scrublet_doublet"}
    )


def run_scdblfinder(adata: AnnData, rscript: str, r_script: Path, workdir: Path, batch_key: str,
                    threads: int = 8, reuse: bool = True) -> pd.DataFrame:
    """Export counts to HDF5, run workflow/R/scdblfinder.R, read the calls back.

    With ``reuse=True`` an existing calls file covering exactly these barcodes is
    read instead of re-running the (slow) R step.
    """
    workdir.mkdir(parents=True, exist_ok=True)
    out = workdir / "scdblfinder_calls.csv"
    if reuse and out.exists():
        res = pd.read_csv(out, index_col="barcode")
        if len(res) == adata.n_obs and res.index.isin(adata.obs_names).all():
            log.info(f"re-using scDblFinder calls from {out}")
            res["scDblFinder_doublet"] = res["scDblFinder_class"].eq("doublet")
            return res.reindex(adata.obs_names)[["scDblFinder_score", "scDblFinder_class", "scDblFinder_doublet"]]
    h5 = write_sparse_h5(adata, workdir / "counts_for_scdblfinder.h5", layer="counts", obs_cols=[batch_key])
    run_r(rscript, r_script, h5, out, str(threads), log_file=workdir / "scdblfinder.log")
    res = pd.read_csv(out, index_col="barcode")
    res = res.reindex(adata.obs_names)
    res["scDblFinder_doublet"] = res["scDblFinder_class"].eq("doublet")
    h5.unlink(missing_ok=True)  # large temporary file
    return res[["scDblFinder_score", "scDblFinder_class", "scDblFinder_doublet"]]


def concordance(obs: pd.DataFrame) -> pd.DataFrame:
    """2x2 agreement table between scDblFinder and Scrublet calls."""
    tab = pd.crosstab(obs["scDblFinder_doublet"], obs["scrublet_doublet"],
                      rownames=["scDblFinder"], colnames=["Scrublet"])
    jac = (obs["scDblFinder_doublet"] & obs["scrublet_doublet"]).sum() / max(
        (obs["scDblFinder_doublet"] | obs["scrublet_doublet"]).sum(), 1)
    log.info(f"doublet call Jaccard (scDblFinder vs Scrublet): {jac:.2f}")
    return tab


def expected_doublet_rate(n_cells_loaded: np.ndarray) -> np.ndarray:
    """10x guideline: ~0.8 % doublets per 1,000 recovered cells."""
    return 0.008 * np.asarray(n_cells_loaded) / 1000
