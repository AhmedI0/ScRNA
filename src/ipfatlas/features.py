"""Highly variable gene selection."""
from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

from .utils import log


def select_hvgs(adata: AnnData, cfg_features: dict, batch_key: str, layer: str = "counts",
                min_cells_per_batch: int = 500) -> pd.Index:
    """Batch-aware Seurat-v3 HVGs on raw counts, excluding nuisance genes.

    * mitochondrial, ribosomal, haemoglobin and sex-chromosome genes are excluded
      so that cells do not cluster by stress, blood contamination or donor sex;
    * the per-batch loess fit of the Seurat v3 method is numerically unstable for
      very small libraries, so only libraries with >= ``min_cells_per_batch``
      cells vote on the HVG set (the selected genes are then used for all cells);
    * if the loess fit still fails, the span is widened (0.3 -> 0.5 -> 0.8).

    Sets ``adata.var['highly_variable']`` (+ rank and n_batches) and returns the HVGs.
    """
    f = cfg_features
    excl = adata.var_names.str.contains("|".join(f["exclude_patterns"])) | adata.var_names.isin(f["exclude_genes"])
    n_per = adata.obs[batch_key].value_counts()
    keep_batches = n_per[n_per >= min_cells_per_batch].index
    cand = adata[adata.obs[batch_key].isin(keep_batches), ~excl].copy()
    cand.obs[batch_key] = cand.obs[batch_key].astype(str).astype("category")
    # genes must be detected in the voting cells, otherwise loess sees constant columns
    cand = cand[:, np.asarray((cand.layers[layer] > 0).sum(axis=0)).ravel() >= 3].copy()
    log.info(f"HVG selection: {excl.sum()} nuisance genes excluded; {len(keep_batches)}/{len(n_per)} "
             f"batches with >= {min_cells_per_batch} cells vote")
    for span in (0.3, 0.5, 0.8):
        try:
            sc.pp.highly_variable_genes(cand, flavor=f["hvg_flavor"], n_top_genes=f["n_top_genes"], layer=layer,
                                        batch_key=batch_key, span=span)
            break
        except ValueError as e:
            log.warning(f"loess failed with span={span} ({e}); retrying with a wider span")
    else:
        raise RuntimeError("HVG selection failed for all spans")
    hvg = cand.var_names[cand.var["highly_variable"]]
    adata.var["highly_variable"] = adata.var_names.isin(hvg)
    adata.var["highly_variable_rank"] = cand.var["highly_variable_rank"].reindex(adata.var_names)
    adata.var["highly_variable_nbatches"] = cand.var["highly_variable_nbatches"].reindex(adata.var_names).fillna(0)
    return pd.Index(hvg)
