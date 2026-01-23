"""Batch integration: Harmony and scVI, compared with scib-metrics.

Why integrate? Libraries differ by donor, centre, sequencing run and handling;
without correction cells cluster partly by library rather than by cell type.
Integration must remove *technical* variation while preserving *biological*
variation (including disease-associated states). We therefore benchmark
methods with scib-metrics (Luecken et al. 2022), which scores both
batch-mixing (iLISI, kBET, graph connectivity, PCR comparison) and
bio-conservation (NMI/ARI, cell-type ASW, cLISI, isolated labels) against
independent labels, and pick the best overall model.

Batch key = library (10x channel), the finest technical unit. Disease status
is never used as a batch covariate, otherwise the disease signal would be
regressed out.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

from .utils import log, timer


def run_harmony(adata: AnnData, batch_key: str, theta: float = 2.0, basis: str = "X_pca",
                seed: int = 0) -> np.ndarray:
    import harmonypy
    from threadpoolctl import threadpool_limits

    # Bit-for-bit reproducible: CPU instead of CUDA, and one BLAS/OpenMP thread for the
    # k-means steps (multithreaded reductions change the result in the 3rd decimal,
    # enough to move Leiden boundaries of rare states). Costs seconds, not minutes.
    with timer("Harmony"), threadpool_limits(limits=1):
        # contiguous copy: PCA can return negatively strided views, which torch rejects
        ho = harmonypy.run_harmony(np.ascontiguousarray(adata.obsm[basis]), adata.obs[[batch_key]], batch_key,
                                   theta=theta, random_state=seed, max_iter_harmony=20, device="cpu")
    z = np.asarray(ho.Z_corr)
    return z.T if z.shape[0] != adata.n_obs else z


def train_scvi(adata: AnnData, batch_key: str, params: dict, layer: str = "counts",
               save_dir=None, categorical_covariates: list[str] | None = None):
    """Train an scVI model on raw counts of the highly variable genes.

    scVI (Lopez et al. 2018) is a variational autoencoder with a negative
    binomial likelihood; the batch is given to the decoder, so the latent space
    ``z`` captures biological variation conditional on batch.
    """
    import scvi
    import torch

    torch.set_float32_matmul_precision("high")
    scvi.settings.seed = 0
    scvi.model.SCVI.setup_anndata(adata, layer=layer, batch_key=batch_key,
                                  categorical_covariate_keys=categorical_covariates)
    model = scvi.model.SCVI(
        adata,
        n_latent=params["n_latent"],
        n_layers=params["n_layers"],
        n_hidden=params["n_hidden"],
        dropout_rate=params["dropout_rate"],
        gene_likelihood=params["gene_likelihood"],
    )
    with timer("scVI training"):
        model.train(
            max_epochs=params["max_epochs"],
            batch_size=params["batch_size"],
            early_stopping=params["early_stopping"],
            early_stopping_patience=15,
            check_val_every_n_epoch=1,
            accelerator="gpu" if torch.cuda.is_available() else "cpu",
        )
    if save_dir is not None:
        model.save(str(save_dir), overwrite=True)
    return model


def benchmark(adata: AnnData, batch_key: str, label_key: str, embeddings: list[str],
              n_jobs: int = 8, subsample: int | None = None, seed: int = 0) -> pd.DataFrame:
    """Score embeddings with scib-metrics; returns the (unscaled) results table."""
    from scib_metrics.benchmark import BatchCorrection, Benchmarker, BioConservation

    ad = adata
    if subsample and adata.n_obs > subsample:
        idx = np.random.default_rng(seed).choice(adata.n_obs, subsample, replace=False)
        ad = adata[np.sort(idx)].copy()
    ad.obs[label_key] = ad.obs[label_key].astype(str).astype("category")
    bm = Benchmarker(
        ad,
        batch_key=batch_key,
        label_key=label_key,
        embedding_obsm_keys=embeddings,
        pre_integrated_embedding_obsm_key="X_pca",
        bio_conservation_metrics=BioConservation(),
        batch_correction_metrics=BatchCorrection(),
        n_jobs=n_jobs,
    )
    with timer(f"scib-metrics benchmark on {ad.n_obs:,} cells"):
        bm.benchmark()
    res = bm.get_results(min_max_scale=False)
    return res


def neighbors_umap(adata: AnnData, use_rep: str, n_neighbors: int = 15, min_dist: float = 0.3,
                   key: str | None = None, seed: int = 0) -> None:
    """kNN graph + UMAP on a given representation (stored under ``X_umap_<key>`` if key given)."""
    sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=n_neighbors, random_state=seed,
                    key_added=key)
    sc.tl.umap(adata, min_dist=min_dist, random_state=seed, neighbors_key=key)
    if key:
        adata.obsm[f"X_umap_{key}"] = adata.obsm["X_umap"].copy()
    log.info(f"UMAP from {use_rep} done")
