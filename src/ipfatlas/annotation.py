"""Cell-type annotation helpers.

Annotation follows a two-level, evidence-combining strategy:

1. **Reference-based**: CellTypist (Domínguez Conde *et al.* 2022) with the
   Human Lung Cell Atlas model (HLCA; Sikkema *et al.* 2023) predicts a label for
   every cell; labels are then smoothed by majority voting within
   over-clustered Leiden communities.
2. **Marker-based**: curated canonical markers (``config/markers.yaml``) are
   scored per cell and averaged per cluster. The HLCA is built from *healthy*
   lungs, so disease-emergent states (aberrant basaloid, KRT8hi transitional
   AT2, SPP1hi macrophages, CTHRC1+ myofibroblasts) are absent from the
   reference and must be identified from markers.

Final labels are assigned per cluster by an explicit, reviewable mapping in
the notebook - never by an opaque automatic rule.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc
import yaml
from anndata import AnnData

from .utils import log


def load_markers(path) -> dict:
    with open(path, encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def flat_markers(markers: dict, sections=("epithelial", "mesenchymal", "endothelial", "immune")) -> dict[str, list]:
    out = {}
    for s in sections:
        out.update(markers[s])
    return out


def run_celltypist(adata: AnnData, model: str, over_clustering: str | None = None, layer: str = "counts"):
    """CellTypist on log1p(CP10k) expression, as required by its models."""
    import celltypist

    tmp = AnnData(X=adata.layers[layer].copy(), obs=adata.obs[[]].copy(), var=adata.var[[]].copy())
    sc.pp.normalize_total(tmp, target_sum=1e4)
    sc.pp.log1p(tmp)
    if over_clustering:
        tmp.obs[over_clustering] = adata.obs[over_clustering].values
    pred = celltypist.annotate(tmp, model=model, majority_voting=over_clustering is not None,
                               over_clustering=over_clustering)
    res = pred.predicted_labels.copy()
    res["conf_score"] = pred.probability_matrix.max(axis=1)
    return res


def overclustering_resolution(n_cells: int) -> int:
    """CellTypist's default Leiden resolution for majority voting on ``n_cells`` cells.

    Majority voting needs clusters much smaller than cell types: in coarse clusters
    a rare population is outvoted by an abundant neighbour.
    """
    for limit, resolution in [(5_000, 5), (20_000, 10), (40_000, 15), (100_000, 20), (200_000, 25)]:
        if n_cells < limit:
            return resolution
    return 30


def score_marker_sets(adata: AnnData, sets: dict[str, list], prefix: str = "score_", seed: int = 0) -> list[str]:
    """Scanpy ``score_genes`` (Tirosh-style, control-gene corrected) for every marker set."""
    cols = []
    for name, genes in sets.items():
        g = [x for x in genes if x in adata.var_names]
        if len(g) < 2:
            log.warning(f"{name}: <2 markers present, skipped")
            continue
        col = prefix + name
        sc.tl.score_genes(adata, g, score_name=col, random_state=seed, use_raw=False)
        cols.append(col)
    return cols


def cluster_marker_table(adata: AnnData, groupby: str, score_cols: list[str]) -> pd.DataFrame:
    """Mean marker-set score per cluster (rows) and set (columns)."""
    t = adata.obs.groupby(groupby, observed=True)[score_cols].mean()
    t.columns = [c.replace("score_", "") for c in t.columns]
    return t


def lineage_coexpression(adata: AnnData, groupby: str, pairs: list[list[str]], thr: float = 0.0) -> pd.DataFrame:
    """Fraction of cells per cluster expressing BOTH genes of each lineage pair.

    Clusters with high co-expression of two lineage programmes (e.g. PTPRC and
    EPCAM) are residual heterotypic doublets that escaped scDblFinder.
    """
    out = {}
    for a, b in pairs:
        if a not in adata.var_names or b not in adata.var_names:
            continue
        xa = np.asarray(adata[:, a].X.todense()).ravel() > thr if hasattr(adata.X, "todense") else adata[:, a].X.ravel() > thr
        xb = np.asarray(adata[:, b].X.todense()).ravel() > thr if hasattr(adata.X, "todense") else adata[:, b].X.ravel() > thr
        out[f"{a}+{b}+"] = pd.Series(xa & xb, index=adata.obs_names).groupby(adata.obs[groupby], observed=True).mean()
    return pd.DataFrame(out)


def majority_label(adata: AnnData, groupby: str, label_key: str) -> pd.DataFrame:
    """Most frequent label per cluster and its purity."""
    ct = pd.crosstab(adata.obs[groupby], adata.obs[label_key])
    top = ct.idxmax(axis=1)
    purity = ct.max(axis=1) / ct.sum(axis=1)
    second = ct.apply(lambda r: r.sort_values(ascending=False).index[1] if (r > 0).sum() > 1 else "", axis=1)
    return pd.DataFrame({"top_label": top, "purity": purity.round(3), "second_label": second,
                         "n_cells": ct.sum(axis=1)})
