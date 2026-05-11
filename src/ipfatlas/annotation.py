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


# Fine discovery labels -> coarser labels shared with the Adams et al. annotation
# (io.ADAMS_TO_PROJECT), used for every cross-cohort comparison.
FINE_TO_HARMONISED = {
    "AT1": "AT1", "AT2": "AT2", "Transitional AT2": "AT2", "Aberrant basaloid": "Aberrant basaloid",
    "Basal": "Basal", "Secretory (club/goblet)": "Secretory (club/goblet)",
    "SCGB3A2+ secretory": "Secretory (club/goblet)", "Ciliated": "Ciliated", "PNEC": "PNEC", "Ionocyte": "Ionocyte",
    "Mesothelial": "Mesothelial",
    "Alveolar fibroblasts": "Fibroblasts", "Adventitial fibroblasts": "Fibroblasts", "HAS1hi fibroblasts": "Fibroblasts",
    "PLIN2+ fibroblasts": "Fibroblasts", "Myofibroblasts": "Myofibroblasts", "Smooth muscle": "Smooth muscle",
    "Pericytes": "Pericytes",
    "Aerocyte capillary EC": "Capillary EC", "General capillary EC": "Capillary EC", "Arterial EC": "Arterial EC",
    "Venous EC": "Venous EC", "Systemic venous EC (COL15A1+)": "Venous EC", "Lymphatic EC": "Lymphatic EC",
    "Alveolar macrophages": "Alveolar macrophages", "SPP1hi macrophages": "Monocyte-derived macrophages",
    "Interstitial macrophages": "Monocyte-derived macrophages",
    "Proliferating macrophages": "Monocyte-derived macrophages",
    "Classical monocytes": "Classical monocytes", "Non-classical monocytes": "Non-classical monocytes",
    "cDC": "cDC", "pDC": "pDC", "Mast cells": "Mast cells", "CD4 T cells": "CD4 T cells", "Treg": "CD4 T cells",
    "CD8 T cells": "CD8 T cells", "NK cells": "NK cells", "B cells": "B cells", "Plasma cells": "Plasma cells",
    "Monocyte-derived macrophages": "Monocyte-derived macrophages",
}


def harmonise(labels: pd.Series) -> pd.Series:
    """Map fine labels to the cross-cohort level; unknown labels become 'Other'."""
    return labels.astype(str).map(FINE_TO_HARMONISED).fillna("Other").astype("category")


# HLCA (CellTypist Human_Lung_Atlas) labels -> labels of this atlas. Used only when the
# marker call of a cluster is ambiguous. Labels without an unambiguous counterpart
# (e.g. AT0, peribronchial fibroblasts) are deliberately left out.
CELLTYPIST_TO_LABEL = {
    "AT1": "AT1", "AT2": "AT2", "AT2 proliferating": "AT2",
    "pre-TB secretory": "SCGB3A2+ secretory", "Club (non-nasal)": "Secretory (club/goblet)",
    "Club (nasal)": "Secretory (club/goblet)", "Goblet (bronchial)": "Secretory (club/goblet)",
    "Goblet (nasal)": "Secretory (club/goblet)", "Goblet (subsegmental)": "Secretory (club/goblet)",
    "Multiciliated (non-nasal)": "Ciliated", "Multiciliated (nasal)": "Ciliated", "Deuterosomal": "Ciliated",
    "Basal resting": "Basal", "Suprabasal": "Basal", "Hillock-like": "Basal",
    "Neuroendocrine": "PNEC", "Ionocyte": "Ionocyte", "Mesothelium": "Mesothelial",
    "Alveolar fibroblasts": "Alveolar fibroblasts", "Adventitial fibroblasts": "Adventitial fibroblasts",
    "Myofibroblasts": "Myofibroblasts", "Smooth muscle": "Smooth muscle",
    "Smooth muscle FAM83D+": "Smooth muscle", "SM activated stress response": "Smooth muscle",
    "Pericytes": "Pericytes",
    "EC aerocyte capillary": "Aerocyte capillary EC", "EC general capillary": "General capillary EC",
    "EC arterial": "Arterial EC", "EC venous pulmonary": "Venous EC",
    "EC venous systemic": "Systemic venous EC (COL15A1+)", "Lymphatic EC mature": "Lymphatic EC",
    "Lymphatic EC differentiating": "Lymphatic EC", "Lymphatic EC proliferating": "Lymphatic EC",
    "Alveolar macrophages": "Alveolar macrophages", "Alveolar Mph CCL3+": "Alveolar macrophages",
    "Alveolar Mph MT-positive": "Alveolar macrophages", "Alveolar Mph proliferating": "Proliferating macrophages",
    "Monocyte-derived Mph": "Monocyte-derived macrophages", "Interstitial Mph perivascular": "Interstitial macrophages",
    "Classical monocytes": "Classical monocytes", "Non-classical monocytes": "Non-classical monocytes",
    "DC1": "cDC", "DC2": "cDC", "Migratory DCs": "cDC", "Plasmacytoid DCs": "pDC", "Mast cells": "Mast cells",
    "CD4 T cells": "CD4 T cells", "CD8 T cells": "CD8 T cells", "NK cells": "NK cells", "B cells": "B cells",
    "Plasma cells": "Plasma cells",
}


def subcluster(adata: AnnData, mask: np.ndarray, use_rep: str, resolution: float, prefix: str,
               n_neighbors: int = 15, seed: int = 0) -> pd.Series:
    """Leiden clustering of a subset on its own kNN graph (exact neighbours, reproducible).

    Only the embedding is copied, not the expression matrix. Returns labels
    ``f"{prefix}{cluster}"`` for the masked cells.
    """
    sub = AnnData(obs=adata.obs.loc[mask, []].copy())
    sub.obsm[use_rep] = np.asarray(adata.obsm[use_rep][mask])
    return (prefix + leiden_on(sub, use_rep, resolution, n_neighbors, seed)).rename("subcluster")


def leiden_on(adata: AnnData, use_rep: str, resolution: float, n_neighbors: int = 15, seed: int = 0) -> pd.Series:
    """Leiden on an exact kNN graph of ``use_rep`` (reproducible for a given embedding)."""
    from sklearn.neighbors import KNeighborsTransformer

    sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=n_neighbors, random_state=seed, key_added="_sub",
                    transformer=KNeighborsTransformer(n_neighbors=n_neighbors, n_jobs=-1))
    sc.tl.leiden(adata, resolution=resolution, flavor="igraph", n_iterations=2, random_state=seed,
                 neighbors_key="_sub", key_added="_leiden")
    return adata.obs.pop("_leiden").astype(str)


def reembed(adata: AnnData, cfg_features: dict, batch_key: str, n_top_genes: int = 2000, n_pcs: int = 30,
            min_cells_per_batch: int = 100, theta: float = 2.0, key: str = "X_sub", seed: int = 0) -> None:
    """New HVGs, PCA and Harmony for a subset of cells (in place, stored in ``obsm[key]``).

    In a global PCA the axes are dominated by differences between lineages, so
    variation that separates states within one lineage (AT1 vs transitional AT2
    vs aberrant basaloid) is diluted. Re-embedding each compartment recovers it,
    as in the original study and the HLCA.
    """
    from .features import select_hvgs
    from .integration import run_harmony

    select_hvgs(adata, {**cfg_features, "n_top_genes": n_top_genes}, batch_key,
                min_cells_per_batch=min_cells_per_batch)
    sc.pp.pca(adata, n_comps=n_pcs, mask_var="highly_variable", random_state=seed)
    adata.obsm[key] = run_harmony(adata, batch_key, theta=theta, seed=seed)


def detection_scores(adata: AnnData, groupby: str, sets: dict[str, list]) -> pd.DataFrame:
    """Per cluster: fraction of cells detecting each marker, averaged over the set.

    Unlike control-gene-subtracted scores (``scanpy.tl.score_genes``), this is not
    distorted by total RNA content: RNA-rich alveolar macrophages get a near-zero
    or negative immune *score* because their control genes are also high, although
    almost all of them detect TYROBP, LAPTM5 and FCER1G.
    """
    genes = sorted({g for gs in sets.values() for g in gs if g in adata.var_names})
    X = adata[:, genes].X > 0
    det = pd.DataFrame(X.toarray() if hasattr(X, "toarray") else np.asarray(X), index=adata.obs_names, columns=genes)
    frac = det.groupby(adata.obs[groupby].to_numpy()).mean()
    return pd.DataFrame({k: frac[[g for g in gs if g in genes]].mean(axis=1) for k, gs in sets.items()})


def assign_by_scores(cluster_scores: pd.DataFrame, min_margin: float = 0.05) -> pd.DataFrame:
    """Best-scoring marker set per cluster, with the runner-up and the margin between them."""
    s = cluster_scores
    order = np.argsort(-s.to_numpy(), axis=1)
    top = s.columns[order[:, 0]]
    second = s.columns[order[:, 1]]
    v = np.sort(s.to_numpy(), axis=1)[:, ::-1]
    out = pd.DataFrame({"marker_label": top, "runner_up": second, "top_score": v[:, 0].round(3),
                        "margin": (v[:, 0] - v[:, 1]).round(3)}, index=s.index)
    out["confident"] = out["margin"] >= min_margin
    return out


def top_markers(adata: AnnData, groupby: str, n: int = 8, max_cells: int = 300, seed: int = 0) -> pd.Series:
    """Top up-regulated genes per group (t-test on a per-group subsample; descriptive only)."""
    rng = np.random.default_rng(seed)
    idx = np.concatenate([rng.choice(np.where(adata.obs[groupby] == g)[0], min(max_cells, (adata.obs[groupby] == g).sum()),
                                     replace=False) for g in adata.obs[groupby].cat.categories])
    sub = adata[np.sort(idx)].copy()
    sc.tl.rank_genes_groups(sub, groupby, method="t-test_overestim_var", n_genes=n)
    names = pd.DataFrame(sub.uns["rank_genes_groups"]["names"])
    return names.apply(lambda c: ", ".join(c.astype(str))).rename("top_genes")
