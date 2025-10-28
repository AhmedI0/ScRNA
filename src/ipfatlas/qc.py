"""Cell- and gene-level quality control.

Strategy (Luecken & Theis 2019; Heumos et al. 2023)
---------------------------------------------------
QC metrics differ systematically between libraries (sequencing depth, tissue
quality, chemistry), so fixed global thresholds remove too many cells from some
libraries and too few from others. We therefore flag outliers *within each
library* using the median absolute deviation (MAD):

    outlier  <=>  |x - median(x)| > n * MAD(x)

on log-transformed library size and detected genes (both tails) and on the
fraction of counts in the top-20 genes (upper tail: low-complexity cells).
Because lung libraries mix very large (macrophages) and very small
(lymphocytes) cells, MAD bounds are wide; hard floors on genes and UMIs
therefore remove ambient-rich, low-capture barcodes. Mitochondrial
fraction uses a one-sided rule ``median + 3 MAD`` bounded to
[``mito_floor``, ``max_pct_mito``]: the floor avoids discarding metabolically
active cells in libraries with very low baseline mito content, the ceiling
matches the 25 % cut-off of the original study.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from scipy.stats import median_abs_deviation


def annotate_gene_groups(adata: AnnData) -> None:
    """Flag mitochondrial, ribosomal-protein and haemoglobin genes in ``adata.var``."""
    names = adata.var_names.str.upper()
    adata.var["mt"] = names.str.startswith("MT-")
    adata.var["ribo"] = names.str.match(r"^RP[SL]\d")
    adata.var["hb"] = names.str.match(r"^HB[ABDEGMQZ]\d?$")


def compute_qc_metrics(adata: AnnData, layer: str | None = "counts") -> None:
    annotate_gene_groups(adata)
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], percent_top=[20], log1p=True, inplace=True, layer=layer
    )


def is_outlier(x: pd.Series, nmads: float, side: str = "both") -> pd.Series:
    med = np.median(x)
    mad = median_abs_deviation(x)
    lo, hi = med - nmads * mad, med + nmads * mad
    if side == "upper":
        return x > hi
    if side == "lower":
        return x < lo
    return (x < lo) | (x > hi)


def flag_outliers(adata: AnnData, cfg_qc: dict, groupby: str = "library_id") -> pd.DataFrame:
    """Add boolean QC flags to ``adata.obs`` and return per-library thresholds."""
    obs = adata.obs
    flags = pd.DataFrame(index=obs.index)
    thresholds = []
    for lib, idx in obs.groupby(groupby, observed=True).groups.items():
        o = obs.loc[idx]
        f_counts = is_outlier(o["log1p_total_counts"], cfg_qc["nmads_counts"])
        f_genes = is_outlier(o["log1p_n_genes_by_counts"], cfg_qc["nmads_genes"])
        f_top = is_outlier(o["pct_counts_in_top_20_genes"], cfg_qc["nmads_top20"], side="upper")
        mito = o["pct_counts_mt"]
        mito_thr = float(np.clip(np.median(mito) + cfg_qc["nmads_mito"] * median_abs_deviation(mito),
                                 cfg_qc.get("mito_floor", 10), cfg_qc["max_pct_mito"]))
        f_mito = mito > mito_thr
        f_hard = (o["n_genes_by_counts"] < cfg_qc["min_genes"]) | (o["total_counts"] < cfg_qc["min_counts"])
        flags.loc[idx, "outlier_counts"] = f_counts.values
        flags.loc[idx, "outlier_genes"] = f_genes.values
        flags.loc[idx, "outlier_top20"] = f_top.values
        flags.loc[idx, "outlier_mito"] = f_mito.values
        flags.loc[idx, "outlier_hard"] = f_hard.values
        med_c, mad_c = np.median(o["log1p_total_counts"]), median_abs_deviation(o["log1p_total_counts"])
        med_g, mad_g = np.median(o["log1p_n_genes_by_counts"]), median_abs_deviation(o["log1p_n_genes_by_counts"])
        thresholds.append({
            groupby: lib,
            "n_barcodes": len(o),
            "counts_low": np.expm1(med_c - cfg_qc["nmads_counts"] * mad_c),
            "counts_high": np.expm1(med_c + cfg_qc["nmads_counts"] * mad_c),
            "genes_low": np.expm1(med_g - cfg_qc["nmads_genes"] * mad_g),
            "genes_high": np.expm1(med_g + cfg_qc["nmads_genes"] * mad_g),
            "mito_max_pct": mito_thr,
        })
    flags = flags.astype(bool)
    for c in flags.columns:
        obs[c] = flags[c]
    obs["qc_fail"] = flags.any(axis=1)
    thr = pd.DataFrame(thresholds).set_index(groupby)
    thr["n_fail"] = obs.groupby(groupby, observed=True)["qc_fail"].sum()
    thr["pct_fail"] = 100 * thr["n_fail"] / thr["n_barcodes"]
    return thr


def qc_summary(adata: AnnData, groupby: str = "library_id") -> pd.DataFrame:
    """Median QC metrics per library (reported as a supplementary table)."""
    cols = ["total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo", "pct_counts_hb",
            "pct_counts_in_top_20_genes"]
    out = adata.obs.groupby(groupby, observed=True)[cols].median().round(2)
    out.insert(0, "n_cells", adata.obs.groupby(groupby, observed=True).size())
    return out


def infer_sex(adata: AnnData, donor_key: str = "donor_id", layer: str = "counts") -> pd.DataFrame:
    """Infer donor genetic sex from XIST (female) vs chrY genes (male) expression.

    Sex is rarely provided in GEO metadata but is an important covariate for
    differential expression. Per donor we compute the fraction of counts from
    XIST and from a panel of broadly expressed chrY genes.
    """
    y_genes = [g for g in ["RPS4Y1", "DDX3Y", "UTY", "KDM5D", "EIF1AY", "USP9Y", "ZFY"] if g in adata.var_names]
    x = adata.layers[layer] if layer else adata.X
    total = np.asarray(x.sum(axis=1)).ravel()
    xist = np.asarray(x[:, adata.var_names.get_loc("XIST")].sum(axis=1)).ravel()
    y_idx = [adata.var_names.get_loc(g) for g in y_genes]
    ysum = np.asarray(x[:, y_idx].sum(axis=1)).ravel()
    df = pd.DataFrame({donor_key: adata.obs[donor_key].values, "total": total, "XIST": xist, "chrY": ysum})
    agg = df.groupby(donor_key, observed=True)[["total", "XIST", "chrY"]].sum()
    agg["XIST_per_10k"] = 1e4 * agg["XIST"] / agg["total"]
    agg["chrY_per_10k"] = 1e4 * agg["chrY"] / agg["total"]
    agg["inferred_sex"] = np.where(agg["chrY_per_10k"] > agg["XIST_per_10k"], "male", "female")
    return agg
