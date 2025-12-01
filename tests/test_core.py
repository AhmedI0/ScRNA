"""Unit tests for the analysis helpers (run with `pytest -q`)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from anndata import AnnData

from ipfatlas import composition, pseudobulk, qc
from ipfatlas.plotting import save_figure


@pytest.fixture
def toy() -> AnnData:
    rng = np.random.default_rng(0)
    n_cells, n_genes = 120, 30
    x = sp.csr_matrix(rng.poisson(2, size=(n_cells, n_genes)).astype(np.float32))
    obs = pd.DataFrame({
        "donor_id": np.repeat([f"d{i}" for i in range(6)], 20),
        "condition": np.repeat(["Control"] * 3 + ["IPF"] * 3, 20),
        "celltype": np.tile(["A", "B"], 60),
        "sex": np.repeat(["male", "female"] * 3, 20),
    }, index=[f"c{i}" for i in range(n_cells)])
    for c in obs:
        obs[c] = obs[c].astype("category")
    var = pd.DataFrame(index=[f"MT-{i}" if i < 3 else f"G{i}" for i in range(n_genes)])
    ad = AnnData(X=x, obs=obs, var=var)
    ad.layers["counts"] = ad.X.copy()
    return ad


def test_is_outlier_flags_extremes():
    x = pd.Series(np.r_[np.ones(50), 100.0])
    flags = qc.is_outlier(x + np.random.default_rng(0).normal(0, 0.01, 51), nmads=5)
    assert flags.iloc[-1] and flags.iloc[:-1].sum() == 0


def test_is_outlier_one_sided():
    x = pd.Series(np.r_[np.linspace(1, 2, 50), -100.0, 100.0])
    assert qc.is_outlier(x, 3, side="upper").iloc[-1]
    assert not qc.is_outlier(x, 3, side="upper").iloc[-2]


def test_gene_groups(toy):
    qc.annotate_gene_groups(toy)
    assert toy.var["mt"].sum() == 3


def test_qc_fail_uses_only_configured_criteria(toy):
    qc.annotate_gene_groups(toy)
    qc.compute_qc_metrics(toy, layer="counts")
    cfg = {"nmads_counts": 4, "nmads_genes": 4, "nmads_top20": 0.01, "nmads_mito": 3, "mito_floor": 15,
           "max_pct_mito": 100, "min_genes": 0, "min_counts": 0}
    qc.flag_outliers(toy, {**cfg, "filter_on": ["counts", "genes", "mito", "hard"]}, groupby="donor_id")
    used = toy.obs[["outlier_counts", "outlier_genes", "outlier_mito", "outlier_hard"]].any(axis=1)
    assert toy.obs["qc_fail"].equals(used)
    only_top20 = toy.obs["outlier_top20"] & ~used   # flagged by the top-20 rule alone ...
    assert only_top20.any() and not toy.obs.loc[only_top20, "qc_fail"].any()   # ... are kept
    qc.flag_outliers(toy, {**cfg, "filter_on": ["top20"]}, groupby="donor_id")
    assert toy.obs["qc_fail"].equals(toy.obs["outlier_top20"])
