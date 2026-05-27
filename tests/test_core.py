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


def test_pseudobulk_sums_counts(toy):
    pb = pseudobulk.aggregate(toy, "celltype", "donor_id", "condition", covariates=["sex"], min_cells=5)
    counts, meta = pb["A"]
    mask = (toy.obs["celltype"] == "A") & (toy.obs["donor_id"] == "d0")
    expected = np.asarray(toy.layers["counts"][mask.to_numpy()].sum(axis=0)).ravel()
    np.testing.assert_array_equal(counts["d0"].to_numpy(), expected.round().astype(int))
    assert set(meta.columns) >= {"condition", "sex", "n_cells"}
    assert meta.loc["d0", "n_cells"] == 10


def test_pseudobulk_min_cells(toy):
    pb = pseudobulk.aggregate(toy, "celltype", "donor_id", "condition", min_cells=11)
    assert pb == {}


def test_safe_key():
    assert pseudobulk.safe_key("KRT5-/KRT17+ cells") == "KRT5_KRT17_cells"


def test_clr_rows_sum_to_zero():
    counts = pd.DataFrame([[10, 5, 1], [3, 3, 3]], columns=list("abc"))
    np.testing.assert_allclose(composition.clr(counts).sum(axis=1), 0, atol=1e-12)


def test_save_figure_writes_all_formats(tmp_path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    paths = save_figure(fig, "t", tmp_path, close=True)
    assert {p.suffix for p in paths} == {".svg", ".jpg", ".png"}
    from PIL import Image

    with Image.open(tmp_path / "t.jpg") as im:
        assert round(im.info["dpi"][0]) == 300


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


def test_assign_by_scores_picks_best_and_margin():
    from ipfatlas.annotation import assign_by_scores

    s = pd.DataFrame({"AT2": [0.9, 0.2], "AT1": [0.1, 0.25]}, index=["c0", "c1"])
    out = assign_by_scores(s, min_margin=0.1)
    assert out.loc["c0", "marker_label"] == "AT2" and out.loc["c0", "confident"]
    assert out.loc["c1", "marker_label"] == "AT1" and not out.loc["c1", "confident"]


def test_harmonise_maps_fine_labels():
    from ipfatlas.annotation import harmonise

    h = harmonise(pd.Series(["Transitional AT2", "Aerocyte capillary EC", "SPP1hi macrophages", "unknown"]))
    assert list(h.astype(str)) == ["AT2", "Capillary EC", "Monocyte-derived macrophages", "Other"]


def test_compartment_palette_groups_hues():
    from ipfatlas.plotting import compartment_palette

    pal = compartment_palette({"AT1": "Epithelial", "AT2": "Epithelial", "B cells": "Immune"})
    assert len(set(pal.values())) == 3 and pal["AT1"] != pal["AT2"]


def test_detection_scores_ignore_expression_level():
    from ipfatlas.annotation import detection_scores

    X = sp.csr_matrix(np.array([[50, 0, 0], [80, 3, 0], [0, 0, 1], [0, 0, 2]], dtype=float))
    ad = AnnData(X, obs=pd.DataFrame({"cl": ["a", "a", "b", "b"]}, index=[f"c{i}" for i in range(4)]),
                 var=pd.DataFrame(index=["G1", "G2", "G3"]))
    d = detection_scores(ad, "cl", {"S1": ["G1", "G2"], "S2": ["G3", "missing"]})
    assert d.loc["a", "S1"] == pytest.approx(0.75) and d.loc["b", "S1"] == 0
    assert d.loc["b", "S2"] == 1 and d.loc["a", "S2"] == 0


def test_celltypist_map_targets_known_labels():
    from ipfatlas.annotation import CELLTYPIST_TO_LABEL, FINE_TO_HARMONISED

    assert set(CELLTYPIST_TO_LABEL.values()) <= set(FINE_TO_HARMONISED)


def test_write_table_gzips_large_tables(tmp_path):
    from ipfatlas.utils import write_table

    assert write_table(pd.DataFrame({"a": range(10)}), tmp_path / "small", index=False).name == "small.csv"
    big = pd.DataFrame({"a": range(50_001)})
    out = write_table(big, tmp_path / "big", index=False)
    assert out.name == "big.csv.gz" and not (tmp_path / "big.xlsx").exists()
    assert pd.read_csv(out)["a"].sum() == big["a"].sum()
    first = out.read_bytes()
    write_table(big, tmp_path / "big", index=False)
    assert out.read_bytes() == first


def test_overclustering_resolution_follows_celltypist():
    from ipfatlas.annotation import overclustering_resolution

    assert [overclustering_resolution(n) for n in (1_000, 30_000, 150_000, 227_434)] == [5, 15, 25, 30]
