# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#   kernelspec:
#     display_name: Python (scrna-ipf)
#     language: python
#     name: scrna-ipf
# ---

# %% [markdown]
# # 02 · Quality control: low-quality cells, doublets and genes
#
# A droplet-based library contains, besides intact single cells,
# **(i)** damaged/dying cells (few genes, high mitochondrial fraction because
# cytoplasmic mRNA leaks out while mitochondrial transcripts are retained),
# **(ii)** droplets with ambient RNA only, and **(iii)** doublets (two cells in
# one droplet). Each artefact can create spurious clusters, so it is removed
# before any biological interpretation.
#
# **Approach**
# 1. per-cell metrics: library size, detected genes, % mitochondrial, % ribosomal,
#    % haemoglobin, % counts in the top-20 genes;
# 2. per-library outlier detection with median ± *n*·MAD (robust to the
#    library-to-library differences in depth seen in Notebook 01);
# 3. doublet detection with **scDblFinder** (primary) and **Scrublet** (secondary);
# 4. removal of genes detected in fewer than 10 cells.

# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from upsetplot import UpSet, from_indicators

from ipfatlas import doublets, qc
from ipfatlas.config import rscript, setup
from ipfatlas.plotting import CONDITION_COLORS, figsize, panel_label, save_figure, set_style
from ipfatlas.utils import log, write_table

cfg, P = setup()
set_style(cfg["figures"]["font_size"])
FIG = P.figures / "02_qc"
TAB = P.tables / "02_qc"
QC = cfg["qc"]

adata = sc.read_h5ad(P.interim / "discovery_raw.h5ad")
adata

# %% [markdown]
# ## 1. QC metrics

# %%
qc.compute_qc_metrics(adata, layer=None)
write_table(qc.qc_summary(adata), TAB / "qc_metrics_per_library_before")
adata.obs[["total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo", "pct_counts_hb",
           "pct_counts_in_top_20_genes"]].describe().round(1)

# %% [markdown]
# ## 2. Distributions per library
# Libraries are ordered by diagnosis; dashed lines mark the hard floors.

# %%
order = (adata.obs.groupby("library_id", observed=True)["diagnosis"].first()
         .sort_values().index.tolist())
pal = {lib: CONDITION_COLORS.get(dx, "grey")
       for lib, dx in adata.obs.groupby("library_id", observed=True)["diagnosis"].first().items()}

fig, axes = plt.subplots(3, 1, figsize=figsize(183, 120), sharex=True)
for ax, col, lab, log_y in [
    (axes[0], "total_counts", "UMI counts", True),
    (axes[1], "n_genes_by_counts", "detected genes", True),
    (axes[2], "pct_counts_mt", "% mitochondrial", False),
]:
    sns.violinplot(data=adata.obs, x="library_id", y=col, order=order, hue="library_id", palette=pal,
                   inner=None, linewidth=0.3, cut=0, density_norm="width", ax=ax, legend=False)
    ax.set_ylabel(lab)
    ax.set_xlabel("")
    if log_y:
        ax.set_yscale("log")
axes[0].axhline(QC["min_counts"], ls="--", lw=0.5, color="black")
axes[1].axhline(QC["min_genes"], ls="--", lw=0.5, color="black")
axes[2].axhline(QC["max_pct_mito"], ls="--", lw=0.5, color="black")
axes[2].tick_params(axis="x", rotation=90, labelsize=4)
for a, l in zip(axes, "abc"):
    panel_label(a, l, x=-0.06)
fig.tight_layout()
save_figure(fig, "S2_qc_violins_per_library", FIG)
plt.show()

# %% [markdown]
# ## 3. Per-library MAD outlier detection

# %%
thr = qc.flag_outliers(adata, QC, groupby="library_id")
write_table(thr.round(2), TAB / "qc_thresholds_per_library")
flag_cols = ["outlier_counts", "outlier_genes", "outlier_top20", "outlier_mito", "outlier_hard"]
print(adata.obs[flag_cols + ["qc_fail"]].sum())
thr.sort_values("pct_fail", ascending=False).head(10)

# %%
fig = plt.figure(figsize=figsize(120, 70))
UpSet(from_indicators(flag_cols, data=adata.obs[flag_cols]), subset_size="count", show_counts=True,
      sort_by="cardinality", min_subset_size=200, element_size=None).plot(fig=fig)
fig.suptitle("Overlap of QC failure criteria (subsets ≥200 cells)", fontsize=7)
save_figure(fig, "S2_qc_failure_upset", FIG)
plt.show()

# %%
fig, axes = plt.subplots(1, 2, figsize=figsize(150, 62))
o = adata.obs.sample(frac=1, random_state=0)
ax = axes[0]
sca = ax.scatter(o["total_counts"], o["n_genes_by_counts"], c=o["pct_counts_mt"], s=0.3, cmap="magma_r",
                 vmax=30, rasterized=True)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("UMI counts")
ax.set_ylabel("detected genes")
plt.colorbar(sca, ax=ax, label="% mitochondrial", shrink=0.8)
panel_label(ax, "a")
ax = axes[1]
ax.scatter(o["total_counts"], o["n_genes_by_counts"], c=np.where(o["qc_fail"], "#D62728", "#BBBBBB"),
           s=0.3, rasterized=True)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("UMI counts")
ax.set_ylabel("detected genes")
ax.set_title(f"removed (red): {int(o['qc_fail'].sum()):,} / {len(o):,} barcodes")
panel_label(ax, "b")
fig.tight_layout()
save_figure(fig, "S2_qc_scatter_counts_genes_mito", FIG)
plt.show()

# %% [markdown]
# ## 4. Filter cells and libraries

# %%
n0 = adata.n_obs
adata = adata[~adata.obs["qc_fail"]].copy()
per_lib = adata.obs["library_id"].value_counts()
small = per_lib[per_lib < QC["min_cells_per_library"]].index.tolist()
if small:
    log.warning(f"libraries below {QC['min_cells_per_library']} cells after QC are removed: {small}")
    adata = adata[~adata.obs["library_id"].isin(small)].copy()
adata.obs["library_id"] = adata.obs["library_id"].cat.remove_unused_categories()
print(f"{n0:,} -> {adata.n_obs:,} cells after QC ({100 * adata.n_obs / n0:.1f} % kept)")
adata.layers["counts"] = adata.X.copy()

# %% [markdown]
# ## 5. Doublet detection
# **scDblFinder** generates artificial doublets by summing random pairs of cells
# within each library and trains a gradient-boosted tree classifier to separate
# them from real droplets using kNN-based features; it was among the most
# accurate methods in an independent benchmark (Xi & Li 2021, *Cell Systems*).
# **Scrublet** uses a similar simulation idea with a kNN doublet score. Doublets
# can only form *within* a 10x channel, so both run per library.

# %%
dbl_r = doublets.run_scdblfinder(adata, rscript(cfg), P.r_scripts / "scdblfinder.R", P.interim / "scdblfinder",
                                 batch_key="library_id", threads=8)
dbl_py = doublets.run_scrublet(adata, batch_key="library_id",
                               expected_rate=cfg["doublets"]["scrublet_expected_rate"])
adata.obs = adata.obs.join(dbl_r).join(dbl_py)
adata.obs["scrublet_doublet"] = adata.obs["scrublet_doublet"].fillna(False).astype(bool)
tab = doublets.concordance(adata.obs)
write_table(tab, TAB / "doublet_concordance")
tab

# %%
rate = adata.obs.groupby("library_id", observed=True).agg(
    n_cells=("library_id", "size"),
    scDblFinder=("scDblFinder_doublet", "mean"),
    Scrublet=("scrublet_doublet", "mean"))
rate[["scDblFinder", "Scrublet"]] *= 100
write_table(rate.round(2), TAB / "doublet_rate_per_library")

fig, axes = plt.subplots(1, 2, figsize=figsize(150, 60))
ax = axes[0]
ax.scatter(rate["n_cells"], rate["scDblFinder"], s=10, label="scDblFinder", color="#4C78A8")
ax.scatter(rate["n_cells"], rate["Scrublet"], s=10, label="Scrublet", color="#F58518", marker="^")
xs = np.linspace(0, rate["n_cells"].max(), 50)
ax.plot(xs, 0.8 * xs / 1000, ls="--", color="grey", lw=0.6, label="10x expectation (0.8 %/1,000 cells)")
ax.set_xlabel("cells per library (after QC)")
ax.set_ylabel("% called doublets")
ax.legend(fontsize=5)
panel_label(ax, "a")
ax = axes[1]
ax.hist(adata.obs["scDblFinder_score"], bins=60, color="#4C78A8")
ax.set_yscale("log")
ax.set_xlabel("scDblFinder score")
ax.set_ylabel("cells")
panel_label(ax, "b")
fig.tight_layout()
save_figure(fig, "S2_doublet_rates", FIG)
plt.show()

# %% [markdown]
# **Decision:** remove scDblFinder doublets (primary method). Scrublet is reported
# for transparency; a final, biology-based doublet check (clusters co-expressing
# two lineage programmes) is done after clustering in Notebook 05.

# %%
n1 = adata.n_obs
adata = adata[~adata.obs["scDblFinder_doublet"]].copy()
print(f"{n1:,} -> {adata.n_obs:,} cells after doublet removal")

# %% [markdown]
# ## 6. Gene filter

# %%
sc.pp.filter_genes(adata, min_cells=QC["min_cells_per_gene"])
print(adata)

# %% [markdown]
# ## 7. Comparison with the authors' QC
# How similar is our cell set to the 114,396 cells retained by Habermann *et al.*?

# %%
kept = pd.Series(True, index=adata.obs_names)
raw_obs = sc.read_h5ad(P.interim / "discovery_raw.h5ad", backed="r").obs
cmp = pd.DataFrame({"ours": raw_obs.index.isin(kept.index), "authors": raw_obs["author_qc_pass"].values},
                   index=raw_obs.index)
ct = pd.crosstab(cmp["ours"], cmp["authors"], rownames=["kept (this study)"], colnames=["kept (authors)"])
write_table(ct, TAB / "qc_overlap_with_authors")
ct

# %%
lost = raw_obs.loc[cmp["authors"] & ~cmp["ours"], "author_celltype"].value_counts()
frac = (lost / raw_obs.loc[cmp["authors"], "author_celltype"].value_counts()).dropna().sort_values(ascending=False) * 100
fig, ax = plt.subplots(figsize=figsize(90, 70))
frac.head(20).iloc[::-1].plot.barh(ax=ax, color="#7F7F7F")
ax.set_xlabel("% of author-retained cells removed by our QC")
ax.set_title("Which author cell types does our stricter QC remove?")
save_figure(fig, "S2_qc_removed_by_author_celltype", FIG)
plt.show()
frac.round(1).head(15)

# %% [markdown]
# ## 8. Final QC summary and save

# %%
write_table(qc.qc_summary(adata), TAB / "qc_metrics_per_library_after")
summary = pd.Series({
    "barcodes (Cell Ranger)": n0,
    "after MAD/hard QC": n1,
    "after doublet removal": adata.n_obs,
    "genes kept": adata.n_vars,
    "libraries": adata.obs.library_id.nunique(),
    "donors": adata.obs.donor_id.nunique(),
})
write_table(summary.to_frame("value"), TAB / "qc_cell_flow")
adata.write_h5ad(P.interim / "discovery_qc.h5ad", compression="lzf")
summary
