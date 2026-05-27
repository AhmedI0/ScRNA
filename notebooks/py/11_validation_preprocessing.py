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
# # 11 · Validation cohort (GSE136831): loading, QC, normalisation and integration
#
# A finding from one cohort may reflect that cohort's donors, tissue handling,
# chemistry or analysis choices. The strongest test is **replication in an
# independent cohort**. Adams *et al.* 2020 profiled 28 control and 32 IPF lungs
# at Yale/Brigham and Women's with 10x 3′ v2 chemistry - a different centre,
# chemistry (3′ vs 5′) and processing pipeline - which makes concordant results
# hard to explain by technical artefacts.
#
# The same processing logic as the discovery cohort is applied. The GEO matrix
# was already filtered by the authors (≥1,000 UMIs per cell, multiplets
# annotated), so we remove the authors' multiplets, exclude COPD donors (not part
# of our question) and apply our mitochondrial and gene/UMI floors.

# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from ipfatlas import features, integration, io, qc
from ipfatlas.config import setup
from ipfatlas.plotting import (
    CONDITION_COLORS,
    embedding_scatter,
    figsize,
    panel_label,
    save_figure,
    set_style,
)
from ipfatlas.utils import write_table

cfg, P = setup()
set_style(cfg["figures"]["font_size"])
FIG = P.figures / "11_validation"
TAB = P.tables / "11_validation"
QC = cfg["qc"]

# %%
adata = io.read_adams(P.raw)
adata = adata[adata.obs["diagnosis"].isin(["Control", "IPF"]) & (adata.obs["author_celltype"] != "Multiplet")].copy()
for c in adata.obs.select_dtypes("category"):
    adata.obs[c] = adata.obs[c].cat.remove_unused_categories()
adata

# %% [markdown]
# ## 1. QC

# %%
qc.compute_qc_metrics(adata, layer=None)
thr = qc.flag_outliers(adata, QC, groupby="library_id")
write_table(thr.round(2), TAB / "validation_qc_thresholds_per_library")
print(adata.obs[["outlier_counts", "outlier_genes", "outlier_top20", "outlier_mito", "outlier_hard", "qc_fail"]].sum())

# %%
fig, axes = plt.subplots(1, 3, figsize=figsize(183, 50))
for ax, col, lab in zip(axes, ["total_counts", "n_genes_by_counts", "pct_counts_mt"],
                        ["UMI counts", "detected genes", "% mitochondrial"]):
    for dx in ["Control", "IPF"]:
        v = adata.obs.loc[adata.obs["diagnosis"] == dx, col]
        ax.hist(np.log10(v + 1) if col != "pct_counts_mt" else v, bins=80, histtype="step",
                color=CONDITION_COLORS[dx], label=dx, density=True)
    ax.set_xlabel(("log10 " if col != "pct_counts_mt" else "") + lab)
    ax.set_ylabel("density")
axes[0].legend()
for a, l in zip(axes, "abc"):
    panel_label(a, l)
fig.tight_layout()
save_figure(fig, "S11_validation_qc_distributions", FIG)
plt.show()

# %%
n0 = adata.n_obs
adata = adata[~adata.obs["qc_fail"]].copy()
sc.pp.filter_genes(adata, min_cells=QC["min_cells_per_gene"])
adata.layers["counts"] = adata.X.copy()
print(f"{n0:,} -> {adata.n_obs:,} cells; {adata.n_vars:,} genes")

# %% [markdown]
# ## 2. Donor sex

# %%
sex = qc.infer_sex(adata, donor_key="donor_id", layer="counts")
sex["diagnosis"] = adata.obs.groupby("donor_id", observed=True)["diagnosis"].first()
write_table(sex, TAB / "validation_inferred_sex")
adata.obs["sex"] = adata.obs["donor_id"].map(sex["inferred_sex"]).astype("category")
pd.crosstab(sex["diagnosis"], sex["inferred_sex"])

# %% [markdown]
# ## 3. Normalisation, HVGs and scVI integration (batch = library)

# %%
sc.pp.normalize_total(adata, target_sum=None)
sc.pp.log1p(adata)
adata.layers["log1p_norm"] = adata.X.copy()
hvg = features.select_hvgs(adata, cfg["features"], batch_key="library_id")

hv = adata[:, adata.var["highly_variable"]].copy()
model = integration.train_scvi(hv, "library_id", cfg["integration"]["scvi"], layer="counts",
                               save_dir=P.processed / "models" / "scvi_validation")
adata.obsm["X_scVI"] = model.get_latent_representation()
del hv
sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=15, random_state=0)
sc.tl.umap(adata, min_dist=0.3, random_state=0)

# %%
fig, axes = plt.subplots(1, 3, figsize=figsize(183, 62))
xy = adata.obsm["X_umap"]
embedding_scatter(axes[0], xy, adata.obs["library_id"], s=0.2)
axes[0].set_title("library")
embedding_scatter(axes[1], xy, adata.obs["diagnosis"], palette=CONDITION_COLORS, s=0.2, legend=True)
axes[1].set_title("diagnosis")
embedding_scatter(axes[2], xy, adata.obs["harmonised_celltype"], s=0.2)
axes[2].set_title("author cell type (harmonised)")
for a, l in zip(axes, "abc"):
    panel_label(a, l, x=0.0)
fig.tight_layout()
save_figure(fig, "S11_validation_umap_overview", FIG)
plt.show()

# %%
summary = adata.obs.groupby("diagnosis", observed=True).agg(donors=("donor_id", "nunique"),
                                                           libraries=("library_id", "nunique"),
                                                           cells=("donor_id", "size"))
write_table(summary, TAB / "validation_cohort_summary")
adata.X = adata.layers["log1p_norm"]
del adata.layers["log1p_norm"]  # .X holds it; keeps the file smaller
adata.write_h5ad(P.processed / "validation_processed.h5ad", compression="lzf")
summary
