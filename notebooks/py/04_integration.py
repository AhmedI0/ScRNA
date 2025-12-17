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
# # 04 · Data integration and benchmarking
#
# Three representations are compared:
#
# | method | idea | input |
# |---|---|---|
# | **PCA** (unintegrated) | linear, no correction | scaled log-normalised HVGs |
# | **Harmony** (Korsunsky *et al.* 2019) | iterative soft clustering in PC space; removes library-specific centroid shifts | PCA |
# | **scVI** (Lopez *et al.* 2018) | variational autoencoder with negative-binomial likelihood; batch is a decoder covariate, so the latent space is batch-conditional | raw counts of HVGs |
#
# They are scored with **scib-metrics** (Luecken *et al.* 2022, *Nat Methods*):
# *batch correction* (iLISI, kBET, graph connectivity, PCR comparison, batch
# silhouette) and *bio-conservation* (isolated labels, Leiden NMI/ARI,
# cell-type silhouette, cLISI). The authors' published labels serve as an
# **independent** ground truth (our own labels are not yet defined, and using
# them would be circular). The batch unit is the 10x library.

# %%
import matplotlib.pyplot as plt
import scanpy as sc

from ipfatlas import integration
from ipfatlas.config import setup
from ipfatlas.plotting import (
    CONDITION_COLORS,
    embedding_scatter,
    figsize,
    palette_for,
    panel_label,
    save_figure,
    set_style,
)
from ipfatlas.utils import write_table

cfg, P = setup()
set_style(cfg["figures"]["font_size"])
FIG = P.figures / "04_integration"
TAB = P.tables / "04_integration"
I = cfg["integration"]
BATCH = I["batch_key"]
MODELS = P.processed / "models"

adata = sc.read_h5ad(P.interim / "discovery_norm.h5ad")
import torch

print("CUDA available:", torch.cuda.is_available(), torch.cuda.get_device_name(0) if torch.cuda.is_available() else "")
adata

# %% [markdown]
# ## 1. Harmony

# %%
adata.obsm["X_harmony"] = integration.run_harmony(adata, BATCH, theta=I["harmony"]["theta"])

# %% [markdown]
# ## 2. scVI
# The model sees raw counts of the 3,000 HVGs. Training stops early when the
# validation ELBO stops improving.

# %%
hv = adata[:, adata.var["highly_variable"]].copy()
model = integration.train_scvi(hv, BATCH, I["scvi"], layer="counts", save_dir=MODELS / "scvi_discovery")
adata.obsm["X_scVI"] = model.get_latent_representation()

hist = model.history
fig, ax = plt.subplots(figsize=figsize(70, 50))
ax.plot(hist["elbo_train"].index, hist["elbo_train"].values.ravel(), label="train")
if "elbo_validation" in hist:
    ax.plot(hist["elbo_validation"].index, hist["elbo_validation"].values.ravel(), label="validation")
ax.set_xlabel("epoch")
ax.set_ylabel("negative ELBO")
ax.set_title("scVI training curve")
ax.legend()
save_figure(fig, "S4_scvi_training_curve", FIG)
plt.show()

# %% [markdown]
# ## 3. Benchmark with scib-metrics
# Run on a random 40,000-cell subsample of cells that carry an author label.

# %%
lab = adata[adata.obs["author_celltype"].notna()].copy()
res = integration.benchmark(lab, BATCH, "author_celltype", ["X_pca", "X_harmony", "X_scVI"],
                            subsample=I["benchmark_subsample"], n_jobs=8)
write_table(res, TAB / "scib_metrics_benchmark")
res

# %%
scores = res.drop(index="Metric Type", errors="ignore").astype(float)
mtype = res.loc["Metric Type"] if "Metric Type" in res.index else None
metrics = [c for c in scores.columns if c not in ("Batch correction", "Bio conservation", "Total")]
name_map = {"X_pca": "Unintegrated (PCA)", "X_harmony": "Harmony", "X_scVI": "scVI"}
scores.index = [name_map.get(i, i) for i in scores.index]

fig, axes = plt.subplots(1, 2, figsize=figsize(183, 60), gridspec_kw={"width_ratios": [3, 1]})
ax = axes[0]
im = ax.imshow(scores[metrics].values, cmap="viridis", vmin=0, vmax=1, aspect="auto")
ax.set_xticks(range(len(metrics)))
ax.set_xticklabels(metrics, rotation=45, ha="right")
ax.set_yticks(range(len(scores)))
ax.set_yticklabels(scores.index)
for i in range(scores.shape[0]):
    for j, m in enumerate(metrics):
        ax.text(j, i, f"{scores.iloc[i][m]:.2f}", ha="center", va="center", fontsize=5,
                color="white" if scores.iloc[i][m] < 0.6 else "black")
plt.colorbar(im, ax=ax, shrink=0.8, label="score")
ax.set_title("scib-metrics (unscaled)")
panel_label(ax, "a", x=-0.15)
ax = axes[1]
agg = scores[["Batch correction", "Bio conservation", "Total"]]
agg.plot.barh(ax=ax, color=["#4C78A8", "#54A24B", "#333333"], width=0.8)
ax.set_xlim(0, 1)
ax.set_xlabel("aggregate score")
ax.legend(fontsize=5, loc="lower right")
panel_label(ax, "b", x=-0.6)
fig.tight_layout()
save_figure(fig, "S4_integration_benchmark", FIG)
plt.show()

best_key = {v: k for k, v in name_map.items()}[scores["Total"].idxmax()]
print("best overall representation:", best_key)

# %% [markdown]
# ## 4. Visual comparison

# %%
emb = {}
for key in ["X_pca", "X_harmony", "X_scVI"]:
    tmp = sc.AnnData(obs=adata.obs[[]].copy())
    tmp.obsm[key] = adata.obsm[key]
    sc.pp.neighbors(tmp, use_rep=key, n_neighbors=15, random_state=0)
    sc.tl.umap(tmp, min_dist=0.3, random_state=0)
    emb[key] = tmp.obsm["X_umap"]
    adata.obsm[f"X_umap_{key.replace('X_', '')}"] = tmp.obsm["X_umap"]

ac = adata.obs["author_celltype"].cat.add_categories("not in author set").fillna("not in author set")
ac_pal = palette_for(ac.cat.categories, {"not in author set": "#DDDDDD"})
fig, axes = plt.subplots(3, 3, figsize=figsize(183, 175))
for r, key in enumerate(["X_pca", "X_harmony", "X_scVI"]):
    embedding_scatter(axes[r, 0], emb[key], adata.obs[BATCH], s=0.2)
    embedding_scatter(axes[r, 1], emb[key], adata.obs["diagnosis"], palette=CONDITION_COLORS, s=0.2,
                      legend=(r == 0))
    embedding_scatter(axes[r, 2], emb[key], ac, palette=ac_pal, s=0.2)
    axes[r, 0].set_ylabel(name_map[key], fontsize=7)
    axes[r, 0].yaxis.set_visible(True)
    axes[r, 0].set_yticks([])
for c, t in enumerate(["library", "diagnosis", "author cell type"]):
    axes[0, c].set_title(t)
fig.tight_layout()
save_figure(fig, "S4_umap_methods_grid", FIG)
plt.show()

# %% [markdown]
# **Choice.** The representation with the best overall scib score is used for all
# downstream neighbourhood-based steps (clustering, UMAP, trajectories). Note that
# integrated embeddings are used *only* for these graph-based steps; differential
# expression is always computed on raw counts with library/donor-aware models,
# never on batch-corrected values.

# %%
adata.obsm["X_integrated"] = adata.obsm[best_key]
adata.uns["integration"] = {"method": best_key, "batch_key": BATCH}
adata.obsm["X_umap"] = adata.obsm[f"X_umap_{best_key.replace('X_', '')}"]
adata.write_h5ad(P.interim / "discovery_integrated.h5ad", compression="lzf")
adata
