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
# # 07 · Cell-type-specific differential expression (pseudobulk)
#
# **Why not a per-cell Wilcoxon test?** Cells from the same donor share genotype,
# age, sampling and handling, so they are *not* independent replicates. Treating
# thousands of cells as replicates makes tiny donor-specific differences
# "significant" (pseudoreplication). Squair *et al.* (2021, *Nat Commun*) showed
# that single-cell methods produce hundreds of false discoveries in data with no
# true differences, whereas **pseudobulk** methods - summing counts per donor and
# cell type and analysing them with bulk RNA-seq GLMs - control the FDR.
#
# **Models** (per cell type, donors with ≥10 cells of that type, ≥3 donors per group):
# * **edgeR quasi-likelihood** (`glmQLFit`, robust dispersion; primary)
# * **DESeq2** Wald test with **apeglm** log-fold-change shrinkage (cross-check)
# * design `~ sex + diagnosis`; genes pre-filtered with `filterByExpr`; TMM normalisation
#
# **Sensitivity analyses** for the site-disease confounding described in Notebook 01:
# `~ site + sex + diagnosis` and a Vanderbilt-only analysis.

# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from adjustText import adjust_text

from ipfatlas import pseudobulk
from ipfatlas.config import rscript, setup
from ipfatlas.plotting import COMPARTMENT_COLORS, figsize, panel_label, save_figure, set_style
from ipfatlas.utils import write_table

cfg, P = setup()
set_style(cfg["figures"]["font_size"])
FIG = P.figures / "07_de"
TAB = P.tables / "07_de"
C, D = cfg["contrast"], cfg["de"]
CASE, CTRL = C["case"], C["control"]
R = rscript(cfg)
WORK = P.interim / "pseudobulk"

adata = sc.read_h5ad(P.processed / "discovery_annotated.h5ad")
sub = adata[adata.obs["diagnosis"].isin([CASE, CTRL])].copy()
ct_comp = adata.obs.groupby("cell_type", observed=True)["compartment"].first()
del adata
sub

# %% [markdown]
# ## 1. Pseudobulk aggregation

# %%
pb = pseudobulk.aggregate(sub, "cell_type", "donor_id", "diagnosis", covariates=["sex", "site"],
                          min_cells=D["min_cells_per_pseudobulk"])
index = pseudobulk.write_inputs(pb, WORK / "main_in", CASE, CTRL, min_donors=D["min_donors_per_group"])
write_table(index, TAB / "pseudobulk_samples_per_celltype", index=False)
index

# %% [markdown]
# ## 2. edgeR-QL and DESeq2: `~ sex + diagnosis`

# %%
summary = pseudobulk.run_de(R, P.r_scripts / "pseudobulk_de.R", WORK / "main_in", WORK / "main_out",
                            covariates=["sex"], case=CASE, control=CTRL)
edger = pseudobulk.read_results(WORK / "main_out", "edgeR")
deseq = pseudobulk.read_results(WORK / "main_out", "DESeq2")
edger["significant"] = edger["FDR"] < D["fdr"]
edger["strong"] = edger["significant"] & (edger["logFC"].abs() > D["lfc"])
edger["compartment"] = edger["celltype"].map(ct_comp)
write_table(edger, TAB / "edgeR_IPF_vs_control_all_celltypes", index=False)
write_table(deseq, TAB / "DESeq2_IPF_vs_control_all_celltypes", index=False)
summary

# %% [markdown]
# ## 3. Agreement between edgeR and DESeq2

# %%
m = edger.merge(deseq, on=["celltype", "gene"], suffixes=("_edgeR", "_DESeq2"))
agree = m.groupby("celltype").apply(lambda d: pd.Series({
    "pearson_logFC": np.corrcoef(d["logFC"], d["log2FC"])[0, 1],
    "n_edgeR": int((d["FDR"] < 0.05).sum()),
    "n_DESeq2": int((d["padj"] < 0.05).sum()),
    "n_both": int(((d["FDR"] < 0.05) & (d["padj"] < 0.05)).sum()),
}), include_groups=False)
agree["jaccard"] = agree["n_both"] / (agree["n_edgeR"] + agree["n_DESeq2"] - agree["n_both"]).replace(0, np.nan)
write_table(agree.round(3), TAB / "edgeR_vs_DESeq2_agreement")
agree.round(3)

# %% [markdown]
# ## 4. Overview: DEGs per cell type

# %%
cnt = edger[edger["significant"]].groupby("celltype").agg(up=("logFC", lambda x: (x > 0).sum()),
                                                          down=("logFC", lambda x: (x < 0).sum()))
cnt = cnt.reindex(summary["celltype"]).fillna(0).astype(int)
cnt["total"] = cnt.sum(axis=1)
cnt = cnt.sort_values("total")
n_donors = summary.set_index("celltype")[["n_case", "n_control"]]

fig, ax = plt.subplots(figsize=figsize(95, 4 + 3.2 * len(cnt)))
y = np.arange(len(cnt))
ax.barh(y, cnt["up"], color="#D62728", label=f"up in {CASE}")
ax.barh(y, -cnt["down"], color="#1F77B4", label=f"down in {CASE}")
ax.set_yticks(y)
ax.set_yticklabels([f"{c}  ({n_donors.loc[c, 'n_case']}/{n_donors.loc[c, 'n_control']})" for c in cnt.index])
for t, c in zip(ax.get_yticklabels(), cnt.index):
    t.set_color(COMPARTMENT_COLORS.get(ct_comp.get(c), "black"))
ax.axvline(0, color="black", lw=0.5)
ax.set_xlabel(f"DEGs (edgeR-QL, FDR < {D['fdr']})")
ax.set_title(f"{CASE} vs {CTRL}: pseudobulk DEGs\n(donors {CASE}/{CTRL} in brackets)")
ax.legend(loc="lower right")
xmax = max(cnt["up"].max(), cnt["down"].max()) * 1.1
ax.set_xlim(-xmax, xmax)
save_figure(fig, "Fig3_deg_counts_per_celltype", FIG)
plt.show()

# %% [markdown]
# ## 5. Volcano plots for the cell types with most DEGs

# %%
top_ct = cnt.sort_values("total", ascending=False).index[:6].tolist()
fig, axes = plt.subplots(2, 3, figsize=figsize(183, 115))
for ax, ct in zip(axes.ravel(), top_ct):
    d = edger[edger["celltype"] == ct].copy()
    d["nlp"] = -np.log10(d["PValue"].clip(lower=1e-300))
    col = np.where(d["strong"] & (d["logFC"] > 0), "#D62728",
                   np.where(d["strong"] & (d["logFC"] < 0), "#1F77B4",
                            np.where(d["significant"], "#999999", "#DDDDDD")))
    ax.scatter(d["logFC"], d["nlp"], s=2, c=col, lw=0, rasterized=True)
    lab = pd.concat([d[d["strong"] & (d["logFC"] > 0)].nsmallest(8, "PValue"),
                     d[d["strong"] & (d["logFC"] < 0)].nsmallest(6, "PValue")])
    texts = [ax.text(r.logFC, r.nlp, r.gene, fontsize=4.5, style="italic") for r in lab.itertuples()]
    adjust_text(texts, ax=ax, arrowprops={"arrowstyle": "-", "lw": 0.3, "color": "grey"})
    ax.set_title(f"{ct} ({int(d['significant'].sum())} DEGs)")
    ax.set_xlabel("log2 fold change (IPF / control)")
    ax.set_ylabel("-log10 P")
for a, l in zip(axes.ravel(), "abcdef"):
    panel_label(a, l)
fig.tight_layout()
save_figure(fig, "Fig3_volcano_top_celltypes", FIG)
plt.show()

# %% [markdown]
# ## 6. Shared vs cell-type-specific responses
# Genes that are DE in many cell types point to tissue-wide programmes (e.g.
# interferon, hypoxia) or residual ambient RNA; cell-type-specific DEGs are more
# likely to reflect cell-intrinsic disease biology.

# %%
sig = edger[edger["strong"]]
n_ct = sig.groupby("gene")["celltype"].nunique().sort_values(ascending=False)
write_table(n_ct.rename("n_celltypes").to_frame(), TAB / "deg_recurrence_across_celltypes")
shared = n_ct.index[:40]
lfc = edger.pivot_table(index="gene", columns="celltype", values="logFC").reindex(shared)
fdr = edger.pivot_table(index="gene", columns="celltype", values="FDR").reindex(shared)
lfc = lfc.loc[:, lfc.notna().sum() > 0]
g = sns.clustermap(lfc.fillna(0), cmap="RdBu_r", center=0, vmin=-3, vmax=3, figsize=(7.2, 6.5),
                   linewidths=0, cbar_kws={"label": "log2 FC"}, xticklabels=True, yticklabels=True,
                   dendrogram_ratio=0.08)
g.ax_heatmap.tick_params(labelsize=5)
g.ax_heatmap.set_xlabel("")
g.ax_heatmap.set_ylabel("")
g.fig.suptitle("Genes differentially expressed in the most cell types", y=1.01, fontsize=8)
save_figure(g.fig, "S7_shared_degs_heatmap", FIG)
plt.show()

# %% [markdown]
# ## 7. Sensitivity analyses for site confounding

# %%
pseudobulk.run_de(R, P.r_scripts / "pseudobulk_de.R", WORK / "main_in", WORK / "site_out",
                  covariates=["site", "sex"], case=CASE, control=CTRL, methods="edgeR")
site = pseudobulk.read_results(WORK / "site_out", "edgeR")

vu = sub[sub.obs["site"] == "Vanderbilt"]
pb_vu = pseudobulk.aggregate(vu, "cell_type", "donor_id", "diagnosis", covariates=["sex"], min_cells=D["min_cells_per_pseudobulk"])
pseudobulk.write_inputs(pb_vu, WORK / "vu_in", CASE, CTRL, min_donors=D["min_donors_per_group"])
pseudobulk.run_de(R, P.r_scripts / "pseudobulk_de.R", WORK / "vu_in", WORK / "vu_out",
                  covariates=["sex"], case=CASE, control=CTRL, methods="edgeR")
vuo = pseudobulk.read_results(WORK / "vu_out", "edgeR")


def concordance(a, b, name):
    mm = a.merge(b, on=["celltype", "gene"], suffixes=("", "_s"))
    rows = []
    for ct, d in mm.groupby("celltype"):
        s = d[d["FDR"] < 0.05]
        rows.append({"celltype": ct, "analysis": name, "n_main_DEGs": len(s),
                     "logFC_pearson_all": np.corrcoef(d["logFC"], d["logFC_s"])[0, 1],
                     "sign_agreement_main_DEGs": (np.sign(s["logFC"]) == np.sign(s["logFC_s"])).mean() if len(s) else np.nan,
                     "replicated_P05": ((s["PValue_s"] < 0.05) & (np.sign(s["logFC"]) == np.sign(s["logFC_s"]))).mean() if len(s) else np.nan})
    return pd.DataFrame(rows)


sens = pd.concat([concordance(edger, site, "~ site + sex + diagnosis"),
                  concordance(edger, vuo, "Vanderbilt only")])
write_table(sens.round(3), TAB / "de_sensitivity_site", index=False)
sens.pivot(index="celltype", columns="analysis", values="sign_agreement_main_DEGs").round(3)

# %%
fig, axes = plt.subplots(1, 2, figsize=figsize(150, 60), sharey=True)
for ax, (name, d) in zip(axes, sens.groupby("analysis")):
    d = d.dropna(subset=["sign_agreement_main_DEGs"]).sort_values("n_main_DEGs")
    d = d[d["n_main_DEGs"] >= 10]
    ax.barh(d["celltype"], d["sign_agreement_main_DEGs"] * 100, color="#7F7F7F", label="same sign")
    ax.barh(d["celltype"], d["replicated_P05"] * 100, color="#2CA02C", label="same sign & P < 0.05")
    ax.set_xlim(0, 100)
    ax.set_xlabel("% of main-analysis DEGs")
    ax.set_title(name)
axes[0].legend(fontsize=5, loc="lower left")
fig.tight_layout()
save_figure(fig, "S7_de_site_sensitivity", FIG)
plt.show()
