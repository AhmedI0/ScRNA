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
# # 06 · Differential cell-type abundance in IPF
#
# **Question:** which cell types/states expand or contract in IPF lungs?
#
# **Statistical considerations**
# * The replicate is the **donor** (n = 10 control vs 12 IPF), not the cell.
# * Proportions are **compositional**: an expansion of one population necessarily
#   lowers the proportions of the others. scCODA models this explicitly
#   (Dirichlet-multinomial with a reference cell type); propeller uses
#   logit-transformed proportions with moderated t-tests.
# * Dissociation and capture efficiency differ between compartments (e.g. immune
#   cells dissociate easily; fibroblasts are under-represented in droplet data).
#   Abundance is therefore also tested **within each compartment** (e.g. AT2 as a
#   fraction of all epithelial cells), which is far less sensitive to such biases.
# * **Milo** (Dann *et al.* 2022) tests differential abundance on overlapping kNN
#   neighbourhoods and does not depend on discrete cluster boundaries - useful for
#   continuous transitions such as AT2 → transitional → aberrant basaloid.

# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

from ipfatlas import composition
from ipfatlas.config import rscript, setup
from ipfatlas.plotting import (
    COMPARTMENT_COLORS,
    CONDITION_COLORS,
    figsize,
    panel_label,
    save_figure,
    set_style,
)
from ipfatlas.utils import log, write_table

cfg, P = setup()
set_style(cfg["figures"]["font_size"])
FIG = P.figures / "06_composition"
TAB = P.tables / "06_composition"
C = cfg["contrast"]
CASE, CTRL = C["case"], C["control"]
R = rscript(cfg)

adata = sc.read_h5ad(P.processed / "discovery_annotated.h5ad")
ct_order = list(adata.obs["cell_type"].cat.categories)
ct_comp = adata.obs.groupby("cell_type", observed=True)["compartment"].first()
ct_pal = dict(zip(ct_order, adata.uns["cell_type_colors"]))
adata

# %% [markdown]
# ## 1. Composition of every donor (all diagnoses)

# %%
counts_all, meta_all = composition.donor_counts(adata, "cell_type", "donor_id", "diagnosis", extra=["site", "sex"])
prop_all = counts_all.div(counts_all.sum(axis=1), axis=0)
dx_order = ["Control", "IPF", "NSIP", "cHP", "Sarcoidosis", "ILD"]
donor_order = meta_all.assign(dx=pd.Categorical(meta_all["diagnosis"], dx_order)).sort_values(["dx"]).index
comp_prop = (counts_all.T.groupby(ct_comp).sum().T).div(counts_all.sum(axis=1), axis=0)

fig, axes = plt.subplots(2, 1, figsize=figsize(183, 110), sharex=True, gridspec_kw={"height_ratios": [1, 2.2]})
bottom = np.zeros(len(donor_order))
for comp in ["Epithelial", "Immune", "Endothelial", "Mesenchymal"]:
    v = comp_prop.loc[donor_order, comp].to_numpy()
    axes[0].bar(range(len(donor_order)), v, bottom=bottom, color=COMPARTMENT_COLORS[comp], label=comp, width=0.85)
    bottom += v
axes[0].set_ylabel("fraction of cells")
axes[0].legend(ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.25))
bottom = np.zeros(len(donor_order))
for ct in ct_order:
    v = prop_all.loc[donor_order, ct].to_numpy()
    axes[1].bar(range(len(donor_order)), v, bottom=bottom, color=ct_pal[ct], label=ct, width=0.85)
    bottom += v
axes[1].set_ylabel("fraction of cells")
axes[1].legend(ncol=1, fontsize=4.5, loc="center left", bbox_to_anchor=(1.0, 0.5))
axes[1].set_xticks(range(len(donor_order)))
axes[1].set_xticklabels([f"{d} ({meta_all.loc[d, 'diagnosis']})" for d in donor_order], rotation=90, fontsize=4.5)
for t, d in zip(axes[1].get_xticklabels(), donor_order):
    t.set_color(CONDITION_COLORS.get(meta_all.loc[d, "diagnosis"], "black"))
panel_label(axes[0], "a", x=-0.06)
panel_label(axes[1], "b", x=-0.06)
fig.tight_layout()
save_figure(fig, "S6_composition_per_donor", FIG)
plt.show()
write_table(prop_all.round(4), TAB / "cell_type_proportions_per_donor")

# %% [markdown]
# ## 2. propeller — all cells and within compartments (IPF vs control)

# %%
sub = adata[adata.obs["diagnosis"].isin([CASE, CTRL])].copy()
for c in ["cell_type", "donor_id", "diagnosis"]:
    sub.obs[c] = sub.obs[c].cat.remove_unused_categories()
print(sub.obs.groupby("diagnosis", observed=True)["donor_id"].nunique())

res = {}
res["all cells"] = composition.run_propeller(sub, R, P.r_scripts / "propeller.R", P.interim / "propeller",
                                             "cell_type", "donor_id", "diagnosis", CASE, CTRL,
                                             covariates=["sex"], tag="_all")
for comp in ["Epithelial", "Immune", "Endothelial", "Mesenchymal"]:
    s = sub[sub.obs["compartment"] == comp].copy()
    s.obs["cell_type"] = s.obs["cell_type"].cat.remove_unused_categories()
    res[comp] = composition.run_propeller(s, R, P.r_scripts / "propeller.R", P.interim / "propeller",
                                          "cell_type", "donor_id", "diagnosis", CASE, CTRL,
                                          covariates=["sex"], tag=f"_{comp}")
prop_tab = []
for scope, r in res.items():
    r = r.copy()
    r["scope"] = scope
    prop_tab.append(r)
prop_tab = pd.concat(prop_tab, ignore_index=True)
prop_tab["log2_prop_ratio"] = np.log2(prop_tab["PropRatio"])
write_table(prop_tab, TAB / "propeller_results", index=False)
within = prop_tab[prop_tab["scope"] != "all cells"].set_index("celltype")
within[["log2_prop_ratio", "P.Value", "FDR"]].sort_values("FDR").round(4)

# %% [markdown]
# ## 3. scCODA — Bayesian compositional model within compartments
# The reference cell type is chosen automatically as the type with the least
# variable relative abundance across donors (present in ≥95 % of donors). An
# effect is "credible" if its spike-and-slab inclusion probability exceeds the
# threshold controlling the expected false-discovery rate at 10 %.

# %%
sccoda = []
for comp in ["Epithelial", "Immune", "Endothelial", "Mesenchymal"]:
    s = sub[sub.obs["compartment"] == comp]
    cnt, met = composition.donor_counts(s, "cell_type", "donor_id", "diagnosis")
    cnt = cnt.loc[:, cnt.sum() > 0]
    try:
        eff = composition.run_sccoda(cnt, met, "diagnosis", CASE, CTRL, fdr=cfg["composition"]["sccoda_fdr"])
        eff["scope"] = comp
        sccoda.append(eff)
    except Exception as e:  # report and continue: scCODA can fail for degenerate compartments
        log.warning(f"scCODA failed for {comp}: {e}")
sccoda = pd.concat(sccoda, ignore_index=True)
sccoda["credible"] = sccoda["Final Parameter"] != 0
write_table(sccoda, TAB / "sccoda_results", index=False)
sccoda[["scope", "celltype", "Final Parameter", "Inclusion probability", "log2-fold change", "credible",
        "reference_cell_type"]].round(3)

# %% [markdown]
# ## 4. Consensus table

# %%
cons = within[["log2_prop_ratio", "FDR"]].rename(columns={"FDR": "propeller_FDR_within"})
cons["propeller_FDR_all"] = prop_tab[prop_tab["scope"] == "all cells"].set_index("celltype")["FDR"]
sc_idx = sccoda.set_index("celltype")
cons["scCODA_log2FC"] = sc_idx["log2-fold change"]
cons["scCODA_credible"] = sc_idx["credible"]
cons["compartment"] = ct_comp
cons["robust"] = (cons["propeller_FDR_within"] < 0.05) & cons["scCODA_credible"].fillna(False).astype(bool)
cons["direction"] = np.where(cons["log2_prop_ratio"] > 0, f"up in {CASE}", f"down in {CASE}")
cons = cons.sort_values(["compartment", "log2_prop_ratio"])
write_table(cons, TAB / "composition_consensus")
cons.round(3)

# %%
plot_types = cons.sort_values("propeller_FDR_within").index[:12].tolist()
long = composition.proportions_long(*composition.donor_counts(sub, "cell_type", "donor_id", "diagnosis"), "diagnosis")
# proportions within compartment for plotting
cnt_s, met_s = composition.donor_counts(sub, "cell_type", "donor_id", "diagnosis")
comp_tot = cnt_s.T.groupby(ct_comp).sum().T
within_prop = cnt_s.copy().astype(float)
for ct in cnt_s.columns:
    within_prop[ct] = cnt_s[ct] / comp_tot[ct_comp[ct]].replace(0, np.nan)
wl = within_prop.reset_index().melt(id_vars="donor_id", var_name="cell_type", value_name="fraction")
wl["diagnosis"] = wl["donor_id"].map(met_s["diagnosis"])

ncol = 6
nrow = int(np.ceil(len(plot_types) / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=figsize(183, 45 * nrow))
for ax, ct in zip(axes.ravel(), plot_types):
    d = wl[wl["cell_type"] == ct]
    sns.boxplot(data=d, x="diagnosis", y="fraction", order=[CTRL, CASE], hue="diagnosis", palette=CONDITION_COLORS,
                width=0.6, fliersize=0, linewidth=0.6, ax=ax, legend=False)
    sns.stripplot(data=d, x="diagnosis", y="fraction", order=[CTRL, CASE], color="black", size=2, ax=ax, jitter=0.15)
    fdr = cons.loc[ct, "propeller_FDR_within"]
    star = "✓ scCODA" if bool(cons.loc[ct, "scCODA_credible"]) else ""
    ax.set_title(f"{ct}\nFDR={fdr:.1e} {star}", fontsize=6)
    ax.set_xlabel("")
    ax.set_ylabel(f"fraction of {ct_comp[ct].lower()}" if ax in axes[:, 0] else "")
for ax in axes.ravel()[len(plot_types):]:
    ax.axis("off")
fig.tight_layout()
save_figure(fig, "Fig2_composition_top_celltypes", FIG)
plt.show()

# %% [markdown]
# ## 5. Milo — neighbourhood-level differential abundance

# %%
import pertpy as pt

milo = pt.tl.Milo()
msub = sub.copy()
sc.pp.neighbors(msub, use_rep="X_integrated", n_neighbors=30, random_state=0)
mdata = milo.load(msub)
milo.make_nhoods(mdata["rna"], prop=0.05, seed=0)
mdata = milo.count_nhoods(mdata, sample_col="donor_id")
mdata["rna"].obs["diagnosis"] = pd.Categorical(mdata["rna"].obs["diagnosis"], categories=[CTRL, CASE])
milo.da_nhoods(mdata, design="~sex+diagnosis", model_contrasts=f"diagnosis{CASE}", solver="pydeseq2")
milo.annotate_nhoods(mdata, anno_col="cell_type")
nh = mdata["milo"].var.copy()
nh.loc[nh["nhood_annotation_frac"] < 0.6, "nhood_annotation"] = "Mixed"
mdata["milo"].var["nhood_annotation"] = nh["nhood_annotation"]
write_table(nh, TAB / "milo_neighbourhoods")
print(f"{len(nh)} neighbourhoods, {(nh['SpatialFDR'] < 0.1).sum()} with SpatialFDR < 0.1")

# %%
d = nh[nh["nhood_annotation"] != "Mixed"].copy()
order = d.groupby("nhood_annotation")["logFC"].median().sort_values().index
fig, ax = plt.subplots(figsize=figsize(110, 120))
rng = np.random.default_rng(0)
for i, ct in enumerate(order):
    v = d[d["nhood_annotation"] == ct]
    y = i + rng.uniform(-0.3, 0.3, len(v))
    sig = v["SpatialFDR"] < 0.1
    ax.scatter(v.loc[~sig, "logFC"], y[~sig.to_numpy()], s=3, color="#CCCCCC", lw=0)
    ax.scatter(v.loc[sig, "logFC"], y[sig.to_numpy()], s=3, c=v.loc[sig, "logFC"], cmap="RdBu_r", vmin=-4, vmax=4, lw=0)
ax.axvline(0, color="black", lw=0.5)
ax.set_yticks(range(len(order)))
ax.set_yticklabels(order)
ax.set_xlabel(f"log fold change ({CASE} vs {CTRL})")
ax.set_title("Milo neighbourhoods (coloured: SpatialFDR < 0.1)")
save_figure(fig, "Fig2_milo_beeswarm", FIG)
plt.show()

# %% [markdown]
# ## 6. Sensitivity: Vanderbilt-only donors
# Vanderbilt is the only site that contributed both controls and IPF; repeating
# the within-compartment test there removes the site-disease confounding.

# %%
vu = sub[sub.obs["site"] == "Vanderbilt"].copy()
print(vu.obs.groupby("diagnosis", observed=True)["donor_id"].nunique())
sens = []
for comp in ["Epithelial", "Immune", "Endothelial", "Mesenchymal"]:
    s = vu[vu.obs["compartment"] == comp].copy()
    s.obs["cell_type"] = s.obs["cell_type"].cat.remove_unused_categories()
    r = composition.run_propeller(s, R, P.r_scripts / "propeller.R", P.interim / "propeller", "cell_type",
                                  "donor_id", "diagnosis", CASE, CTRL, covariates=None, tag=f"_VU_{comp}")
    sens.append(r)
sens = pd.concat(sens).set_index("celltype")
cmp = pd.DataFrame({"all_sites": cons["log2_prop_ratio"], "vanderbilt_only": np.log2(sens["PropRatio"])}).dropna()
write_table(cmp, TAB / "composition_sensitivity_vanderbilt")
r = cmp.corr(method="spearman").iloc[0, 1]
fig, ax = plt.subplots(figsize=figsize(65, 60))
ax.scatter(cmp["all_sites"], cmp["vanderbilt_only"], s=8, c=[COMPARTMENT_COLORS[ct_comp[c]] for c in cmp.index])
lim = np.nanmax(np.abs(cmp.to_numpy())) * 1.05
ax.plot([-lim, lim], [-lim, lim], ls="--", color="grey", lw=0.5)
for c in cmp.reindex(plot_types).dropna().index[:8]:
    ax.annotate(c, (cmp.loc[c, "all_sites"], cmp.loc[c, "vanderbilt_only"]), fontsize=4.5)
ax.set_xlabel("log2 proportion ratio, all sites")
ax.set_ylabel("log2 proportion ratio, Vanderbilt only")
ax.set_title(f"Site sensitivity (Spearman ρ = {r:.2f})")
save_figure(fig, "S6_composition_site_sensitivity", FIG)
plt.show()
