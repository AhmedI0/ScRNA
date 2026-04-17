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
# # 09 · Alveolar epithelial trajectories: AT2 → transitional → AT1 or aberrant basaloid
#
# **Biological question.** After injury, AT2 cells act as alveolar stem cells:
# they proliferate and differentiate into AT1 cells through a KRT8-high
# *transitional* state (mouse: Strunz *et al.* 2020; Kobayashi *et al.* 2020;
# Choi *et al.* 2020). In IPF lungs, Habermann *et al.* and Adams *et al.* (2020)
# described **aberrant basaloid (KRT5⁻/KRT17⁺) cells** that line fibroblastic foci
# and appear to be a *pathological end point* of failed AT2-to-AT1
# differentiation. We test whether the data support this lineage structure and
# which genes track the commitment towards the aberrant state.
#
# **Methods and caveats**
# * **PAGA** (Wolf *et al.* 2019): graph abstraction showing which cell types are
#   connected at all in the kNN graph.
# * **Diffusion pseudotime** (DPT; Haghverdi *et al.* 2016) from a mature AT2 root.
# * **CellRank 2** (Weiler *et al.* 2024) with the **CytoTRACE kernel** - direction
#   is inferred from *transcriptional diversity* (number of expressed genes,
#   Gulati *et al.* 2020), which is high in less-differentiated cells. The 5′ data
#   have no reliable spliced/unspliced information for RNA velocity
#   (5′ chemistry captures transcript starts, not intronic reads near 3′ ends).
# * **Snapshot caveat.** scRNA-seq of human tissue is a snapshot: a trajectory is a
#   *hypothesis about ordering* consistent with transcriptional similarity, not a
#   demonstration of lineage (which needs lineage tracing, possible only in mice).
#
# GPCCA (pyGPCCA) needs a Schur decomposition of the transition matrix. Without
# PETSc/SLEPc (not available on Windows) the dense method is used, so cells are
# subsampled to ≤1,800 per cell type, which keeps rare states fully represented.

# %%
import anndata as ad
import cellrank as cr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv

from ipfatlas import annotation, pseudobulk
from ipfatlas.config import rscript, setup
from ipfatlas.plotting import (
    CONDITION_COLORS,
    embedding_scatter,
    figsize,
    label_centroids,
    panel_label,
    save_figure,
    set_style,
)
from ipfatlas.utils import log, write_table

cfg, P = setup()
set_style(cfg["figures"]["font_size"])
FIG = P.figures / "09_trajectory"
TAB = P.tables / "09_trajectory"
LINEAGE = cfg["trajectory"]["lineage"]
CAP = 1800

full = sc.read_h5ad(P.processed / "discovery_annotated.h5ad")
colors = dict(zip(full.obs["cell_type"].cat.categories, full.uns["cell_type_colors"]))
pd.crosstab(full.obs.loc[full.obs.cell_type.isin(LINEAGE), "cell_type"].astype(str),
            full.obs.loc[full.obs.cell_type.isin(LINEAGE), "diagnosis"])

# %% [markdown]
# ## 1. Lineage subset (stratified subsample)

# %%
rng = np.random.default_rng(0)
idx = []
for ct in LINEAGE:
    cells = np.where(full.obs["cell_type"] == ct)[0]
    idx.extend(rng.choice(cells, min(CAP, len(cells)), replace=False))
adata = full[np.sort(idx)].copy()
adata.obs["cell_type"] = adata.obs["cell_type"].cat.remove_unused_categories()
adata.uns["cell_type_colors"] = [colors[c] for c in adata.obs["cell_type"].cat.categories]
print(adata.obs["cell_type"].value_counts())

sc.pp.neighbors(adata, use_rep="X_integrated", n_neighbors=30, random_state=0)
sc.tl.paga(adata, groups="cell_type")
sc.pl.paga(adata, show=False, threshold=0.1)
plt.close()
sc.tl.umap(adata, init_pos="paga", min_dist=0.3, random_state=0)
sc.tl.diffmap(adata, n_comps=15)

# %% [markdown]
# ## 2. PAGA connectivity
# Edge weight = how much more connected two cell types are in the kNN graph than
# expected by chance. A path AT2 → transitional AT2 → AT1 and a branch from the
# transitional state to aberrant basaloid cells is the pattern predicted by the
# mouse injury models.

# %%
con = pd.DataFrame(adata.uns["paga"]["connectivities"].toarray(), index=adata.obs["cell_type"].cat.categories,
                   columns=adata.obs["cell_type"].cat.categories)
write_table(con.round(3), TAB / "paga_connectivities")

fig, axes = plt.subplots(1, 2, figsize=figsize(150, 65))
embedding_scatter(axes[0], adata.obsm["X_umap"], adata.obs["cell_type"], palette=colors, s=1.5)
label_centroids(axes[0], adata.obsm["X_umap"], adata.obs["cell_type"], fontsize=5)
axes[0].set_title("alveolar epithelial lineage")
sc.pl.paga(adata, ax=axes[1], show=False, threshold=0.1, fontsize=5, node_size_scale=0.6, edge_width_scale=0.6,
           frameon=False, color="cell_type")
axes[1].set_title("PAGA (edges ≥ 0.1)")
for a, l in zip(axes, "ab"):
    panel_label(a, l, x=0.0)
fig.tight_layout()
save_figure(fig, "Fig4_lineage_umap_paga", FIG)
plt.show()
con.round(2)

# %% [markdown]
# ## 3. Diffusion pseudotime
# Root: the AT2 cell with the highest AT2 marker score (most mature AT2 cell).

# %%
markers = annotation.load_markers(P.root / "config" / "markers.yaml")
annotation.score_marker_sets(adata, {"AT2": markers["epithelial"]["AT2"]}, prefix="mk_")
at2 = np.where(adata.obs["cell_type"] == cfg["trajectory"]["root_cell_type"])[0]
adata.uns["iroot"] = int(at2[np.argmax(adata.obs["mk_AT2"].to_numpy()[at2])])
sc.tl.dpt(adata)

# %% [markdown]
# ## 4. CellRank 2: CytoTRACE kernel + GPCCA
# The CytoTRACE kernel needs smoothed expression ("Ms"), computed with scVelo's
# kNN moments on log-normalised expression (spliced = unspliced = X is the
# documented way to use moments without velocity).

# %%
adata.layers["spliced"] = adata.X
adata.layers["unspliced"] = adata.X
scv.pp.moments(adata, n_pcs=None, n_neighbors=30, use_rep="X_integrated")
ctk = cr.kernels.CytoTRACEKernel(adata).compute_cytotrace()
ctk.compute_transition_matrix(threshold_scheme="soft", nu=0.5)
ck = cr.kernels.ConnectivityKernel(adata).compute_transition_matrix()
kernel = 0.8 * ctk + 0.2 * ck  # small diffusion term for robustness
del adata.layers["spliced"], adata.layers["unspliced"]

pot = adata.obs.groupby("cell_type", observed=True)["ct_pseudotime"].agg(["median", "mean"]).sort_values("median")
write_table(pot.round(3), TAB / "cytotrace_pseudotime_per_celltype")
pot

# %%
g = cr.estimators.GPCCA(kernel)
g.compute_schur(n_components=20, method="brandts")
g.compute_macrostates(n_states=6, cluster_key="cell_type")
ms_comp = pd.crosstab(g.macrostates.astype(str), adata.obs["cell_type"].astype(str))
write_table(ms_comp, TAB / "gpcca_macrostates_composition")
ms_comp

# %% [markdown]
# Terminal states are set to the hypothesised end points (AT1, aberrant
# basaloid, basal). Where GPCCA found a macrostate dominated by that cell type,
# its cells are used; otherwise the 30 most differentiated cells (lowest
# CytoTRACE potency) of that type define the state. The macrostate table above
# shows how well the unsupervised states match.

# %%
TERMINAL = [t for t in ["AT1", "Aberrant basaloid", "Basal"] if t in adata.obs["cell_type"].cat.categories]
terminal = {}
for t in TERMINAL:
    hits = [m for m in g.macrostates.cat.categories if m == t or m.startswith(f"{t}_")]
    if hits:
        terminal[t] = adata.obs_names[g.macrostates == hits[0]].tolist()
        src = f"macrostate {hits[0]}"
    else:
        cells = adata.obs.index[adata.obs["cell_type"] == t]
        terminal[t] = adata.obs.loc[cells, "ct_pseudotime"].nlargest(30).index.tolist()
        src = "30 least-potent cells"
    log.info(f"terminal state {t}: {len(terminal[t])} cells ({src})")
g.set_terminal_states(states=terminal)
g.compute_fate_probabilities(use_petsc=False, show_progress_bar=False)
fates = pd.DataFrame(np.asarray(g.fate_probabilities), index=adata.obs_names, columns=g.fate_probabilities.names)
for c in fates:
    adata.obs[f"fate_{c}"] = fates[c]
fate_ct = fates.groupby(adata.obs["cell_type"].values).mean()
write_table(fate_ct.round(3), TAB / "fate_probabilities_mean_per_celltype")
fate_ct.round(2)

# %%
fig, axes = plt.subplots(1, 4, figsize=figsize(183, 50))
xy = adata.obsm["X_umap"]
for ax, (col, title, cmap) in zip(axes, [("ct_pseudotime", "CytoTRACE pseudotime\n(0 = high potency)", "viridis"),
                                         ("dpt_pseudotime", "diffusion pseudotime\n(root: mature AT2)", "viridis"),
                                         ("fate_AT1", "fate probability: AT1", "magma_r"),
                                         ("fate_Aberrant basaloid", "fate probability: aberrant basaloid",
                                          "magma_r")]):
    sca = embedding_scatter(ax, xy, adata.obs[col], cmap=cmap, s=1.2, order="value")
    plt.colorbar(sca, ax=ax, shrink=0.55)
    ax.set_title(title)
for a, l in zip(axes, "abcd"):
    panel_label(a, l, x=0.0)
fig.tight_layout()
save_figure(fig, "Fig4_pseudotime_fates", FIG)
plt.show()

# %% [markdown]
# ### Do transitional cells from IPF lungs lean towards the aberrant fate?
# Compared per **donor** (mean over cells), not per cell.

# %%
tr = adata.obs[adata.obs["cell_type"] == "Transitional AT2"]
don = tr.groupby("donor_id", observed=True).agg(diagnosis=("diagnosis", "first"), n=("diagnosis", "size"),
                                                fate_basaloid=("fate_Aberrant basaloid", "mean"),
                                                fate_AT1=("fate_AT1", "mean"))
don = don[(don["n"] >= 10) & don["diagnosis"].isin(["Control", "IPF"])]
write_table(don.round(3), TAB / "transitional_fate_per_donor")
from scipy.stats import mannwhitneyu  # noqa: E402

p = mannwhitneyu(don.loc[don.diagnosis == "IPF", "fate_basaloid"], don.loc[don.diagnosis == "Control", "fate_basaloid"]).pvalue
fig, ax = plt.subplots(figsize=figsize(45, 50))
for i, dx in enumerate(["Control", "IPF"]):
    v = don.loc[don.diagnosis == dx, "fate_basaloid"]
    ax.boxplot(v, positions=[i], widths=0.5, showfliers=False, medianprops={"color": "black"})
    ax.scatter(np.full(len(v), i) + rng.uniform(-0.12, 0.12, len(v)), v, s=8, color=CONDITION_COLORS[dx], zorder=3)
ax.set_xticks([0, 1], ["Control", "IPF"])
ax.set_ylabel("mean fate probability\n→ aberrant basaloid")
ax.set_title(f"Transitional AT2 cells\n(per donor, Mann-Whitney P = {p:.2g})")
save_figure(fig, "Fig4_transitional_fate_by_donor", FIG)
plt.show()

# %% [markdown]
# ## 5. Genes tracking commitment to the aberrant basaloid fate
# CellRank lineage drivers: correlation of each gene's expression with the fate
# probability, among AT2, transitional and aberrant basaloid cells.

# %%
drivers = g.compute_lineage_drivers(lineages=["Aberrant basaloid"], cluster_key="cell_type",
                                    clusters=["AT2", "Transitional AT2", "Aberrant basaloid"])
write_table(drivers.head(300), TAB / "lineage_drivers_aberrant_basaloid")
top_drivers = drivers.index[:30].tolist()
print(", ".join(top_drivers))

# %%
# expression trends along the AT2 -> transitional -> aberrant basaloid path (cells ordered by
# fate-weighted DPT, binned)
path = adata[adata.obs["cell_type"].isin(["AT2", "Transitional AT2", "Aberrant basaloid"])
             & (adata.obs["fate_Aberrant basaloid"] >= adata.obs["fate_AT1"])].copy()
order = path.obs["dpt_pseudotime"].rank(method="first")
bins = pd.qcut(order, 40, labels=False)
expr = pd.DataFrame(path[:, top_drivers].X.toarray(), index=path.obs_names, columns=top_drivers)
trend = expr.groupby(bins.values).mean()
trend = (trend - trend.mean()) / trend.std()
comp = pd.crosstab(bins.values, path.obs["cell_type"].astype(str), normalize="index")

fig = plt.figure(figsize=figsize(100, 110))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 14], width_ratios=[30, 1], hspace=0.05, wspace=0.05)
axb = fig.add_subplot(gs[0, 0])
bottom = np.zeros(len(comp))
for ct in comp.columns:
    axb.bar(np.arange(len(comp)), comp[ct], bottom=bottom, color=colors[ct], width=1.0, label=ct)
    bottom += comp[ct].to_numpy()
axb.set_xlim(-0.5, len(comp) - 0.5)
axb.axis("off")
axb.legend(ncol=3, fontsize=5, loc="lower center", bbox_to_anchor=(0.5, 1.0), frameon=False)
axh = fig.add_subplot(gs[1, 0])
peak = trend.idxmax().sort_values().index
im = axh.imshow(trend[peak].T, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5, interpolation="nearest")
axh.set_yticks(range(len(peak)), peak, fontsize=5)
axh.set_xticks([0, len(trend) - 1], ["AT2", "aberrant basaloid"])
axh.set_xlabel("pseudotime bins (DPT)")
cax = fig.add_subplot(gs[1, 1])
fig.colorbar(im, cax=cax, label="z-score")
save_figure(fig, "Fig4_driver_gene_trends", FIG)
plt.show()

# %% [markdown]
# ## 6. Aberrant basaloid signature by **paired pseudobulk** (within IPF donors)
# A signature for the validation cohort is defined by comparing aberrant
# basaloid cells with the other alveolar/basal epithelial cells *of the same
# donor* (edgeR, design `~ donor + group`). Pairing removes donor effects; using
# the donor-level sums avoids the pseudoreplication of a per-cell Wilcoxon test.

# %%
ref_types = ["AT2", "Transitional AT2", "AT1", "Basal"]
mask = (full.obs["cell_type"].isin(["Aberrant basaloid", *ref_types]) & (full.obs["diagnosis"] == "IPF")).to_numpy()
epi = ad.AnnData(X=full.layers["counts"][mask], obs=full.obs.loc[mask, ["donor_id", "cell_type"]].copy(),
                 var=full.var[[]].copy())
epi.obs["group"] = np.where(epi.obs["cell_type"] == "Aberrant basaloid", "Basaloid", "OtherAlveolar")
epi.obs["pb_sample"] = epi.obs["donor_id"].astype(str) + "__" + epi.obs["group"]
epi.obs["signature"] = pd.Categorical(["aberrant_basaloid"] * epi.n_obs)
pb = pseudobulk.aggregate(epi, "signature", "pb_sample", "group", covariates=["donor_id"], layer=None, min_cells=10)
counts, meta = pb["aberrant_basaloid"]
both = meta.groupby("donor_id")["condition"].nunique()
keep = meta.index[meta["donor_id"].isin(both[both == 2].index)]
print(f"{len(keep) // 2} IPF donors with >= 10 aberrant basaloid and >= 10 reference cells")
sig_in, sig_out = P.processed / "pseudobulk" / "signature_in", P.processed / "pseudobulk" / "signature_out"
pseudobulk.write_inputs({"aberrant_basaloid": (counts[keep], meta.loc[keep])}, sig_in, case="Basaloid",
                        control="OtherAlveolar", min_donors=3)
pseudobulk.run_de(rscript(cfg), P.r_scripts / "pseudobulk_de.R", sig_in, sig_out, ["donor_id"], "Basaloid",
                  "OtherAlveolar", methods="edgeR")
sig_res = pseudobulk.read_results(sig_out, "edgeR")
signature = (sig_res[(sig_res.FDR < 0.01) & (sig_res.logFC > 2) & (sig_res.logCPM > 3)]
             .sort_values("logFC", ascending=False).head(50))
write_table(sig_res, TAB / "aberrant_basaloid_vs_alveolar_paired_edgeR", index=False)
write_table(signature[["gene", "logFC", "logCPM", "FDR"]], TAB / "aberrant_basaloid_signature", index=False)
print(len(signature), "signature genes:", ", ".join(signature.gene))

# %%
fig, ax = plt.subplots(figsize=figsize(80, 70))
sig_res["mlog10"] = -np.log10(sig_res["FDR"].clip(lower=1e-300))
ax.scatter(sig_res.logFC, sig_res.mlog10, s=1, color="#BBBBBB", rasterized=True)
s = sig_res.gene.isin(signature.gene)
ax.scatter(sig_res.loc[s, "logFC"], sig_res.loc[s, "mlog10"], s=3, color=colors["Aberrant basaloid"])
for _, r in sig_res[s].nlargest(12, "mlog10").iterrows():
    ax.annotate(r.gene, (r.logFC, r.mlog10), fontsize=5)
ax.set_xlabel("log2 fold change (aberrant basaloid / other alveolar epithelium)")
ax.set_ylabel("-log10 FDR")
ax.set_title("Paired pseudobulk, IPF donors")
save_figure(fig, "S9_basaloid_signature_volcano", FIG)
plt.show()

# %% [markdown]
# ## 7. Export input for pySCENIC (Linux workstation)
# Gene-regulatory network inference (GRNBoost2 + cisTarget motif pruning +
# AUCell; Aibar *et al.* 2017; Van de Sande *et al.* 2020) is run on the Linux
# machine (`linux_upstream/scripts/08_pyscenic.sh`). We export raw counts of the
# subsampled lineage cells, with genes detected in ≥1 % of them.

# %%
import loompy  # noqa: E402

det = np.asarray((adata.layers["counts"] > 0).mean(axis=0)).ravel()
genes = adata.var_names[det >= 0.01]
loom_path = P.root / "linux_upstream" / "input" / "scenic_input.loom"
loom_path.parent.mkdir(parents=True, exist_ok=True)
loom_path.unlink(missing_ok=True)
mat = adata[:, genes].layers["counts"].T.tocsr().astype(np.float32)
loompy.create(str(loom_path), mat.toarray(),
              row_attrs={"Gene": np.array(genes)},
              col_attrs={"CellID": np.array(adata.obs_names), "cell_type": adata.obs["cell_type"].astype(str).values,
                         "diagnosis": adata.obs["diagnosis"].astype(str).values,
                         "donor_id": adata.obs["donor_id"].astype(str).values})
print(f"wrote {loom_path} ({len(genes):,} genes x {adata.n_obs:,} cells)")

# %%
adata.obsm["fate_probabilities"] = fates.to_numpy()
adata.uns["fate_names"] = list(fates.columns)
for k in ["Ms", "Mu"]:
    adata.layers.pop(k, None)
# CellRank stores macrostate memberships as its own Lineage array type, which h5ad cannot write
for k in list(adata.obsm):
    if type(adata.obsm[k]).__module__.startswith("cellrank"):
        adata.obsm[k] = np.asarray(adata.obsm[k])
adata.write_h5ad(P.processed / "epithelial_trajectory.h5ad", compression="lzf")
