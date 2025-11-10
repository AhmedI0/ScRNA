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
# # 01 · Data acquisition, cohort design and metadata harmonisation
#
# **Question addressed by the project.** Idiopathic pulmonary fibrosis (IPF) is a
# progressive, fatal interstitial lung disease in which repetitive injury of the
# alveolar epithelium is followed by aberrant repair, fibroblast activation and
# extracellular-matrix (ECM) deposition. Single-cell RNA-seq (scRNA-seq) studies
# in 2019-2020 described disease-emergent cell states - most notably
# **KRT5⁻/KRT17⁺ "aberrant basaloid" epithelial cells**, **KRT8^hi transitional
# alveolar cells**, **CTHRC1⁺/POSTN⁺ myofibroblasts** and **SPP1^hi macrophages**.
#
# This project re-analyses two independent public cohorts from raw count matrices
# with a modern, statistically rigorous workflow and asks:
#
# 1. Which cell types and states change in abundance in IPF (donor-level statistics)?
# 2. Which cell-type-specific transcriptional programmes change (pseudobulk DE)?
# 3. What transcriptional trajectory links AT2 cells to aberrant basaloid cells?
# 4. Which ligand-receptor circuits connect aberrant epithelium, fibroblasts and macrophages?
# 5. Do these findings **replicate in an independent cohort**?
#
# | role | GEO | study | technology | donors |
# |---|---|---|---|---|
# | discovery | GSE135893 | Habermann *et al.* 2020, *Sci Adv* | 10x 5′ v1 | 10 control, 12 IPF, 10 other ILD |
# | validation | GSE136831 | Adams *et al.* 2020, *Sci Adv* | 10x 3′ v2 | 28 control, 32 IPF, 18 COPD |
#
# Files were downloaded with `workflow/scripts/download_data.py`, which verifies the
# remote file size and records a SHA-256 checksum for each file.

# %%
import matplotlib.pyplot as plt
import pandas as pd

from ipfatlas import io, qc
from ipfatlas.config import setup
from ipfatlas.plotting import CONDITION_COLORS, SITE_COLORS, figsize, panel_label, save_figure, set_style
from ipfatlas.utils import package_versions, write_table

cfg, P = setup()
set_style(cfg["figures"]["font_size"])
FIG = P.figures / "01_data"
TAB = P.tables / "01_data"
package_versions()

# %% [markdown]
# ## 1. Provenance of the raw files

# %%
manifest = pd.read_csv(P.raw / "download_manifest.tsv", sep="\t")
manifest["MB"] = (manifest["bytes"] / 1e6).round(1)
manifest[["cohort", "accession", "file", "MB", "sha256", "retrieved"]]

# %% [markdown]
# ## 2. Discovery cohort (GSE135893)
#
# The GEO matrix contains **all Cell Ranger-called barcodes** (220,213) of 43
# libraries, i.e. *before* the authors' QC (≥1,000 genes, <25 % mitochondrial
# reads), which retained 114,396 cells. Starting from this matrix lets us perform
# QC ourselves rather than inherit the authors' choices. The authors' labels are
# kept only for benchmarking.

# %%
adata = io.read_habermann(P.raw)
adata

# %%
lib = (adata.obs.groupby("library_id", observed=True)
       .agg(donor_id=("donor_id", "first"), diagnosis=("diagnosis", "first"), status=("status", "first"),
            site=("site", "first"), n_barcodes=("library_id", "size"),
            n_author_cells=("author_qc_pass", "sum"))
       .sort_values(["status", "diagnosis", "donor_id"]))
write_table(lib, TAB / "discovery_libraries")
lib

# %%
donors = lib.groupby("donor_id").agg(diagnosis=("diagnosis", "first"), site=("site", "first"),
                                     n_libraries=("n_barcodes", "size"), n_barcodes=("n_barcodes", "sum"))
pd.crosstab(donors["diagnosis"], donors["site"], margins=True)

# %% [markdown]
# **Design caveat — site is partly confounded with disease.** All donors from the
# Donor Network of Arizona (DNA) are controls and all from NTI are fibrotic; only
# Vanderbilt contributed both. A naive IPF-vs-control contrast therefore mixes
# disease and site effects. We address this later with (i) a site-adjusted model,
# (ii) a Vanderbilt-only sensitivity analysis and, most importantly, (iii)
# replication in an independent cohort from a different centre.
#
# Several donors contributed two libraries (e.g. upper and lower lobe, `VUILD59-1/-2`).
# These are **technical/regional replicates of one biological replicate**: all
# statistics below use the *donor* as the unit of replication.

# %% [markdown]
# ## 3. Genetic sex inference
# Sex is not recorded in GEO, but sex-linked expression is an important covariate.
# Females express *XIST*; males express broadly transcribed chrY genes (*RPS4Y1*,
# *DDX3Y*, *UTY*, *KDM5D*, *EIF1AY*...). Both are measured from the raw counts.

# %%
sex = qc.infer_sex(adata, donor_key="donor_id", layer=None)
sex = sex.join(donors[["diagnosis", "site"]])
write_table(sex, TAB / "discovery_inferred_sex")
sex.sort_values("diagnosis")

# %%
fig, ax = plt.subplots(figsize=figsize(60, 55))
for dx, d in sex.groupby("diagnosis"):
    ax.scatter(d["XIST_per_10k"], d["chrY_per_10k"], s=14, color=CONDITION_COLORS.get(dx, "grey"),
               edgecolor="black", lw=0.3, label=dx)
ax.axline((0, 0), slope=1, color="grey", lw=0.5, ls="--")
ax.set_xlabel("XIST (counts per 10k)")
ax.set_ylabel("chrY genes (counts per 10k)")
ax.set_title("Donor genetic sex (discovery)")
ax.legend(loc="upper right", fontsize=5, handletextpad=0.1)
save_figure(fig, "S1a_inferred_sex_discovery", FIG)
plt.show()

# %%
adata.obs["sex"] = adata.obs["donor_id"].map(sex["inferred_sex"]).astype("category")
pd.crosstab(sex["diagnosis"], sex["inferred_sex"])

# %% [markdown]
# ## 4. Validation cohort (GSE136831) — metadata only
# The 2 GB validation matrix is loaded in notebook 11; here we only describe its design.

# %%
vmeta = pd.read_csv(P.raw / "GSE136831" / "GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz", sep="\t")
vdon = vmeta.groupby("Subject_Identity").agg(diagnosis=("Disease_Identity", "first"),
                                              n_libraries=("Library_Identity", "nunique"),
                                              n_cells=("CellBarcode_Identity", "size"))
write_table(vdon, TAB / "validation_donors")
vdon["diagnosis"].value_counts()

# %% [markdown]
# ## 5. Figure: study design overview

# %%
fig, axes = plt.subplots(1, 3, figsize=figsize(183, 55), gridspec_kw={"width_ratios": [1, 1, 2.2]})

ax = axes[0]
ct = pd.crosstab(donors["diagnosis"], donors["site"]).loc[["Control", "IPF", "NSIP", "cHP", "Sarcoidosis", "ILD"]]
ct.plot.barh(stacked=True, ax=ax, color=[SITE_COLORS[c] for c in ct.columns], width=0.8)
ax.invert_yaxis()
ax.set_xlabel("donors")
ax.set_ylabel("")
ax.set_title("Discovery (GSE135893)")
ax.legend(title="site", fontsize=5, title_fontsize=6)
panel_label(ax, "a", x=-0.35)

ax = axes[1]
vc = vdon["diagnosis"].value_counts().loc[["Control", "IPF", "COPD"]]
ax.barh(vc.index, vc.values, color=[CONDITION_COLORS[c] for c in vc.index])
ax.invert_yaxis()
ax.set_xlabel("donors")
ax.set_title("Validation (GSE136831)")
panel_label(ax, "b", x=-0.3)

ax = axes[2]
order = lib.sort_values(["diagnosis", "n_barcodes"]).index
colors = [CONDITION_COLORS.get(d, "grey") for d in lib.loc[order, "diagnosis"]]
ax.bar(range(len(order)), lib.loc[order, "n_barcodes"], color=colors, width=0.8, alpha=0.35, label="Cell Ranger barcodes")
ax.bar(range(len(order)), lib.loc[order, "n_author_cells"], color=colors, width=0.8, label="cells passing author QC")
ax.set_xticks(range(len(order)))
ax.set_xticklabels(order, rotation=90, fontsize=4)
ax.set_ylabel("cells per library")
ax.set_title("Discovery libraries: barcodes before vs after author QC")
handles = [plt.Rectangle((0, 0), 1, 1, color=CONDITION_COLORS[d]) for d in ["Control", "IPF", "NSIP", "cHP", "Sarcoidosis", "ILD"]]
ax.legend(handles, ["Control", "IPF", "NSIP", "cHP", "Sarcoidosis", "ILD"], ncol=3, fontsize=5)
panel_label(ax, "c", x=-0.08)
fig.tight_layout()
save_figure(fig, "Fig1_study_design_cohorts", FIG)
plt.show()

# %% [markdown]
# ## 6. Save
# Raw UMI counts are in `.X`. From notebook 02 on, they are kept unchanged in
# `layers["counts"]` and every transformation writes to `.X` or new layers.

# %%
out = P.interim / "discovery_raw.h5ad"
adata.write_h5ad(out, compression="lzf")
print(out, f"{out.stat().st_size / 1e9:.2f} GB")
