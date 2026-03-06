"""Donor-level pseudobulk aggregation and orchestration of the R DE models."""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData

from .utils import log, run_r


def safe_key(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", str(name)).strip("_")


def aggregate(
    adata: AnnData,
    celltype_key: str,
    donor_key: str,
    condition_key: str,
    covariates: list[str] | None = None,
    layer: str = "counts",
    min_cells: int = 10,
    celltypes: list[str] | None = None,
) -> dict[str, tuple[pd.DataFrame, pd.DataFrame]]:
    """Sum raw counts per donor within each cell type.

    Returns ``{celltype: (counts genes x samples, sample metadata)}``. Donors
    with fewer than ``min_cells`` cells of a type are dropped for that type,
    because sums over a handful of cells are dominated by sampling noise.
    """
    covariates = covariates or []
    x = adata.layers[layer] if layer else adata.X
    x = sp.csr_matrix(x)
    obs = adata.obs
    out = {}
    cts = celltypes or list(obs[celltype_key].cat.categories)
    for ct in cts:
        mask = (obs[celltype_key] == ct).to_numpy()
        if mask.sum() == 0:
            continue
        sub = obs.loc[mask]
        donors = sub[donor_key].astype(str)
        ncell = donors.value_counts()
        keep_donors = ncell[ncell >= min_cells].index
        if len(keep_donors) == 0:
            continue
        d_codes = pd.Categorical(donors, categories=keep_donors)
        valid = d_codes.codes >= 0
        rows = np.where(mask)[0][valid]
        ind = sp.csr_matrix((np.ones(valid.sum()), (d_codes.codes[valid], np.arange(valid.sum()))),
                            shape=(len(keep_donors), valid.sum()))
        summed = ind @ x[rows]
        counts = pd.DataFrame(np.asarray(summed.todense()).T.round().astype(np.int64),
                              index=adata.var_names, columns=list(keep_donors))
        meta = (sub.loc[sub[donor_key].astype(str).isin(keep_donors), [donor_key, condition_key, *covariates]]
                .astype(str).drop_duplicates(donor_key).set_index(donor_key).loc[list(keep_donors)])
        meta = meta.rename(columns={condition_key: "condition"})
        meta["n_cells"] = ncell.loc[list(keep_donors)].values
        out[ct] = (counts, meta)
    return out


def write_inputs(pb: dict, outdir: Path, case: str, control: str, min_donors: int = 3) -> pd.DataFrame:
    """Write pseudobulk matrices for the R script; skip under-powered cell types."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rows = []
    for ct, (counts, meta) in pb.items():
        n_case = (meta["condition"] == case).sum()
        n_ctrl = (meta["condition"] == control).sum()
        status = "tested" if (n_case >= min_donors and n_ctrl >= min_donors) else "skipped (too few donors)"
        rows.append({"celltype": ct, "key": safe_key(ct), "n_case": n_case, "n_control": n_ctrl,
                     "median_cells_per_donor": int(meta["n_cells"].median()), "status": status})
        if status != "tested":
            continue
        counts = counts.loc[counts.sum(axis=1) > 0]
        counts.to_csv(outdir / f"{safe_key(ct)}__counts.csv.gz", compression="gzip")
        meta.to_csv(outdir / f"{safe_key(ct)}__meta.csv")
    idx = pd.DataFrame(rows)
    idx[idx.status == "tested"][["key", "celltype"]].to_csv(outdir / "index.csv", index=False)
    return idx


def run_de(rscript: str, r_script: Path, indir: Path, outdir: Path, covariates: list[str], case: str,
           control: str, methods: str = "edgeR,DESeq2") -> pd.DataFrame:
    cov = ",".join(covariates) if covariates else "none"
    run_r(rscript, r_script, indir, outdir, cov, case, control, methods, log_file=Path(outdir).with_suffix(".log"))
    return pd.read_csv(Path(outdir) / "summary.csv")


def read_results(outdir: Path, method: str = "edgeR") -> pd.DataFrame:
    """Concatenate per-cell-type result files into one long table."""
    outdir = Path(outdir)
    frames = []
    for f in sorted(outdir.glob(f"*__{method}.csv")):
        key = f.name.split("__")[0]
        df = pd.read_csv(f)
        df.insert(0, "key", key)
        frames.append(df)
    res = pd.concat(frames, ignore_index=True)
    summary = pd.read_csv(outdir / "summary.csv")[["key", "celltype"]]
    res = res.merge(summary, on="key", how="left")
    log.info(f"{method}: {res.celltype.nunique()} cell types, {len(res):,} gene tests")
    return res


def stat_matrix(res: pd.DataFrame, stat: str = "stat", min_celltypes: int = 1) -> pd.DataFrame:
    """Wide genes x cell types matrix of a DE statistic (input for decoupler)."""
    mat = res.pivot_table(index="celltype", columns="gene", values=stat)
    mat = mat.loc[:, mat.notna().sum(axis=0) >= min_celltypes]
    return mat.fillna(0.0)
