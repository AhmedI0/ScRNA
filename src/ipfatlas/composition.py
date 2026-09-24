"""Differential cell-type abundance between conditions.

Cell-type proportions are *compositional*: if one population expands, the
proportions of all others shrink even if their absolute numbers are unchanged.
Two complementary donor-level methods are used:

* **propeller** (speckle; Phipson *et al.* 2022): logit-transformed proportions,
  empirical-Bayes moderated t-tests (limma). Fast, frequentist, supports covariates.
* **scCODA** (Büttner *et al.* 2021): Bayesian Dirichlet-multinomial model with
  an explicit reference cell type, i.e. it models compositionality directly.

Only cell-type changes supported by both are called "robust".
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from anndata import AnnData

from .utils import log, run_r


def donor_counts(adata: AnnData, celltype_key: str, donor_key: str, condition_key: str,
                 extra: list[str] | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Cells per donor x cell type and the donor metadata table."""
    counts = pd.crosstab(adata.obs[donor_key], adata.obs[celltype_key])
    cols = [condition_key, *(extra or [])]
    meta = adata.obs.groupby(donor_key, observed=True)[cols].first().loc[counts.index]
    return counts, meta


def proportions_long(counts: pd.DataFrame, meta: pd.DataFrame, condition_key: str) -> pd.DataFrame:
    prop = counts.div(counts.sum(axis=1), axis=0)
    long = prop.reset_index().melt(id_vars=prop.index.name or "index", var_name="celltype", value_name="proportion")
    long = long.rename(columns={prop.index.name or "index": "donor"})
    long[condition_key] = long["donor"].map(meta[condition_key])
    return long


def run_propeller(adata: AnnData, rscript: str, r_script: Path, workdir: Path, celltype_key: str,
                  donor_key: str, condition_key: str, case: str, control: str,
                  covariates: list[str] | None = None, tag: str = "") -> pd.DataFrame:
    workdir.mkdir(parents=True, exist_ok=True)
    cols = [celltype_key, donor_key, condition_key, *(covariates or [])]
    cells = adata.obs[cols].astype(str).rename(columns={celltype_key: "celltype", donor_key: "donor_id",
                                                         condition_key: "condition"})
    inp = workdir / f"cells{tag}.csv"
    out = workdir / f"propeller{tag}.csv"
    cells.to_csv(inp, index=False)
    run_r(rscript, r_script, inp, out, case, control, ",".join(covariates) if covariates else "none",
          log_file=workdir / f"propeller{tag}.log")
    return pd.read_csv(out)


def run_sccoda(counts: pd.DataFrame, meta: pd.DataFrame, condition_key: str, case: str, control: str,
               reference: str | None = "automatic", fdr: float = 0.1, seed: int = 0,
               num_samples: int = 20000, num_warmup: int = 5000, max_tries: int = 3,
               accept_range: tuple[float, float] = (0.6, 0.95)) -> pd.DataFrame:
    """scCODA via pertpy (NUTS sampler). Returns per-cell-type effect table.

    A chain whose mean acceptance rate is outside ``accept_range`` (pertpy's own
    warning bounds) has not explored the posterior; it is restarted with the next
    seed, up to ``max_tries`` times, and a RuntimeError is raised if none converges.
    """
    import anndata as ad
    import numpyro
    import pertpy as pt

    sample = ad.AnnData(X=counts.values.astype(float), obs=meta.copy(),
                        var=pd.DataFrame(index=counts.columns.astype(str)))
    sample.obs[condition_key] = pd.Categorical(sample.obs[condition_key], categories=[control, case])
    model = pt.tl.Sccoda()
    rates = []
    for key in range(seed, seed + max_tries):
        mdata = model.load(sample.copy(), type="sample_level", generate_sample_level=False)
        mdata = model.prepare(mdata, modality_key="coda", formula=condition_key, reference_cell_type=reference)
        # scCODA replaces zero counts by 0.5, so per-sample totals are not integers. numpyro >= 0.22
        # validates distribution arguments by default and rejects that; the Dirichlet-multinomial
        # likelihood (gamma functions) is defined for real-valued totals, as in the original scCODA.
        # (set_fdr re-runs the model for the posterior summary, so it needs the same setting)
        with numpyro.validation_enabled(False):
            model.run_nuts(mdata, modality_key="coda", rng_key=key, num_samples=num_samples, num_warmup=num_warmup)
            rate = float(mdata["coda"].uns["scCODA_params"]["mcmc"]["acceptance_rate"])
            rates.append(rate)
            if not accept_range[0] <= rate <= accept_range[1]:
                log.warning(f"scCODA chain with rng_key={key} did not converge (acceptance rate {rate:.2f})")
                continue
            model.set_fdr(mdata, modality_key="coda", est_fdr=fdr)
            eff = model.get_effect_df(mdata, modality_key="coda").reset_index()
        break
    else:
        raise RuntimeError(f"no converged chain in {max_tries} tries (acceptance rates {[round(r, 2) for r in rates]})")
    eff = eff.rename(columns={"Cell Type": "celltype"})
    ref = mdata["coda"].uns["scCODA_params"]["reference_cell_type"]
    eff["reference_cell_type"] = ref
    eff["acceptance_rate"] = rate
    eff["rng_key"] = key
    log.info(f"scCODA reference cell type: {ref}; acceptance rate {rate:.2f} (rng_key={key}); "
             f"credible effects: {(eff['Final Parameter'] != 0).sum()}")
    return eff


def clr(counts: pd.DataFrame, pseudocount: float = 0.5) -> pd.DataFrame:
    """Centred log-ratio transform of a donors x cell types count table."""
    x = np.log(counts + pseudocount)
    return x.sub(x.mean(axis=1), axis=0)
