"""Pathway and transcription-factor activity inference (decoupler 2.x).

Rather than testing over-representation of thresholded gene lists, footprint
methods use the *full* signed DE statistic of every gene:

* **PROGENy** (Schubert *et al.* 2018): 14 signalling pathways, each defined by
  genes consistently *responding* to pathway perturbation (a footprint, not the
  pathway members) - captures e.g. TGF-β, TNF-α, hypoxia, p53 activity.
* **CollecTRI** (Müller-Dott *et al.* 2023): signed TF-target regulons curated
  from 12 resources; TF activity is inferred from its targets' changes.
* **MSigDB Hallmarks** (Liberzon *et al.* 2015): 50 coherent biological programmes.

Activities are estimated with a univariate linear model (ULM): for each source,
the DE statistics of all genes are regressed on the source's target weights;
the t-value of the slope is the activity score (Badia-i-Mompel *et al.* 2022).
"""
from __future__ import annotations

import decoupler as dc
import numpy as np
import pandas as pd

_RES: dict = {}


def resource(name: str, **kw) -> pd.DataFrame:
    """Cached OmniPath resources (downloaded once per session)."""
    key = (name, tuple(sorted(kw.items())))
    if key not in _RES:
        if name == "progeny":
            _RES[key] = dc.op.progeny(organism="human", **kw)
        elif name == "collectri":
            _RES[key] = dc.op.collectri(organism="human", **kw)
        elif name == "hallmark":
            _RES[key] = dc.op.hallmark(organism="human", **kw)
        else:
            raise ValueError(name)
    return _RES[key].copy()


def activity_from_stats(stat_mat: pd.DataFrame, net: pd.DataFrame, tmin: int = 5) -> tuple[pd.DataFrame, pd.DataFrame]:
    """ULM activities for a contrasts x genes matrix of DE statistics."""
    scores, padj = dc.mt.ulm(stat_mat, net, tmin=tmin)
    return scores, padj


def ulm_cells(adata, net: pd.DataFrame, key: str, tmin: int = 5) -> pd.DataFrame:
    """Per-cell activities (stored in ``adata.obsm[f'score_{key}']``)."""
    dc.mt.ulm(adata, net, tmin=tmin)
    adata.obsm[f"score_{key}"] = adata.obsm.pop("score_ulm")
    adata.obsm[f"padj_{key}"] = adata.obsm.pop("padj_ulm")
    return adata.obsm[f"score_{key}"]


def top_sources(scores: pd.DataFrame, padj: pd.DataFrame, n: int = 15, alpha: float = 0.05) -> list[str]:
    """Union of the n most significant sources per contrast."""
    keep = set()
    for c in scores.index:
        sig = padj.loc[c][padj.loc[c] < alpha].index
        keep |= set(scores.loc[c, sig].abs().sort_values(ascending=False).head(n).index)
    return sorted(keep)


def fisher_overlap(a: set, b: set, universe: set) -> tuple[float, float]:
    """Odds ratio and one-sided P of the overlap between two gene sets."""
    from scipy.stats import fisher_exact

    a, b = a & universe, b & universe
    tab = np.array([[len(a & b), len(a - b)], [len(b - a), len(universe - a - b)]])
    odds, p = fisher_exact(tab, alternative="greater")
    return odds, p
