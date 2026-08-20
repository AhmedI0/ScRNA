#!/usr/bin/env python
"""Combine pySCENIC runs into robust consensus regulons.

A regulon (TF + motif-supported targets) is kept if the TF is recovered in at
least ``--min-fraction`` of runs; its consensus targets are the genes present in
that same fraction of runs. Per-cell AUCell activities of the robust regulons
are averaged over runs.

Outputs (in <outdir>):
  regulon_recovery.tsv     TF, number of runs, n consensus targets
  consensus_regulons.tsv   TF, target gene, number of runs supporting the edge
  auc_consensus.tsv.gz     cells x regulons mean AUCell matrix
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from pathlib import Path

import loompy
import numpy as np
import pandas as pd
from pyscenic.utils import load_motifs


def regulon_targets(regulons_csv: Path) -> dict[str, set[str]]:
    """Union of TargetGenes over all enriched motifs of each TF (as pySCENIC does)."""
    df = load_motifs(str(regulons_csv))
    out = defaultdict(set)
    for (tf, _motif), row in df.iterrows():
        for gene, _weight in row[("Enrichment", "TargetGenes")]:
            out[tf].add(gene)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("runs_dir", type=Path)
    ap.add_argument("outdir", type=Path)
    ap.add_argument("--min-fraction", type=float, default=0.66)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    runs = sorted(p for p in args.runs_dir.glob("run*") if (p / "regulons.csv").exists())
    n = len(runs)
    need = int(np.ceil(args.min_fraction * n))
    tf_count, edge_count = Counter(), Counter()
    for r in runs:
        reg = regulon_targets(r / "regulons.csv")
        tf_count.update(reg.keys())
        for tf, genes in reg.items():
            edge_count.update((tf, g) for g in genes)

    robust = sorted(tf for tf, c in tf_count.items() if c >= need)
    edges = [(tf, g, c) for (tf, g), c in edge_count.items() if tf in robust and c >= need]
    edges = pd.DataFrame(edges, columns=["TF", "target", "n_runs"])
    edges.sort_values(["TF", "n_runs"], ascending=[True, False]).to_csv(
        args.outdir / "consensus_regulons.tsv", sep="\t", index=False)
    rec = pd.DataFrame({"TF": list(tf_count), "n_runs": list(tf_count.values())})
    rec["robust"] = rec["n_runs"] >= need
    rec["n_consensus_targets"] = rec["TF"].map(edges.groupby("TF").size()).fillna(0).astype(int)
    rec.sort_values(["robust", "n_runs"], ascending=False).to_csv(args.outdir / "regulon_recovery.tsv", sep="\t",
                                                                  index=False)

    aucs = []
    for r in runs:
        with loompy.connect(str(r / "auc.loom"), mode="r", validate=False) as ds:
            auc = pd.DataFrame(ds.ca["RegulonsAUC"], index=ds.ca["CellID"])
            auc.columns = [c.replace("(+)", "") for c in auc.columns]
        aucs.append(auc[[c for c in auc.columns if c in robust]])
    mean_auc = pd.concat(aucs).groupby(level=0).mean()
    mean_auc.to_csv(args.outdir / "auc_consensus.tsv.gz", sep="\t", compression="gzip")
    print(f"{n} runs; {len(robust)} robust regulons (>= {need} runs); {len(edges)} consensus edges")


if __name__ == "__main__":
    main()
