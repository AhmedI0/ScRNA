#!/usr/bin/env python
"""GRNBoost2 over a fork-shared expression matrix.

Same regressions as pySCENIC's ``arboreto_with_multiprocessing.py`` (arboreto's
``infer_partial_network`` with the GRNBoost2 settings and the same per-target seed),
so the adjacencies are identical for a given seed. That script passes the whole
expression and TF matrices to the pool with every target gene; with ~6,000 cells x
13,000 genes the parent process spends its time pickling and the workers idle
(about one gene per second). Here the matrices are module globals set before the
pool forks, so workers read them from shared memory and only gene indices travel.

    grnboost2_fork.py <loom> <tfs.txt> --output adjacencies.tsv --num_workers 16 --seed 1
"""
from __future__ import annotations

import argparse
import multiprocessing as mp
import time
from pathlib import PurePath

import pandas as pd
import tqdm
from arboreto.algo import _prepare_input
from arboreto.core import (
    EARLY_STOP_WINDOW_LENGTH,
    SGBM_KWARGS,
    infer_partial_network,
    target_gene_indices,
    to_tf_matrix,
)
from arboreto.utils import load_tf_names
from pyscenic.cli.utils import load_exp_matrix, suffixes_to_separator

_SHARED: dict = {}  # filled in the parent before fork; read-only in the workers


def _infer(target: int) -> pd.DataFrame:
    s = _SHARED
    return infer_partial_network(
        regressor_type="GBM", regressor_kwargs=SGBM_KWARGS,
        tf_matrix=s["tf_matrix"], tf_matrix_gene_names=s["tf_names"],
        target_gene_name=s["gene_names"][target], target_gene_expression=s["ex"][:, target],
        include_meta=False, early_stop_window_length=EARLY_STOP_WINDOW_LENGTH, seed=s["seed"],
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("expression_mtx_fname")
    ap.add_argument("tfs_fname")
    ap.add_argument("--output", required=True)
    ap.add_argument("--num_workers", type=int, default=mp.cpu_count())
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--chunksize", type=int, default=4)
    args = ap.parse_args()

    t0 = time.time()
    ex = load_exp_matrix(args.expression_mtx_fname, False, False, "CellID", "Gene")
    ex, gene_names, tf_names = _prepare_input(ex, ex.columns, load_tf_names(args.tfs_fname))
    tf_matrix, tf_matrix_names = to_tf_matrix(ex, gene_names, tf_names)
    _SHARED.update(ex=ex, gene_names=gene_names, tf_matrix=tf_matrix, tf_names=tf_matrix_names,
                   seed=args.seed)
    print(f"{ex.shape[0]} cells, {ex.shape[1]} genes, {len(tf_matrix_names)} TFs "
          f"(loaded in {time.time() - t0:.0f} s); {args.num_workers} workers", flush=True)

    targets = list(target_gene_indices(gene_names, target_genes="all"))
    t0 = time.time()
    with mp.get_context("fork").Pool(args.num_workers) as pool:
        adjs = list(tqdm.tqdm(pool.imap(_infer, targets, chunksize=args.chunksize),
                              total=len(targets), mininterval=30))
    adj = pd.concat(adjs).sort_values(by="importance", ascending=False)
    print(f"GRNBoost2 done in {(time.time() - t0) / 60:.1f} min: {len(adj):,} TF-target links", flush=True)
    adj.to_csv(args.output, index=False, sep=suffixes_to_separator(PurePath(args.output).suffixes))


if __name__ == "__main__":
    main()
