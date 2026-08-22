#!/usr/bin/env python
"""Convert STARsolo output of one library into a single AnnData (.h5ad).

  X                   Gene (exonic) UMI counts of EmptyDrops-called cells
  layers["genefull"]  GeneFull (exon + intron) counts
  layers["spliced"], layers["unspliced"], layers["ambiguous"]   Velocyto counts
Cell barcodes are prefixed with the library id ("F01157_AAACCTG...") to match
the GEO matrix of the original study.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import pandas as pd
import scipy.io
import scipy.sparse as sp


def read_mtx_dir(d: Path, matrix: str = "matrix.mtx") -> tuple[sp.csr_matrix, list[str], pd.DataFrame]:
    def pick(name):
        for cand in (d / name, d / f"{name}.gz"):
            if cand.exists():
                return cand
        raise FileNotFoundError(d / name)

    m = sp.csr_matrix(scipy.io.mmread(pick(matrix)).T)  # cells x genes
    bc = pd.read_csv(pick("barcodes.tsv"), header=None)[0].tolist()
    feat = pd.read_csv(pick("features.tsv"), header=None, sep="\t")
    return m, bc, feat


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("solo_dir", type=Path, help="<starsolo>/<library>/Solo.out")
    ap.add_argument("library_id")
    ap.add_argument("out", type=Path)
    args = ap.parse_args()

    x, bc, feat = read_mtx_dir(args.solo_dir / "Gene" / "filtered")
    var = pd.DataFrame({"gene_ids": feat[0].values}, index=feat[1].astype(str).values)
    adata = ad.AnnData(X=x.astype("float32"), obs=pd.DataFrame(index=bc), var=var)
    adata.var_names_make_unique()

    gf, gf_bc, _ = read_mtx_dir(args.solo_dir / "GeneFull" / "filtered")
    adata.layers["genefull"] = sp.csr_matrix(gf[pd.Index(gf_bc).get_indexer(bc)]).astype("float32") \
        if set(bc) <= set(gf_bc) else sp.csr_matrix(adata.shape, dtype="float32")

    vdir = args.solo_dir / "Velocyto" / "raw"
    if vdir.exists():
        idx = None
        for layer in ("spliced", "unspliced", "ambiguous"):
            m, vbc, _ = read_mtx_dir(vdir, f"{layer}.mtx")
            if idx is None:
                idx = pd.Index(vbc).get_indexer(bc)
            adata.layers[layer] = m[idx].astype("float32")

    adata.obs_names = [f"{args.library_id}_{b}" for b in adata.obs_names]
    adata.obs["library_id"] = args.library_id
    adata.obs["n_umi_gene"] = adata.X.sum(axis=1).A1
    if "spliced" in adata.layers:
        s = adata.layers["spliced"].sum(axis=1).A1
        u = adata.layers["unspliced"].sum(axis=1).A1
        adata.obs["frac_unspliced"] = u / (s + u + 1e-9)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.out, compression="gzip")
    print(f"{args.library_id}: {adata.n_obs} cells x {adata.n_vars} genes -> {args.out}")


if __name__ == "__main__":
    main()
