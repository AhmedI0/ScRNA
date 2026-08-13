#!/usr/bin/env python
"""Identify the barcode read and the cDNA read of every 10x FASTQ pair.

SRA/ENA do not label files as "R1/R2": depending on the submission, ``_1`` may be
the 26-bp barcode read or a 150-bp read whose first 26 bases are the barcode.
For each run we sample reads from both files and measure
  * the read length distribution, and
  * the fraction of reads whose first 16 bases are in the 10x whitelist.
The barcode read has a whitelist hit rate of ~90 %; a cDNA read has ~0 %.

Output: <fastq_dir>/<library>/read_layout.tsv with columns
run, barcode_fastq, cdna_fastq, barcode_len, cdna_len, wl_rate_barcode, wl_rate_cdna
"""
from __future__ import annotations

import argparse
import gzip
import itertools
import statistics
import sys
from pathlib import Path


def sample_reads(path: Path, n: int = 20000) -> list[str]:
    seqs = []
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(itertools.islice(fh, 4 * n)):
            if i % 4 == 1:
                seqs.append(line.strip())
    return seqs


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("library_dir", type=Path)
    ap.add_argument("whitelist", type=Path)
    ap.add_argument("--cb-len", type=int, default=16)
    args = ap.parse_args()

    wl = set(args.whitelist.read_text().split())
    rows = []
    for f1 in sorted(args.library_dir.glob("*_1.fastq.gz")):
        f2 = f1.with_name(f1.name.replace("_1.fastq.gz", "_2.fastq.gz"))
        if not f2.exists():
            sys.exit(f"missing mate for {f1}")
        stats = {}
        for f in (f1, f2):
            s = sample_reads(f)
            hit = sum(r[: args.cb_len] in wl for r in s) / max(len(s), 1)
            stats[f] = (statistics.median(len(r) for r in s), hit)
        bc, cdna = (f1, f2) if stats[f1][1] >= stats[f2][1] else (f2, f1)
        if stats[bc][1] < 0.5:
            sys.exit(f"{f1.name}: no read matches the whitelist well ({stats[bc][1]:.2f}) - wrong chemistry?")
        rows.append([f1.name.split("_")[0], bc.name, cdna.name, int(stats[bc][0]), int(stats[cdna][0]),
                     round(stats[bc][1], 3), round(stats[cdna][1], 3)])
        print(f"{f1.name.split('_')[0]}: barcode={bc.name} ({stats[bc][0]:.0f} bp, WL {stats[bc][1]:.1%})  "
              f"cDNA={cdna.name} ({stats[cdna][0]:.0f} bp)")
    out = args.library_dir / "read_layout.tsv"
    with open(out, "w") as fh:
        fh.write("run\tbarcode_fastq\tcdna_fastq\tbarcode_len\tcdna_len\twl_rate_barcode\twl_rate_cdna\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")


if __name__ == "__main__":
    main()
