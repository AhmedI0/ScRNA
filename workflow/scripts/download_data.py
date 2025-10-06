#!/usr/bin/env python
"""Download the public GEO datasets listed in config/datasets.yaml.

Features
--------
* resumable downloads (HTTP Range requests), so an interrupted 2 GB transfer
  continues where it stopped;
* the remote Content-Length is checked against the local file size;
* SHA-256 checksums and timestamps are written to
  ``data/raw/download_manifest.tsv`` so every input file is traceable.

Only the Python standard library (plus PyYAML) is used, so the script runs
before the analysis environment exists.

Usage
-----
    python workflow/scripts/download_data.py                 # all cohorts
    python workflow/scripts/download_data.py --cohort discovery
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import sys
import time
import urllib.request
from pathlib import Path

try:
    import yaml
except ImportError:  # pragma: no cover - fallback for a bare interpreter
    yaml = None

ROOT = Path(__file__).resolve().parents[2]
CHUNK = 1 << 20  # 1 MiB


def load_config(path: Path) -> dict:
    if yaml is None:
        sys.exit("PyYAML is required: pip install pyyaml")
    with open(path, encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def remote_size(url: str) -> int | None:
    req = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(req, timeout=60) as resp:
        size = resp.headers.get("Content-Length")
    return int(size) if size is not None else None


def download(url: str, dest: Path, retries: int = 5) -> None:
    """Download ``url`` to ``dest``, resuming a partial ``.part`` file."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    expected = remote_size(url)
    if dest.exists() and (expected is None or dest.stat().st_size == expected):
        print(f"[skip] {dest.name} already complete")
        return

    part = dest.with_suffix(dest.suffix + ".part")
    for attempt in range(1, retries + 1):
        offset = part.stat().st_size if part.exists() else 0
        req = urllib.request.Request(url)
        if offset:
            req.add_header("Range", f"bytes={offset}-")
        try:
            with urllib.request.urlopen(req, timeout=120) as resp, open(part, "ab") as fh:
                done, t0, last = offset, time.time(), 0.0
                while chunk := resp.read(CHUNK):
                    fh.write(chunk)
                    done += len(chunk)
                    if time.time() - last > 15:
                        rate = (done - offset) / max(time.time() - t0, 1e-6) / 1e6
                        pct = f"{100 * done / expected:5.1f}%" if expected else "?"
                        print(f"  {dest.name}: {pct} ({done / 1e9:.2f} GB, {rate:.1f} MB/s)", flush=True)
                        last = time.time()
            break
        except Exception as exc:  # network hiccup -> resume
            print(f"  attempt {attempt}/{retries} failed: {exc}; retrying in 10 s", flush=True)
            time.sleep(10)
    else:
        raise RuntimeError(f"could not download {url}")

    if expected is not None and part.stat().st_size != expected:
        raise RuntimeError(f"size mismatch for {dest.name}: {part.stat().st_size} != {expected}")
    part.replace(dest)
    print(f"[done] {dest.name}")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while chunk := fh.read(CHUNK * 8):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=ROOT / "config" / "datasets.yaml", type=Path)
    ap.add_argument("--outdir", default=ROOT / "data" / "raw", type=Path)
    ap.add_argument("--cohort", choices=["discovery", "validation"], action="append")
    args = ap.parse_args()

    cfg = load_config(args.config)
    cohorts = args.cohort or list(cfg)
    manifest = args.outdir / "download_manifest.tsv"
    rows = []
    for cohort in cohorts:
        ds = cfg[cohort]
        outdir = args.outdir / ds["accession"]
        jobs = [(ds["base_url"] + f, outdir / f) for f in ds["files"]]
        jobs += [(e["url"], outdir / e["file"]) for e in ds.get("extra", [])]
        for url, dest in jobs:
            print(f"[{cohort}] {url}", flush=True)
            download(url, dest)
            rows.append((cohort, ds["accession"], dest.name, str(dest.stat().st_size),
                         sha256(dest), url, dt.datetime.now().isoformat(timespec="seconds")))

    header = "cohort\taccession\tfile\tbytes\tsha256\turl\tretrieved\n"
    existing = {}
    if manifest.exists():
        for line in manifest.read_text(encoding="utf-8").splitlines()[1:]:
            existing[line.split("\t")[2]] = line
    for r in rows:
        existing[r[2]] = "\t".join(r)
    manifest.write_text(header + "\n".join(existing[k] for k in sorted(existing)) + "\n", encoding="utf-8")
    print(f"manifest written: {manifest}")


if __name__ == "__main__":
    main()
