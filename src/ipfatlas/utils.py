"""Small utilities: logging, timing, environment report and R calls."""
from __future__ import annotations

import logging
import platform
import subprocess
import sys
import time
from contextlib import contextmanager
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import pandas as pd

log = logging.getLogger("ipfatlas")
if not log.handlers:
    _h = logging.StreamHandler(sys.stdout)
    _h.setFormatter(logging.Formatter("%(asctime)s | %(levelname)s | %(message)s", "%H:%M:%S"))
    log.addHandler(_h)
    log.setLevel(logging.INFO)


@contextmanager
def timer(msg: str):
    t0 = time.time()
    log.info(f"{msg} ...")
    yield
    log.info(f"{msg} done in {time.time() - t0:.1f} s")


def package_versions(packages: list[str] | None = None) -> pd.DataFrame:
    """Versions of the key packages (printed at the top of each notebook for reproducibility)."""
    packages = packages or [
        "python", "numpy", "scipy", "pandas", "anndata", "scanpy", "scvi-tools", "torch", "harmonypy",
        "scib-metrics", "celltypist", "decoupler", "liana", "cellrank", "pydeseq2", "gseapy", "matplotlib",
    ]
    rows = []
    for pkg in packages:
        if pkg == "python":
            rows.append(("python", platform.python_version()))
            continue
        try:
            rows.append((pkg, version(pkg)))
        except PackageNotFoundError:
            rows.append((pkg, "not installed"))
    return pd.DataFrame(rows, columns=["package", "version"]).set_index("package")


def run_r(rscript: str, script: str | Path, *args: str, log_file: str | Path | None = None) -> None:
    """Run an R script with ``Rscript`` and stream its output; raise if it fails.

    R is called as a separate process (rather than through rpy2) because that is
    robust on Windows and keeps the R and Python environments independent.
    """
    cmd = [rscript, "--vanilla", str(script), *map(str, args)]
    log.info("R: " + " ".join(cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True)
    out = proc.stdout + ("\n" + proc.stderr if proc.stderr else "")
    if log_file:
        Path(log_file).write_text(out, encoding="utf-8")
    tail = "\n".join(out.strip().splitlines()[-15:])
    if proc.returncode != 0:
        raise RuntimeError(f"R script failed ({script}):\n{tail}")
    log.info(f"R finished: {Path(script).name}\n{tail}")


def write_table(df: pd.DataFrame, path: str | Path, index: bool = True) -> Path:
    """Write a results table as CSV (and xlsx for small tables) creating parent dirs.

    Tables with more than 50,000 rows are written gzip-compressed (``.csv.gz``, which
    ``pd.read_csv`` reads directly) and without an xlsx copy; the fixed gzip timestamp
    keeps the file byte-identical between runs. Returns the path written.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if len(df) > 50_000:
        out = path.with_suffix(".csv.gz")
        df.to_csv(out, index=index, compression={"method": "gzip", "mtime": 0})
        return out
    out = path.with_suffix(".csv")
    df.to_csv(out, index=index)
    try:
        df.to_excel(path.with_suffix(".xlsx"), index=index)
    except Exception:  # openpyxl missing or sheet too wide
        pass
    return out
