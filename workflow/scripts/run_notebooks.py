#!/usr/bin/env python
"""Convert the Jupytext sources in notebooks/py/ to .ipynb and execute them in order.

The ``.py`` (percent format) files are the version-controlled source; the
executed ``.ipynb`` files (with outputs and figures) are produced by this
script with papermill, so every notebook in the repository is reproducible
from a clean environment.

Usage
-----
    python workflow/scripts/run_notebooks.py                # all notebooks in order
    python workflow/scripts/run_notebooks.py 03 04          # selected notebooks
    python workflow/scripts/run_notebooks.py --from 05      # 05 and everything after it
"""
from __future__ import annotations

import argparse
import json
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "notebooks" / "py"
OUT = ROOT / "notebooks"
LOGS = ROOT / "logs" / "notebooks"
# notebooks that analyse results returned from the Linux upstream package
NEEDS_LINUX = {"14", "15"}
LINUX_RESULTS = ROOT / "data" / "linux_results" / "results_to_return"


def activate_env() -> None:
    """Emulate `conda activate` for this process and the kernels it starts.

    On Windows, MKL loads some of its DLLs lazily through PATH. Without the
    environment's Library\\bin first on PATH, an older MKL from the base
    installation can be picked up and the kernel crashes (0xc06d007f) at the
    first BLAS call. Prepending the env folders avoids that.
    """
    prefix = Path(sys.prefix)
    if platform.system() == "Windows":
        dirs = [prefix, prefix / "Library" / "mingw-w64" / "bin", prefix / "Library" / "usr" / "bin",
                prefix / "Library" / "bin", prefix / "Scripts", prefix / "bin"]
    else:
        dirs = [prefix / "bin"]
    os.environ["PATH"] = os.pathsep.join([str(d) for d in dirs] + [os.environ.get("PATH", "")])
    os.environ["CONDA_PREFIX"] = str(prefix)
    os.environ.setdefault("PYTHONUTF8", "1")
    tmp = ROOT / ".cache" / "tmp"
    tmp.mkdir(parents=True, exist_ok=True)
    for var in ("TMP", "TEMP", "TMPDIR"):
        os.environ[var] = str(tmp)
    # kernel history and connection files, which IPython/Jupyter otherwise keep in the user profile
    for var, sub in (("IPYTHONDIR", "ipython"), ("JUPYTER_RUNTIME_DIR", "jupyter_runtime")):
        (ROOT / ".cache" / sub).mkdir(parents=True, exist_ok=True)
        os.environ[var] = str(ROOT / ".cache" / sub)


def run(nb_py: Path, kernel: str) -> float:
    ipynb = OUT / f"{nb_py.stem}.ipynb"
    subprocess.run([sys.executable, "-m", "jupytext", "--to", "ipynb", "--update-metadata",
                    json.dumps({"kernelspec": {"name": kernel, "display_name": "Python (scrna-ipf)", "language": "python"}}),
                    str(nb_py), "-o", str(ipynb)], check=True)
    t0 = time.time()
    LOGS.mkdir(parents=True, exist_ok=True)
    with open(LOGS / f"{nb_py.stem}.log", "w", encoding="utf-8") as log:
        proc = subprocess.run(
            [sys.executable, "-m", "papermill", str(ipynb), str(ipynb), "-k", kernel, "--cwd", str(OUT),
             "--log-output", "--request-save-on-cell-execute", "--autosave-cell-every", "60"],
            stdout=log, stderr=subprocess.STDOUT,
        )
    dt = time.time() - t0
    status = "OK" if proc.returncode == 0 else f"FAILED (see logs/notebooks/{nb_py.stem}.log)"
    print(f"{nb_py.stem:45s} {status:10s} {dt / 60:6.1f} min", flush=True)
    if proc.returncode != 0:
        sys.exit(proc.returncode)
    return dt


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("select", nargs="*", help="notebook number prefixes, e.g. 01 02")
    ap.add_argument("--from", dest="start", help="run this notebook and all later ones")
    ap.add_argument("--kernel", default="scrna-ipf")
    args = ap.parse_args()
    activate_env()

    nbs = sorted(SRC.glob("[0-9][0-9]_*.py"))
    if args.select:
        nbs = [nb for nb in nbs if nb.name[:2] in args.select]
    if args.start:
        nbs = [nb for nb in nbs if nb.name[:2] >= args.start]
    if not LINUX_RESULTS.exists():
        skipped = [nb.stem for nb in nbs if nb.name[:2] in NEEDS_LINUX]
        nbs = [nb for nb in nbs if nb.name[:2] not in NEEDS_LINUX]
        for name in skipped:
            print(f"{name:45s} SKIPPED (no Linux results in {LINUX_RESULTS.relative_to(ROOT)})", flush=True)
    total = sum(run(nb, args.kernel) for nb in nbs)
    print(f"all done in {total / 60:.1f} min")


if __name__ == "__main__":
    main()
