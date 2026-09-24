#!/usr/bin/env python
"""CellBender remove-background on the raw STARsolo matrices returned from the Linux package.

CellBender fits its model in minutes on a CUDA GPU and in hours on a CPU. When the
Linux workstation has no NVIDIA GPU, step 07 of ``linux_upstream/`` is skipped
there (``USE_GPU=no``) and step 10 returns the raw droplet matrices instead; this
script then runs CellBender on the analysis PC with the parameters of
``linux_upstream/scripts/07_cellbender.sh`` and writes the output where
Notebook 14 reads it:

    data/linux_results/results_to_return/starsolo/<library>/Gene/raw/   input
    data/linux_results/results_to_return/cellbender/<library>/           output

Usage
-----
    python workflow/scripts/run_cellbender.py              # every library without output
    python workflow/scripts/run_cellbender.py F01214       # selected libraries
"""
from __future__ import annotations

import argparse
import os
import platform
import shutil
import subprocess
import sys
import time
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "data" / "linux_results" / "results_to_return"
LOGS = ROOT / "logs" / "cellbender"


def expected_cells(summary: Path) -> int:
    """STARsolo's EmptyDrops_CR cell count, as in step 07."""
    for line in summary.read_text().splitlines():
        if line.startswith("Estimated Number of Cells"):
            return int(float(line.split(",")[1]))
    raise SystemExit(f"{summary}: no 'Estimated Number of Cells' line")


def env_command(env: str) -> tuple[list[str], dict]:
    """CellBender (through cellbender_compat.py) in ``env``, with an activated PATH.

    An empty ``env`` means the Python running this script, which must have CellBender.
    """
    compat = str(Path(__file__).with_name("cellbender_compat.py"))
    cache = ROOT / ".cache"  # matplotlib/torch caches stay in the project, as in config.set_cache_dirs
    run_env = {"MPLCONFIGDIR": str(cache / "matplotlib"), "TORCH_HOME": str(cache / "torch"),
               "XDG_CACHE_HOME": str(cache / "xdg"), **os.environ}
    if not env:
        return [sys.executable, compat], run_env
    prefix = Path(env)
    if platform.system() == "Windows":
        # PyTables and PyTorch load DLLs from Library\bin, as after `conda activate`
        dirs = [prefix, prefix / "Library" / "bin", prefix / "Scripts"]
        python = prefix / "python.exe"
    else:
        dirs, python = [prefix / "bin"], prefix / "bin" / "python"
    run_env["PATH"] = os.pathsep.join([*map(str, dirs), run_env.get("PATH", "")])
    if not python.exists():
        raise SystemExit(f"no Python at {python}; create the environment from envs/cellbender.yml")
    return [str(python), compat], run_env


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("libraries", nargs="*", help="library ids (default: all with a raw matrix)")
    ap.add_argument("--cpu", action="store_true", help="run without --cuda (hours per library)")
    args = ap.parse_args()

    cfg = yaml.safe_load((ROOT / "config" / "config.yaml").read_text(encoding="utf-8"))["cellbender"]
    cmd0, run_env = env_command(cfg["env_windows" if platform.system() == "Windows" else "env_linux"])
    solo = RESULTS / "starsolo"
    libs = args.libraries or sorted(p.name for p in solo.iterdir() if (p / "Gene" / "raw").is_dir())
    if not libs:
        raise SystemExit(f"no raw matrices under {solo} (run linux_upstream step 10 with USE_GPU=no)")
    LOGS.mkdir(parents=True, exist_ok=True)

    failed = []
    for lib in libs:
        raw, out = solo / lib / "Gene" / "raw", RESULTS / "cellbender" / lib
        if (out / "cellbender_filtered.h5").exists():
            print(f"[skip] {lib}: output exists")
            continue
        n = expected_cells(solo / lib / "Gene" / "Summary.csv")
        ndrop = max(cfg["droplets_per_cell"] * n, cfg["min_droplets"])
        out.mkdir(parents=True, exist_ok=True)
        cmd = [*cmd0, "remove-background", "--input", str(raw), "--output", str(out / "cellbender.h5"),
               "--expected-cells", str(n), "--total-droplets-included", str(ndrop),
               "--fpr", str(cfg["fpr"]), "--epochs", str(cfg["epochs"])]
        if not args.cpu:
            cmd.append("--cuda")
        print(f"{lib}: expected cells {n:,}, droplets included {ndrop:,}", flush=True)
        t0, log_path = time.time(), LOGS / f"{lib}.log"
        tmp = out / "_tmp"  # CellBender's temporary checkpoint files stay next to the output
        tmp.mkdir(exist_ok=True)
        lib_env = {**run_env, "TMP": str(tmp), "TEMP": str(tmp), "TMPDIR": str(tmp)}
        with open(log_path, "w", encoding="utf-8") as fh:  # checkpoints go to the cwd
            code = subprocess.run(cmd, cwd=out, env=lib_env, stdout=fh, stderr=subprocess.STDOUT).returncode
        shutil.rmtree(tmp, ignore_errors=True)
        minutes = (time.time() - t0) / 60
        if code == 0 and (out / "cellbender_filtered.h5").exists():
            print(f"{lib}: done in {minutes:.1f} min", flush=True)
        else:
            print(f"{lib}: FAILED (exit {code}) after {minutes:.1f} min, see {log_path}", flush=True)
            failed.append(lib)
    if failed:
        raise SystemExit(f"CellBender failed for: {' '.join(failed)}")


if __name__ == "__main__":
    main()
