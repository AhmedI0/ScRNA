"""Project configuration, paths and runtime environment.

Every notebook starts with::

    from ipfatlas.config import setup
    cfg, P = setup()

which loads ``config/config.yaml``, fixes random seeds, redirects library
caches into ``.cache/`` (so nothing large lands on the system drive) and
returns a :class:`Paths` object with absolute project paths.
"""
from __future__ import annotations

import os
import platform
import random
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[2]


@dataclass(frozen=True)
class Paths:
    root: Path
    raw: Path
    interim: Path
    processed: Path
    external: Path
    linux_results: Path
    figures: Path
    tables: Path
    r_scripts: Path
    cache: Path

    def mkdirs(self) -> None:
        for p in (self.interim, self.processed, self.external, self.figures, self.tables, self.cache):
            p.mkdir(parents=True, exist_ok=True)


def load_config(path: str | Path | None = None) -> dict:
    path = Path(path) if path else ROOT / "config" / "config.yaml"
    with open(path, encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def get_paths(cfg: dict) -> Paths:
    p = cfg["paths"]
    return Paths(
        root=ROOT,
        raw=ROOT / p["raw"],
        interim=ROOT / p["interim"],
        processed=ROOT / p["processed"],
        external=ROOT / p["external"],
        linux_results=ROOT / p["linux_results"],
        figures=ROOT / p["figures"],
        tables=ROOT / p["tables"],
        r_scripts=ROOT / p["r_scripts"],
        cache=ROOT / ".cache",
    )


def set_cache_dirs(cache: Path) -> None:
    """Point third-party caches (matplotlib, numba, CellTypist, ...) at ``cache``."""
    env = {
        "MPLCONFIGDIR": cache / "matplotlib",
        "NUMBA_CACHE_DIR": cache / "numba",
        "CELLTYPIST_FOLDER": cache / "celltypist",
        "SCVI_CACHE": cache / "scvi",
        "XDG_CACHE_HOME": cache / "xdg",
        "PYSTOW_HOME": cache / "pystow",
        "OMNIPATH_CACHE": cache / "omnipath",
    }
    for key, val in env.items():
        val.mkdir(parents=True, exist_ok=True)
        os.environ.setdefault(key, str(val))


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    os.environ["PYTHONHASHSEED"] = str(seed)
    try:
        import torch

        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)
    except ImportError:
        pass


def rscript(cfg: dict) -> str:
    """Path to the Rscript executable for this operating system."""
    key = "rscript_windows" if platform.system() == "Windows" else "rscript_linux"
    return cfg["r"][key]


def setup(config_path: str | Path | None = None, seed: int | None = None) -> tuple[dict, Paths]:
    cfg = load_config(config_path)
    paths = get_paths(cfg)
    paths.mkdirs()
    set_cache_dirs(paths.cache)
    set_seed(cfg["project"]["seed"] if seed is None else seed)
    return cfg, paths
