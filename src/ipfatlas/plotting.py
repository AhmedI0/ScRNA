"""Publication figure style and a single save helper.

All figures go through :func:`save_figure`, which writes SVG (vector, editable
text), JPG and PNG at 300 dpi to ``results/figures/<section>/``. Using one
function guarantees identical format, resolution and naming for every panel.
"""
from __future__ import annotations

from collections.abc import Iterable, Sequence
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

MM = 1 / 25.4  # inches per millimetre

# Colour-blind-safe palettes ---------------------------------------------------
CONDITION_COLORS = {
    "Control": "#4C78A8",
    "IPF": "#E45756",
    "NSIP": "#F2CF5B",
    "cHP": "#54A24B",
    "Sarcoidosis": "#B279A2",
    "ILD": "#9D755D",
    "COPD": "#72B7B2",
}

# alveolar epithelial lineage (Fig. 4 and the replication figures): distinct, colour-blind safe
# Okabe-Ito colours, because the atlas palette gives all epithelial types similar oranges
LINEAGE_COLORS = {
    "AT2": "#0072B2",
    "Transitional AT2": "#CC79A7",
    "AT1": "#009E73",
    "Aberrant basaloid": "#D55E00",
    "Basal": "#E69F00",
}
COMPARTMENT_COLORS = {
    "Epithelial": "#E15759",
    "Immune": "#4E79A7",
    "Endothelial": "#F28E2B",
    "Mesenchymal": "#59A14F",
}
SITE_COLORS = {"Vanderbilt": "#1b9e77", "NTI": "#d95f02", "DNA": "#7570b3"}
# 32-colour qualitative palette (Kelly + Glasbey style), readable on white
CATEGORICAL = [
    "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", "#E377C2", "#7F7F7F",
    "#BCBD22", "#17BECF", "#AEC7E8", "#FFBB78", "#98DF8A", "#FF9896", "#C5B0D5", "#C49C94",
    "#F7B6D2", "#C7C7C7", "#DBDB8D", "#9EDAE5", "#393B79", "#637939", "#8C6D31", "#843C39",
    "#7B4173", "#3182BD", "#E6550D", "#31A354", "#756BB1", "#636363", "#6BAED6", "#FD8D3C",
    "#74C476", "#9E9AC8", "#969696", "#A55194",
]
DIVERGING_CMAP = "RdBu_r"
SEQUENTIAL_CMAP = "viridis"


def set_style(font_size: float = 7, font_family: str = "Arial") -> None:
    """Journal-style matplotlib defaults (small fonts, thin lines, editable SVG text)."""
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": [font_family, "Helvetica", "DejaVu Sans"],
            "font.size": font_size,
            "axes.titlesize": font_size + 1,
            "axes.labelsize": font_size,
            "xtick.labelsize": font_size - 1,
            "ytick.labelsize": font_size - 1,
            "legend.fontsize": font_size - 1,
            "legend.title_fontsize": font_size,
            "figure.titlesize": font_size + 2,
            "axes.linewidth": 0.6,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "xtick.major.size": 2.5,
            "ytick.major.size": 2.5,
            "lines.linewidth": 1.0,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "savefig.dpi": 300,
            "figure.dpi": 110,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.02,
            "svg.fonttype": "none",  # keep text as text in SVG (editable in Illustrator/Inkscape)
            "pdf.fonttype": 42,
            "image.cmap": SEQUENTIAL_CMAP,
        }
    )
    try:
        import scanpy as sc

        sc.settings.set_figure_params(dpi=110, dpi_save=300, frameon=False, vector_friendly=True, fontsize=font_size)
        # scanpy's style turns on grid lines and resets fonts; restore the journal style
        mpl.rcParams.update({"font.sans-serif": [font_family, "Helvetica", "DejaVu Sans"], "svg.fonttype": "none",
                             "axes.grid": False, "font.size": font_size, "axes.labelsize": font_size,
                             "xtick.labelsize": font_size - 1, "ytick.labelsize": font_size - 1,
                             "legend.fontsize": font_size - 1, "axes.titlesize": font_size + 1})
    except ImportError:
        pass


def save_figure(
    fig: plt.Figure,
    name: str,
    outdir: str | Path,
    formats: Iterable[str] = ("svg", "jpg", "png"),
    dpi: int = 300,
    close: bool = False,
) -> list[Path]:
    """Save ``fig`` as ``<outdir>/<name>.<fmt>`` for every format at ``dpi``.

    JPG is written with quality 95 and a white background (JPG has no alpha).
    Returns the list of written paths.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    written = []
    for fmt in formats:
        path = outdir / f"{name}.{fmt}"
        kwargs = {"dpi": dpi, "bbox_inches": "tight", "facecolor": "white"}
        if fmt in ("jpg", "jpeg"):
            kwargs["pil_kwargs"] = {"quality": 95, "dpi": (dpi, dpi)}
        fig.savefig(path, format="jpeg" if fmt == "jpg" else fmt, **kwargs)
        written.append(path)
    if close:
        plt.close(fig)
    return written


def figsize(width_mm: float, height_mm: float) -> tuple[float, float]:
    return width_mm * MM, height_mm * MM


def palette_for(categories: Sequence[str], base: dict | None = None) -> dict[str, str]:
    """Stable colour mapping: use ``base`` where defined, otherwise the categorical list."""
    base = base or {}
    out, i = {}, 0
    for c in categories:
        if c in base:
            out[c] = base[c]
        else:
            out[c] = CATEGORICAL[i % len(CATEGORICAL)]
            i += 1
    return out


# colour families per compartment (light -> dark); cell types of one compartment share a hue
_COMPARTMENT_RAMPS = {
    "Epithelial": ["#FDD0A2", "#FDAE6B", "#FD8D3C", "#F16913", "#D94801", "#A63603", "#E7298A", "#C51B7D",
                   "#8E0152", "#F4A582", "#B2182B", "#67001F", "#FEC44F"],
    "Immune": ["#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#08519C", "#08306B", "#BCBDDC",
               "#9E9AC8", "#807DBA", "#6A51A3", "#54278F", "#3F007D", "#80CDC1", "#35978F", "#01665E"],
    "Endothelial": ["#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45", "#005A32"],
    "Mesenchymal": ["#F6E8C3", "#DFC27D", "#BF812D", "#8C510A", "#543005", "#C2A5CF", "#9970AB", "#762A83"],
}


def compartment_palette(celltype_to_compartment: dict[str, str]) -> dict[str, str]:
    """Colours that group cell types by compartment (warm epithelium, blue/purple immune,
    green endothelium, brown/purple mesenchyme), in the order the cell types are given."""
    used: dict[str, int] = {}
    out = {}
    for ct, comp in celltype_to_compartment.items():
        ramp = _COMPARTMENT_RAMPS.get(comp, CATEGORICAL)
        i = used.get(comp, 0)
        out[ct] = ramp[i % len(ramp)]
        used[comp] = i + 1
    return out


def panel_label(ax: plt.Axes, label: str, x: float = -0.12, y: float = 1.06, size: float = 9) -> None:
    """Bold lowercase panel letter (Nature style) in axes coordinates."""
    ax.text(x, y, label, transform=ax.transAxes, fontsize=size, fontweight="bold", va="bottom", ha="left")


def panel_labels(fig: plt.Figure, axes: Sequence[plt.Axes], letters: str | Sequence[str], size: float = 9,
                 pad: float = 2) -> None:
    """Panel letters left of each axes' tick labels, level with the top of its title.

    Call after the final layout. Unlike :func:`panel_label` the position comes from
    the drawn extent of the axes, so letters cannot collide with long y tick labels
    or a centred multi-line title; ``pad`` (points) separates letter and content.
    """
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    to_fig = fig.transFigure.inverted()
    for ax, letter in zip(axes, letters):
        box = ax.get_tightbbox(renderer)
        x0, y1 = to_fig.transform((box.x0 - pad * fig.dpi / 72, box.y1))
        fig.text(x0, y1, letter, fontsize=size, fontweight="bold", ha="right", va="top")


def embedding_scatter(
    ax: plt.Axes,
    coords: np.ndarray,
    values: pd.Series | np.ndarray,
    palette: dict | None = None,
    cmap: str = SEQUENTIAL_CMAP,
    s: float = 0.6,
    rasterized: bool = True,
    order: str = "random",
    legend: bool = False,
    vmin: float | None = None,
    vmax: float | None = None,
    seed: int = 0,
):
    """Fast UMAP scatter with rasterised points (small SVG files, crisp text).

    Categorical ``values`` use ``palette``; numeric values use ``cmap``.
    Points are drawn in random order so no category hides another.
    """
    values = pd.Series(np.asarray(values))
    idx = np.arange(len(values))
    if order == "random":
        idx = np.random.default_rng(seed).permutation(len(values))
    elif order == "value" and pd.api.types.is_numeric_dtype(values):
        idx = np.argsort(values.to_numpy())
    x, y = coords[idx, 0], coords[idx, 1]
    v = values.iloc[idx]
    if pd.api.types.is_numeric_dtype(v) and not isinstance(v.dtype, pd.CategoricalDtype):
        sca = ax.scatter(x, y, c=v.to_numpy(), s=s, cmap=cmap, vmin=vmin, vmax=vmax, lw=0, rasterized=rasterized)
    else:
        cats = list(pd.Categorical(values).categories)
        palette = palette or palette_for(cats)
        sca = ax.scatter(x, y, c=v.map(palette).to_numpy(), s=s, lw=0, rasterized=rasterized)
        if legend:
            handles = [
                plt.Line2D([], [], marker="o", ls="", color=palette[c], markersize=3, label=c) for c in cats
            ]
            ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.0, 0.5), handletextpad=0.1, borderaxespad=0)
    ax.set_xticks([])
    ax.set_yticks([])
    for side in ("left", "bottom", "top", "right"):
        ax.spines[side].set_visible(False)
    ax.set_aspect("equal", adjustable="datalim")
    return sca


def umap_axis_arrows(ax: plt.Axes, label: str = "UMAP", length: float = 0.12) -> None:
    """Small L-shaped axis indicator in the lower-left corner instead of full axes."""
    kw = {"xycoords": "axes fraction", "textcoords": "axes fraction",
          "arrowprops": {"arrowstyle": "-|>", "lw": 0.6, "color": "black", "shrinkA": 0, "shrinkB": 0}}
    ax.annotate("", xy=(length, 0.0), xytext=(0, 0), **kw)
    ax.annotate("", xy=(0, length), xytext=(0, 0), **kw)
    ax.text(length / 2, -0.035, f"{label}1", transform=ax.transAxes, ha="center", va="top", fontsize=5)
    ax.text(-0.035, length / 2, f"{label}2", transform=ax.transAxes, ha="right", va="center", rotation=90, fontsize=5)


def label_centroids(ax: plt.Axes, coords: np.ndarray, labels: pd.Series, fontsize: float = 5, min_cells: int = 30,
                    numbered: bool = False, order: Sequence[str] | None = None) -> dict:
    """Write each category at the densest point of its cells (with a white halo).

    The densest bin of a 2-D histogram, not the median: a few scattered cells
    would pull the median off the cluster body. With ``numbered=True`` a number
    is written instead of the name (for dense atlases; pair it with a legend);
    categories are numbered in ``order`` (default: alphabetical).
    Returns ``{category: text written}``.
    """
    import matplotlib.patheffects as pe

    labels = pd.Series(np.asarray(labels))
    extent = [[coords[:, 0].min(), coords[:, 0].max()], [coords[:, 1].min(), coords[:, 1].max()]]
    written = {}
    for i, cat in enumerate(order if order is not None else pd.Categorical(labels).categories, 1):
        mask = (labels == cat).to_numpy()
        if mask.sum() < min_cells:
            continue
        h, xe, ye = np.histogram2d(coords[mask, 0], coords[mask, 1], bins=60, range=extent)
        bx, by = np.unravel_index(h.argmax(), h.shape)
        cx, cy = (xe[bx] + xe[bx + 1]) / 2, (ye[by] + ye[by + 1]) / 2
        written[cat] = str(i) if numbered else cat
        ax.text(
            cx, cy, written[cat], fontsize=fontsize, ha="center", va="center",
            fontweight="bold" if numbered else "normal",
            path_effects=[pe.withStroke(linewidth=1.5, foreground="white")],
        )
    return written
