"""Nature journal style for matplotlib figures.

Usage:
    from nature_style import apply_rcparams, save_fig

    apply_rcparams()  # call once at the start of main()
    fig, ax = plt.subplots(figsize=(3.5, 2.5))
    # ... plotting ...
    save_fig(fig, "plots/fig1_coverage.png")  # writes both .png and .pdf

The styling spec follows Nature family journal requirements:
- Arial/Helvetica, 7pt axis/tick labels, 8pt panel titles
- 1-column = 3.5" (89 mm), 2-column = 7.2" (183 mm)
- Top and right spines removed
- 0.75pt axis and tick line width
- 300 dpi for export, 150 dpi for display
- Legend without frame, 6pt font
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt


# Standard figure widths in inches (Nature spec)
COL1_INCHES = 3.5  # ~89 mm
COL2_INCHES = 7.2  # ~183 mm

# Panel label position in axes coordinates (top-left, slightly outside)
PANEL_LABEL_XY = (-0.12, 1.05)


NATURE_RCPARAMS = {
    # Fonts
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "axes.labelsize": 7,
    "axes.titlesize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 6,
    "figure.titlesize": 8,
    # Spines and lines
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.75,
    "xtick.major.width": 0.75,
    "ytick.major.width": 0.75,
    "xtick.minor.width": 0.5,
    "ytick.minor.width": 0.5,
    "lines.linewidth": 1.0,
    # Resolution
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
    # Legend
    "legend.frameon": False,
    # PDF/PS: keep text as text rather than outlines, so editors can edit it
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
}


def apply_rcparams(extra: dict | None = None) -> None:
    """Apply the Nature-style rcParams globally.

    Call once at the start of `main()`. If you need to override or add
    entries, pass them via `extra` rather than mutating rcParams elsewhere.
    """
    matplotlib.rcParams.update(NATURE_RCPARAMS)
    if extra:
        matplotlib.rcParams.update(extra)


def save_fig(fig, path, also_pdf: bool = True, also_eps: bool = False) -> Path:
    """Save a figure as PNG (for the report) and PDF (for Adobe Illustrator).

    PDF is preferred over EPS as the editable vector format: it preserves
    fonts as text (with ``pdf.fonttype = 42`` set in rcParams), supports
    transparency, and is the modern default for figure submission and
    Illustrator workflows.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The figure to save.
    path : str or pathlib.Path
        Output path. If the suffix is not ``.png``, it is replaced with ``.png``.
    also_pdf : bool, default True
        Also write a sibling ``.pdf`` file. Set to False for quick exploratory
        figures where you don't need the editable vector version.
    also_eps : bool, default False
        Also write a sibling ``.eps`` file. Off by default; enable only if a
        downstream tool specifically requires EPS.

    Returns
    -------
    pathlib.Path
        The path to the written PNG.
    """
    path = Path(path).with_suffix(".png")
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, format="png")
    if also_pdf:
        fig.savefig(path.with_suffix(".pdf"), format="pdf")
    if also_eps:
        fig.savefig(path.with_suffix(".eps"), format="eps")
    return path


def add_panel_label(ax, letter: str, xy: tuple[float, float] = PANEL_LABEL_XY) -> None:
    """Add a bold lowercase panel label (a, b, c, ...) to an axes.

    Places the label at the standard Nature position (top-left, slightly
    outside the axes box) using axes coordinates.
    """
    ax.text(
        xy[0],
        xy[1],
        letter,
        transform=ax.transAxes,
        fontsize=8,
        fontweight="bold",
        va="bottom",
        ha="left",
    )


def figsize(columns: int = 1, aspect: float = 0.75) -> tuple[float, float]:
    """Return a (width, height) tuple for a Nature-spec figure.

    Parameters
    ----------
    columns : {1, 2}
        1 for single-column (3.5 in), 2 for double-column (7.2 in).
    aspect : float
        Height as a fraction of width. Default 0.75 (4:3-ish).
    """
    width = COL1_INCHES if columns == 1 else COL2_INCHES
    return (width, width * aspect)


# ---------------------------------------------------------------------------
# Schematic / cartoon helpers
#
# For data-driven methods schematics (boxes + arrows that should regenerate
# with the analysis, e.g. a pipeline with real counts at each step). For purely
# conceptual cartoons that never change with the data, hand-author an .svg
# instead (see SKILL.md "Schematics and cartoons").
# ---------------------------------------------------------------------------

# Muted palette that reads on white and in grayscale print.
SCHEMATIC_FILL = "#e8eef5"
SCHEMATIC_EDGE = "#33475b"
SCHEMATIC_ACCENT = "#c44e52"


def new_schematic(width_columns: int = 2, aspect: float = 0.42):
    """Create a blank axes sized for a schematic, with no spines/ticks.

    Returns ``(fig, ax)`` with a 0–100 × 0–100 coordinate space so boxes and
    arrows can be placed in intuitive percentages. Call ``save_fig`` as usual.
    """
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=figsize(width_columns, aspect))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.set_aspect("auto")
    ax.axis("off")
    return fig, ax


def add_box(
    ax,
    x: float,
    y: float,
    text: str,
    width: float = 22,
    height: float = 16,
    fill: str = SCHEMATIC_FILL,
    edge: str = SCHEMATIC_EDGE,
    text_color: str = "#1a1a1a",
    fontsize: float = 7,
):
    """Draw a rounded labeled box centered at ``(x, y)`` in 0–100 coordinates.

    Returns ``(x, y, width, height)`` so callers can anchor arrows to edges.
    """
    from matplotlib.patches import FancyBboxPatch

    box = FancyBboxPatch(
        (x - width / 2, y - height / 2),
        width,
        height,
        boxstyle="round,pad=0.6,rounding_size=2",
        linewidth=0.75,
        edgecolor=edge,
        facecolor=fill,
    )
    ax.add_patch(box)
    ax.text(x, y, text, ha="center", va="center", fontsize=fontsize,
            color=text_color, wrap=True)
    return (x, y, width, height)


def add_arrow(ax, start, end, color: str = SCHEMATIC_EDGE, label: str | None = None):
    """Draw an arrow between two ``(x, y)`` points in 0–100 coordinates.

    ``start``/``end`` may be points or ``(x, y, w, h)`` box tuples from
    :func:`add_box` (the arrow then anchors to the facing box edges).
    """
    sx, sy = start[0], start[1]
    ex, ey = end[0], end[1]
    if len(start) == 4 and len(end) == 4:  # box -> box: anchor to facing edges
        sx = start[0] + (start[2] / 2 if end[0] > start[0] else -start[2] / 2)
        ex = end[0] + (-end[2] / 2 if end[0] > start[0] else end[2] / 2)
    ax.annotate(
        "", xy=(ex, ey), xytext=(sx, sy),
        arrowprops=dict(arrowstyle="-|>", color=color, linewidth=0.9,
                        shrinkA=0, shrinkB=0),
    )
    if label:
        ax.text((sx + ex) / 2, (sy + ey) / 2 + 4, label, ha="center",
                va="bottom", fontsize=6, color=color)


__all__ = [
    "NATURE_RCPARAMS",
    "COL1_INCHES",
    "COL2_INCHES",
    "PANEL_LABEL_XY",
    "SCHEMATIC_FILL",
    "SCHEMATIC_EDGE",
    "SCHEMATIC_ACCENT",
    "apply_rcparams",
    "save_fig",
    "add_panel_label",
    "figsize",
    "new_schematic",
    "add_box",
    "add_arrow",
]
