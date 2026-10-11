"""Shared figure style (journal-like: 7.2 in wide, 600 dpi PNG + vector PDF) for the real-data analysis."""
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "results" / "figures"
TAB = ROOT / "results" / "tables"
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

UP, DOWN, NS = "#C0392B", "#2E6DA4", "#BDBDBD"
GROUP_ORDER = ["WT_control", "WT_LFS", "delta_hfq_control", "delta_hfq_LFS"]
GROUP_COL = {"WT_control": "#0072B2", "WT_LFS": "#56B4E9", "delta_hfq_control": "#D55E00", "delta_hfq_LFS": "#E69F00"}
GROUP_LAB = {"WT_control": "WT, control", "WT_LFS": "WT, LFS", "delta_hfq_control": "\u0394hfq, control", "delta_hfq_LFS": "\u0394hfq, LFS"}
GROUP_MARK = {"WT_control": "o", "WT_LFS": "^", "delta_hfq_control": "o", "delta_hfq_LFS": "^"}
DATA_NOTE = "GEO GSE146347, Salmonella-only TPM table \u2013 exploratory TPM-based analysis (not raw counts, not DESeq2)"

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 7.5, "axes.labelsize": 8, "axes.titlesize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "legend.fontsize": 6.8, "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "xtick.major.size": 3, "ytick.major.size": 3, "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42, "ps.fonttype": 42, "figure.dpi": 150, "savefig.dpi": 600,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.05,
})


def panel(ax, letter, dx=-0.16, dy=1.04):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="left")


def save(fig, name, note=True, text=None):
    msg = text or DATA_NOTE
    if note and msg:
        fig.canvas.draw()
        bb = fig.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans + fig.transFigure.inverted())
        fig.text(bb.x0 + 0.004, bb.y0 - 0.006, msg, fontsize=5.3, color="#555", ha="left", va="top")
    fig.savefig(FIG / f"{name}.png")
    fig.savefig(FIG / f"{name}.pdf")
    plt.close(fig)
    print(f"  saved {name}")
