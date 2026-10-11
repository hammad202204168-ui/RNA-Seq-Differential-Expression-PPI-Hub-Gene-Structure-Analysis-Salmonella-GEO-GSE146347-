"""Annotated 'Volcano Plot Interpretation' chart (real contrast) and a schematic reading guide (simulated points)."""
import matplotlib.pyplot as plt
import numpy as np
from adjustText import adjust_text
from matplotlib.lines import Line2D

from plotstyle import DOWN, NS, UP, save

GUIDE_NEG = "Negative change in gene\nexpression compared to control"
GUIDE_POS = "Positive change in gene\nexpression compared to control"


def _title(fig, text="Volcano Plot Interpretation"):
    fig.add_artist(plt.Rectangle((0.012, 0.915), 0.007, 0.065, transform=fig.transFigure, color="black", lw=0))
    fig.text(0.03, 0.947, text, fontsize=13, fontweight="bold", va="center", ha="left")


def draw(fig, ax, x, y, status, xlim, ylim, thr_y, hi_y, thr_label, hi_label, xlabel, ylabel,
         text_xy, text_ha, vlines=(), labels=None, neg_hdr=GUIDE_NEG, pos_hdr=GUIDE_POS, size=14):
    ax.set_xlim(*xlim); ax.set_ylim(*ylim)
    ax.axhspan(0, thr_y, color="#F2F2F2", lw=0, zorder=0)
    x, y, status = np.asarray(x), np.asarray(y), np.asarray(status)
    for st, c, z in (("NS", NS, 2), ("Down", DOWN, 3), ("Up", UP, 3)):
        m = status == st
        ax.scatter(x[m], y[m], s=size, c=c, lw=0, alpha=0.85 if st != "NS" else 0.6, zorder=z, rasterized=True)
    ax.axvline(0, color="#777", lw=1.0, zorder=1)
    for v in vlines:
        ax.axvline(v, color="#999", lw=0.7, ls=(0, (3, 3)), zorder=1)
    for yy in (thr_y, hi_y):
        ax.axhline(yy, color="#888", lw=0.8, ls=(0, (3, 3)), zorder=1)
    xr = xlim[1] - xlim[0]
    ax.text(xlim[1] - 0.01 * xr, thr_y, thr_label, fontsize=6.3, color="#555", ha="right", va="bottom")
    ax.text(xlim[1] - 0.01 * xr, hi_y, hi_label, fontsize=6.3, color="#555", ha="right", va="bottom")
    tx, ty_ins, ty_sig, ty_hi = text_xy
    ax.text(tx, ty_hi, "Most significant change\nin differential gene\nexpression", fontsize=6.5, ha=text_ha, va="center")
    ax.text(tx, ty_sig, "Statistically significant\nchange in differential\ngene expression", fontsize=6.5, ha=text_ha, va="center")
    ax.text(tx, ty_ins, "Statistically insignificant\nchange in differential\ngene expression", fontsize=6.5, ha=text_ha, va="center")
    tr = ax.get_xaxis_transform()
    ax.text(0, 0.965, "Zero point", transform=tr, ha="center", va="center", fontsize=6.8,
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#999", lw=0.6), zorder=5)
    for sgn in (-1, 1):
        ax.annotate("", xy=(sgn * 0.30 * xr / 2, 0.965), xytext=(sgn * 0.07 * xr / 2, 0.965), xycoords=tr, textcoords=tr,
                    arrowprops=dict(arrowstyle="-|>", color="#999", lw=1.0), zorder=4)
    ax.text(xlim[0] + 0.25 * xr, 1.04, neg_hdr, transform=tr, ha="center", va="bottom", fontsize=7.2, color=DOWN,
            fontweight="bold", linespacing=1.25)
    ax.text(xlim[0] + 0.75 * xr, 1.04, pos_hdr, transform=tr, ha="center", va="bottom", fontsize=7.2, color=UP,
            fontweight="bold", linespacing=1.25)
    ax.legend(handles=[Line2D([], [], marker="o", ls="", c=DOWN, ms=5, label="Negative"),
                       Line2D([], [], marker="o", ls="", c=UP, ms=5, label="Positive"),
                       Line2D([], [], marker="o", ls="", c=NS, ms=5, label="Insignificant")],
              frameon=False, loc="upper left", bbox_to_anchor=(0.0, 0.90), handletextpad=0.2, borderaxespad=0.2)
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    if labels:
        texts = [ax.text(xv, yv, g, fontsize=6.3, style="italic", zorder=6) for g, xv, yv in labels]
        adjust_text(texts, ax=ax, expand=(1.6, 2.0), force_text=(0.6, 1.0),
                    arrowprops=dict(arrowstyle="-", lw=0.4, color="#555"))
    _title(fig)


def interpretation(res, contrast_label, stem, thr_fdr=0.05, hi_fdr=0.001, lfc_cut=1.0, n_label=10, ctrl_word="control"):
    r = res.dropna(subset=["FDR"]).copy()
    r["nlp"] = -np.log10(r.FDR.clip(lower=1e-300))
    xmax = float(np.ceil(np.abs(r.log2FC).max() + 0.5))
    ymax = float(np.ceil(r.nlp.max() * 1.12))
    side = -1 if (r.log2FC < -0.5 * xmax).sum() <= (r.log2FC > 0.5 * xmax).sum() else 1
    tx = side * 0.42 * xmax if side > 0 else -0.97 * xmax
    fig = plt.figure(figsize=(7.2, 4.3))
    ax = fig.add_axes([0.09, 0.14, 0.88, 0.62])
    lab = r[r.status != "NS"].sort_values("FDR").head(n_label)
    draw(fig, ax, r.log2FC, r.nlp, r.status, (-xmax, xmax), (0, ymax), -np.log10(thr_fdr), -np.log10(hi_fdr),
         f"Threshold (FDR = {thr_fdr})", f"FDR = {hi_fdr}", contrast_label, "\u2212log$_{10}$ (FDR)",
         (tx, -np.log10(thr_fdr) * 0.5, (-np.log10(thr_fdr) + -np.log10(hi_fdr)) / 2, (-np.log10(hi_fdr) + ymax) / 2 + 0.1 * ymax),
         "left", vlines=(-lfc_cut, lfc_cut), labels=[(row.label if "label" in lab.columns else g, row.log2FC, row.nlp) for g, row in lab.iterrows()], size=6)
    save(fig, stem)


def schematic(stem="Fig_volcano_reading_guide"):
    rng = np.random.default_rng(42)
    n = 420
    x = rng.normal(0, 0.45, n)
    y = np.minimum(np.abs(x) ** 0.9 * rng.gamma(3.0, 1.9, n) + rng.uniform(0, 0.4, n), 6.8)
    st = np.where((y > 2) & (x < -0.12), "Down", np.where((y > 2) & (x > 0.12), "Up", "NS"))
    fig = plt.figure(figsize=(7.2, 4.3))
    ax = fig.add_axes([0.09, 0.14, 0.88, 0.62])
    draw(fig, ax, x, y, st, (-2.5, 2.5), (0, 7.2), 2.0, 4.9, "Threshold", "", "Fold change (log$_2$)",
         "\u2212log$_{10}$ (P-value)", (-2.4, 1.0, 3.4, 6.0), "left", size=14)
    ax.set_xticks(np.arange(-2.5, 2.6, 0.5))
    fig.text(0.97, 0.025, "Schematic with simulated points \u2013 not study data", fontsize=6, color="#555", ha="right")
    save(fig, stem, note=False)
