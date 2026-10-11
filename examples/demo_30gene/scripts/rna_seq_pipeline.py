#!/usr/bin/env python3
"""
RNA-Seq differential expression, functional enrichment and hub-gene pipeline
Salmonella enterica serovar Typhimurium: wild type vs. delta-hfq (plus acid / heat stress)

Workflow (mirrors the reference-based workflow of Batut et al., 2021, Methods Mol Biol 2284:367-392):
    count matrix -> QC -> DESeq2 normalisation + Wald test (PyDESeq2) -> volcano / MA / heatmap
    -> functional enrichment -> STRING network hub analysis

Usage (from anywhere):
    python scripts/rna_seq_pipeline.py
Outputs:
    results/figures/*.png (600 dpi) and *.pdf (vector)
    results/tables/*.csv
"""
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from adjustText import adjust_text
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Patch
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage
from scipy.spatial.distance import pdist, squareform
from scipy.stats import fisher_exact
from sklearn.decomposition import PCA

if __package__:
    from .structure_utils import is_hfq_pdb, pdb_header_text
else:
    from structure_utils import is_hfq_pdb, pdb_header_text

ROOT = Path(__file__).resolve().parents[1]
DATA, FIG, TAB = ROOT / "data", ROOT / "results" / "figures", ROOT / "results" / "tables"
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

PADJ_CUTOFF = 0.05
LFC_CUTOFF = 1.0          # |log2FC| >= 1  <=>  fold change >= 2 (same rule as the reference chapter)
SEED = 42

# ----------------------------------------------------------------------------- style
COND_COL = {"Control": "#0072B2", "Mutant": "#D55E00", "AcidStress": "#009E73", "HeatShock": "#CC79A7"}
COND_LAB = {"Control": "Wild type", "Mutant": "\u0394hfq", "AcidStress": "Acid stress", "HeatShock": "Heat shock"}
UP, DOWN, NS = "#C0392B", "#2E6DA4", "#BDBDBD"

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 7.5, "axes.labelsize": 8, "axes.titlesize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 6.8,
    "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "xtick.major.size": 3, "ytick.major.size": 3,
    "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "figure.dpi": 150, "savefig.dpi": 600, "savefig.bbox": "tight", "savefig.pad_inches": 0.05,
})


def panel(ax, letter, dx=-0.16, dy=1.04):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="left")


DEMO_NOTE = "Demonstration data \u2013 provenance unverified; not a biological result"   # set to None for verified data


def save(fig, name, note=True, text=None):
    """Save PNG (600 dpi) + vector PDF; adds a small provenance footer unless note=False."""
    msg = text or DEMO_NOTE
    if note and msg:
        fig.canvas.draw()
        bb = fig.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans + fig.transFigure.inverted())
        fig.text(bb.x0 + 0.004, bb.y0 - 0.006, msg, fontsize=5.5, color="#8A3B12", ha="left", va="top")
    fig.savefig(FIG / f"{name}.png")
    fig.savefig(FIG / f"{name}.pdf")
    plt.close(fig)
    print(f"  saved {name}")


# ----------------------------------------------------------------------------- data
def validate_inputs(counts, meta, edges):
    if (counts.empty or not counts.index.is_unique or not counts.columns.is_unique
            or counts.index.isna().any() or counts.columns.isna().any()):
        raise ValueError("Count matrix must contain data with unique gene and sample identifiers.")
    if not np.issubdtype(counts.to_numpy().dtype, np.number):
        raise ValueError("Count matrix must contain numeric raw counts.")
    values = counts.to_numpy(dtype=float)
    if not np.isfinite(values).all() or (values < 0).any() or not np.equal(values, np.floor(values)).all():
        raise ValueError("Count matrix must contain finite, non-negative integer raw counts.")
    if not meta.index.is_unique or set(meta.index) != set(counts.columns):
        raise ValueError("Sample metadata IDs must uniquely match the count-matrix columns.")
    if "Condition" not in meta.columns or meta["Condition"].isna().any():
        raise ValueError("Sample metadata must include a non-empty Condition for every sample.")
    unknown_conditions = set(meta["Condition"]) - set(COND_COL)
    if unknown_conditions:
        raise ValueError(f"Unknown sample conditions: {sorted(unknown_conditions)}")
    if not {"source", "target", "score"}.issubset(edges.columns):
        raise ValueError("STRING edge table must include source, target, and score columns.")
    if edges.empty:
        raise ValueError("STRING edge table must contain at least one interaction.")
    edges["score"] = pd.to_numeric(edges["score"], errors="coerce")
    if edges[["source", "target", "score"]].isna().any().any():
        raise ValueError("STRING edges must have non-empty source, target, and numeric score values.")
    if (edges["source"] == edges["target"]).any():
        raise ValueError("STRING edge table must not contain self-interactions.")
    if not np.isfinite(edges["score"]).all() or not edges["score"].between(0, 1, inclusive="right").all():
        raise ValueError("STRING edge scores must be greater than 0 and no greater than 1.")
    meta = meta.loc[counts.columns]
    counts = counts.astype(np.int64)
    return counts, meta, edges


def load():
    counts = pd.read_csv(DATA / "count_matrix.csv", index_col="GeneID")
    meta = pd.read_csv(DATA / "sample_metadata.csv", index_col="SampleID")
    edges = pd.read_csv(DATA / "string_ppi_network_edges.csv")
    return validate_inputs(counts, meta, edges)


def run_contrast(counts, meta, group, ref="Control"):
    """DESeq2 Wald test for `group` vs `ref` (PyDESeq2, same model as DESeq2: NB GLM, BH-adjusted p)."""
    sub = meta[meta.Condition.isin([ref, group])]
    if group == ref or set(sub.Condition) != {ref, group}:
        raise ValueError(f"Both distinct conditions must be present for contrast {group!r} vs {ref!r}.")
    dds = DeseqDataSet(counts=counts[sub.index].T.astype(int), metadata=sub,
                       design="~Condition", quiet=True)
    dds.deseq2()
    stats = DeseqStats(dds, contrast=["Condition", group, ref], quiet=True)
    stats.summary()
    res = stats.results_df.copy()
    res.index.name = "GeneID"
    return res, dds


def classify(res):
    sig = (res.padj < PADJ_CUTOFF) & (res.log2FoldChange.abs() >= LFC_CUTOFF)
    res["status"] = np.where(sig & (res.log2FoldChange > 0), "Up",
                             np.where(sig & (res.log2FoldChange < 0), "Down", "NS"))
    return res


def zscore_rows(df):
    return df.sub(df.mean(axis=1), axis=0).div(df.std(axis=1, ddof=0).replace(0, 1), axis=0)


# ----------------------------------------------------------------------------- Fig 1: workflow
def fig_workflow():
    bands = [
        ("Example upstream\n(not established or run)", ["FASTQ", "QC / trimming", "Alignment / counting"], True),
        ("Input", ["count_matrix.csv", "sample_metadata.csv"], False),
        ("Quality control", ["Library sizes", "Count distributions", "PCA / sample\ndistances"], False),
        ("Differential\nexpression", ["Median-of-ratios\nnormalisation", "Dispersion\nshrinkage", "Wald test +\nBH FDR"], False),
        ("Visualisation", ["Volcano plot", "MA plot", "Expression heatmap"], False),
        ("Functional\nenrichment", ["Curated functional\ngene sets", "Fisher exact + BH"], False),
        ("Network &\nhub genes", ["STRING PPI edges", "Degree, betweenness,\ncloseness", "Consensus hub rank"], False),
        ("Structure", ["Hfq 3D structure\n(PDB)"], False),
    ]
    palette = ["#EEEEEE", "#FFF4D6", "#E3F0FA", "#E8F5E9", "#FDE8E0", "#F3E5F5", "#E0F2F1", "#ECEFF1"]
    fig, ax = plt.subplots(figsize=(7.2, 8.0))
    ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")
    n = len(bands); h = 100 / n
    centers = []
    for i, (lab, boxes, dashed) in enumerate(bands):
        y0 = 100 - (i + 1) * h
        ax.add_patch(plt.Rectangle((0, y0), 100, h, color=palette[i], lw=0, zorder=0))
        ax.text(2, y0 + h / 2, lab, fontsize=9, va="center", ha="left", fontweight="bold", color="#333")
        k = len(boxes); x_start, x_end = 30, 98; w_tot = x_end - x_start
        bw = min(w_tot / k - 3, 21)
        row = []
        for j, b in enumerate(boxes):
            cx = x_start + (j + 0.5) * w_tot / k
            ax.add_patch(FancyBboxPatch((cx - bw / 2, y0 + h * 0.2), bw, h * 0.6,
                                        boxstyle="round,pad=0.25,rounding_size=1.2",
                                        fc="white", ec="#555", lw=0.9, ls="--" if dashed else "-", zorder=2))
            ax.text(cx, y0 + h / 2, b, fontsize=7.2, ha="center", va="center",
                    color="#777" if dashed else "#111", zorder=3)
            row.append(cx)
        centers.append((y0, row))
    for i in range(n - 1):
        y_from = centers[i][0] + h * 0.2 - 0.6
        y_to = centers[i + 1][0] + h * 0.8 + 0.6
        ax.annotate("", xy=(64, y_to), xytext=(64, y_from),
                    arrowprops=dict(arrowstyle="-|>", lw=1.2, color="#444"), zorder=1)
    save(fig, "Fig1_workflow_overview", note=False)


# ----------------------------------------------------------------------------- Fig 2: QC
def fig_qc(counts, meta):
    order = meta.index.tolist()
    cols = [COND_COL[c] for c in meta.Condition]
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.6), gridspec_kw={"wspace": 0.55})
    ax = axs[0]
    ax.bar(range(len(order)), counts.sum().values / 1e3, color=cols, edgecolor="k", lw=0.5)
    ax.set_xticks(range(len(order))); ax.set_xticklabels(order, rotation=90)
    ax.set_ylabel("Total counts (\u00d710\u00b3)"); panel(ax, "a")
    ax = axs[1]
    lg = np.log2(counts + 1)
    bp = ax.boxplot([lg[s] for s in order], patch_artist=True, widths=0.65, showfliers=False,
                    medianprops=dict(color="k", lw=1), whiskerprops=dict(lw=0.7), capprops=dict(lw=0.7),
                    boxprops=dict(lw=0.7))
    for p, c in zip(bp["boxes"], cols):
        p.set_facecolor(c); p.set_alpha(0.85)
    ax.set_xticks(range(1, len(order) + 1)); ax.set_xticklabels(order, rotation=90)
    ax.set_ylabel("log$_2$(count + 1)"); panel(ax, "b")
    ax = axs[2]
    mu, var = counts.mean(axis=1), counts.var(axis=1)
    ax.scatter(mu, var, s=14, c="#444", alpha=0.8, lw=0)
    lim = [mu.min() * 0.6, mu.max() * 1.6]
    ax.plot(lim, lim, ls="--", c="#C0392B", lw=1, label="Poisson (var = mean)")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Mean count per gene"); ax.set_ylabel("Variance")
    ax.xaxis.set_major_locator(mpl.ticker.LogLocator(subs=(1.0, 3.0))); ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter()); ax.legend(frameon=False, loc="lower right", fontsize=6)
    panel(ax, "c")
    handles = [Patch(fc=COND_COL[k], ec="k", lw=0.5, label=COND_LAB[k]) for k in COND_COL]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.30))
    save(fig, "Fig2_quality_control")


# ----------------------------------------------------------------------------- Fig 3: PCA + distances
def fig_pca_distance(norm, meta):
    lg = np.log2(norm + 1)
    pca = PCA(n_components=2).fit(lg.T.values)
    sc = pca.transform(lg.T.values)
    ev = pca.explained_variance_ratio_ * 100
    fig = plt.figure(figsize=(7.2, 3.5))
    ax = fig.add_axes([0.07, 0.14, 0.36, 0.78]); panel(ax, "a", -0.17)
    for cnd in COND_COL:
        idx = np.where(meta.Condition.values == cnd)[0]
        ax.scatter(sc[idx, 0], sc[idx, 1], s=55, c=COND_COL[cnd], ec="k", lw=0.6, label=COND_LAB[cnd], zorder=3)
    texts = [ax.text(sc[i, 0], sc[i, 1], s, fontsize=6.3) for i, s in enumerate(meta.index)]
    adjust_text(texts, ax=ax, expand=(1.3, 1.6), arrowprops=dict(arrowstyle="-", lw=0.4, color="#666"))
    ax.axhline(0, c="#ddd", lw=0.6, zorder=0); ax.axvline(0, c="#ddd", lw=0.6, zorder=0)
    ax.set_xlabel(f"PC1: {ev[0]:.0f}% variance"); ax.set_ylabel(f"PC2: {ev[1]:.0f}% variance")
    ax.legend(frameon=False, loc="best", title="Condition", title_fontsize=7)
    # distance heatmap with dendrogram
    d = pdist(lg.T.values, "euclidean"); Z = linkage(d, "average")
    order = leaves_list(Z); D = squareform(d)[np.ix_(order, order)]
    names = meta.index[order]
    axd = fig.add_axes([0.54, 0.76, 0.30, 0.16]); axd.axis("off")
    dendrogram(Z, ax=axd, no_labels=True, color_threshold=0, above_threshold_color="#555", link_color_func=lambda k: "#555")
    axh = fig.add_axes([0.54, 0.14, 0.30, 0.60])
    cax = fig.add_axes([0.89, 0.14, 0.015, 0.60])
    im = axh.imshow(D, cmap="Blues_r", aspect="auto")
    axh.set_xticks(range(len(names))); axh.set_yticks(range(len(names)))
    axh.set_xticklabels(names, rotation=90); axh.set_yticklabels(names)
    for lab_ in axh.get_xticklabels() + axh.get_yticklabels():
        lab_.set_color(COND_COL[meta.loc[lab_.get_text(), "Condition"]])
    for s in axh.spines.values(): s.set_visible(False)
    axh.tick_params(length=0)
    cb = fig.colorbar(im, cax=cax); cb.set_label("Euclidean distance", fontsize=7); cb.outline.set_visible(False)
    axh.text(-0.22, 1.30, "b", transform=axh.transAxes, fontsize=11, fontweight="bold")
    save(fig, "Fig3_PCA_and_sample_distance")
    pd.DataFrame(sc, index=meta.index, columns=["PC1", "PC2"]).assign(Condition=meta.Condition).to_csv(TAB / "pca_coordinates.csv")


# ----------------------------------------------------------------------------- Fig 4: DESeq2 diagnostics
def fig_diagnostics(res, dds):
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.5), gridspec_kw={"wspace": 0.45})
    ax = axs[0]
    v = dds.var
    gw, fit, fin = v["genewise_dispersions"], v["fitted_dispersions"], v["dispersions"]
    mean = dds.layers["normed_counts"].mean(axis=0)
    ax.scatter(mean, gw, s=9, c="k", label="Gene-wise", lw=0)
    ax.scatter(mean, fit, s=9, c="#E74C3C", label="Fitted", lw=0)
    ax.scatter(mean, fin, s=9, c="#2E86DE", label="Final (MAP)", lw=0)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_ylim(1e-4, 1)  # genes with near-zero dispersion fall below the axis
    ax.legend(frameon=False, markerscale=1.4, loc="lower left")
    ax.set_xlabel("Mean of normalised counts"); ax.set_ylabel("Dispersion"); panel(ax, "a")
    ax = axs[1]
    ax.hist(res.pvalue.dropna(), bins=np.linspace(0, 1, 21), color="#7F8C8D", ec="white", lw=0.5)
    ax.set_xlabel("p-value"); ax.set_ylabel("Number of genes"); panel(ax, "b")
    ax = axs[2]
    r = res.dropna(subset=["padj"])
    colors = np.where(r.padj < PADJ_CUTOFF, UP, NS)
    ax.scatter(r.baseMean, r.log2FoldChange, c=colors, s=16, lw=0, alpha=0.9)
    ax.axhline(0, c="k", lw=0.6)
    ax.set_xscale("log"); ax.set_xlabel("Mean of normalised counts"); ax.set_ylabel("log$_2$ fold change")
    ax.legend(handles=[Line2D([], [], marker="o", ls="", c=UP, label=f"padj < {PADJ_CUTOFF}"),
                       Line2D([], [], marker="o", ls="", c=NS, label="n.s.")], frameon=False, loc="lower left")
    panel(ax, "c")
    save(fig, "Fig4_DESeq2_diagnostics")


# ----------------------------------------------------------------------------- Fig 5: heatmap + volcano
def fig_heatmap_volcano(res, norm, meta):
    sig = res[res.status != "NS"].sort_values("padj")
    fig = plt.figure(figsize=(7.2, 3.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.15], wspace=0.5)
    ax = fig.add_subplot(gs[0]); panel(ax, "a", -0.32)
    lg = np.log2(norm.loc[sig.index] + 1)
    z = zscore_rows(lg)
    col_order = []
    for cnd in COND_COL:
        col_order += meta.index[meta.Condition == cnd].tolist()
    if len(z) > 1:
        z = z.iloc[leaves_list(linkage(z.values, "average", metric="correlation"))]
    z = z[col_order]
    sns.heatmap(z, ax=ax, cmap="RdBu_r", center=0, vmin=-2.5, vmax=2.5, linewidths=0.4, linecolor="white",
                cbar_kws=dict(label="Row z-score of log$_2$(norm. count + 1)", shrink=0.7, pad=0.03),
                yticklabels=True, xticklabels=True)
    ax.set_xticklabels(col_order, rotation=90); ax.set_yticklabels(z.index, rotation=0, style="italic")
    ax.set_xlabel(""); ax.set_ylabel("")
    for i, s in enumerate(col_order):
        ax.add_patch(plt.Rectangle((i, -0.9), 1, 0.7, color=COND_COL[meta.loc[s, "Condition"]], clip_on=False, lw=0))
    ax.legend(handles=[Patch(fc=COND_COL[k], ec="none", label=COND_LAB[k]) for k in COND_COL], ncol=4, frameon=False,
              loc="lower center", bbox_to_anchor=(0.45, 1.16), fontsize=6, columnspacing=0.8, handlelength=1.0, handletextpad=0.3)
    ax = fig.add_subplot(gs[1]); panel(ax, "b", -0.18)
    r = res.dropna(subset=["padj"]).copy()
    floor = r.padj[r.padj > 0].min() / 10 if (r.padj == 0).any() else 0
    r["nlp"] = -np.log10(r.padj.clip(lower=max(floor, 1e-300)))
    for st, c in [("NS", NS), ("Down", DOWN), ("Up", UP)]:
        s = r[r.status == st]
        ax.scatter(s.log2FoldChange, s.nlp, s=22 if st != "NS" else 14, c=c, ec="k" if st != "NS" else "none",
                   lw=0.35, alpha=0.95 if st != "NS" else 0.75, zorder=3 if st != "NS" else 2)
    ax.axhline(-np.log10(PADJ_CUTOFF), ls="--", c="#777", lw=0.7)
    for x in (-LFC_CUTOFF, LFC_CUTOFF): ax.axvline(x, ls="--", c="#777", lw=0.7)
    top = r[r.status != "NS"].sort_values("padj").head(10)
    texts = [ax.text(row.log2FoldChange, row.nlp, g, fontsize=6.8, style="italic") for g, row in top.iterrows()]
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", lw=0.4, color="#555"),
                expand=(1.4, 1.8), force_text=(0.4, 0.6))
    ax.set_xlabel(r"log$_2$ fold change ($\Delta$$\it{hfq}$ / wild type)"); ax.set_ylabel("\u2212log$_{10}$ adjusted p-value")
    ax.legend(handles=[Line2D([], [], marker="o", ls="", c=UP, mec="k", mew=0.4, label=f"Up ({(r.status=='Up').sum()})"),
                       Line2D([], [], marker="o", ls="", c=DOWN, mec="k", mew=0.4, label=f"Down ({(r.status=='Down').sum()})"),
                       Line2D([], [], marker="o", ls="", c=NS, label="Not significant")],
              frameon=False, loc="center right", bbox_to_anchor=(1.0, 0.55))
    save(fig, "Fig5_heatmap_and_volcano")


# ----------------------------------------------------------------------------- Fig 6: DEG counts
def fig_gene_counts(res, norm, meta):
    genes = res[res.status != "NS"].sort_values("padj").index[:6]
    n = len(genes); nc = 3; nr = int(np.ceil(n / nc))
    fig, axs = plt.subplots(nr, nc, figsize=(7.2, 2.2 * nr), gridspec_kw={"wspace": 0.5, "hspace": 0.7})
    axs = np.atleast_1d(axs).ravel()
    conds = list(COND_COL)
    rng = np.random.default_rng(SEED)
    for ax, g in zip(axs, genes):
        for i, cnd in enumerate(conds):
            v = norm.loc[g, meta.index[meta.Condition == cnd]].values
            ax.scatter(i + rng.uniform(-0.12, 0.12, len(v)), v, c=COND_COL[cnd], ec="k", lw=0.5, s=26, zorder=3)
            ax.hlines(v.mean(), i - 0.28, i + 0.28, color="k", lw=1.2, zorder=4)
        ax.set_xticks(range(len(conds))); ax.set_xticklabels([COND_LAB[c] for c in conds], rotation=45, ha="right")
        ax.set_title(g, style="italic", fontsize=8); ax.set_ylabel("Normalised count")
        ax.set_ylim(bottom=-0.04 * norm.loc[g].max())
    for ax in axs[n:]: ax.axis("off")
    save(fig, "Fig6_top_DEG_normalised_counts")


# ----------------------------------------------------------------------------- Fig 7: condition comparison
def fig_contrasts(results):
    lfc = pd.DataFrame({k: v.log2FoldChange for k, v in results.items()})
    padj = pd.DataFrame({k: v.padj for k, v in results.items()})
    order = lfc.index[leaves_list(linkage(lfc.fillna(0).values, "average"))]
    lfc, padj = lfc.loc[order], padj.loc[order]
    fig, ax = plt.subplots(figsize=(3.6, 6.4))
    im = ax.imshow(lfc.clip(-6, 6).values, cmap="RdBu_r", vmin=-6, vmax=6, aspect="auto")
    ax.set_xticks(range(lfc.shape[1])); ax.set_xticklabels([COND_LAB[c] + "\nvs wild type" for c in lfc.columns])
    ax.set_yticks(range(len(lfc))); ax.set_yticklabels(lfc.index, style="italic")
    for i in range(lfc.shape[0]):
        for j in range(lfc.shape[1]):
            p = padj.iloc[i, j]
            if pd.notna(p) and p < PADJ_CUTOFF:
                ax.text(j, i, "*" if p >= 0.001 else "**", ha="center", va="center", fontsize=9)
    for s in ax.spines.values(): s.set_visible(False)
    ax.tick_params(length=0)
    cb = fig.colorbar(im, ax=ax, shrink=0.35, pad=0.04); cb.set_label("log$_2$ fold change (clipped at \u00b16)"); cb.outline.set_visible(False)
    ax.set_title("* padj < 0.05    ** padj < 0.001", fontsize=7, pad=6)
    save(fig, "Fig7_log2FC_across_conditions")
    return lfc, padj


# ----------------------------------------------------------------------------- enrichment (curated sets)
GENE_SETS = {
    # Curated functional groups; membership rationale and references are in docs/METHODS.md.
    "SPI-1-associated invasion genes": ["hilA", "hilC", "hilD", "invA", "invF", "sipA", "sipB", "sipC", "sipD"],
    "SPI-1 transcriptional regulators": ["hilA", "hilC", "hilD", "invF"],
    "SPI-2 (intracellular survival) genes": ["ssrA", "ssrB", "ssaG", "ssaH", "ssaE", "sseB", "sseC", "sseD"],
    "T3SS effectors and translocon proteins": ["sipA", "sipB", "sipC", "sipD", "sopB", "sopE", "sseB", "sseC", "sseD"],
    "Two-component systems": ["phoP", "phoQ", "ssrA", "ssrB", "ompR", "envZ"],
    "DNA maintenance (gyrase and repair)": ["gyrA", "gyrB", "recA"],
    "Global regulators": ["hfq", "fnr", "slyA", "rpoS"],
}


def enrichment(res, background):
    deg = set(res.index[res.status != "NS"])
    bg = set(background)
    rows = []
    for name, genes in GENE_SETS.items():
        gs = set(genes) & bg
        a = len(deg & gs); b = len(deg - gs); c = len(gs - deg); d = len(bg - deg - gs)
        odds, p = fisher_exact([[a, b], [c, d]], alternative="greater")
        rows.append(dict(gene_set=name, set_size=len(gs), n_DEG_in_set=a, n_DEG=len(deg), background=len(bg),
                         gene_ratio=a / max(len(deg), 1), pvalue=p, DEGs=",".join(sorted(deg & gs))))
    df = pd.DataFrame(rows)
    df["padj"] = adjust_pvalues_bh(df.pvalue)
    return df.sort_values("pvalue").reset_index(drop=True)


def adjust_pvalues_bh(pvalues):
    """Benjamini-Hochberg adjusted p-values, preserving NaNs and tied p-values."""
    values = pvalues.to_numpy(dtype=float)
    valid = np.isfinite(values)
    adjusted = np.full(values.shape, np.nan)
    valid_values = values[valid]
    if valid_values.size == 0:
        return pd.Series(adjusted, index=pvalues.index)
    order = np.argsort(valid_values, kind="stable")
    ranked = valid_values[order] * valid_values.size / np.arange(1, valid_values.size + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1].clip(max=1)
    sorted_adjusted = np.empty(valid_values.size)
    sorted_adjusted[order] = ranked
    adjusted[valid] = sorted_adjusted
    return pd.Series(adjusted, index=pvalues.index)


def fig_enrichment(enr):
    e = enr[enr.n_DEG_in_set > 0].sort_values("padj", ascending=False)
    fig, ax = plt.subplots(figsize=(5.2, 0.55 * len(e) + 1.1))
    colors = [UP if p < PADJ_CUTOFF else NS for p in e.padj]
    ax.barh(e.gene_set, -np.log10(e.padj), color=colors, ec="k", lw=0.5, height=0.62)
    ax.axvline(-np.log10(PADJ_CUTOFF), ls="--", c="#555", lw=0.8)
    for i, (_, r) in enumerate(e.iterrows()):
        ax.text(-np.log10(r.padj) + 0.04, i, f"{r.n_DEG_in_set}/{r.set_size} genes", va="center", fontsize=6.8)
    ax.set_xlabel("\u2212log$_{10}$ adjusted p-value (Fisher exact, BH)")
    ax.set_xlim(0, max(-np.log10(e.padj).max() * 1.5, 2.5))
    save(fig, "Fig8_functional_enrichment")


# ----------------------------------------------------------------------------- network / hubs
def hub_analysis(edges, res):
    G = nx.Graph()
    for _, r in edges.iterrows():
        G.add_edge(r.source, r.target, weight=r.score, dist=1 / r.score)
    deg = dict(G.degree()); strength = dict(G.degree(weight="weight"))
    btw = nx.betweenness_centrality(G, weight="dist", normalized=True)
    clo = nx.closeness_centrality(G, distance="dist")
    df = pd.DataFrame({"degree": deg, "weighted_degree": strength, "betweenness": btw, "closeness": clo})
    df["mean_rank"] = df.rank(ascending=False, method="average").mean(axis=1)
    df = df.join(res[["log2FoldChange", "padj", "status"]]).sort_values("mean_rank")
    df["is_DEG"] = df.status.isin(["Up", "Down"])
    df["hub"] = False; df.iloc[:5, df.columns.get_loc("hub")] = True
    df.index.name = "GeneID"
    return G, df


# ----------------------------------------------------------------------------- data diagnostics
def data_diagnostics(counts, meta):
    """Describe the supplied count matrix without treating it as genome-wide validation."""
    wt = counts[meta.index[meta.Condition == "Control"]]
    mu, var = wt.mean(axis=1), wt.var(axis=1, ddof=1)
    phi = ((var - mu) / mu ** 2).clip(lower=0)
    d = {
        "n_genes": counts.shape[0], "n_samples": counts.shape[1],
        "median_mean_count_WT": float(mu.median()),
        "dynamic_range_fold_(max/min_WT_mean)": float(mu.max() / max(mu.min(), 1e-9)),
        "median_variance_to_mean_WT": float((var / mu).median()),
        "median_NB_dispersion_moment_WT": float(phi.median()),
    }
    out = pd.DataFrame(d, index=["value"]).T
    out.to_csv(TAB / "data_diagnostics.csv")
    note = "descriptive statistics for the supplied panel only; raw-data provenance is unverified"
    print("  data diagnostics:", {k: round(v, 3) for k, v in d.items()}, "->", note)
    return d


# ----------------------------------------------------------------------------- Hfq structure (only if a validated file exists)
# ----------------------------------------------------------------------------- graphical abstract

# ----------------------------------------------------------------------------- Fig 11/12: annotated volcano (interpretation chart)
GUIDE_NEG = "Negative change in gene\nexpression compared to control"
GUIDE_POS = "Positive change in gene\nexpression compared to control"


def _guide_title(fig, text="Volcano Plot Interpretation"):
    fig.add_artist(plt.Rectangle((0.012, 0.915), 0.007, 0.065, transform=fig.transFigure, color="black", lw=0))
    fig.text(0.03, 0.947, text, fontsize=13, fontweight="bold", va="center", ha="left")


def _volcano_guide(fig, ax, x, y, status, xlim, ylim, thr_y, hi_y, thr_label, hi_label,
                   xlabel, ylabel, text_xy, text_ha, vlines=(), labels=None, symlog=False):
    """Draw an annotated volcano in the style of a teaching 'interpretation' chart."""
    ax.set_xlim(*xlim)
    if symlog:
        ax.set_yscale("symlog", linthresh=3, linscale=1.2)
        ax.set_yticks([0, 1, 2, 3, 10, 30, 100])
        ax.get_yaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    ax.set_ylim(*ylim)
    ax.axhspan(0, thr_y, color="#F2F2F2", lw=0, zorder=0)
    for st, c, z in (("NS", NS, 2), ("Down", DOWN, 3), ("Up", UP, 3)):
        m = np.asarray(status) == st
        ax.scatter(np.asarray(x)[m], np.asarray(y)[m], s=16 if len(x) > 100 else 34, c=c, lw=0.3 if len(x) <= 100 else 0,
                   ec="k" if len(x) <= 100 and st != "NS" else "none", alpha=0.9, zorder=z)
    ax.axvline(0, color="#777", lw=1.0, zorder=1)
    for v in vlines:
        ax.axvline(v, color="#999", lw=0.7, ls=(0, (3, 3)), zorder=1)
    ax.axhline(thr_y, color="#888", lw=0.8, ls=(0, (3, 3)), zorder=1)
    ax.axhline(hi_y, color="#888", lw=0.8, ls=(0, (3, 3)), zorder=1)
    xr = xlim[1] - xlim[0]
    ax.text(xlim[1] - 0.01 * xr, thr_y, thr_label, fontsize=6.3, color="#555", ha="right", va="bottom")
    ax.text(xlim[1] - 0.01 * xr, hi_y, hi_label, fontsize=6.3, color="#555", ha="right", va="bottom")
    # explanatory region notes
    (tx, ty_ins, ty_sig, ty_hi) = text_xy
    ax.text(tx, ty_hi, "Most significant change\nin differential gene\nexpression", fontsize=6.5, ha=text_ha, va="center")
    ax.text(tx, ty_sig, "Statistically significant\nchange in differential\ngene expression", fontsize=6.5, ha=text_ha, va="center")
    ax.text(tx, ty_ins, "Statistically insignificant\nchange in differential\ngene expression", fontsize=6.5, ha=text_ha, va="center")
    # zero-point marker
    tr = ax.get_xaxis_transform()
    ax.text(0, 0.965, "Zero point", transform=tr, ha="center", va="center", fontsize=6.8,
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#999", lw=0.6), zorder=5)
    for sgn in (-1, 1):
        ax.annotate("", xy=(sgn * 0.30 * xr / 2, 0.965), xytext=(sgn * 0.07 * xr / 2, 0.965), xycoords=tr, textcoords=tr,
                    arrowprops=dict(arrowstyle="-|>", color="#999", lw=1.0), zorder=4)
    # column headers
    ax.text(xlim[0] + 0.25 * xr, 1.04, GUIDE_NEG, transform=tr, ha="center", va="bottom", fontsize=7.2,
            color=DOWN, fontweight="bold", linespacing=1.25)
    ax.text(xlim[0] + 0.75 * xr, 1.04, GUIDE_POS, transform=tr, ha="center", va="bottom", fontsize=7.2,
            color=UP, fontweight="bold", linespacing=1.25)
    ax.legend(handles=[Line2D([], [], marker="o", ls="", c=DOWN, ms=5, label="Negative"),
                       Line2D([], [], marker="o", ls="", c=UP, ms=5, label="Positive"),
                       Line2D([], [], marker="o", ls="", c=NS, ms=5, label="Insignificant")],
              frameon=False, loc="upper right", bbox_to_anchor=(1.0, 0.93), handletextpad=0.2, borderaxespad=0.2)
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    if labels is not None and len(labels):
        texts = [ax.text(xv, yv, g, fontsize=6.6, style="italic", zorder=6) for g, xv, yv in labels]
        adjust_text(texts, ax=ax, expand=(1.8, 2.4), force_text=(0.8, 1.2), force_static=(0.5, 0.8),
                    arrowprops=dict(arrowstyle="-", lw=0.4, color="#555"))
    _guide_title(fig)


def fig_volcano_interpretation(res):
    """Annotated volcano of the actual mutant-vs-control contrast."""
    r = res.dropna(subset=["padj"]).copy()
    r["nlp"] = -np.log10(r.padj.clip(lower=1e-300))
    fig = plt.figure(figsize=(7.2, 4.3))
    ax = fig.add_axes([0.09, 0.14, 0.88, 0.62])
    sig = r[r.status != "NS"]
    labels = [(g, row.log2FoldChange, row.nlp) for g, row in sig.iterrows()]
    _volcano_guide(fig, ax, r.log2FoldChange, r.nlp, r.status, xlim=(-11, 11), ylim=(0, 150),
                   thr_y=-np.log10(PADJ_CUTOFF), hi_y=3.0,
                   thr_label="Threshold (adj. p = 0.05)", hi_label="adj. p = 0.001",
                   xlabel=r"log$_2$ fold change ($\Delta$$\it{hfq}$ / wild type)",
                   ylabel="\u2212log$_{10}$ (adjusted p-value)",
                   text_xy=(2.0, 0.55, 2.1, 22), text_ha="left", vlines=(-LFC_CUTOFF, LFC_CUTOFF),
                   labels=labels, symlog=True)
    save(fig, "Fig11_volcano_interpretation")


def fig_volcano_guide():
    """Schematic reading guide: SIMULATED points, not study data."""
    rng = np.random.default_rng(SEED)
    n = 420
    x = rng.normal(0, 0.45, n)
    y = np.minimum(np.abs(x) ** 0.9 * rng.gamma(3.0, 1.9, n) + rng.uniform(0, 0.4, n), 6.8)
    status = np.where((y > 2) & (x < -0.12), "Down", np.where((y > 2) & (x > 0.12), "Up", "NS"))
    fig = plt.figure(figsize=(7.2, 4.3))
    ax = fig.add_axes([0.09, 0.14, 0.88, 0.62])
    ax.set_xticks(np.arange(-2.5, 2.6, 0.5))
    _volcano_guide(fig, ax, x, y, status, xlim=(-2.5, 2.5), ylim=(0, 7.2), thr_y=2.0, hi_y=4.9,
                   thr_label="Threshold", hi_label="", xlabel="Fold change (log$_2$)",
                   ylabel="\u2212log$_{10}$ (P-value)", text_xy=(-2.4, 1.0, 3.4, 6.0), text_ha="left")
    fig.text(0.97, 0.025, "Schematic with simulated points \u2013 not study data", fontsize=6, color="#555", ha="right")
    save(fig, "Fig12_volcano_reading_guide", note=False)


# ----------------------------------------------------------------------------- Fig 0: compact graphical summary
def fig_summary(results, enr):
    res = results["Mutant"].dropna(subset=["padj"])
    fig = plt.figure(figsize=(7.4, 3.6), layout="constrained")
    grid = fig.add_gridspec(2, 3, height_ratios=[0.16, 1], width_ratios=[0.95, 1.05, 1.5])
    t = fig.add_subplot(grid[0, :]); t.set_axis_off()
    t.text(0.5, 0.80, "Salmonella RNA-seq workflow output", ha="center", va="center", fontsize=11,
           fontweight="bold", color="#17324D")
    t.text(0.5, 0.12, "DEMONSTRATION ONLY  |  count-matrix provenance and sample labels are unverified", ha="center",
           va="center", fontsize=6.8, fontweight="bold", color="#8A3B12",
           bbox=dict(boxstyle="round,pad=0.3", fc="#FFF1D6", ec="#E7B45A", lw=0.6))
    ax = fig.add_subplot(grid[1, 0]); panel(ax, "a", -0.2, 1.03)
    y = -np.log10(res.padj.clip(lower=1e-300))
    for st, c, lab in (("NS", NS, "Not significant"), ("Down", DOWN, "Down"), ("Up", UP, "Up")):
        m = res.status == st
        ax.scatter(res.loc[m, "log2FoldChange"], y[m], c=c, s=18, lw=0, label=lab, zorder=3 if st != "NS" else 2)
    ax.axhline(-np.log10(PADJ_CUTOFF), color="#555", ls="--", lw=0.7)
    for v in (-LFC_CUTOFF, LFC_CUTOFF):
        ax.axvline(v, color="#555", ls="--", lw=0.7)
    sig = res[res.status != "NS"]
    texts = [ax.text(row.log2FoldChange, -np.log10(max(row.padj, 1e-300)), g, fontsize=6.3, style="italic")
             for g, row in sig.iterrows()]
    if texts:
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="#555", lw=0.4))
    ax.set_xlabel("log$_2$ fold change"); ax.set_ylabel("\u2212log$_{10}$ adjusted p-value")
    ax.set_title("Mutant vs control", loc="left", fontsize=7.5, fontweight="bold")
    ax.legend(frameon=False, fontsize=5.8, loc="center", bbox_to_anchor=(0.62, 0.42), handletextpad=0.1)

    ax = fig.add_subplot(grid[1, 1]); panel(ax, "b", -0.2, 1.03)
    spi1 = ["hilA", "hilC", "hilD", "invF", "sipA"]
    for off, (key, lab) in zip((-0.2, 0.0, 0.2), (("Mutant", "mutant"), ("AcidStress", "acid"), ("HeatShock", "heat"))):
        c = results[key].loc[spi1]
        ax.errorbar(c.log2FoldChange, np.arange(len(spi1)) + off, xerr=1.96 * c.lfcSE, fmt="o", ms=3.6, capsize=1.8,
                    elinewidth=0.9, color=COND_COL[key], label=f"{lab} vs control")
    ax.axvline(0, color="#555", lw=0.7)
    ax.set_yticks(np.arange(len(spi1)), labels=spi1, style="italic"); ax.invert_yaxis()
    ax.set_xlabel("log$_2$ fold change (\u00b1 1.96 SE)")
    ax.set_title("Selected invasion-associated genes", loc="left", fontsize=7.5, fontweight="bold")
    ax.legend(frameon=False, fontsize=6, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=3, columnspacing=0.8,
              handletextpad=0.2)

    ax = fig.add_subplot(grid[1, 2]); panel(ax, "c", -0.12, 1.03)
    sizes = enr.set_size.replace(0, np.nan)
    frac = (enr.n_DEG_in_set / sizes).fillna(0)
    order = np.argsort(frac.to_numpy(), kind="stable")
    pl, pf = enr.iloc[order], frac.iloc[order]
    yy = np.arange(len(pl))
    ax.barh(yy, pf, color=[UP if p < PADJ_CUTOFF else "#8A9BA8" for p in pl.padj], height=0.68, ec="white", lw=0.5)
    short = {"SPI-1 transcriptional regulators": "SPI-1 regulators", "SPI-1-associated invasion genes": "SPI-1 invasion genes",
             "T3SS effectors and translocon proteins": "T3SS effectors/translocon",
             "SPI-2 (intracellular survival) genes": "SPI-2 genes", "DNA maintenance (gyrase and repair)": "DNA maintenance"}
    ax.set_yticks(yy, labels=[short.get(g, g) for g in pl.gene_set], fontsize=6)
    for i, (_, row) in enumerate(pl.iterrows()):
        ax.text(min(pf.iloc[i] + 0.02, 0.86), i, f"{int(row.n_DEG_in_set)}/{int(row.set_size)}", va="center", fontsize=5.8)
    ax.set_xlim(0, 1.05); ax.set_xlabel("Fraction of set that are DEGs")
    ax.set_title("Curated-set overlap", loc="left", fontsize=7.5, fontweight="bold")
    fig.text(0.5, -0.025, "Stress groups n = 2 each (exploratory). Enrichment uses a 30-gene background; not genome-wide GO.",
             ha="center", va="top", fontsize=5.8, color="#555")
    save(fig, "Fig0_graphical_summary", note=False)


# ----------------------------------------------------------------------------- Fig 9: network (readable layout)
def fig_network(G, hubs):
    fig = plt.figure(figsize=(7.2, 4.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.75, 1], wspace=0.28)
    ax = fig.add_subplot(gs[0]); panel(ax, "a", -0.02, 1.0); ax.axis("off")
    pos = nx.kamada_kawai_layout(G, weight="dist")
    strong = [(u, v) for u, v in G.edges() if G[u][v]["weight"] >= 0.7]
    weak = [(u, v) for u, v in G.edges() if G[u][v]["weight"] < 0.7]
    nx.draw_networkx_edges(G, pos, ax=ax, edgelist=weak, width=0.3, edge_color="#BDBDBD", alpha=0.35, style="dotted")
    nx.draw_networkx_edges(G, pos, ax=ax, edgelist=strong,
                           width=[0.4 + 1.6 * (G[u][v]["weight"] - 0.7) / 0.3 for u, v in strong],
                           edge_color="#8E8E8E", alpha=0.5)
    meas = [n for n in G.nodes if pd.notna(hubs.loc[n, "log2FoldChange"])]
    unmeas = [n for n in G.nodes if n not in meas]
    size = lambda ns: [200 + 22 * hubs.loc[n, "degree"] for n in ns]
    ring = lambda ns: ["#000" if hubs.loc[n, "hub"] else "#777" for n in ns]
    lw = lambda ns: [2.0 if hubs.loc[n, "hub"] else 0.6 for n in ns]
    if unmeas:
        nx.draw_networkx_nodes(G, pos, ax=ax, nodelist=unmeas, node_size=size(unmeas), node_color="#E8E8E8",
                               node_shape="o", edgecolors=ring(unmeas), linewidths=lw(unmeas))
    nodes = None
    if meas:
        nodes = nx.draw_networkx_nodes(G, pos, ax=ax, nodelist=meas, node_size=size(meas),
                                       node_color=hubs.loc[meas, "log2FoldChange"].clip(-6, 6).values,
                                       cmap="RdBu_r", vmin=-6, vmax=6, edgecolors=ring(meas), linewidths=lw(meas))
    for n_, (x_, y_) in pos.items():
        dark = n_ in meas and abs(hubs.loc[n_, "log2FoldChange"]) > 3.5
        ax.text(x_, y_, n_, fontsize=5.8, style="italic", ha="center", va="center", zorder=5,
                color="white" if dark else "black")
    ax.margins(0.08)
    if nodes is not None:
        cb = fig.colorbar(nodes, ax=ax, shrink=0.5, pad=0.0, aspect=18)
        cb.set_label("log$_2$FC (\u0394hfq / WT)", fontsize=6.5); cb.outline.set_visible(False)
    handles = [Line2D([], [], marker="o", ls="", mfc="w", mec="k", mew=2, ms=6, label="Top-5 hub"),
               Line2D([], [], marker="o", ls="", mfc="#E8E8E8", mec="#777", ms=6, label="Gene not in count panel"),
               Line2D([], [], color="#8E8E8E", lw=1.3, label="STRING score \u2265 0.70"),
               Line2D([], [], color="#BDBDBD", lw=0.8, ls=":", label="score 0.40\u20130.70")]
    ax.legend(handles=handles, frameon=False, loc="lower left", fontsize=6, bbox_to_anchor=(-0.02, -0.04))
    ax = fig.add_subplot(gs[1]); panel(ax, "b", -0.38, 1.0)
    top = hubs.head(8).iloc[::-1]
    ax.barh(np.arange(len(top)), top.weighted_degree, color=["#6C3483" if h else "#B0BEC5" for h in top.hub],
            ec="k", lw=0.5, height=0.65)
    ax.set_yticks(np.arange(len(top)), labels=top.index, style="italic")
    ax.set_xlabel("Weighted degree")
    ax.set_title("Top 8 by consensus rank", loc="left", fontsize=7.5)
    ax.legend(handles=[Patch(fc="#6C3483", ec="k", lw=0.5, label="Top-5 hub"), Patch(fc="#B0BEC5", ec="k", lw=0.5, label="Other")],
              frameon=False, loc="upper center", bbox_to_anchor=(0.45, -0.18), ncol=2, fontsize=6)
    save(fig, "Fig9_PPI_network_and_hubs")


# ----------------------------------------------------------------------------- Fig 10: Hfq structure (hexamer ring, two views)
def fig_structure():
    files = sorted((ROOT / "structure").glob("*.pdb"))
    if not files:
        print("  structure: no PDB in structure/ -> skipped (run scripts/fetch_hfq_structure.py)")
        return
    pdb = None
    lines = None
    for candidate in files:
        candidate_lines = candidate.read_text(encoding="ascii", errors="replace").splitlines()
        candidate_text = "\n".join(candidate_lines)
        if is_hfq_pdb(candidate_text):
            pdb, lines = candidate, candidate_lines
            break
        print(f"  structure: {candidate.name} does not identify Hfq "
              f"(header: '{pdb_header_text(candidate_text)[:70]}') -> skipped")
    if pdb is None:
        print("  structure: no validated Hfq PDB -> Fig10 skipped")
        return
    chains = {}
    organism = None
    for l in lines:
        if "ORGANISM_SCIENTIFIC:" in l and organism is None:
            organism = l.split("ORGANISM_SCIENTIFIC:")[1].strip(" ;\r").capitalize()
        if l.startswith("ENDMDL"):
            break                                             # first model only (NMR ensembles)
        if l.startswith("ATOM") and l[12:16].strip() == "CA":
            chains.setdefault(l[21], []).append((int(l[22:26]), float(l[30:38]), float(l[38:46]), float(l[46:54])))
    if not chains:
        print(f"  structure: {pdb.name} has no C-alpha coordinates -> Fig10 skipped")
        return
    allxyz = np.vstack([[r[1:] for r in rows] for rows in chains.values()])
    centre = allxyz.mean(axis=0)
    _, _, vt = np.linalg.svd(allxyz - centre, full_matrices=False)   # principal axes; last = ring axis for a toroid
    rot = vt[[0, 1, 2]].T
    palette = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#7F7F7F", "#000000"]
    fig = plt.figure(figsize=(7.2, 3.5))
    lim = np.abs((allxyz - centre) @ rot).max() * 0.62
    for k, (elev, azim, ttl, letter) in enumerate(((90, -90, "Top view (ring axis towards viewer)", "a"), (8, -90, "Side view", "b"))):
        ax = fig.add_subplot(1, 2, k + 1, projection="3d")
        for i, (ch, rows) in enumerate(sorted(chains.items())):
            a = (np.array([r[1:] for r in rows]) - centre) @ rot
            ax.plot(a[:, 0], a[:, 1], a[:, 2], color=palette[i % len(palette)], lw=1.7, solid_capstyle="round", label=f"Chain {ch}")
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, lim)
        ax.set_box_aspect((1, 1, 1)); ax.set_axis_off(); ax.view_init(elev=elev, azim=azim)
        ax.set_title(ttl, fontsize=7.5, pad=-4)
        ax.text2D(0.02, 0.95, letter, transform=ax.transAxes, fontsize=11, fontweight="bold")
    fig.legend(*fig.axes[0].get_legend_handles_labels(), loc="lower center", ncol=min(len(chains), 6), frameon=False,
               fontsize=6.5, bbox_to_anchor=(0.5, -0.03))
    src = f"{organism} " if organism else ""
    fig.suptitle(f"Hfq C\u03b1 backbone ({pdb.stem}; {src}homolog structure, not a Salmonella-specific model)",
                 fontsize=8, y=0.99)
    save(fig, "Fig10_Hfq_structure", note=False)


# ----------------------------------------------------------------------------- main
def main():
    print("Loading data ...")
    counts, meta, edges = load()
    print(f"  {counts.shape[0]} genes x {counts.shape[1]} samples; conditions: {meta.Condition.value_counts().to_dict()}")

    data_diagnostics(counts, meta)
    print("Differential expression (PyDESeq2) ...")
    results, dds_main = {}, None
    for grp in ["Mutant", "AcidStress", "HeatShock"]:
        r, d = run_contrast(counts, meta, grp)
        results[grp] = classify(r)
        if grp == "Mutant":
            dds_main = d
        results[grp].sort_values("padj").to_csv(TAB / f"de_{grp}_vs_Control.csv")
    res = results["Mutant"]

    # size-factor normalisation across all ten samples (used for QC, PCA, heatmaps)
    dds_all = DeseqDataSet(counts=counts.T.astype(int), metadata=meta, design="~Condition", quiet=True)
    dds_all.fit_size_factors()
    norm = pd.DataFrame(dds_all.layers["normed_counts"], index=counts.columns, columns=counts.index).T
    norm.to_csv(TAB / "normalized_counts.csv")

    sig = res[res.status != "NS"].sort_values("padj")
    sig.assign(Mean_WT=norm.loc[sig.index, meta.index[meta.Condition == "Control"]].mean(axis=1),
               Mean_Mutant=norm.loc[sig.index, meta.index[meta.Condition == "Mutant"]].mean(axis=1)
               ).to_csv(TAB / "significant_degs.csv")
    print(f"  significant DEGs (padj<{PADJ_CUTOFF}, |log2FC|>={LFC_CUTOFF}): "
          f"{(res.status=='Up').sum()} up, {(res.status=='Down').sum()} down -> {', '.join(sig.index)}")

    print("Figures ...")
    fig_workflow()
    fig_qc(counts, meta)
    fig_pca_distance(norm, meta)
    fig_diagnostics(res, dds_main)
    fig_heatmap_volcano(res, norm, meta)
    fig_volcano_interpretation(res)
    fig_volcano_guide()
    fig_gene_counts(res, norm, meta)
    fig_contrasts(results)

    print("Functional enrichment ...")
    enr = enrichment(res, counts.index)
    enr.to_csv(TAB / "functional_enrichment.csv", index=False)
    fig_enrichment(enr)
    fig_summary(results, enr)
    print(enr[["gene_set", "n_DEG_in_set", "set_size", "padj"]].round(4).to_string(index=False))

    print("Network hub analysis ...")
    G, hubs = hub_analysis(edges, res)
    hubs.to_csv(TAB / "hub_genes.csv")
    hubs.reset_index().to_csv(TAB / "cytoscape_node_table.csv", index=False)   # import in Cytoscape: File > Import > Table
    fig_network(G, hubs)
    top = hubs.head(5)
    print("  top hubs:", ", ".join(top.index), f"| {int(top.is_DEG.sum())}/5 are DEGs in this analysis")
    if top.is_DEG.sum() == 0:
        print("  WARNING: no hub gene is a DEG; the network does not match the current DEG list "
              "(run scripts/update_with_internet.sh).")
    fig_structure()
    print("Done.")


if __name__ == "__main__":
    main()
