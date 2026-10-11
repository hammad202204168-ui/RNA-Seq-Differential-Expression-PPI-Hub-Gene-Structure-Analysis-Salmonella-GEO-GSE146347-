#!/usr/bin/env python3
"""
Exploratory TPM-based analysis of the Salmonella-only table of GEO GSE146347
(Barrila et al. 2022, doi:10.3389/fcimb.2022.705647): S. Typhimurium wild type vs delta-hfq,
grown in a rotating-wall-vessel bioreactor in control or low-fluid-shear (LFS) orientation; n = 3 per group.

IMPORTANT: only TPM values are available here (no raw counts), so DESeq2/edgeR cannot be used. Statistics are
ordinary pooled-variance t-tests on log2(TPM + 1), Benjamini-Hochberg adjusted. Treat results as exploratory.

Run:  python scripts/prepare_geo_tpm.py && python scripts/real_data_analysis.py
"""
import warnings

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from adjustText import adjust_text
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Patch
from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage
from scipy.spatial.distance import pdist, squareform
from scipy.stats import fisher_exact, spearmanr, t as tdist
from sklearn.decomposition import PCA

import robustness as rb
import volcano_guide
from plotstyle import (DOWN, GROUP_COL, GROUP_LAB, GROUP_MARK, GROUP_ORDER, NS, ROOT, TAB, UP, panel, save)

warnings.filterwarnings("ignore")
PROC = ROOT / "data" / "processed"
FDR_CUT, LFC_CUT, MIN_TPM, MIN_SAMPLES = 0.05, 1.0, 1.0, 3

CONTRASTS = {   # name: (numerator group, denominator group)
    "hfq_effect_control": ("delta_hfq_control", "WT_control"),
    "hfq_effect_LFS": ("delta_hfq_LFS", "WT_LFS"),
    "LFS_effect_WT": ("WT_LFS", "WT_control"),
    "LFS_effect_hfq": ("delta_hfq_LFS", "delta_hfq_control"),
}
CONTRAST_LAB = {
    "hfq_effect_control": r"$\Delta$$\it{hfq}$ vs WT (control)", "hfq_effect_LFS": r"$\Delta$$\it{hfq}$ vs WT (LFS)",
    "LFS_effect_WT": "LFS vs control (WT)", "LFS_effect_hfq": r"LFS vs control ($\Delta$$\it{hfq}$)",
}
PRIMARY = "hfq_effect_control"


# ----------------------------------------------------------------------------- data + statistics
def load():
    tpm = pd.read_csv(PROC / "tpm_matrix.csv", index_col="GeneID")
    ann = pd.read_csv(PROC / "gene_annotation.csv", index_col="GeneID")
    meta = pd.read_csv(PROC / "sample_metadata.csv").set_index("sample").loc[tpm.columns]
    return tpm, ann, meta


STRUCTURAL = ("rRNA", "tRNA", "ncRNA")


def renormalise(tpm_raw, ann):
    """Drop rRNA/tRNA/ncRNA genes and rescale each sample to 1e6, so that differing rRNA-depletion efficiency
    between libraries does not distort the TPM of the remaining genes."""
    keep = ~ann.loc[tpm_raw.index, "type"].isin(STRUCTURAL)
    t = tpm_raw[keep]
    return t / t.sum() * 1e6


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n); out[o] = np.minimum(q, 1)
    return out


class Model:
    """One-way (4 cell) model on log2(TPM+1): group means, pooled within-group variance, df = N - 4."""

    def __init__(self, Y, meta):
        self.groups = GROUP_ORDER
        self.n = meta.group.value_counts().reindex(self.groups).values
        cols = {g: np.where(meta.group.values == g)[0] for g in self.groups}
        self.M = pd.DataFrame({g: Y.values[:, i].mean(axis=1) for g, i in cols.items()}, index=Y.index)
        ssw = sum(((Y.values[:, i] - self.M[g].values[:, None]) ** 2).sum(axis=1) for g, i in cols.items())
        self.df = Y.shape[1] - len(self.groups)
        self.s2 = np.maximum(ssw / self.df, 1e-8)

    def test(self, weights):
        c = pd.Series(weights, index=self.groups).reindex(self.groups).fillna(0.0).values
        est = self.M.values @ c
        se = np.sqrt(self.s2 * (c ** 2 / self.n).sum())
        t = est / se
        p = 2 * tdist.sf(np.abs(t), self.df)
        return pd.DataFrame({"log2FC": est, "se": se, "t": t, "pvalue": p, "FDR": bh(p)}, index=self.M.index)


def annotate(res, tpm, ann, meta, a=None, b=None):
    out = res.join(ann[["label", "symbol", "description", "type", "replicon"]])
    for g in GROUP_ORDER:
        out[f"meanTPM_{g}"] = tpm.loc[out.index, meta.index[meta.group == g]].mean(axis=1)
    sig = (out.FDR < FDR_CUT) & (out.log2FC.abs() >= LFC_CUT)
    out["status"] = np.where(sig & (out.log2FC > 0), "Up", np.where(sig & (out.log2FC < 0), "Down", "NS"))
    out.index.name = "GeneID"
    return out.sort_values("FDR")


# ----------------------------------------------------------------------------- gene sets (rule-based, from the supplied annotation)
def build_sets(ann):
    sym = ann.symbol.fillna(ann.name).astype(str)
    desc = ann.description.fillna("").astype(str)
    rx = lambda pat: sym.str.match(pat)
    kw = lambda pat: desc.str.contains(pat, case=False, regex=True)
    s = {
        "SPI-1 T3SS, invasion & effectors": rx(r"^(inv|spa|prg|org|sic|iag|sip|hil)[A-Z]$") | rx(r"^(sopB|sopE2?|sopA|sopD2?)$"),
        "SPI-2 T3SS & effectors": rx(r"^(ssa|sse|ssr|ssc)[A-Z]$") | rx(r"^(sifA|sifB|sseL|sspH2|spiC|pipB2?)$"),
        "Flagellar assembly & motility": rx(r"^(fli|flg|flh|mot)[A-Z]$") | kw("flagell"),
        "Chemotaxis": rx(r"^(che[A-Z]|tar|tsr|trg|tap|aer)$") | kw("chemotaxis|methyl-accepting"),
        "Ribosomal proteins": rx(r"^(rpl|rps|rpm)[A-Z]$") | kw("ribosomal protein"),
        "ABC transporters": kw("ABC transporter|ATP-binding component|permease component"),
        "Fimbriae & pili": kw("fimbri|pilus|pilin"),
        "Prophage genes": kw("phage"),
        "Transposases / IS elements": kw("transposase|insertion sequence"),
        "Respiration & ATP synthase": kw("ATP synthase|NADH|cytochrome|oxidase subunit|succinate dehydrogenase|fumarate reductase"),
        "PTS sugar uptake": kw("phosphotransferase|PTS system"),
        "Iron / siderophore / heme": kw("iron|siderophore|enterobactin|ferric|ferrous|heme|hemin"),
        "Transcriptional regulators": kw("transcriptional regulator|transcriptional repressor|transcriptional activator"),
        "Two-component systems": kw("response regulator|sensor histidine kinase|sensor kinase|two-component"),
        "Toxin-antitoxin systems": kw("toxin|antitoxin"),
        "tRNA genes": ann.type == "tRNA",
    }
    return {k: v.fillna(False) for k, v in s.items()}


def enrichment(deg_ids, background, sets, direction, min_size=5):
    rows = []
    deg = set(deg_ids) & set(background); bg = set(background)
    for name, mask in sets.items():
        gs = set(mask.index[mask]) & bg
        if len(gs) < min_size:
            continue
        a, b = len(deg & gs), len(deg - gs); c, d = len(gs - deg), len(bg - deg - gs)
        odds, p = fisher_exact([[a, b], [c, d]], alternative="greater")
        rows.append(dict(direction=direction, gene_set=name, set_size=len(gs), n_DEG_in_set=a, n_DEG=len(deg),
                         background=len(bg), odds_ratio=odds, pvalue=p))
    df = pd.DataFrame(rows)
    if len(df):
        df["FDR"] = bh(df.pvalue.values)
    return df


# ----------------------------------------------------------------------------- figures
def fig_workflow():
    bands = [
        ("Source data", ["GEO GSE146347\nSalmonella-only TPM table", "SRA run table\n(design check)"], False),
        ("Preparation", ["Keep TPM columns\n(12 bacterial samples)", "Gene annotation +\nsample groups"], False),
        ("Quality control", ["Detected genes,\ndistributions", "PCA and sample\ndistances", "hfq positive\ncontrol"], False),
        ("Statistics", ["Drop rRNA/tRNA/ncRNA;\nrenormalise TPM", "Filter TPM \u2265 1 in\n\u2265 3 samples; log2(TPM+1)", "Pooled-variance t +\nBH FDR (4 contrasts)"], False),
        ("Robustness\nchecks", ["Flag group-exclusive\ngenes", "Refit with technical\ncovariates", "High-confidence\nDEG set"], False),
        ("Visualisation", ["Volcano and\ninterpretation chart", "DEG heatmap", "Contrast\ncomparison"], False),
        ("Functional &\nnetwork", ["Rule-based gene-set\nenrichment", "Co-expression\nnetwork + hubs"], False),
        ("Not done\n(needs raw reads)", ["SRA download +\nalignment", "Raw gene counts", "DESeq2"], True),
    ]
    palette = ["#EEEEEE", "#FFF4D6", "#E3F0FA", "#E8F5E9", "#EFE6F7", "#FDE8E0", "#E0F2F1", "#F4F4F4"]
    fig, ax = plt.subplots(figsize=(7.2, 7.8)); ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")
    n = len(bands); h = 100 / n; ys = []
    for i, (lab, boxes, dashed) in enumerate(bands):
        y0 = 100 - (i + 1) * h; ys.append(y0)
        ax.add_patch(plt.Rectangle((0, y0), 100, h, color=palette[i], lw=0, zorder=0))
        ax.text(2, y0 + h / 2, lab, fontsize=9, va="center", fontweight="bold", color="#333")
        k = len(boxes); xs, xe = 30, 98; bw = min((xe - xs) / k - 2.5, 22)
        for j, b in enumerate(boxes):
            cx = xs + (j + 0.5) * (xe - xs) / k
            ax.add_patch(FancyBboxPatch((cx - bw / 2, y0 + h * 0.2), bw, h * 0.6, boxstyle="round,pad=0.25,rounding_size=1.2",
                                        fc="white", ec="#555", lw=0.9, ls="--" if dashed else "-", zorder=2))
            ax.text(cx, y0 + h / 2, b, fontsize=7, ha="center", va="center", color="#777" if dashed else "#111", zorder=3)
    for i in range(n - 2):
        ax.annotate("", xy=(64, ys[i + 1] + h * 0.8 + 0.6), xytext=(64, ys[i] + h * 0.2 - 0.6),
                    arrowprops=dict(arrowstyle="-|>", lw=1.2, color="#444"), zorder=1)
    save(fig, "Fig1_workflow_overview", note=False)


def fig_qc(tpm_raw, tpm, meta, ann, mean_log, kept):
    Y = np.log2(tpm + 1)
    order = meta.index.tolist(); cols = [GROUP_COL[g] for g in meta.group]
    sra = pd.read_csv(ROOT / "data" / "sra" / "SraRunTable_GSE146347.csv")
    sra = sra[sra.Organism.str.startswith("Salmonella")].copy()
    sra["group"] = sra.bacterial_genotype.map({"wild type": "WT", "delta_hfq": "delta_hfq"}) + "_" + sra.bacterial_growth_orientation
    fig, axs = plt.subplots(2, 3, figsize=(7.2, 5.6), gridspec_kw={"wspace": 0.5, "hspace": 0.9})
    ax = axs[0, 0]
    bp = ax.boxplot([Y[s_] for s_ in order], patch_artist=True, widths=0.65, showfliers=False,
                    medianprops=dict(color="k", lw=1), whiskerprops=dict(lw=0.7), capprops=dict(lw=0.7), boxprops=dict(lw=0.7))
    for p, c in zip(bp["boxes"], cols):
        p.set_facecolor(c); p.set_alpha(0.9)
    ax.set_xticks(range(1, 13)); ax.set_xticklabels(order, rotation=90, fontsize=6); ax.set_ylabel("log$_2$(TPM + 1)"); panel(ax, "a")
    ax = axs[0, 1]
    ax.bar(range(12), (tpm >= MIN_TPM).sum().values, color=cols, ec="k", lw=0.5)
    ax.set_xticks(range(12)); ax.set_xticklabels(order, rotation=90, fontsize=6); ax.set_ylabel(f"Genes with TPM \u2265 {MIN_TPM:g}"); panel(ax, "b")
    ax = axs[0, 2]
    ax.hist(mean_log, bins=60, color="#7F8C8D", ec="white", lw=0.3)
    ax.hist(mean_log[kept], bins=60, color="#2E4A6B", ec="white", lw=0.3, label=f"tested ({kept.sum():,})")
    ax.set_xlabel("Mean log$_2$(TPM + 1)"); ax.set_ylabel("Number of genes"); ax.legend(frameon=False, loc="upper right", fontsize=6); panel(ax, "c")
    ax = axs[1, 0]
    typ = ann.loc[tpm_raw.index, "type"]
    share = pd.DataFrame({"rRNA": tpm_raw[typ == "rRNA"].sum(), "tRNA + ncRNA": tpm_raw[typ.isin(["tRNA", "ncRNA"])].sum()}) / 1e4
    ax.bar(range(12), share["rRNA"], color="#8E44AD", ec="k", lw=0.4, label="rRNA")
    ax.bar(range(12), share["tRNA + ncRNA"], bottom=share["rRNA"], color="#BDC3C7", ec="k", lw=0.4, label="tRNA + ncRNA")
    ax.set_xticks(range(12)); ax.set_xticklabels(order, rotation=90, fontsize=6); ax.set_ylabel("% of original TPM")
    ax.legend(frameon=False, fontsize=6, loc="upper left"); panel(ax, "d")
    ax = axs[1, 1]
    rng = np.random.default_rng(3)
    for i, g in enumerate(GROUP_ORDER):
        v = sra.loc[sra.group == g, "Bases"].values / 1e6
        ax.scatter(i + rng.uniform(-0.12, 0.12, len(v)), v, c=GROUP_COL[g], marker=GROUP_MARK[g], ec="k", lw=0.5, s=28, zorder=3)
        ax.hlines(np.median(v), i - 0.28, i + 0.28, color="k", lw=1.1, zorder=4)
    ax.set_yscale("log"); ax.set_ylabel("Sequenced bases (Mb)")
    ax.set_xticks(range(4)); ax.set_xticklabels([GROUP_LAB[g] for g in GROUP_ORDER], rotation=45, ha="right", fontsize=6)
    ax.set_title("from SRA run table (group level)", fontsize=6, loc="left"); panel(ax, "e")
    ax = axs[1, 2]; ax.axis("off")
    ax.legend(handles=[Patch(fc=GROUP_COL[g], ec="k", lw=0.5, label=GROUP_LAB[g]) for g in GROUP_ORDER], loc="center", frameon=False, fontsize=7)
    save(fig, "Fig2_quality_control")


def fig_pca_distance(Yk, meta):
    pca = PCA(n_components=3).fit(Yk.T.values)
    sc = pca.transform(Yk.T.values); ev = pca.explained_variance_ratio_ * 100
    fig = plt.figure(figsize=(7.2, 3.5))
    ax = fig.add_axes([0.07, 0.14, 0.36, 0.78]); panel(ax, "a", -0.17)
    for g in GROUP_ORDER:
        i = np.where(meta.group.values == g)[0]
        ax.scatter(sc[i, 0], sc[i, 1], s=60, c=GROUP_COL[g], marker=GROUP_MARK[g], ec="k", lw=0.6, label=GROUP_LAB[g], zorder=3)
    texts = [ax.text(sc[i, 0], sc[i, 1], s, fontsize=6) for i, s in enumerate(meta.index)]
    adjust_text(texts, ax=ax, expand=(1.4, 1.8), arrowprops=dict(arrowstyle="-", lw=0.4, color="#666"))
    ax.axhline(0, c="#ddd", lw=0.6, zorder=0); ax.axvline(0, c="#ddd", lw=0.6, zorder=0)
    ax.set_xlabel(f"PC1: {ev[0]:.0f}% variance"); ax.set_ylabel(f"PC2: {ev[1]:.0f}% variance")
    ax.legend(frameon=False, loc="best", fontsize=6)
    d = pdist(Yk.T.values, "euclidean"); Z = linkage(d, "average"); o = leaves_list(Z)
    D = squareform(d)[np.ix_(o, o)]; names = meta.index[o]
    axd = fig.add_axes([0.54, 0.76, 0.30, 0.16]); axd.axis("off")
    dendrogram(Z, ax=axd, no_labels=True, color_threshold=0, above_threshold_color="#555", link_color_func=lambda k: "#555")
    axh = fig.add_axes([0.54, 0.14, 0.30, 0.60]); cax = fig.add_axes([0.89, 0.14, 0.015, 0.60])
    im = axh.imshow(D, cmap="Blues_r", aspect="auto")
    axh.set_xticks(range(12)); axh.set_yticks(range(12)); axh.set_xticklabels(names, rotation=90); axh.set_yticklabels(names)
    for lb in axh.get_xticklabels() + axh.get_yticklabels():
        lb.set_color(GROUP_COL[meta.loc[lb.get_text(), "group"]])
    for sp in axh.spines.values(): sp.set_visible(False)
    axh.tick_params(length=0)
    cb = fig.colorbar(im, cax=cax); cb.set_label("Euclidean distance", fontsize=7); cb.outline.set_visible(False)
    axh.text(-0.22, 1.30, "b", transform=axh.transAxes, fontsize=11, fontweight="bold")
    save(fig, "Fig3_PCA_and_sample_distance")
    scores = pd.DataFrame(sc, index=meta.index, columns=["PC1", "PC2", "PC3"])
    scores.assign(group=meta.group).to_csv(TAB / "pca_coordinates.csv")
    return ev, scores


def fig_top_genes(res, tpm, meta, ann):
    hfq = ann.index[ann.name == "hfq"][0]
    mt = res[[c for c in res.columns if c.startswith("meanTPM")]].max(axis=1)
    cand = res[(res.status != "NS") & (res.robust_DEG) & (mt >= 10) & (res.index != hfq)]
    genes = [hfq] + cand[cand.status == "Down"].index[:3].tolist() + cand[cand.status == "Up"].index[:2].tolist()
    fig, axs = plt.subplots(2, 3, figsize=(7.2, 4.6), gridspec_kw={"wspace": 0.5, "hspace": 0.85})
    rng = np.random.default_rng(1)
    for ax, g in zip(axs.ravel(), genes):
        for i, grp in enumerate(GROUP_ORDER):
            v = np.log2(tpm.loc[g, meta.index[meta.group == grp]].values + 1)
            ax.scatter(i + rng.uniform(-0.12, 0.12, len(v)), v, c=GROUP_COL[grp], marker=GROUP_MARK[grp], ec="k", lw=0.5, s=30, zorder=3)
            ax.hlines(v.mean(), i - 0.28, i + 0.28, color="k", lw=1.2, zorder=4)
        ax.set_xticks(range(4)); ax.set_xticklabels([GROUP_LAB[g_] for g_ in GROUP_ORDER], rotation=45, ha="right")
        ax.set_title(ann.loc[g, "label"] + ("  (positive control)" if g == hfq else ""), style="italic", fontsize=7.5)
        ax.set_ylabel("log$_2$(TPM + 1)")
    save(fig, "Fig4_hfq_control_and_top_genes")


def fig_heatmap_volcano(res, Yk, meta, ann):
    sig = res[res.status != "NS"].sort_values("FDR")
    sel = sig.head(40).index
    fig = plt.figure(figsize=(7.2, 5.0))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.15], wspace=0.55)
    ax = fig.add_subplot(gs[0]); panel(ax, "a", -0.34)
    z = Yk.loc[sel]; z = z.sub(z.mean(axis=1), axis=0).div(z.std(axis=1, ddof=0).replace(0, 1), axis=0)
    z = z.iloc[leaves_list(linkage(z.values, "average", metric="correlation"))]
    cols = [s for g in GROUP_ORDER for s in meta.index[meta.group == g]]
    sns.heatmap(z[cols], ax=ax, cmap="RdBu_r", center=0, vmin=-2.5, vmax=2.5, linewidths=0.2, linecolor="white",
                cbar_kws=dict(label="Row z-score of log$_2$(TPM + 1)", shrink=0.5, pad=0.03), yticklabels=True, xticklabels=True)
    ax.set_xticklabels(cols, rotation=90, fontsize=6); ax.set_yticklabels([ann.loc[g, "label"] + (" *" if res.loc[g, "robust_DEG"] else "") + (" \u2020" if res.loc[g, "group_exclusive"] else "")
                        for g in z.index], rotation=0, style="italic", fontsize=5.6)
    ax.text(0.0, -0.17, "* high-confidence DEG    \u2020 expressed in one group only", transform=ax.transAxes, fontsize=5.6, color="#444")
    ax.set_xlabel(""); ax.set_ylabel("")
    for i, s in enumerate(cols):
        ax.add_patch(plt.Rectangle((i, -1.6), 1, 1.2, color=GROUP_COL[meta.loc[s, "group"]], clip_on=False, lw=0))
    ax.legend(handles=[Patch(fc=GROUP_COL[g], label=GROUP_LAB[g]) for g in GROUP_ORDER], ncol=2, frameon=False, loc="lower center",
              bbox_to_anchor=(0.45, 1.06), fontsize=6, columnspacing=0.8, handlelength=1.0, handletextpad=0.3)
    ax2 = fig.add_subplot(gs[1]); panel(ax2, "b", -0.18)
    volcano(ax2, res, ann, n_label=12)
    save(fig, "Fig5_heatmap_and_volcano")


def volcano(ax, res, ann, n_label=10, title=None, small=False, legend_loc="upper right"):
    r = res.dropna(subset=["FDR"]).copy(); r["nlp"] = -np.log10(r.FDR.clip(lower=1e-300))
    for st, c, z in (("NS", NS, 1), ("Down", DOWN, 2), ("Up", UP, 2)):
        s = r[r.status == st]
        ax.scatter(s.log2FC, s.nlp, s=(3 if small else 6) if st == "NS" else (6 if small else 10), c=c, lw=0, alpha=0.7 if st == "NS" else 0.9, zorder=z, rasterized=True)
    if "robust_DEG" in r.columns and r.robust_DEG.any():
        rr = r[r.robust_DEG]
        ax.scatter(rr.log2FC, rr.nlp, s=(14 if small else 26), facecolors="none", edgecolors="k", lw=0.6, zorder=4, rasterized=True)
    ax.axhline(-np.log10(FDR_CUT), ls="--", c="#777", lw=0.7)
    for x in (-LFC_CUT, LFC_CUT): ax.axvline(x, ls="--", c="#777", lw=0.7)
    if n_label:
        lab = r[r.status != "NS"].sort_values("FDR").head(n_label)
        texts = [ax.text(row.log2FC, row.nlp, row.label, fontsize=6, style="italic") for _, row in lab.iterrows()]
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", lw=0.4, color="#555"), expand=(1.3, 1.6))
    ax.set_xlabel("log$_2$ fold change"); ax.set_ylabel("\u2212log$_{10}$ FDR")
    n_up, n_dn = (r.status == "Up").sum(), (r.status == "Down").sum()
    if title: ax.set_title(f"{title}\n\u2191 {n_up}   \u2193 {n_dn}", fontsize=7, loc="left")
    elif legend_loc is None: pass
    else: ax.legend(handles=[Line2D([], [], marker="o", ls="", c=UP, label=f"Up ({n_up})"), Line2D([], [], marker="o", ls="", c=DOWN, label=f"Down ({n_dn})"),
                             Line2D([], [], marker="o", ls="", c=NS, label="Not significant")]
                             + ([Line2D([], [], marker="o", ls="", mfc="none", mec="k", label="High-confidence")] if "robust_DEG" in r.columns else []),
                             frameon=False, loc=legend_loc)


def fig_contrast_grid(results, ann):
    fig, axs = plt.subplots(2, 2, figsize=(7.2, 6.0), gridspec_kw={"wspace": 0.3, "hspace": 0.5})
    for ax, (k, r), l in zip(axs.ravel(), results.items(), "abcd"):
        volcano(ax, r, ann, n_label=6, title=CONTRAST_LAB[k], small=True); panel(ax, l, -0.2, 1.1)
    save(fig, "Fig7_all_contrast_volcanoes")


def fig_diagnostics(results):
    fig, axs = plt.subplots(2, 4, figsize=(7.2, 3.8), gridspec_kw={"wspace": 0.55, "hspace": 0.65})
    for j, (k, r) in enumerate(results.items()):
        ax = axs[0, j]; ax.hist(r.pvalue.dropna(), bins=np.linspace(0, 1, 21), color="#7F8C8D", ec="white", lw=0.4)
        ax.set_title(CONTRAST_LAB[k], fontsize=6, loc="left"); ax.set_xlabel("p-value"); ax.set_ylabel("Genes" if j == 0 else "")
        ax = axs[1, j]; mean_log = np.log2((r[[c for c in r.columns if c.startswith("meanTPM")]].mean(axis=1)) + 1)
        ax.scatter(mean_log, r.log2FC, s=2, c=np.where(r.FDR < FDR_CUT, UP, NS), lw=0, rasterized=True)
        ax.axhline(0, c="k", lw=0.5); ax.set_xlabel("Mean log$_2$(TPM + 1)"); ax.set_ylabel("log$_2$FC" if j == 0 else "")
    panel(axs[0, 0], "a", -0.45, 1.15); panel(axs[1, 0], "b", -0.45, 1.15)
    save(fig, "Fig8_pvalue_histograms_and_MA")


def fig_interaction(wt, hfq, inter, ann):
    d = pd.DataFrame({"wt": wt.log2FC, "hfq": hfq.log2FC}).join(inter[["FDR", "log2FC"]].rename(columns={"FDR": "ifdr", "log2FC": "idelta"})).dropna()
    d = d.join(ann[["label"]])
    sig = d.ifdr < FDR_CUT
    fig, ax = plt.subplots(figsize=(4.6, 4.4))
    ax.scatter(d.wt[~sig], d.hfq[~sig], s=4, c=NS, lw=0, alpha=0.6, rasterized=True)
    ax.scatter(d.wt[sig], d.hfq[sig], s=10, c="#6C3483", lw=0, label=f"genotype \u00d7 culture FDR < {FDR_CUT} ({sig.sum()})", rasterized=True)
    lim = float(np.ceil(max(d.wt.abs().max(), d.hfq.abs().max())))
    ax.plot([-lim, lim], [-lim, lim], c="k", lw=0.7, ls="--"); ax.axhline(0, c="#aaa", lw=0.5); ax.axvline(0, c="#aaa", lw=0.5)
    lab = d[sig].assign(a=lambda x: x.idelta.abs()).sort_values("a", ascending=False).head(10)
    texts = [ax.text(r.wt, r.hfq, r.label, fontsize=6, style="italic") for _, r in lab.iterrows()]
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", lw=0.4, color="#555"))
    ax.set_xlabel("LFS effect in WT (log$_2$FC)"); ax.set_ylabel(r"LFS effect in $\Delta$$\it{hfq}$ (log$_2$FC)")
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.legend(frameon=False, loc="upper left", fontsize=6)
    ax.text(0.98, 0.03, "dashed line: equal effect in both strains", transform=ax.transAxes, ha="right", fontsize=5.8, color="#555")
    save(fig, "Fig9_LFS_effect_by_genotype")
    return int(sig.sum())


def fig_enrichment(enr, enr_rob):
    """Heat-dot summary of rule-based gene-set enrichment for all four contrasts and both directions."""
    cols = [(k, d) for k in CONTRASTS for d in ("Down", "Up")]
    e = enr.assign(col=list(zip(enr.contrast, enr.direction)))
    keep = e[e.FDR < 0.25].gene_set.unique()
    if len(keep) == 0:
        keep = e.sort_values("FDR").gene_set.unique()[:6]
    order = e[e.gene_set.isin(keep)].groupby("gene_set").FDR.min().sort_values().index.tolist()[::-1]
    fig, ax = plt.subplots(figsize=(7.2, 0.45 * len(order) + 2.6))
    for i, gset in enumerate(order):
        for j, c in enumerate(cols):
            r = e[(e.gene_set == gset) & (e.col == c)]
            if r.empty:
                continue
            r = r.iloc[0]; v = -np.log10(max(r.FDR, 1e-300))
            col = (DOWN if c[1] == "Down" else UP) if r.FDR < FDR_CUT else "#D5D8DC"
            rb_ok = ((enr_rob.gene_set == gset) & (enr_rob.contrast == c[0]) & (enr_rob.direction == c[1]) & (enr_rob.FDR < FDR_CUT)).any()
            ax.scatter(j, i, s=30 + 70 * min(v, 6) / 6 if r.n_DEG_in_set else 8, c=col, ec="k" if r.FDR < FDR_CUT else "#999",
                       lw=2.0 if rb_ok else 0.6, zorder=3)
            if r.n_DEG_in_set:
                ax.text(j, i - 0.38, f"{int(r.n_DEG_in_set)}/{int(r.set_size)}", ha="center", va="top", fontsize=5.4, color="#333")
    ax.set_yticks(range(len(order))); ax.set_yticklabels(order); ax.set_ylim(-0.9, len(order) - 0.4)
    ax.set_xticks(range(len(cols))); ax.set_xticklabels([f"{CONTRAST_LAB[k]}\n{d}-regulated".replace("-regulated", "") for k, d in cols], rotation=40, ha="right", fontsize=6.6)
    ax.set_xlim(-0.6, len(cols) - 0.4); ax.grid(axis="both", lw=0.3, color="#eee", zorder=0)
    for sp in ("left", "bottom"): ax.spines[sp].set_visible(False)
    ax.tick_params(length=0)
    ax.legend(handles=[Line2D([], [], marker="o", ls="", mfc=DOWN, mec="k", ms=6, label="FDR < 0.05, lower"),
                       Line2D([], [], marker="o", ls="", mfc=UP, mec="k", ms=6, label="FDR < 0.05, higher"),
                       Line2D([], [], marker="o", ls="", mfc="#D5D8DC", mec="#999", ms=6, label="FDR \u2265 0.05"),
                       Line2D([], [], marker="o", ls="", mfc="w", mec="k", mew=2, ms=6, label="still FDR < 0.05 with\nhigh-confidence DEGs only")],
              frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=6.2, title="DEGs in set / set size\n(label under dot); size = \u2212log$_{10}$FDR", title_fontsize=5.8)
    save(fig, "Fig10_functional_enrichment")


def coexpression(Yk, res, ann, n=60, rho_cut=0.9):
    """Spearman co-expression among the n most significant DEGs of one contrast (12 samples); isolated genes are dropped."""
    genes = res[res.status != "NS"].sort_values("FDR").head(n).index
    X = Yk.loc[genes].T
    rho = pd.DataFrame(spearmanr(X.values).statistic, index=genes, columns=genes)
    G = nx.Graph(); G.add_nodes_from(genes)
    for i, a in enumerate(genes):
        for b in genes[i + 1:]:
            if abs(rho.loc[a, b]) >= rho_cut:
                G.add_edge(a, b, weight=abs(rho.loc[a, b]), sign=np.sign(rho.loc[a, b]), dist=1 - abs(rho.loc[a, b]) + 1e-3)
    n_iso = nx.number_of_isolates(G)
    G.remove_nodes_from(list(nx.isolates(G)))
    deg = dict(G.degree()); stre = dict(G.degree(weight="weight"))
    btw = nx.betweenness_centrality(G, weight="dist"); clo = nx.closeness_centrality(G, distance="dist")
    t = pd.DataFrame({"degree": deg, "weighted_degree": stre, "betweenness": btw, "closeness": clo})
    t["mean_rank"] = t.rank(ascending=False, method="average").mean(axis=1)
    t = t.join(ann[["label", "description"]]).join(res[["log2FC", "FDR", "status"]]).sort_values("mean_rank")
    t["hub"] = False; t.iloc[:5, t.columns.get_loc("hub")] = True
    t.index.name = "GeneID"
    return G, t, n_iso


def fig_network(G, hubs, ann, n_iso, contrast_label):
    fig = plt.figure(figsize=(7.2, 4.4))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.8, 1], wspace=0.3)
    ax = fig.add_subplot(gs[0]); panel(ax, "a", -0.02, 1.0); ax.axis("off")
    pos = nx.spring_layout(G, seed=42, k=2.2 / np.sqrt(len(G)), iterations=600, weight="weight")
    nx.draw_networkx_edges(G, pos, ax=ax, width=0.35, edge_color=["#C0392B" if G[u][v]["sign"] < 0 else "#8E8E8E" for u, v in G.edges()], alpha=0.35)
    nodes = nx.draw_networkx_nodes(G, pos, ax=ax, node_size=[25 + 7 * hubs.loc[n, "degree"] for n in G.nodes],
                                   node_color=hubs.loc[list(G.nodes), "log2FC"].clip(-6, 6).values, cmap="RdBu_r", vmin=-6, vmax=6,
                                   edgecolors=["#000" if hubs.loc[n, "hub"] else "#777" for n in G.nodes],
                                   linewidths=[1.8 if hubs.loc[n, "hub"] else 0.4 for n in G.nodes])
    show = set(hubs.sort_values("degree", ascending=False).index[:18]) | set(hubs.index[hubs.hub])
    texts = [ax.text(x_, y_, ann.loc[n_, "label"], fontsize=5.8, style="italic", fontweight="bold" if hubs.loc[n_, "hub"] else "normal", zorder=5)
             for n_, (x_, y_) in pos.items() if n_ in show]
    adjust_text(texts, ax=ax, expand=(1.5, 1.8), arrowprops=dict(arrowstyle="-", lw=0.3, color="#777"))
    ax.text(0.99, 0.99, f"{n_iso} further DEGs with no |rho| \u2265 0.9 partner not shown; labels: top-5 hubs and 18 highest-degree genes", transform=ax.transAxes, ha="right", va="top", fontsize=5.4, color="#555")
    cb = fig.colorbar(nodes, ax=ax, shrink=0.45, pad=0.0, aspect=18); cb.set_label("log$_2$FC: " + contrast_label, fontsize=6.0); cb.outline.set_visible(False)
    ax.legend(handles=[Line2D([], [], marker="o", ls="", mfc="w", mec="k", mew=1.8, ms=6, label="Top-5 hub"),
                       Line2D([], [], color="#8E8E8E", lw=1.2, label="positive correlation"), Line2D([], [], color="#C0392B", lw=1.2, label="negative correlation")],
              frameon=False, loc="lower left", fontsize=6, bbox_to_anchor=(-0.02, -0.04))
    ax = fig.add_subplot(gs[1]); panel(ax, "b", -0.5, 1.0)
    top = hubs.head(8).iloc[::-1]
    ax.barh(np.arange(len(top)), top.weighted_degree, color=["#6C3483" if h else "#B0BEC5" for h in top.hub], ec="k", lw=0.5, height=0.65)
    ax.set_yticks(np.arange(len(top)), labels=top.label, style="italic"); ax.set_xlabel("Weighted degree")
    ax.set_title("Top 8 by consensus rank", loc="left", fontsize=7.5)
    save(fig, "Fig11_coexpression_network_and_hubs")


def fig_robustness(summ, pcs, meta):
    fig = plt.figure(figsize=(7.2, 3.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.35, 0.8], wspace=0.75)
    ks = list(CONTRASTS)
    ax = fig.add_subplot(gs[0]); panel(ax, "a", -0.45, 1.04)
    b = summ[summ.analysis == "baseline"].set_index("contrast").loc[ks]; r_ = summ[summ.analysis == "robust"].set_index("contrast").loc[ks]
    y = np.arange(len(ks))
    ax.barh(y + 0.19, b.n_up + b.n_down, height=0.36, color="#B0BEC5", ec="k", lw=0.4, label="all DEGs")
    ax.barh(y - 0.19, r_.n_up + r_.n_down, height=0.36, color="#6C3483", ec="k", lw=0.4, label="high-confidence")
    for i, k in enumerate(ks):
        ax.text(b.loc[k, "n_up"] + b.loc[k, "n_down"] + 3, i + 0.19, str(int(b.loc[k, "n_up"] + b.loc[k, "n_down"])), va="center", fontsize=6)
        ax.text(r_.loc[k, "n_up"] + r_.loc[k, "n_down"] + 3, i - 0.19, str(int(r_.loc[k, "n_up"] + r_.loc[k, "n_down"])), va="center", fontsize=6)
    ax.set_yticks(y, labels=["\u0394hfq vs WT\n(control)", "\u0394hfq vs WT\n(LFS)", "LFS vs ctrl\n(WT)", "LFS vs ctrl\n(\u0394hfq)"], fontsize=6)
    ax.invert_yaxis(); ax.set_xlabel("Number of DEGs"); ax.legend(frameon=False, fontsize=5.8, loc="upper center", bbox_to_anchor=(0.4, -0.2), ncol=1)
    ax = fig.add_subplot(gs[1]); panel(ax, "b", -0.3, 1.04)
    an = ["no_exclusive", "detect_adj", "rrna_adj", "batch2_adj", "consistent", "robust"]
    M = np.array([[summ[(summ.contrast == k) & (summ.analysis == a)].pct_retained.iloc[0] for k in ks] for a in an])
    im = ax.imshow(M, cmap="Purples", vmin=0, vmax=100, aspect="auto")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            ax.text(j, i, f"{M[i, j]:.0f}", ha="center", va="center", fontsize=6.3, color="white" if M[i, j] > 60 else "k")
    ax.set_xticks(range(4), labels=["\u0394hfq vs WT\n(control)", "\u0394hfq vs WT\n(LFS)", "LFS vs ctrl\n(WT)", "LFS vs ctrl\n(\u0394hfq)"], fontsize=5.8, rotation=40, ha="right")
    ax.set_yticks(range(len(an)), labels=["no group-exclusive\ngenes", "+ detected-genes\ncovariate", "+ rRNA-share\ncovariate", "+ '_2' suffix\ncovariate", "consistent\ndetection", "high-confidence"], fontsize=6)
    ax.tick_params(length=0)
    for sp in ax.spines.values(): sp.set_visible(False)
    cb = fig.colorbar(im, ax=ax, shrink=0.7, pad=0.03); cb.set_label("% of baseline DEGs retained", fontsize=6.3); cb.outline.set_visible(False)
    ax = fig.add_subplot(gs[2]); panel(ax, "c", -0.55, 1.04)
    p = pcs[pcs.PC.isin(["PC1", "PC2"])]
    w = 0.26
    for i, (f, c) in enumerate((("genotype", "#D55E00"), ("culture", "#56B4E9"), ("suffix_2", "#7F7F7F"))):
        v = p[p.factor == f].set_index("PC").loc[["PC1", "PC2"], "eta2"]
        ax.bar(np.arange(2) + (i - 1) * w, v, w, color=c, ec="k", lw=0.4, label={"suffix_2": "name ends '_2'"}.get(f, f))
    ax.set_xticks(range(2), labels=["PC1", "PC2"]); ax.set_ylabel("Variance explained ($\\eta^2$)"); ax.set_ylim(0, 1.05)
    ax.legend(frameon=False, fontsize=5.8, loc="upper right")
    save(fig, "Fig13_robustness_checks")


def fig_summary(results, enr_all, ev):
    fig = plt.figure(figsize=(7.4, 3.5), layout="constrained")
    grid = fig.add_gridspec(2, 3, height_ratios=[0.14, 1], width_ratios=[1.05, 0.9, 1.35])
    t = fig.add_subplot(grid[0, :]); t.set_axis_off()
    t.text(0.5, 0.78, "Salmonella \u0394hfq vs wild type: exploratory TPM analysis of GSE146347", ha="center", va="center", fontsize=10.5, fontweight="bold", color="#17324D")
    t.text(0.5, 0.1, "12 bacterial samples, 2 x 2 design (genotype x culture), n = 3 per group; TPM only, no raw counts", ha="center", va="center", fontsize=6.6, color="#555")
    ax = fig.add_subplot(grid[1, 0]); panel(ax, "a", -0.2, 1.03)
    r = results[PRIMARY]; volcano(ax, r, ANN, n_label=4, title=CONTRAST_LAB[PRIMARY], small=True)
    ax = fig.add_subplot(grid[1, 1]); panel(ax, "b", -0.25, 1.03)
    ks = list(results); up = [(results[k].status == "Up").sum() for k in ks]; dn = [(results[k].status == "Down").sum() for k in ks]
    rup = [((results[k].status == "Up") & results[k].robust_DEG).sum() for k in ks]; rdn = [((results[k].status == "Down") & results[k].robust_DEG).sum() for k in ks]
    y = np.arange(len(ks))
    ax.barh(y, up, color=UP, alpha=0.3, ec="k", lw=0.4, height=0.7); ax.barh(y, [-d for d in dn], color=DOWN, alpha=0.3, ec="k", lw=0.4, height=0.7)
    ax.barh(y, rup, color=UP, ec="k", lw=0.4, height=0.35); ax.barh(y, [-d for d in rdn], color=DOWN, ec="k", lw=0.4, height=0.35)
    ax.set_yticks(y, labels=["\u0394hfq vs WT\n(control)", "\u0394hfq vs WT\n(LFS)", "LFS vs ctrl\n(WT)", "LFS vs ctrl\n(\u0394hfq)"], fontsize=6)
    for i, (u, d) in enumerate(zip(up, dn)):
        ax.text(u + 3, i, str(u), va="center", fontsize=6); ax.text(-d - 3, i, str(d), va="center", ha="right", fontsize=6)
    ax.axvline(0, c="k", lw=0.6); ax.invert_yaxis(); ax.set_xlabel("Genes lower (blue) / higher (red) in first-named group\npale = all DEGs, dark = high-confidence", fontsize=6.3)
    m = max(max(up), max(dn)) * 1.25; ax.set_xlim(-m, m)
    ax = fig.add_subplot(grid[1, 2]); panel(ax, "c", -0.12, 1.03)
    short = {"hfq_effect_control": "\u0394hfq vs WT, ctrl", "hfq_effect_LFS": "\u0394hfq vs WT, LFS", "LFS_effect_WT": "LFS vs ctrl, WT", "LFS_effect_hfq": "LFS vs ctrl, \u0394hfq"}
    e = enr_all[enr_all.FDR < FDR_CUT].sort_values("FDR").head(7).iloc[::-1]
    if len(e):
        ax.barh([f"{g}\n({short[c]}, {d.lower()})" for g, c, d in zip(e.gene_set, e.contrast, e.direction)], -np.log10(e.FDR.clip(lower=1e-300)),
                color=[UP if d == "Up" else DOWN for d in e.direction], ec="k", lw=0.4, height=0.65)
        ax.axvline(-np.log10(FDR_CUT), ls="--", c="#555", lw=0.7); ax.tick_params(axis="y", labelsize=5.4)
    else:
        ax.text(0.5, 0.5, "no gene set at FDR < 0.05", ha="center", transform=ax.transAxes)
    ax.set_xlabel("\u2212log$_{10}$ FDR"); ax.set_title("Gene sets enriched at FDR < 0.05", fontsize=7, loc="left")
    save(fig, "Fig0_graphical_summary", note=False)


# ----------------------------------------------------------------------------- main
ANN = None


def main():
    global ANN
    tpm_raw, ann, meta = load(); ANN = ann
    tpm = renormalise(tpm_raw, ann)
    tpm.to_csv(TAB / "tpm_renormalised_no_structural_RNA.csv")
    print(f"{tpm_raw.shape[0]} genes x {tpm_raw.shape[1]} samples; groups {meta.group.value_counts().to_dict()}")
    print(f"excluded {tpm_raw.shape[0] - tpm.shape[0]} rRNA/tRNA/ncRNA genes and renormalised TPM over {tpm.shape[0]} genes")
    Y = np.log2(tpm + 1)
    kept = (tpm >= MIN_TPM).sum(axis=1) >= MIN_SAMPLES
    Yk = Y[kept]; print(f"tested genes (TPM >= {MIN_TPM:g} in >= {MIN_SAMPLES} samples): {kept.sum():,} of {len(kept):,}")
    model = Model(Yk, meta)

    results = {}
    for k, (a, b) in CONTRASTS.items():
        results[k] = annotate(model.test({a: 1.0, b: -1.0}), tpm, ann, meta)
        results[k].to_csv(TAB / f"de_{k}.csv")
        print(f"  {k}: {(results[k].status=='Up').sum()} up, {(results[k].status=='Down').sum()} down")
    inter = annotate(model.test({"delta_hfq_LFS": 1, "delta_hfq_control": -1, "WT_LFS": -1, "WT_control": 1}), tpm, ann, meta)
    inter.to_csv(TAB / "de_interaction_genotype_x_culture.csv")
    results, rob_summary, extra = rb.run(tpm, tpm_raw, ann, meta, Yk, results, CONTRASTS, FDR_CUT, LFC_CUT)
    for k, df in results.items():
        df.to_csv(TAB / f"de_{k}.csv")
        df[df.robust_DEG].to_csv(TAB / f"high_confidence_DEGs_{k}.csv")
    rob_summary.to_csv(TAB / "robustness_summary.csv", index=False)
    for n_, df_ in extra.items():
        df_.to_csv(TAB / f"{n_}.csv")
    print(rob_summary[rob_summary.analysis.isin(["baseline", "robust"])].pivot(index="contrast", columns="analysis", values="n_up").add(
        rob_summary[rob_summary.analysis.isin(["baseline", "robust"])].pivot(index="contrast", columns="analysis", values="n_down")).rename(columns=lambda c: c + " DEGs").to_string())
    hfq_id = ann.index[ann.name == "hfq"][0]
    r = results[PRIMARY].loc[hfq_id]
    print(f"positive control hfq: log2FC {r.log2FC:.2f}, FDR {r.FDR:.2e}, meanTPM WT {r.meanTPM_WT_control:.0f} vs dhfq {r.meanTPM_delta_hfq_control:.1f}")
    results[PRIMARY][results[PRIMARY].status != "NS"].to_csv(TAB / "significant_DEGs_hfq_effect_control.csv")
    print(f"  hfq positive control robust? {bool(r.robust_DEG)}  | group-exclusive genes among tested: {int(results[PRIMARY].group_exclusive.sum())}")

    print("Figures ...")
    fig_workflow()
    fig_qc(tpm_raw, tpm, meta, ann, Y.mean(axis=1), kept)
    ev, scores = fig_pca_distance(Yk, meta)
    pcs = rb.pc_association(scores, meta); pcs.to_csv(TAB / "pc_variance_by_factor.csv", index=False)
    fig_top_genes(results[PRIMARY], tpm, meta, ann)
    fig_robustness(rob_summary, pcs, meta)
    fig_heatmap_volcano(results[PRIMARY], Yk, meta, ann)
    volcano_guide.interpretation(results[PRIMARY].assign(label=results[PRIMARY].label), r"log$_2$ fold change ($\Delta$$\it{hfq}$ / WT, control culture)",
                                 "Fig6_volcano_interpretation", thr_fdr=FDR_CUT, lfc_cut=LFC_CUT, n_label=10)
    volcano_guide.schematic("Fig12_volcano_reading_guide")
    fig_contrast_grid(results, ann)
    fig_diagnostics(results)
    n_int = fig_interaction(results["LFS_effect_WT"], results["LFS_effect_hfq"], inter, ann)
    print(f"  genotype x culture interaction genes (FDR<{FDR_CUT}): {n_int}")

    print("Enrichment ...")
    sets = build_sets(ann)
    pd.DataFrame({k: v.reindex(ann.index) for k, v in sets.items()}).astype(int).to_csv(TAB / "gene_set_membership.csv")
    bg = Yk.index
    enr_all = []
    for k in ("hfq_effect_control", "hfq_effect_LFS", "LFS_effect_WT", "LFS_effect_hfq"):
        for d in ("Up", "Down"):
            e = enrichment(results[k].index[results[k].status == d], bg, sets, d)
            if len(e): enr_all.append(e.assign(contrast=k))
    enr_all = pd.concat(enr_all).sort_values("FDR")
    enr_all.to_csv(TAB / "functional_enrichment.csv", index=False)
    enr_rob = []
    for k in CONTRASTS:
        for d in ("Up", "Down"):
            e = enrichment(results[k].index[(results[k].status == d) & results[k].robust_DEG], bg, sets, d)
            if len(e): enr_rob.append(e.assign(contrast=k))
    enr_rob = pd.concat(enr_rob).sort_values("FDR")
    enr_rob.to_csv(TAB / "functional_enrichment_high_confidence.csv", index=False)
    print("  enriched (FDR<0.05) with high-confidence DEGs only:",
          enr_rob[enr_rob.FDR < FDR_CUT][["contrast", "direction", "gene_set", "n_DEG_in_set", "set_size"]].to_dict("records"))
    fig_enrichment(enr_all, enr_rob)
    print(enr_all[enr_all.contrast == PRIMARY].head(6)[["direction", "gene_set", "n_DEG_in_set", "set_size", "FDR"]].round(4).to_string(index=False))

    print("Co-expression network ...")
    NET = "LFS_effect_hfq"
    net_res = results[NET].copy(); net_res.loc[~net_res.robust_DEG, "status"] = "NS"
    G, hubs, n_iso = coexpression(Yk, net_res, ann)
    hubs.to_csv(TAB / "coexpression_hubs_LFS_effect_hfq.csv")
    fig_network(G, hubs, ann, n_iso, CONTRAST_LAB[NET].replace("$", "").replace("\\Delta", "\u0394").replace("\\it{hfq}", "hfq"))
    print(f"  {G.number_of_nodes()} nodes, {G.number_of_edges()} edges ({n_iso} isolated DEGs dropped); hubs: {', '.join(hubs.label.head(5))}")
    fig_summary(results, enr_all, ev)
    print("Done.")


if __name__ == "__main__":
    main()
