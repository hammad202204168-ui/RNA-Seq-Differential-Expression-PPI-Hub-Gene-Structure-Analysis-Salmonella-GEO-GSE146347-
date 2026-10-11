#!/usr/bin/env python3
"""
IGV-style alignment view (Fig13) and Sashimi-style junction plot (Fig14), after Fig. 4 of Batut et al. (2021).

Both figures are drawn from the SYNTHETIC BAMs made by scripts/simulate_alignments.py; they illustrate the display,
not real sequencing data. Run simulate_alignments.py first.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
from matplotlib.collections import LineCollection
from matplotlib.patches import FancyArrowPatch, Rectangle

from rna_seq_pipeline import COND_COL, ROOT, panel, save

AL = ROOT / "data" / "synthetic_alignments"
SYN_NOTE = "Synthetic reads simulated from the demonstration count matrix \u2013 not sequencing data"
FWD, REV = "#7FA3C7", "#D49A9A"


def coverage(bam, contig, start, end):
    cov = np.zeros(end - start)
    for a in bam.fetch(contig, start, end):
        for bs, be in a.get_blocks():
            lo, hi = max(bs, start), min(be, end)
            if hi > lo:
                cov[lo - start:hi - start] += 1
    return cov


def pack_rows(reads, max_rows):
    """Greedy packing of (start, end, strand) into rows; returns segments (x0, x1, row, strand)."""
    ends, segs = [], []
    for s, e, st in sorted(reads):
        for r, last in enumerate(ends):
            if s > last + 5:
                ends[r] = e; break
        else:
            r = len(ends); ends.append(e)
        if r < max_rows:
            segs.append((s, e, r, st))
    return segs, len(ends)


def gene_track(ax, genes, start, end, scale=1000):
    ax.set_xlim((start) / scale, end / scale); ax.set_ylim(-1, 1)
    for sp in ("left", "top", "right"):
        ax.spines[sp].set_visible(False)
    ax.set_yticks([])
    ax.hlines(0, start / scale, end / scale, color="#333", lw=0.6)
    for g, r in genes.iterrows():
        if r.end < start or r.start > end:
            continue
        x0, x1 = max(r.start, start) / scale, min(r.end, end) / scale
        ax.add_patch(Rectangle((x0, -0.38), x1 - x0, 0.76, fc="#2E4A6B", ec="none"))
        ax.text((x0 + x1) / 2, 0.62, g, ha="center", va="bottom", fontsize=6.3, style="italic")
        tip = 0.06 if r.strand == "+" else -0.06
        xm = (x0 + x1) / 2
        for dx in np.linspace(x0 + 0.1 * (x1 - x0), x1 - 0.1 * (x1 - x0), 4):
            ax.plot([dx - tip, dx, dx - tip], [-0.2, 0, 0.2], color="white", lw=0.6)


def fig_alignment_view():
    genes = pd.read_csv(AL / "synthetic_gene_coordinates.csv", index_col="gene")
    first, last = "hilA", "sipA"
    start = int(genes.loc[first, "start"]) - 300
    end = int(genes.loc[last, "end"]) + 300
    samples = [("WT_Rep1", "Wild type (WT_Rep1)", COND_COL["Control"]), ("Mut_Rep1", "\u0394hfq (Mut_Rep1)", COND_COL["Mutant"])]
    covs, reads = {}, {}
    for sample, _, _ in samples:
        with pysam.AlignmentFile(str(AL / f"{sample}.bam")) as bam:
            covs[sample] = coverage(bam, "synthetic_ref", start, end)
            reads[sample] = [(a.reference_start, a.reference_end, "-" if a.is_reverse else "+")
                             for a in bam.fetch("synthetic_ref", start, end)]
    ymax = max(c.max() for c in covs.values()) * 1.08
    fig = plt.figure(figsize=(7.2, 6.2))
    gs = fig.add_gridspec(5, 1, height_ratios=[1.0, 1.7, 1.0, 1.7, 0.55], hspace=0.12)
    xs = np.arange(start, end) / 1000
    for i, (sample, label, col) in enumerate(samples):
        axc = fig.add_subplot(gs[2 * i]); axr = fig.add_subplot(gs[2 * i + 1], sharex=axc)
        axc.fill_between(xs, covs[sample], step="mid", color=col, alpha=0.85, lw=0)
        axc.set_ylim(0, ymax); axc.set_ylabel("Coverage"); axc.tick_params(labelbottom=False, bottom=False)
        panel(axc, "ab"[i], -0.1, 1.0)
        axc.text(0.995, 0.93, label, transform=axc.transAxes, ha="right", va="top", fontsize=7.2, fontweight="bold", color=col)
        segs, nrows = pack_rows(reads[sample], 45)
        lines = [[(s / 1000, -r), (e / 1000, -r)] for s, e, r, st in segs]
        cols = [FWD if st == "+" else REV for *_, st in segs]
        axr.add_collection(LineCollection(lines, colors=cols, linewidths=0.9))
        axr.set_ylim(-46, 1); axr.set_yticks([]); axr.tick_params(bottom=False, labelbottom=False)
        axr.spines["left"].set_visible(False); axr.spines["bottom"].set_visible(False)
        if nrows > 45:
            axc.text(0.995, 0.60, f"read track: first 45 of {nrows} rows", transform=axc.transAxes, ha="right", va="top",
                     fontsize=5.8, color="#666")
    axg = fig.add_subplot(gs[4], sharex=axc)
    gene_track(axg, genes, start, end)
    axc.set_xlim(start / 1000, end / 1000)
    axg.set_xlabel("Position on synthetic reference (kb)")
    fig.legend(handles=[plt.Line2D([], [], color=FWD, lw=2, label="Forward-strand read"),
                        plt.Line2D([], [], color=REV, lw=2, label="Reverse-strand read")],
               loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 0.0), fontsize=6.5)
    save(fig, "Fig13_alignment_view_IGV_style", text=SYN_NOTE)


def fig_sashimi():
    gtf = pd.read_csv(AL / "spliced_locus_demo.gtf", sep="\t", header=None)
    exons = [(int(r[3]) - 1, int(r[4])) for _, r in gtf[gtf[2] == "exon"].iterrows()]
    start, end = exons[0][0] - 150, exons[-1][1] + 150
    with pysam.AlignmentFile(str(AL / "spliced_locus_demo.bam")) as bam:
        cov = coverage(bam, "spliced_locus_demo", start, end)
        junc = {}
        for a in bam.fetch("spliced_locus_demo"):
            blocks = a.get_blocks()
            for (b0, b1), (c0, c1) in zip(blocks[:-1], blocks[1:]):
                junc[(b1, c0)] = junc.get((b1, c0), 0) + 1
    fig = plt.figure(figsize=(7.2, 3.6))
    gs = fig.add_gridspec(2, 1, height_ratios=[3.2, 0.6], hspace=0.05)
    ax = fig.add_subplot(gs[0]); panel(ax, "a", -0.08, 1.0)
    xs = np.arange(start, end)
    ax.fill_between(xs, cov, step="mid", color="#C0392B", alpha=0.85, lw=0)
    ax.set_xlim(start, end); ax.set_ylim(-cov.max() * 0.75, cov.max() * 1.9)
    ax.set_ylabel("Coverage"); ax.spines["bottom"].set_visible(False); ax.tick_params(labelbottom=False, bottom=False)
    ax.set_yticks([t for t in ax.get_yticks() if 0 <= t <= cov.max()]); ax.spines["left"].set_bounds(0, cov.max())
    top = cov.max()
    for k, ((a, b), n) in enumerate(sorted(junc.items(), key=lambda t: t[0][0] - t[0][1])):
        above = k % 2 == 0
        h = (0.55 + 0.35 * (b - a) / (end - start)) * top
        xsa = np.linspace(a, b, 100); ya = np.sin(np.linspace(0, np.pi, 100)) * h
        ax.plot(xsa, ya if above else -ya * 0.8, color="#7B1F16", lw=0.6 + 3.0 * n / max(junc.values()))
        ax.text((a + b) / 2, (h * 1.0 + 0.03 * top) if above else -h * 0.8 - 0.03 * top, str(n), ha="center",
                va="bottom" if above else "top", fontsize=7, fontweight="bold")
    axg = fig.add_subplot(gs[1], sharex=ax); axg.set_ylim(-1, 1); axg.axis("off")
    axg.hlines(0, exons[0][0], exons[-1][1], color="#333", lw=0.8)
    for s, e in exons:
        axg.add_patch(Rectangle((s, -0.4), e - s, 0.8, fc="#2E4A6B", ec="none"))
    axg.text(end, -0.95, "", fontsize=6)
    fig.text(0.5, 0.01, "Synthetic eukaryote-style locus for illustrating the display. Salmonella mRNAs are not spliced.",
             ha="center", fontsize=6.5, color="#555")
    save(fig, "Fig14_sashimi_style_junctions", text=SYN_NOTE)
    return junc


if __name__ == "__main__":
    fig_alignment_view()
    print("junctions:", fig_sashimi())
