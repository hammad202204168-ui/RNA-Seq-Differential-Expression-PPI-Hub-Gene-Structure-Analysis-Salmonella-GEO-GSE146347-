#!/usr/bin/env python3
"""PPI network (STRING) -> hub genes -> predicted structures (AlphaFold DB) for the high-confidence DEGs.

REQUIRES INTERNET (string-db.org, rest.uniprot.org, alphafold.ebi.ac.uk). NOT executed against the live
servers in the sandbox that wrote it; the parsing, hub, pLDDT and plotting logic is unit-tested offline
(tests/test_ppi_structures.py) with fabricated inputs.

Run:  python scripts/ppi_hubs_structures.py --string_edges <STRING short tsv> --string_annotations <STRING annotations tsv>
          (offline: uses a STRING web export;  without them it calls the STRING API at --score 700)
      add --skip_structures to stop before the UniProt/AlphaFold step

Outputs
  results/tables/ppi_string_ids.csv        input gene -> STRING id (which genes STRING could / could not resolve)
  results/tables/ppi_edges.csv             STRING edges (experimental+database+text-mining combined score)
  results/tables/ppi_hubs.csv              degree / weighted degree / betweenness per gene, hub flag
  results/tables/ppi_vs_coexpression_hubs.csv   overlap with the co-expression hubs already in the repo
  results/tables/hub_structures.csv        UniProt id, AlphaFold entry, mean pLDDT, % residues < 70, file
  results/structures/<gene>_<uniprot>.pdb  AlphaFold models;   results/figures/Fig14_*, Fig15_*

Interpretation limits (also written into the tables' README line below):
  * STRING edges are functional-association evidence, not proof of physical binding.
  * A hub is central in THIS network (partly because motility/SPI-1 genes are well studied and densely
    annotated in STRING); it is not shown to be more important biologically.
  * AlphaFold models are predictions; low pLDDT (< 70) regions should not be interpreted structurally.
"""
import argparse, io, json, os, re, sys, urllib.parse, urllib.request
import numpy as np, pandas as pd
import networkx as nx

SPECIES = 99287            # Salmonella enterica subsp. enterica serovar Typhimurium LT2 (the annotation used here)
UA = {"User-Agent": "hfq-analysis/1.0"}


# ------------------------------------------------------------------------------------------------ network access
def _post(url, params, timeout=120):
    data = urllib.parse.urlencode(params).encode()
    return urllib.request.urlopen(urllib.request.Request(url, data=data, headers=UA), timeout=timeout).read().decode()


def _get(url, timeout=120):
    return urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=timeout).read().decode()


def string_resolve(genes):
    txt = _post("https://string-db.org/api/tsv/get_string_ids",
                dict(identifiers="\r".join(genes), species=SPECIES, limit=1, caller_identity="hfq-analysis"))
    return parse_string_ids(txt, genes)


def string_network(string_ids, score):
    txt = _post("https://string-db.org/api/tsv/network",
                dict(identifiers="\r".join(string_ids), species=SPECIES, required_score=score, caller_identity="hfq-analysis"))
    return parse_string_edges(txt)


def uniprot_lookup(gene):
    """Best UniProt entry for a gene symbol or locus tag in LT2 (reviewed entries preferred)."""
    q = urllib.parse.quote(f"(gene:{gene}) AND (organism_id:{SPECIES})")
    url = f"https://rest.uniprot.org/uniprotkb/search?query={q}&fields=accession,reviewed,gene_names,protein_name&format=tsv&size=5"
    return pick_uniprot(_get(url))


def alphafold_entry(acc):
    j = json.loads(_get(f"https://alphafold.ebi.ac.uk/api/prediction/{acc}"))
    return j[0] if j else None


# ------------------------------------------------------------------------------------------------ pure logic (tested)
def parse_string_ids(tsv, wanted):
    d = pd.read_csv(io.StringIO(tsv), sep="\t") if tsv.strip() else pd.DataFrame(columns=["queryItem", "stringId", "preferredName"])
    d = d.drop_duplicates("queryItem")[["queryItem", "stringId", "preferredName"]]
    out = pd.DataFrame({"gene": list(wanted)}).merge(d, left_on="gene", right_on="queryItem", how="left").drop(columns="queryItem")
    return out


def parse_string_edges(tsv):
    d = pd.read_csv(io.StringIO(tsv), sep="\t")
    d = d[["stringId_A", "stringId_B", "preferredName_A", "preferredName_B", "score"]].copy()
    key = d.apply(lambda r: tuple(sorted((r.stringId_A, r.stringId_B))), axis=1)
    return d.loc[~key.duplicated()].reset_index(drop=True)


def hubs_from_edges(edges, top=10, a="preferredName_A", b="preferredName_B"):
    G = nx.Graph()
    for _, r in edges.iterrows():
        if r[a] != r[b]:
            G.add_edge(r[a], r[b], weight=float(r["score"]))
    if G.number_of_nodes() == 0:
        return pd.DataFrame(columns=["gene", "degree", "weighted_degree", "betweenness", "mean_rank", "hub"]), G
    t = pd.DataFrame({"gene": list(G.nodes),
                      "degree": [G.degree(n) for n in G.nodes],
                      "weighted_degree": [G.degree(n, weight="weight") for n in G.nodes]})
    bc = nx.betweenness_centrality(G)
    t["betweenness"] = t.gene.map(bc)
    t["mean_rank"] = t[["degree", "weighted_degree", "betweenness"]].rank(ascending=False, method="min").mean(axis=1)
    t = t.sort_values(["mean_rank", "degree"], ascending=[True, False]).reset_index(drop=True)
    t["hub"] = t.index < top
    return t, G


def pick_uniprot(tsv):
    if not tsv.strip():
        return None
    d = pd.read_csv(io.StringIO(tsv), sep="\t")
    if d.empty:
        return None
    d["rev"] = d["Reviewed"].astype(str).str.lower().eq("reviewed")
    r = d.sort_values("rev", ascending=False).iloc[0]
    return dict(uniprot=r["Entry"], reviewed=bool(r["rev"]), protein=r.get("Protein names", ""))


def plddt_from_pdb(pdb_text):
    """Per-residue pLDDT (AlphaFold stores it in the B-factor column of each C-alpha) and coordinates."""
    ca = [l for l in pdb_text.splitlines() if l.startswith("ATOM") and l[12:16].strip() == "CA"]
    xyz = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in ca])
    pl = np.array([float(l[60:66]) for l in ca])
    return xyz, pl


def find_structure_file(folder, acc):
    hits = sorted(f for f in os.listdir(folder) if acc in f and f.lower().endswith(".pdb"))
    return os.path.join(folder, hits[-1]) if hits else None      # latest model version if several


def plddt_summary(pl):
    return dict(n_residues=int(len(pl)), mean_pLDDT=float(np.mean(pl)), pct_residues_below_70=float(100 * np.mean(pl < 70)),
                pct_residues_above_90=float(100 * np.mean(pl > 90)))



STRING_CHANNELS = ["neighborhood_on_chromosome", "gene_fusion", "phylogenetic_cooccurrence", "homology", "coexpression",
                   "experimentally_determined_interaction", "database_annotated", "automated_textmining"]


def map_inputs_to_nodes(inputs, ann):
    """Map our gene labels to STRING node names: exact name first, then a UNIQUE alias hit
    (e.g. our 'spaN' is STRING's 'invJ'; our 'groES' is 'mopB'). ann needs columns node, other_names_and_aliases."""
    al = {r.node: set(re.split(r"[,; ]+", str(r.other_names_and_aliases))) for r in ann.itertuples()}
    out = []
    for l in inputs:
        if l in al:
            out.append((l, l, "exact"))
        else:
            hit = [n for n, a in al.items() if l in a]
            out.append((l, hit[0], "alias") if len(hit) == 1 else (l, None, "ambiguous_alias" if hit else "not_found_in_STRING"))
    return pd.DataFrame(out, columns=["input_gene", "string_node", "match"])


def combined_without(row, drop, prior=0.041):
    """STRING combined score recomputed without some evidence channels (STRING's own formula, prior-corrected)."""
    prod = 1.0
    for c in STRING_CHANNELS:
        if c not in drop:
            prod *= 1 - max((row[c] - prior) / (1 - prior), 0)
    s = 1 - prod
    return s + prior * (1 - s)


def hub_robustness(edges, top=10, threshold=0.7):
    """Top-`top` hubs when whole evidence classes are removed. Returns dict variant -> hubs table."""
    variants = {"all_channels": [], "no_textmining": ["automated_textmining"],
                "no_textmining_no_database": ["automated_textmining", "database_annotated"]}
    out = {}
    for k, drop in variants.items():
        sc = edges.apply(lambda r: combined_without(r, drop), axis=1) if drop else edges["combined_score"]
        ek = edges[sc >= threshold].copy(); ek["score"] = sc[sc >= threshold]
        out[k], _ = hubs_from_edges(ek, top, "node1", "node2")
    return out

# ------------------------------------------------------------------------------------------------ figures
def fig_network(G, hubs, path, top, status=None, highlight=(), tm_dependent=None, stable=()):
    """Publication figure (double-column, 7.2 in). Panel A: largest component with its two modularity communities;
    panel B: smaller components. Saves <base>.png (600 dpi), .pdf and .svg (text stays editable).
      status        node -> 'Up'/'Down' in LFS (Okabe-Ito vermilion / blue, colour-blind safe)
      tm_dependent  set of frozenset({a, b}) edges that fall below the threshold once text-mining evidence is removed (light grey)
      stable        genes that stay top-`top` hubs when text-mining evidence is removed (asterisk)"""
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Polygon
    from networkx.algorithms import community as nxc
    from scipy.spatial import ConvexHull
    plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "svg.fonttype": "none", "font.size": 6})
    UP, DOWN = "#D55E00", "#0072B2"
    tm_dependent, stable = tm_dependent or set(), set(stable)
    comps = sorted((c for c in nx.connected_components(G) if len(c) > 1), key=len, reverse=True)
    n_iso = sum(1 for n in G if G.degree(n) == 0)
    hubset = set(hubs.gene[:top])
    colour = lambda n: {"Up": UP, "Down": DOWN}.get((status or {}).get(n), "#999999")
    lab = lambda n: n + ("*" if n in stable else "")

    fig = plt.figure(figsize=(7.2, 5.3))
    a1 = fig.add_axes([0.005, 0.135, 0.60, 0.845]); a2 = fig.add_axes([0.62, 0.135, 0.375, 0.845])

    def draw(ax, H, pos, fs, node_scale):
        e_tm = [(u, v) for u, v in H.edges if frozenset((u, v)) in tm_dependent]
        e_ok = [(u, v) for u, v in H.edges if frozenset((u, v)) not in tm_dependent]
        w = lambda u, v: 0.18 + 0.6 * (H[u][v].get("weight", 0.7) - 0.7) / 0.3
        nx.draw_networkx_edges(H, pos, edgelist=e_tm, width=[w(u, v) for u, v in e_tm], edge_color="#cfcfcf", ax=ax)
        nx.draw_networkx_edges(H, pos, edgelist=e_ok, width=[w(u, v) for u, v in e_ok], edge_color="#5a5a5a" if tm_dependent else "#9a9a9a", ax=ax)
        deg = dict(H.degree)
        nx.draw_networkx_nodes(H, pos, node_size=[node_scale * (1 + 0.55 * deg[n]) for n in H], node_color=[colour(n) for n in H], ax=ax,
                               linewidths=[1.3 if n in hubset else 0.35 for n in H],
                               edgecolors=["black" if n in hubset else "white" for n in H])
        import matplotlib.patheffects as pe
        xs = np.array([pos[n][0] for n in H]); ys = np.array([pos[n][1] for n in H])
        cx, cy = xs.mean(), ys.mean()
        texts = []
        for n in H:
            dx, dy = pos[n][0] - cx, pos[n][1] - cy; d = max(np.hypot(dx, dy), 1e-6)
            off = 0.045 if H.number_of_nodes() > 12 else 0.03
            texts.append(ax.text(pos[n][0] + off * dx / d, pos[n][1] + off * dy / d, lab(n), fontsize=fs, ha="center", va="center",
                                 fontweight="bold" if n in hubset else "normal", zorder=5,
                                 path_effects=[pe.withStroke(linewidth=1.6, foreground="white")]))
        try:
            from adjustText import adjust_text
            adjust_text(texts, x=xs, y=ys, ax=ax, expand=(1.2, 1.4), force_text=(0.4, 0.6), force_static=(0.6, 0.8))
        except ImportError:
            print("adjustText not installed: labels are only offset from node centres (pip install adjustText for better placement)")
        ax.set_axis_off(); ax.margins(0.04)

    big = G.subgraph(comps[0])
    pos = nx.kamada_kawai_layout(big, weight=None)
    draw(a1, big, pos, 5.6, 9)
    mods = sorted(nxc.greedy_modularity_communities(big, weight="weight"), key=lambda c: -np.mean([pos[n][1] for n in c]))
    names = {"inv": "SPI-1 type III secretion", "fli": "Flagella and chemotaxis"}
    for c in mods:
        key = "inv" if sum(g[:3] in ("inv", "prg", "sip", "hil") for g in c) >= 3 else ("fli" if sum(g[:3] in ("fli", "flg", "che") for g in c) >= 3 else None)
        pts = np.array([pos[n] for n in c]); r = 0.09
        ang = np.linspace(0, 2 * np.pi, 16, endpoint=False)
        cloud = np.vstack([p + r * np.c_[np.cos(ang), np.sin(ang)] for p in pts])
        hull = cloud[ConvexHull(cloud).vertices]
        a1.add_patch(Polygon(hull, closed=True, fc="#eaeaea", ec="none", alpha=.55, zorder=0))
        if key:                                  # caption on the outer side of its group (away from the other group)
            upper = pts[:, 1].mean() >= np.mean([pos[n][1] for n in big])
            a1.text(pts[:, 0].min(), pts[:, 1].max() + 0.07 if upper else pts[:, 1].min() - 0.09, names[key],
                    fontsize=7, style="italic", color="#333333", ha="left", va="bottom" if upper else "top")
    a1.text(0.0, 1.0, "A", transform=a1.transAxes, fontsize=10, fontweight="bold", va="top")

    rest = comps[1:]
    if rest:
        pos2, ncol = {}, 3
        for k, c in enumerate(rest):
            lp = nx.spring_layout(G.subgraph(c), seed=2, scale=.36)
            for n, (x, y) in lp.items():
                pos2[n] = (x + (k % ncol) * 1.0, y - (k // ncol) * 0.9)
        draw(a2, G.subgraph(set().union(*rest)), pos2, 5, 7)
    a2.text(0.0, 1.0, "B", transform=a2.transAxes, fontsize=10, fontweight="bold", va="top")
    a2.text(0.12, 0.995, f"Smaller components (n = {len(rest)})", transform=a2.transAxes, fontsize=7, style="italic", color="#333333", va="top")

    handles = [Line2D([], [], marker="o", ls="", mfc=UP, mec="white", ms=6, label="up in LFS"),
               Line2D([], [], marker="o", ls="", mfc=DOWN, mec="white", ms=6, label="down in LFS"),
               Line2D([], [], marker="o", ls="", mfc="#bbbbbb", mec="black", mew=1.3, ms=6, label=f"top-{top} hub"),
               Line2D([], [], color="#5a5a5a", lw=1.2, label="edge, also supported without text-mining"),
               Line2D([], [], color="#cfcfcf", lw=1.2, label="edge, depends on text-mining")]
    fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False, fontsize=6, bbox_to_anchor=(0.5, 0.05), columnspacing=1.2, handletextpad=.4)
    fig.text(0.5, 0.020, f"Node size = degree; edge width = STRING combined score (0.7-1.0); * = stays a top-{top} hub when text-mining evidence is removed.\n"
             f"{n_iso} genes without a partner are not shown. Associations indicate functional linkage, not proof of physical binding.",
             ha="center", va="center", fontsize=5.4, color="#333333", linespacing=1.5)
    base = os.path.splitext(path)[0]
    for ext, kw in (("png", dict(dpi=600)), ("pdf", {}), ("svg", {})):
        fig.savefig(f"{base}.{ext}", **kw)
    plt.close(fig)


def fig_structures(models, path):
    """models: list of (label, xyz, plddt). C-alpha traces coloured by AlphaFold confidence bands."""
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    bands = [(90, "#0053d6"), (70, "#65cbf3"), (50, "#ffdb13"), (-1, "#ff7d45")]
    col = lambda v: next(c for t, c in bands if v > t)
    n = len(models); cols = min(5, n); rows = int(np.ceil(n / cols))
    fig = plt.figure(figsize=(3.2 * cols, 3.4 * rows))
    for i, (label, xyz, pl) in enumerate(models):
        ax = fig.add_subplot(rows, cols, i + 1, projection="3d")
        for j in range(len(xyz) - 1):
            ax.plot(*zip(xyz[j], xyz[j + 1]), color=col(pl[j]), lw=1.4)
        ax.set_title(f"{label}  (mean pLDDT {pl.mean():.0f})\n{100*np.mean(pl<70):.0f}% of residues < 70", fontsize=9, pad=2); ax.set_axis_off()
    fig.text(.5, .01, "AlphaFold C-alpha traces: blue >90, light blue 70-90, yellow 50-70, orange <50 (pLDDT). Predictions, not experimental structures.",
             ha="center", fontsize=8)
    fig.tight_layout(rect=(0, .03, 1, .97)); fig.savefig(path, dpi=300); plt.close(fig)


def fig_plddt_bands(models, path, uniprot=None):
    """Stacked horizontal bars: share of residues in each AlphaFold confidence band, one bar per protein."""
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "svg.fonttype": "none"})
    cut = [(">90 very high", lambda v: v > 90, "#0053D6"), ("70-90 confident", lambda v: (v > 70) & (v <= 90), "#65CBF3"),
           ("50-70 low", lambda v: (v > 50) & (v <= 70), "#FFDB13"), ("<50 very low", lambda v: v <= 50, "#FF7D45")]
    labels = [m[0] + (f" ({uniprot[m[0]]})" if uniprot and m[0] in uniprot else "") for m in models]
    fig, ax = plt.subplots(figsize=(7.2, 0.33 * len(models) + 1.3))
    left = np.zeros(len(models))
    for name, f, c in cut:
        w = np.array([100 * np.mean(f(m[2])) for m in models])
        ax.barh(range(len(models)), w, left=left, color=c, label=name, edgecolor="white", linewidth=.4); left += w
    ax.set_yticks(range(len(models))); ax.set_yticklabels(labels, fontsize=7); ax.invert_yaxis()
    for i, m in enumerate(models):
        ax.text(101, i, f"{m[2].mean():.0f}  ({len(m[2])} aa)", va="center", fontsize=6.5)
    ax.set_xlim(0, 100); ax.set_xlabel("Residues (%) by AlphaFold pLDDT band", fontsize=8)
    ax.text(101, -0.9, "mean pLDDT (length)", fontsize=6.5, style="italic")
    ax.legend(ncol=4, fontsize=6.5, frameon=False, loc="upper center", bbox_to_anchor=(0.45, -0.13))
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    fig.tight_layout()
    base = os.path.splitext(path)[0]
    for ext, kw in (("png", dict(dpi=600)), ("pdf", {}), ("svg", {})):
        fig.savefig(f"{base}.{ext}", **kw)
    plt.close(fig)


# ------------------------------------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", default="results/tables/high_confidence_DEGs_LFS_effect_hfq.csv")
    ap.add_argument("--coexp_hubs", default="results/tables/coexpression_hubs_LFS_effect_hfq.csv")
    ap.add_argument("--string_edges", help="STRING 'short tabular text output' TSV (offline mode; no STRING API call)")
    ap.add_argument("--string_annotations", help="STRING 'protein annotations' TSV (needed with --string_edges for alias matching)")
    ap.add_argument("--score", type=int, default=700, help="STRING combined score threshold (API mode only)")
    ap.add_argument("--top", type=int, default=10)
    ap.add_argument("--uniprot_map", help="UniProt ID-mapping TSV (From=STRING id, Entry=accession); avoids UniProt API calls")
    ap.add_argument("--structures_dir", help="folder with AlphaFold PDB files you downloaded (any name containing the UniProt accession); avoids AlphaFold API calls")
    ap.add_argument("--skip_structures", action="store_true")
    a = ap.parse_args()
    T, F, S = "results/tables", "results/figures", "results/structures"
    for d in (T, F, S):
        os.makedirs(d, exist_ok=True)

    deg = pd.read_csv(a.genes).drop_duplicates("label")
    if a.string_edges:
        e = pd.read_csv(a.string_edges, sep="\t").rename(columns={"#node1": "node1"})
        ann = pd.read_csv(a.string_annotations, sep="\t").rename(columns={"#node": "node"})
        mp = map_inputs_to_nodes(deg.label.astype(str), ann)
        thr = e.combined_score.min()
        print(f"offline mode: {len(e)} edges, minimum combined score in file = {thr:.3f} (STRING threshold actually used)")
    else:
        ids = string_resolve(deg.label.astype(str).tolist())
        ed = string_network(ids.stringId.dropna().tolist(), a.score)
        e = ed.rename(columns={"preferredName_A": "node1", "preferredName_B": "node2", "score": "combined_score"})
        mp = ids.rename(columns={"gene": "input_gene", "preferredName": "string_node"}).assign(match="exact")
        for c in STRING_CHANNELS:
            e[c] = np.nan
    mp = mp.merge(deg[["label", "GeneID", "log2FC", "FDR", "status", "description"]], left_on="input_gene", right_on="label").drop(columns="label")
    mp.to_csv(f"{T}/ppi_string_ids.csv", index=False)
    print(f"STRING matched {mp.string_node.notna().sum()}/{len(mp)} genes; not found: {mp[mp.string_node.isna()].input_gene.tolist()}; "
          f"matched by alias: {mp[mp.match=='alias'][['input_gene','string_node']].values.tolist()}")
    n2 = mp.dropna(subset=["string_node"]).set_index("string_node")
    e.to_csv(f"{T}/ppi_edges.csv", index=False)

    ee = e.assign(preferredName_A=e.node1, preferredName_B=e.node2, score=e.combined_score)
    hubs, G = hubs_from_edges(ee, a.top)
    for n in n2.index:
        G.add_node(n)
    hubs["input_gene"] = hubs.gene.map(n2.input_gene)
    for c in ["log2FC", "FDR", "status", "description"]:
        hubs[c] = hubs.gene.map(n2[c])
    co = pd.read_csv(a.coexp_hubs) if os.path.exists(a.coexp_hubs) else pd.DataFrame(columns=["label", "hub"])
    co_hubs = set(co[co.hub].label) if "hub" in co else set(co.label)
    if a.string_edges:
        rob = hub_robustness(e, a.top)
        for k, v in rob.items():
            hubs["rank_" + k] = hubs.gene.map(pd.Series(range(1, len(v) + 1), index=v.gene))
        rk = [c for c in hubs if c.startswith("rank_")]
        hubs["stable_hub"] = hubs[rk].le(a.top).all(axis=1)
        ch = {c.split("_")[0]: int((e[c] >= .4).sum()) for c in STRING_CHANNELS}
        print("edges supported (channel score >= 0.4):", ch, "of", len(e))
    hubs["coexpression_hub"] = hubs.input_gene.isin(co_hubs)
    hubs.to_csv(f"{T}/ppi_hubs.csv", index=False)
    comps = sorted((len(c) for c in nx.connected_components(G)), reverse=True)
    print(f"network: {G.number_of_nodes()} genes, {G.number_of_edges()} edges, {comps.count(1)} without a partner, largest component {comps[0]}")
    print("top hubs:", hubs.gene[:a.top].tolist())
    tm_dep, stable = set(), set()
    if a.string_edges:
        sc = e.apply(lambda r: combined_without(r, ["automated_textmining"]), axis=1)
        tm_dep = {frozenset((u, v)) for u, v in zip(e.node1[sc < 0.7], e.node2[sc < 0.7])}
        stable = set(hubs.gene[hubs.stable_hub])
        print(f"{len(tm_dep)} of {len(e)} edges fall below 0.7 without text-mining")
    for u, v, sc_ in zip(e.node1, e.node2, e.combined_score):
        G[u][v]["weight"] = sc_
    fig_network(G, hubs, f"{F}/Fig14_PPI_network_STRING.png", a.top, status=n2.status.to_dict(),
                highlight=[g for g in n2.index if n2.input_gene[g] in co_hubs], tm_dependent=tm_dep, stable=stable)

    targets = list(dict.fromkeys(list(hubs.gene[:a.top]) +
                                 ([g for g in hubs.gene if hubs.set_index("gene").stable_hub[g]] if "stable_hub" in hubs else []) +
                                 [n for n in n2.index if n2.input_gene[n] in co_hubs]))
    hubs["structure_target"] = hubs.gene.isin(targets)
    hubs.to_csv(f"{T}/ppi_hubs.csv", index=False)
    print("structure targets:", targets)
    if a.skip_structures:
        return
    rows, models = [], []
    acc_of = {}
    if a.uniprot_map:
        um = pd.read_csv(a.uniprot_map, sep="\t")
        sid2 = um.set_index("From")
        acc_of = {r.node: r.identifier for r in ann.itertuples()} if a.string_edges else {}
        acc_of = {g: (sid2.Entry.get(sid), sid2["Protein names"].get(sid), sid2.Reviewed.get(sid) == "reviewed") for g, sid in acc_of.items() if sid in sid2.index}
    for g in targets:
        if g in acc_of:
            acc, prot, rev = acc_of[g]
            u = dict(uniprot=acc, reviewed=rev, protein=prot)
        else:
            names = [g] + ([n2.input_gene[g]] if n2.input_gene.get(g) not in (None, g) else [])    # STRING name, then our own label
            u = next((x for x in (uniprot_lookup(n) for n in names) if x), None)
        if not u:
            rows.append(dict(gene=g, status="no UniProt entry found - give --uniprot_map or search by locus tag manually")); continue
        if a.structures_dir:
            path = find_structure_file(a.structures_dir, u["uniprot"])
            if not path:
                rows.append(dict(gene=g, uniprot=u["uniprot"], protein=u["protein"], status=f"PDB for {u['uniprot']} not found in {a.structures_dir}")); continue
            pdb, af = open(path).read(), {"entryId": os.path.basename(path)}
        else:
            af = alphafold_entry(u["uniprot"])
            if not af:
                rows.append(dict(gene=g, uniprot=u["uniprot"], status="no AlphaFold DB entry")); continue
            pdb = _get(af["pdbUrl"])
            path = f"{S}/{g}_{u['uniprot']}.pdb"
            open(path, "w").write(pdb)
        xyz, pl = plddt_from_pdb(pdb)
        rows.append(dict(gene=g, uniprot=u["uniprot"], uniprot_reviewed=u["reviewed"], protein=u["protein"], af_model=af.get("entryId"),
                         file=path, status="ok", **plddt_summary(pl)))
        models.append((g, xyz, pl))
    pd.DataFrame(rows).to_csv(f"{T}/hub_structures.csv", index=False)
    if models:
        fig_structures(models, f"{F}/Fig15_hub_gene_structures_AlphaFold.png")
        fig_plddt_bands(models, f"{F}/Fig16_hub_pLDDT_confidence_bands.png", {r["gene"]: r["uniprot"] for r in rows if r.get("status") == "ok"})
    print(pd.DataFrame(rows).drop(columns=["file"], errors="ignore").to_string(index=False))


if __name__ == "__main__":
    main()
