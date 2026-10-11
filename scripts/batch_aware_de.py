#!/usr/bin/env python3
"""Batch-aware differential expression and dropout diagnostics for raw counts.

Inputs
  --counts   genes x samples RAW integer counts (e.g. from featureCounts via rebuild_counts_from_sra.sh)
  --meta     CSV with columns: sample, genotype, culture, batch   (batch = the *documented* batch; see
             data/sra/column_to_run_map.csv / batch_evidence.csv from collect_batch_evidence.py)
  --outdir   output directory

What it does
  1. Design report: rank, estimability and variance inflation of every genotype/culture contrast with
     and without the batch term. Shows how many samples *actually* carry each contrast once batch is modelled.
  2. PyDESeq2 with  ~ group  (M0)  and  ~ batch + group  (M1), four contrasts each.
  3. Concordance table: which M0 hits survive batch adjustment (and keep sign).
  4. Dropout diagnosis (counts only): for every gene x group, is an all-zero observation explained by
     sequencing depth ("depth_limited_uninformative": the gene could be as high as in its best group and
     still give these zeros) or not ("zero_not_explained_by_depth": expression is genuinely lower than in
     the best group - reduced or absent; the call does not distinguish those two)?
     Also reports a 95 % *upper detection limit* of expression for each group.

What it can NOT do: create reads that were never sequenced. A zero in a low-depth sample is turned into an
honest bound ("expression could be as high as X"), not into a measured value.
"""
import argparse, json, os, warnings
import numpy as np, pandas as pd
from scipy import optimize

warnings.filterwarnings("ignore")

CONTRASTS = {
    "hfq_vs_WT_control": ("delta_hfq_control", "WT_control"),
    "hfq_vs_WT_LFS": ("delta_hfq_LFS", "WT_LFS"),
    "LFS_vs_control_WT": ("WT_LFS", "WT_control"),
    "LFS_vs_control_hfq": ("delta_hfq_LFS", "delta_hfq_control"),
}


# --------------------------------------------------------------------------------------------- design
def design_matrix(meta, with_batch):
    groups = ["WT_control", "delta_hfq_control", "WT_LFS", "delta_hfq_LFS"]
    X = pd.DataFrame({g: (meta["group"] == g).astype(float) for g in groups}, index=meta.index)
    if with_batch:
        for b in sorted(meta["batch"].astype(str).unique())[1:]:  # first level = reference
            X["batch_" + b] = (meta["batch"].astype(str) == b).astype(float)
    return X


def design_report(meta):
    rows = []
    for with_batch in (False, True):
        X = design_matrix(meta, with_batch)
        rank, p = np.linalg.matrix_rank(X.values), X.shape[1]
        XtXi = np.linalg.pinv(X.values.T @ X.values)
        for name, (a, b) in CONTRASTS.items():
            c = np.zeros(p)
            c[list(X.columns).index(a)], c[list(X.columns).index(b)] = 1, -1
            estimable = np.linalg.matrix_rank(np.vstack([X.values, c])) == rank
            rows.append(dict(model="with_batch" if with_batch else "no_batch", contrast=name, rank=f"{rank}/{p}",
                             estimable=bool(estimable), var_factor=float(c @ XtXi @ c) if estimable else np.nan))
    d = pd.DataFrame(rows)
    base = d[d.model == "no_batch"].set_index("contrast")["var_factor"]
    d["se_inflation_vs_no_batch"] = d.apply(lambda r: np.sqrt(r.var_factor / base[r.contrast]), axis=1)
    d["effective_n_per_group"] = d.apply(lambda r: 3 / (r.var_factor / base[r.contrast]) if r.estimable else np.nan, axis=1)
    return d


# --------------------------------------------------------------------------------------------- DE
def fit(counts, meta, design):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    dds = DeseqDataSet(counts=counts.T, metadata=meta, design=design, refit_cooks=True, quiet=True)
    dds.deseq2()
    out = {}
    for name, (a, b) in CONTRASTS.items():
        st = DeseqStats(dds, contrast=["group", a, b], quiet=True)
        st.summary()
        r = st.results_df.copy()
        r["contrast"] = name
        out[name] = r
    return dds, out


# --------------------------------------------------------------------------------------------- dropout
def nb_p0(mu, alpha):
    """P(count = 0) for NB with mean mu and dispersion alpha (alpha -> 0 gives Poisson)."""
    alpha = np.maximum(alpha, 1e-8)
    return (1.0 + alpha * mu) ** (-1.0 / alpha)


def upper_limit(sf, alpha, level=0.05):
    """Largest normalised expression lambda still compatible (>= level) with all-zero counts in these samples."""
    f = lambda lam: np.sum(-np.log1p(alpha * lam * sf) / alpha) - np.log(level)
    hi = 1.0
    while f(hi) > 0 and hi < 1e7:
        hi *= 2
    return float(optimize.brentq(f, 0.0, hi))


def dropout_diagnosis(dds, meta, depth_alpha=0.05):
    sf = dds.obs["size_factors"].values
    norm = pd.DataFrame(dds.layers["normed_counts"], index=meta.index, columns=dds.var_names)
    raw = pd.DataFrame(dds.X, index=meta.index, columns=dds.var_names)
    alpha_g = pd.Series(dds.var["dispersions"].values, index=dds.var_names)
    trend = pd.Series(dds.var["fitted_dispersions"].values, index=dds.var_names) if "fitted_dispersions" in dds.var else alpha_g
    alpha_g = alpha_g.fillna(trend).fillna(np.nanmedian(trend.values)).clip(lower=1e-3)
    groups = meta["group"]
    gmeans = norm.groupby(groups).mean()                       # group x gene
    top_group = gmeans.idxmax()                                # group where gene is most expressed
    top_mean = gmeans.max()
    rows = []
    for g in dds.var_names:
        for grp in gmeans.index:
            idx = (groups == grp).values
            counts_g, sf_g, a = raw.loc[idx, g].values, sf[idx], alpha_g[g]
            n0 = int((counts_g == 0).sum())
            # probability of seeing this many zeros if the gene were expressed at its best group's level
            p0 = nb_p0(top_mean[g] * sf_g, a)
            p_all_zero = float(np.prod(p0[counts_g == 0])) * float(np.prod(1 - p0[counts_g > 0])) if n0 else 1.0
            ul = upper_limit(sf_g, a) if counts_g.sum() == 0 else np.nan
            rows.append((g, grp, int(idx.sum()), n0, float(norm.loc[idx, g].mean()), top_group[g], float(top_mean[g]),
                         p_all_zero, ul))
    d = pd.DataFrame(rows, columns=["gene", "group", "n", "n_zero", "mean_norm", "best_group", "best_group_mean",
                                    "p_obs_pattern_if_expressed_like_best", "upper95_norm_expr_if_all_zero"])

    def call(r):
        if r.n_zero == 0:
            return "detected"
        if r.group == r.best_group:
            return "zeros_in_best_group" if r.n_zero < r.n else "all_zero_in_best_group"
        if r.n_zero == r.n:
            return "depth_limited_uninformative" if r.p_obs_pattern_if_expressed_like_best >= 0.05 else "zero_not_explained_by_depth"
        return "partial_detection"
    d["call"] = d.apply(call, axis=1)
    return d


# --------------------------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", required=True)
    ap.add_argument("--meta", required=True)
    ap.add_argument("--outdir", default="results/count_based")
    ap.add_argument("--min_total", type=int, default=10, help="drop genes with fewer total counts")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    counts = pd.read_csv(a.counts, index_col=0)
    meta = pd.read_csv(a.meta).set_index("sample").loc[counts.columns]
    meta["group"] = meta["genotype"].astype(str) + "_" + meta["culture"].astype(str)
    meta["group"] = meta["group"].replace({"WT_control": "WT_control", "WT_LFS": "WT_LFS"})
    meta["batch"] = meta["batch"].astype(str)
    counts = counts.loc[counts.sum(axis=1) >= a.min_total].astype(int)

    rep = design_report(meta)
    rep.to_csv(f"{a.outdir}/design_report.csv", index=False)
    print(rep.round(3).to_string(index=False))

    dds0, res0 = fit(counts, meta, "~ group")
    pd.concat(res0.values()).to_csv(f"{a.outdir}/de_no_batch.csv")
    summary = {"M0_no_batch": {k: int(((v.padj < .05) & (v.log2FoldChange.abs() >= 1)).sum()) for k, v in res0.items()}}

    if meta["batch"].nunique() > 1 and np.linalg.matrix_rank(design_matrix(meta, True).values) == design_matrix(meta, True).shape[1]:
        dds1, res1 = fit(counts, meta, "~ batch + group")
        pd.concat(res1.values()).to_csv(f"{a.outdir}/de_batch_adjusted.csv")
        summary["M1_batch_adjusted"] = {k: int(((v.padj < .05) & (v.log2FoldChange.abs() >= 1)).sum()) for k, v in res1.items()}
        conc = []
        for k in CONTRASTS:
            m = res0[k].join(res1[k], lsuffix="_M0", rsuffix="_M1")
            s0 = (m.padj_M0 < .05) & (m.log2FoldChange_M0.abs() >= 1)
            s1 = (m.padj_M1 < .05) & (m.log2FoldChange_M1.abs() >= 1) & (np.sign(m.log2FoldChange_M0) == np.sign(m.log2FoldChange_M1))
            conc.append(dict(contrast=k, hits_M0=int(s0.sum()), hits_M0_surviving_batch_adjustment=int((s0 & s1).sum())))
        pd.DataFrame(conc).to_csv(f"{a.outdir}/batch_concordance.csv", index=False)
        print(pd.DataFrame(conc).to_string(index=False))
    else:
        print("Batch term not fitted: single batch level or aliased with group (see design_report.csv).")

    d = dropout_diagnosis(dds0, meta)
    d.to_csv(f"{a.outdir}/dropout_diagnosis.csv", index=False)
    print(d["call"].value_counts().to_string())
    json.dump(summary, open(f"{a.outdir}/summary.json", "w"), indent=1)


if __name__ == "__main__":
    main()
