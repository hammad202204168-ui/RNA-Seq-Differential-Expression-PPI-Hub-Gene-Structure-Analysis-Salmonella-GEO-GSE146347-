"""
Robustness safeguards for the TPM-based GSE146347 analysis.

The data have four known problems that cannot be removed by analysis (unequal rRNA depletion, shallow/unequal depth,
a possible library batch confounded with genotype, and genes expressed in only one group). This module therefore does
not "fix" them; it (1) flags affected genes, (2) re-tests every contrast under alternative models/filters, and
(3) defines a stricter, transparent set of high-confidence DEGs.

Analyses compared with the baseline model (pooled-variance t-test on log2(TPM+1), 4 group means):
    no_exclusive      refit after removing genes expressed in one group only
    detect_adj        + covariate: number of detected genes per sample (depth/complexity proxy)
    rrna_adj          + covariate: rRNA share of the original TPM per sample
    batch2_adj        + covariate: sample name ends in "_2" (nearly confounded with genotype; ill-conditioned)
    consistent        baseline call kept only if the higher group has TPM >= 1 in all 3 replicates and mean TPM >= 5
Robust DEG = baseline DEG that is not group-exclusive, passes `consistent`, and keeps the same direction with
FDR < 0.05 and |log2FC| >= 1 in `detect_adj`.
"""
import numpy as np
import pandas as pd
from scipy.stats import t as tdist

GROUPS = ["WT_control", "WT_LFS", "delta_hfq_control", "delta_hfq_LFS"]
ANALYSES = ["baseline", "no_exclusive", "detect_adj", "rrna_adj", "batch2_adj", "consistent", "robust"]


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    out = np.empty(n); out[o] = np.minimum(q, 1)
    return out


class CovModel:
    """Group-means model on log2 expression with optional z-scored sample covariates (OLS)."""

    def __init__(self, Y, meta, covars=None):
        cols = [(meta.group.values == g).astype(float) for g in GROUPS]
        if covars is not None:
            for c in covars.columns:
                v = covars[c].reindex(meta.index).astype(float).values
                cols.append((v - v.mean()) / (v.std() or 1.0))
        self.X = np.column_stack(cols)
        self.p, self.df = self.X.shape[1], len(meta) - self.X.shape[1]
        self.XtXi = np.linalg.inv(self.X.T @ self.X)
        self.B = Y.values @ self.X @ self.XtXi
        resid = Y.values - self.B @ self.X.T
        self.s2 = np.maximum((resid ** 2).sum(axis=1) / self.df, 1e-8)
        self.index = Y.index
        self.cond = np.linalg.cond(self.X)

    def test(self, weights):
        c = np.zeros(self.p)
        for g, w in weights.items():
            c[GROUPS.index(g)] = w
        est = self.B @ c
        se = np.sqrt(self.s2 * (c @ self.XtXi @ c))
        t = est / se
        p = 2 * tdist.sf(np.abs(t), self.df)
        return pd.DataFrame({"log2FC": est, "pvalue": p, "FDR": bh(p)}, index=self.index)


def sign_calls(df, fdr, lfc):
    sig = (df.FDR < fdr) & (df.log2FC.abs() >= lfc)
    return pd.Series(np.where(sig, np.sign(df.log2FC), 0), index=df.index).astype(int)


def exclusive_genes(tpm, meta):
    """Genes with mean TPM >= 5 in one group and TPM < 1 in every sample of the other three groups."""
    out = pd.Series(False, index=tpm.index)
    for g in GROUPS:
        inside = tpm[meta.index[meta.group == g]].mean(axis=1)
        others = tpm[meta.index[meta.group != g]].max(axis=1)
        out |= (inside >= 5) & (others < 1)
    return out


def consistent_detection(tpm, meta, a, b, log2fc):
    """Higher group of the contrast must have TPM >= 1 in all replicates and mean TPM >= 5."""
    hi_a = (tpm[meta.index[meta.group == a]] >= 1).all(axis=1) & (tpm[meta.index[meta.group == a]].mean(axis=1) >= 5)
    hi_b = (tpm[meta.index[meta.group == b]] >= 1).all(axis=1) & (tpm[meta.index[meta.group == b]].mean(axis=1) >= 5)
    return pd.Series(np.where(log2fc.reindex(tpm.index) > 0, hi_a, hi_b), index=tpm.index)


def technical_covariates(tpm, tpm_raw, ann, meta):
    typ = ann.loc[tpm_raw.index, "type"]
    return pd.DataFrame({
        "detected_genes": (tpm >= 1).sum(),
        "rrna_share": tpm_raw[typ == "rRNA"].sum() / 1e6,
        "batch2": [1.0 if s.endswith("_2") else 0.0 for s in meta.index],
    }).loc[meta.index]


def run(tpm, tpm_raw, ann, meta, Yk, results, contrasts, fdr, lfc):
    """Augment each contrast table with robustness columns; return (results, summary, extra_tables)."""
    cov = technical_covariates(tpm, tpm_raw, ann, meta)
    excl = exclusive_genes(tpm.loc[Yk.index], meta)
    models = {
        "detect_adj": CovModel(Yk, meta, cov[["detected_genes"]]),
        "rrna_adj": CovModel(Yk, meta, cov[["rrna_share"]]),
        "batch2_adj": CovModel(Yk, meta, cov[["batch2"]]),
    }
    Yne = Yk.loc[~excl]
    base_ne = CovModel(Yne, meta)
    summary, out = [], {}
    for k, (a, b) in contrasts.items():
        r = results[k].copy()
        w = {a: 1.0, b: -1.0}
        calls = {"baseline": sign_calls(r, fdr, lfc)}
        calls["no_exclusive"] = sign_calls(base_ne.test(w), fdr, lfc).reindex(r.index).fillna(0).astype(int)
        for name, m in models.items():
            calls[name] = sign_calls(m.test(w), fdr, lfc).reindex(r.index).fillna(0).astype(int)
        cons = consistent_detection(tpm.loc[r.index], meta, a, b, r.log2FC)
        calls["consistent"] = calls["baseline"].where(cons, 0)
        base = calls["baseline"]
        same = lambda c: (base != 0) & (c == base)
        robust = same(calls["no_exclusive"]) & same(calls["consistent"]) & same(calls["detect_adj"]) & ~excl.reindex(r.index).fillna(False)
        calls["robust"] = base.where(robust, 0)
        r["group_exclusive"] = excl.reindex(r.index).fillna(False)
        r["consistent_detection"] = cons
        for name in ("no_exclusive", "detect_adj", "rrna_adj", "batch2_adj"):
            r[f"call_{name}"] = calls[name]
        r["robust_DEG"] = calls["robust"] != 0
        out[k] = r
        nb = int((base != 0).sum())
        for name in ANALYSES:
            c = calls[name]
            kept = int(same(c).sum()) if name != "baseline" else nb
            summary.append(dict(contrast=k, analysis=name, n_up=int((c > 0).sum()), n_down=int((c < 0).sum()),
                                baseline_DEGs=nb, retained=kept, pct_retained=100 * kept / nb if nb else np.nan))
    extra = {
        "technical_covariates": cov.assign(group=meta.group),
        "model_conditioning": pd.DataFrame({"model": list(models), "condition_number": [m.cond for m in models.values()],
                                            "residual_df": [m.df for m in models.values()]}),
        "group_exclusive_genes": ann.loc[excl[excl].index, ["label", "replicon", "description"]],
    }
    return out, pd.DataFrame(summary), extra


def pc_association(scores, meta):
    """Share of variance of each PC explained (eta^2) by genotype, culture and the '_2' name suffix."""
    rows = []
    factors = {"genotype": meta.genotype, "culture": meta.culture, "suffix_2": pd.Series(meta.index.str.endswith("_2"), index=meta.index)}
    for pc in scores.columns:
        y = scores[pc].reindex(meta.index)
        tot = ((y - y.mean()) ** 2).sum()
        for f, v in factors.items():
            ssb = sum(len(g) * (g.mean() - y.mean()) ** 2 for _, g in y.groupby(v.values))
            rows.append(dict(PC=pc, factor=f, eta2=ssb / tot))
    return pd.DataFrame(rows)
