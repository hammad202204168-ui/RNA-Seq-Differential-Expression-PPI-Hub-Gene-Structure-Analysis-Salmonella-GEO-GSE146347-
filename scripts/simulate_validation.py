#!/usr/bin/env python3
"""Validate batch_aware_de.py against SIMULATED counts with known ground truth.

The simulation reproduces the real design (12 samples, 3/3/3/3, `_2` batch on 5 of 6 delta-hfq samples,
per-group read depths from data/sra/SraRunTable_GSE146347.csv). Nothing here is real expression data; it
tests what the statistics can and cannot recover, so the conclusions drawn from real counts can be calibrated.

Outputs: results/tables/sim_validation_*.csv
"""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
import batch_aware_de as B

OUT = "results/tables"
os.makedirs(OUT, exist_ok=True)

SAMPLES = ["WT2_1xg", "WT3_1xg", "WT4_1xg", "hfq1_1xg", "hfq3_1xg_2", "hfq4_1xg_2",
           "WT2_LS", "WT3_LS", "WT4_LS", "hfq1_LS_2", "hfq3_LS_2", "hfq4_LS_2"]


def make_meta():
    m = pd.DataFrame({"sample": SAMPLES})
    m["genotype"] = np.where(m["sample"].str.startswith("hfq"), "delta_hfq", "WT")
    m["culture"] = np.where(m["sample"].str.contains("_LS"), "LFS", "control")
    m["batch"] = np.where(m["sample"].str.endswith("_2"), "b2", "b1")
    return m.set_index("sample", drop=False)


def depths(meta, rng):
    """Reads per sample from the real SRA table (Bases / ~73 nt); order within group is arbitrary because the
    TPM columns were never linked to SRR runs."""
    t = pd.read_csv("data/sra/SraRunTable_GSE146347.csv")
    t = t[t.Run.between("SRR11235557", "SRR11235568")]
    t["grp"] = np.where(t.bacterial_genotype == "wild type", "WT", "delta_hfq") + "_" + t.bacterial_growth_orientation.replace({"LFS": "LFS", "control": "control"})
    out = {}
    for g, sub in t.groupby("grp"):
        reads = sub["Bases"].values / 73.0 * 0.7          # ~70 % of reads in genes (rRNA/other lost)
        members = meta.index[(meta.genotype + "_" + meta.culture) == g]
        out.update(dict(zip(members, rng.permutation(reads))))
    return pd.Series(out)[meta.index]


def simulate(seed, G=3000):
    rng = np.random.default_rng(seed)
    meta = make_meta()
    lib = depths(meta, rng)
    p = np.exp(rng.normal(0, 1.9, G)); p /= p.sum()                       # relative abundance
    lfc = lambda frac, sd: np.where(rng.random(G) < frac, rng.normal(0, sd, G) + np.sign(rng.normal(size=G)) * 0.8, 0.0)
    eff_culture, eff_geno, eff_batch = lfc(.10, .5), lfc(.04, .5), lfc(.12, .6)
    absent_in = np.where(rng.random(G) < .03, rng.choice(["WT_control", "WT_LFS"], G), "")   # truly off in one WT group
    alpha = np.exp(rng.normal(np.log(0.12), 0.35, G))
    mu = np.zeros((G, len(meta)))
    for j, s in enumerate(meta.index):
        lf = (meta.loc[s, "culture"] == "LFS") * eff_culture + (meta.loc[s, "genotype"] == "delta_hfq") * eff_geno \
             + (meta.loc[s, "batch"] == "b2") * eff_batch
        m = lib[s] * p * 2.0 ** lf
        grp = meta.loc[s, "genotype"] + "_" + meta.loc[s, "culture"]
        m = np.where(absent_in == grp, m * 1e-4, m)
        mu[:, j] = m
    r = 1 / alpha[:, None]
    counts = rng.negative_binomial(r, r / (r + mu))
    counts = pd.DataFrame(counts, index=[f"g{i}" for i in range(G)], columns=meta.index)
    truth = pd.DataFrame(dict(eff_culture=eff_culture, eff_geno=eff_geno, eff_batch=eff_batch, absent_in=absent_in,
                              p=p, alpha=alpha), index=counts.index)
    return counts, meta, truth, lib


def evaluate(seed):
    counts, meta, truth, lib = simulate(seed)
    meta["group"] = meta.genotype + "_" + meta.culture
    keep = counts.sum(axis=1) >= 10
    counts, truth = counts[keep], truth[keep]
    _, r0 = B.fit(counts, meta, "~ group")
    _, r1 = B.fit(counts, meta, "~ batch + group")
    dds0, _ = B.fit(counts, meta, "~ group")
    rows = []
    for name, res in (("M0_no_batch", r0), ("M1_batch_adjusted", r1)):
        r = res["hfq_vs_WT_control"]
        called = (r.padj < .05) & (r.log2FoldChange.abs() >= 1)
        true_geno = truth.eff_geno.abs() >= 1
        batch_only = (truth.eff_batch.abs() >= 1) & (truth.eff_geno == 0) & (truth.absent_in == "")
        rows.append(dict(seed=seed, model=name, contrast="hfq_vs_WT_control",
                         called=int(called.sum()),
                         sensitivity_true_genotype_genes=float((called & true_geno).sum() / max(true_geno.sum(), 1)),
                         false_calls_from_batch_only_genes=int((called & batch_only).sum()),
                         fp_rate_batch_only=float((called & batch_only).sum() / max(batch_only.sum(), 1)),
                         median_lfc_SE=float(r.lfcSE.median()),
                         lfc_corr_with_truth=float(np.corrcoef(r.log2FoldChange.fillna(0), truth.eff_geno)[0, 1])))
    # --- dropout diagnostic
    dm = B.dropout_diagnosis(dds0, meta)
    dm = dm.merge(truth[["p", "eff_culture", "eff_geno", "eff_batch", "absent_in"]], left_on="gene", right_index=True)
    sf_grp = {g: lib[(meta.group == g).values].mean() for g in meta.group.unique()}
    z = dm[(dm.n_zero == dm.n) & (dm.group != dm.best_group)].copy()          # all-zero cell not in the gene's best group
    z["truly_absent"] = z.apply(lambda r: r.absent_in == r.group, axis=1)
    z["call_absent"] = z.call == "zero_not_explained_by_depth"
    tab = pd.crosstab(z.truly_absent, z.call)
    # what does "zero_not_explained_by_depth" mean in truth? true expression of that cell relative to the gene's best group
    gm0 = np.exp(np.mean(np.log(lib.values)))
    rep = {g: meta.index[(meta.group == g).values][0] for g in meta.group.unique()}
    def tn(gene_row, grp):
        s_ = rep[grp]
        lf = (meta.loc[s_, "culture"] == "LFS") * gene_row.eff_culture + (meta.loc[s_, "genotype"] == "delta_hfq") * gene_row.eff_geno \
             + (meta.loc[s_, "batch"] == "b2") * gene_row.eff_batch
        v = gm0 * gene_row.p * 2.0 ** lf
        return v * (1e-4 if gene_row.absent_in == grp else 1)
    z["true_ratio_to_best"] = [tn(truth.loc[g], grp) / max(tn(truth.loc[g], bg), 1e-12) for g, grp, bg in zip(z.gene, z.group, z.best_group)]
    flagged = z[z.call_absent]
    reduced_frac = float((flagged.true_ratio_to_best < 0.5).mean()) if len(flagged) else np.nan
    # upper-limit coverage for all-zero cells of expressed genes (true normalised expr should be <= limit in >=95 %)
    ze = z[(~z.truly_absent) & z.upper95_norm_expr_if_all_zero.notna()].copy()
    # true expression on the DESeq2 normalised scale ~ p * mean_lib * 2^effects (size factors ~ lib/geomean)
    gm = np.exp(np.mean(np.log(lib.values)))
    def true_norm(r):
        s = meta.index[(meta.group == r.group).values][0]
        lf = (meta.loc[s, "culture"] == "LFS") * r.eff_culture + (meta.loc[s, "genotype"] == "delta_hfq") * r.eff_geno + \
             (meta.loc[s, "batch"] == "b2") * r.eff_batch
        return gm * r.p * 2.0 ** lf
    ze["true_norm"] = ze.apply(true_norm, axis=1)
    cov = float((ze.true_norm <= ze.upper95_norm_expr_if_all_zero).mean()) if len(ze) else np.nan
    drop = dict(seed=seed, all_zero_cells=len(z), truly_absent=int(z.truly_absent.sum()),
                absent_recall=float((z.truly_absent & z.call_absent).sum() / max(z.truly_absent.sum(), 1)),
                expressed_called_not_explained=int(((~z.truly_absent) & z.call_absent).sum()),
                flagged_cells_truly_below_half_of_best=reduced_frac,
                expressed_called_depth_limited=int(((~z.truly_absent) & (z.call == "depth_limited_uninformative")).sum()),
                upper95_coverage_expressed=cov, n_expressed_cells=len(ze))
    return pd.DataFrame(rows), pd.DataFrame([drop]), tab


if __name__ == "__main__":
    de, dr = [], []
    for seed in (1, 2, 3):
        a, b, tab = evaluate(seed)
        de.append(a); dr.append(b)
        print(f"seed {seed}\n{tab}\n")
    de, dr = pd.concat(de), pd.concat(dr)
    de.to_csv(f"{OUT}/sim_validation_de.csv", index=False)
    dr.to_csv(f"{OUT}/sim_validation_dropout.csv", index=False)
    pd.set_option("display.width", 250)
    print(de.groupby("model").mean(numeric_only=True).drop(columns="seed").round(3).T)
    print(dr.mean(numeric_only=True).drop("seed").round(3))
