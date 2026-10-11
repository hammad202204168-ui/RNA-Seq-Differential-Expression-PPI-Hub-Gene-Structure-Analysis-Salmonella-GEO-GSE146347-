"""Regression tests for the GSE146347 TPM analysis (run after scripts/prepare_geo_tpm.py and real_data_analysis.py)."""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
T, P, F = ROOT / "results" / "tables", ROOT / "data" / "processed", ROOT / "results" / "figures"
CONTRASTS = ["hfq_effect_control", "hfq_effect_LFS", "LFS_effect_WT", "LFS_effect_hfq"]


def test_design_is_balanced_2x2_with_three_replicates():
    m = pd.read_csv(P / "sample_metadata.csv")
    assert len(m) == 12 and m.group.value_counts().eq(3).all() and m.group.nunique() == 4
    assert set(m.genotype) == {"WT", "delta_hfq"} and set(m.culture) == {"control", "LFS"}


def test_tpm_columns_sum_to_one_million_and_are_nonnegative():
    t = pd.read_csv(P / "tpm_matrix.csv", index_col=0)
    assert t.shape == (4743, 12) and (t >= 0).all().all()
    assert np.allclose(t.sum().values, 1e6, atol=5)


def test_renormalised_tpm_excludes_structural_rna():
    r = pd.read_csv(T / "tpm_renormalised_no_structural_RNA.csv", index_col=0)
    a = pd.read_csv(P / "gene_annotation.csv", index_col=0)
    assert not a.loc[r.index, "type"].isin(["rRNA", "tRNA", "ncRNA"]).any()
    assert np.allclose(r.sum().values, 1e6, atol=1)


def test_hfq_positive_control():
    d = pd.read_csv(T / "de_hfq_effect_control.csv", index_col=0)
    h = d[d.label == "hfq"].iloc[0]
    assert h.log2FC < -5 and h.FDR < 0.05 and h.status == "Down"          # deleted gene must drop
    assert h.meanTPM_WT_control > 100 and h.meanTPM_delta_hfq_control < 20


def test_contrast_tables_are_consistent():
    for c in CONTRASTS + ["interaction_genotype_x_culture"]:
        d = pd.read_csv(T / f"de_{c}.csv", index_col=0)
        assert len(d) == 4420 and d.FDR.between(0, 1).all() and d.pvalue.between(0, 1).all()
        assert (d.FDR >= d.pvalue - 1e-12).all()                        # BH-adjusted p is never below raw p
    d = pd.read_csv(T / "de_hfq_effect_control.csv")
    sig = d[d.status != "NS"]
    assert ((sig.FDR < 0.05) & (sig.log2FC.abs() >= 1)).all()


def test_figures_and_tables_exist():
    names = ["Fig0_graphical_summary", "Fig1_workflow_overview", "Fig2_quality_control", "Fig3_PCA_and_sample_distance",
             "Fig4_hfq_control_and_top_genes", "Fig5_heatmap_and_volcano", "Fig6_volcano_interpretation",
             "Fig7_all_contrast_volcanoes", "Fig8_pvalue_histograms_and_MA", "Fig9_LFS_effect_by_genotype",
             "Fig10_functional_enrichment", "Fig11_coexpression_network_and_hubs", "Fig12_volcano_reading_guide", "Fig13_robustness_checks"]
    for n in names:
        assert (F / f"{n}.png").stat().st_size > 10_000 and (F / f"{n}.pdf").stat().st_size > 1_000, n
    for t in ("functional_enrichment.csv", "coexpression_hubs_LFS_effect_hfq.csv", "gene_set_membership.csv"):
        assert (T / t).stat().st_size > 100


def test_provenance_file_does_not_overclaim():
    a = pd.read_csv(ROOT / "data" / "source_accessions.csv")
    txt = " ".join(str(v) for v in a.values.ravel()).lower()
    assert "gse146347" in txt and "not linked to an individual" in txt and "raw counts" in txt


# ----------------------------------------------------------------------------- robustness safeguards
import sys  # noqa: E402

sys.path.insert(0, str(ROOT / "scripts"))
import robustness as rb  # noqa: E402


def test_covariate_model_reproduces_baseline_tables():
    tpm = pd.read_csv(T / "tpm_renormalised_no_structural_RNA.csv", index_col=0)
    meta = pd.read_csv(P / "sample_metadata.csv").set_index("sample").loc[tpm.columns]
    d = pd.read_csv(T / "de_hfq_effect_control.csv", index_col=0)
    Y = np.log2(tpm.loc[d.index] + 1)
    t = rb.CovModel(Y, meta).test({"delta_hfq_control": 1, "WT_control": -1})
    assert np.allclose(t.log2FC, d.log2FC) and np.allclose(t.FDR, d.FDR)


def test_high_confidence_is_subset_of_baseline_with_all_checks():
    for c in CONTRASTS:
        d = pd.read_csv(T / f"de_{c}.csv", index_col=0)
        rob = d[d.robust_DEG]
        assert (rob.status != "NS").all() and not rob.group_exclusive.any() and rob.consistent_detection.all()
        assert (np.sign(rob.log2FC) == rob.call_detect_adj).all()
        assert (T / f"high_confidence_DEGs_{c}.csv").exists()


def test_robustness_summary_complete():
    s = pd.read_csv(T / "robustness_summary.csv")
    assert set(s.analysis) == set(rb.ANALYSES) and set(s.contrast) == set(CONTRASTS) and len(s) == 28
    b = s[s.analysis == "baseline"].set_index("contrast")
    r = s[s.analysis == "robust"].set_index("contrast")
    assert ((r.n_up + r.n_down) <= (b.n_up + b.n_down)).all()


def test_group_exclusive_genes_are_wt_control_only():
    g = pd.read_csv(T / "group_exclusive_genes.csv", index_col=0)
    tpm = pd.read_csv(T / "tpm_renormalised_no_structural_RNA.csv", index_col=0)
    meta = pd.read_csv(P / "sample_metadata.csv").set_index("sample")
    assert len(g) == 9
    assert (tpm.loc[g.index, meta.index[meta.group == "WT_control"]].mean(axis=1) >= 5).all()
    assert (tpm.loc[g.index, meta.index[meta.group != "WT_control"]].max(axis=1) < 1).all()


def test_lfs_response_in_hfq_survives_high_confidence_filter():
    e = pd.read_csv(T / "functional_enrichment_high_confidence.csv")
    e = e[(e.contrast == "LFS_effect_hfq") & (e.direction == "Up") & (e.FDR < 0.05)]
    assert {"SPI-1 T3SS, invasion & effectors", "Flagellar assembly & motility", "Chemotaxis"} <= set(e.gene_set)


def test_pc1_tracks_suffix_2_at_least_as_well_as_genotype():
    p = pd.read_csv(T / "pc_variance_by_factor.csv")
    p = p[p.PC == "PC1"].set_index("factor").eta2
    assert p["suffix_2"] >= p["genotype"] > 0.5            # documents the genotype / library-batch confounding
