import pytest
import pandas as pd
from scripts.rna_seq_pipeline import adjust_pvalues_bh, load, validate_inputs
from scripts.structure_utils import is_hfq_pdb
from pathlib import Path
T = Path(__file__).resolve().parents[1] / "results" / "tables"


def test_degs_meet_thresholds():
    d = pd.read_csv(T / "significant_degs.csv")
    assert len(d) > 0
    assert (d.padj < 0.05).all() and (d.log2FoldChange.abs() >= 1).all()


def test_hfq_group_difference_in_demonstration_matrix():
    d = pd.read_csv(T / "de_Mutant_vs_Control.csv", index_col=0)
    assert d.loc["hfq", "log2FoldChange"] < -1 and d.loc["hfq", "padj"] < 0.05


def test_normalised_counts_shape():
    assert pd.read_csv(T / "normalized_counts.csv", index_col=0).shape == (30, 10)


def test_data_diagnostics_are_written():
    d = pd.read_csv(T / "data_diagnostics.csv", index_col=0)["value"]
    assert d["n_genes"] == 30
    assert d["n_samples"] == 10
    assert d["dynamic_range_fold_(max/min_WT_mean)"] > 1


def test_metadata_does_not_claim_unverified_sequencing_details():
    metadata = pd.read_csv(T.parents[1] / "data" / "sample_metadata.csv")
    assert metadata["Sequencing_Type"].eq("Unknown").all()
    assert metadata["Instrument"].eq("Unknown").all()


def test_enrichment_reports_curated_panel_sets():
    d = pd.read_csv(T / "functional_enrichment.csv").set_index("gene_set")["set_size"]
    expected = {
        "SPI-1-associated invasion genes": 9,
        "SPI-1 transcriptional regulators": 4,
        "SPI-2 (intracellular survival) genes": 8,
        "T3SS effectors and translocon proteins": 9,
        "Two-component systems": 6,
        "DNA maintenance (gyrase and repair)": 3,
        "Global regulators": 4,
    }
    assert d.to_dict() == expected


def test_hfq_structure_header_validation():
    assert is_hfq_pdb("HEADER    RNA BINDING PROTEIN\nTITLE     CRYSTAL STRUCTURE OF THE HFQ PROTEIN")
    assert not is_hfq_pdb("HEADER    OXIDOREDUCTASE\nTITLE     ALCOHOL DEHYDROGENASE")


def test_bh_adjustment_handles_ties_and_nan():
    p = pd.Series([0.01, 0.01, 0.04, float("nan")])
    adjusted = adjust_pvalues_bh(p)
    assert adjusted.iloc[:3].tolist() == pytest.approx([0.015, 0.015, 0.04])
    assert pd.isna(adjusted.iloc[3])


def test_input_validation_accepts_project_inputs():
    counts, metadata, edges = load()
    checked_counts, checked_metadata, checked_edges = validate_inputs(counts, metadata, edges)
    assert checked_counts.shape == (30, 10)
    assert checked_metadata.index.tolist() == checked_counts.columns.tolist()
    assert checked_edges.score.between(0, 1, inclusive="right").all()


def test_input_validation_rejects_mismatched_sample_metadata():
    counts = pd.DataFrame([[1, 2]], index=["gene"], columns=["sample1", "sample2"])
    metadata = pd.DataFrame({"Condition": ["Control", "Mutant"]}, index=["sample1", "other"])
    edges = pd.DataFrame({"source": ["a"], "target": ["b"], "score": [0.5]})
    with pytest.raises(ValueError, match="must uniquely match"):
        validate_inputs(counts, metadata, edges)


@pytest.mark.parametrize(
    ("count_values", "score", "message"),
    [
        ([[1.5, 2]], 0.5, "integer raw counts"),
        ([[-1, 2]], 0.5, "non-negative integer raw counts"),
        ([[1, 2]], 0, "edge scores"),
        ([[1, 2]], 1.1, "edge scores"),
    ],
)
def test_input_validation_rejects_invalid_counts_and_edges(count_values, score, message):
    counts = pd.DataFrame(count_values, index=["gene"], columns=["sample1", "sample2"])
    metadata = pd.DataFrame({"Condition": ["Control", "Mutant"]}, index=counts.columns)
    edges = pd.DataFrame({"source": ["a"], "target": ["b"], "score": [score]})
    with pytest.raises(ValueError, match=message):
        validate_inputs(counts, metadata, edges)


def test_all_figures_written_as_png_and_pdf():
    figs = T.parents[0] / "figures"
    names = ["Fig0_graphical_summary", "Fig1_workflow_overview", "Fig2_quality_control", "Fig3_PCA_and_sample_distance",
             "Fig4_DESeq2_diagnostics", "Fig5_heatmap_and_volcano", "Fig6_top_DEG_normalised_counts",
             "Fig7_log2FC_across_conditions", "Fig8_functional_enrichment", "Fig9_PPI_network_and_hubs",
             "Fig10_Hfq_structure", "Fig11_volcano_interpretation", "Fig12_volcano_reading_guide"]
    for n in names:
        assert (figs / f"{n}.png").stat().st_size > 10_000, n
        assert (figs / f"{n}.pdf").stat().st_size > 1_000, n


def test_every_contrast_table_has_all_genes():
    for g in ("Mutant", "AcidStress", "HeatShock"):
        d = pd.read_csv(T / f"de_{g}_vs_Control.csv", index_col=0)
        assert len(d) == 30 and {"log2FoldChange", "padj", "status"} <= set(d.columns)


def test_synthetic_bams_roundtrip_matrix_and_are_labelled():
    import pysam
    root = T.parents[1]
    counts = pd.read_csv(root / "data" / "count_matrix.csv", index_col="GeneID")
    derived = pd.read_csv(T / "bam_derived_counts.csv", index_col="GeneID")
    assert derived.equals(counts)
    for sample in counts.columns:
        with pysam.AlignmentFile(str(root / "data" / "synthetic_alignments" / f"{sample}.bam")) as bam:
            assert any("SYNTHETIC" in c for c in bam.header["CO"])
            assert bam.header["HD"]["SO"] == "coordinate" and bam.mapped == counts[sample].sum()


def test_spliced_demo_junctions():
    import pysam
    root = T.parents[1]
    junc = {}
    with pysam.AlignmentFile(str(root / "data" / "synthetic_alignments" / "spliced_locus_demo.bam")) as bam:
        for a in bam.fetch("spliced_locus_demo"):
            b = a.get_blocks()
            for (_, e), (s, _) in zip(b[:-1], b[1:]):
                junc[(e, s)] = junc.get((e, s), 0) + 1
    assert junc == {(900, 1700): 230, (2000, 3000): 200, (900, 3000): 50}


def test_alignment_figures_written():
    figs = T.parent / "figures"
    for n in ("Fig13_alignment_view_IGV_style", "Fig14_sashimi_style_junctions"):
        assert (figs / f"{n}.png").stat().st_size > 10_000 and (figs / f"{n}.pdf").stat().st_size > 1_000


def test_accessions_do_not_overclaim_provenance():
    a = pd.read_csv(T.parents[1] / "data" / "source_accessions.csv").set_index("Accession")
    assert a.loc["GSE114853", "Verification"] == "Rejected"
    assert a.loc["Not established", "Verification"] == "Unverified"
    assert "Not the origin" in a.loc["GSE146347", "Finding"]
    assert "counts not retrieved" in a.loc["GSE146347", "Verification"]


def test_geo_soft_parser():
    from scripts.geo_helper import parse_soft_samples, series_url
    df = parse_soft_samples("^SAMPLE = GSM1\n!Sample_title = a\n!Sample_characteristics_ch1 = strain: wild type\n^SAMPLE = GSM2\n!Sample_title = b\n")
    assert list(df.accession) == ["GSM1", "GSM2"] and df.loc[0, "char_strain"] == "wild type"
    assert series_url("GSE146347").endswith("/GSE146nnn/GSE146347/")
