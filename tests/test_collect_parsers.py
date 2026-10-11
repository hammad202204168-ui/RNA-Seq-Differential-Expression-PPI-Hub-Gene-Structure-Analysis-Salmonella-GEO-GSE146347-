"""Offline tests for the parsers in scripts/collect_batch_evidence.py (network calls are NOT tested)."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import collect_batch_evidence as C

SOFT = """^SAMPLE = GSM4378554
!Sample_title = hfq1_1xg
!Sample_characteristics_ch1 = strain: chi3339
!Sample_characteristics_ch1 = bacterial_genotype: delta_hfq
!Sample_relation = SRA: https://www.ncbi.nlm.nih.gov/sra?term=SRX7847384
"""

def test_parse_soft():
    r = C.parse_soft(SOFT)
    assert r["title"] == ["hfq1_1xg"] and len(r["characteristics_ch1"]) == 2
    assert C.srr_from_soft(r) == "SRX7847384"

def test_read_name_casava():
    p = C.parse_read_name("@NB501234:45:H7ABCBGX5:1:11101:1234:5678 1:N:0:ATCACG")
    assert p == dict(instrument="NB501234", run="45", flowcell="H7ABCBGX5", lane="1")

def test_read_name_sra_reformatted_has_no_flowcell():
    assert C.parse_read_name("@SRR11235557.1 1 length=73") is None

def test_fastq_qc():
    reads = []
    for i in range(10):
        reads += [f"@NB5:1:FC1:2:1101:{i}:1 1:N:0:X", "ACGT" * 5 + "AGATCGGAAGAGC" if i < 3 else "GGGGCCCCAAAATTTT" + str(i % 2) * 0 + "ACGTAC", "+", "I" * 10]
    q = C.fastq_qc(reads, 100)
    assert q["n_reads_sampled"] == 10 and abs(q["adapter_frac"] - 0.3) < 1e-9
    assert q["flowcells"].startswith("NB5:1:FC1:L2") and q["header_parsed_frac"] == 1.0

def test_link_columns():
    titles = {"GSM1": "WT2_1xg", "GSM2": "hfq3_1xg_2", "GSM3": "hfq1_1xg"}
    m, miss = C.link_columns(["WT2_1xg", "hfq3_1xg_2", "hfq1_1xg", "WT9_LS"], titles)
    assert m == {"WT2_1xg": "GSM1", "hfq3_1xg_2": "GSM2", "hfq1_1xg": "GSM3"} and miss == ["WT9_LS"]
