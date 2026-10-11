import os, sys, textwrap
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import pandas as pd, featurecounts_to_matrix as F

def test_roundtrip(tmp_path):
    fc = tmp_path / "fc.txt"
    fc.write_text("# Program:featureCounts\nGeneid\tChr\tStart\tEnd\tStrand\tLength\tdata/rebuild/bam/SRR2_Aligned.sortedByCoord.out.bam\tdata/rebuild/bam/SRR1_Aligned.sortedByCoord.out.bam\n"
                  "STM0001\tNC\t1\t9\t+\t9\t5\t7\nSTM0002\tNC\t1\t9\t+\t9\t0\t3\n")
    mp = tmp_path / "map.csv"
    pd.DataFrame(dict(tpm_column=["hfq1_1xg", "WT2_1xg"], gsm=["G1", "G2"], srr=["SRR1", "SRR2"],
                      genotype=["delta_hfq", "wild type"], orientation=["control", "control"], bases=[1, 2])).to_csv(mp, index=False)
    c, m = F.run(str(fc), str(mp), None, str(tmp_path / "c.csv"), str(tmp_path / "m.csv"))
    assert list(c.columns) == ["hfq1_1xg", "WT2_1xg"] and c.loc["STM0001", "hfq1_1xg"] == 7   # SRR1 column reordered correctly
    assert m.set_index("sample").loc["WT2_1xg", "genotype"] == "WT"
    assert m.set_index("sample").loc["hfq1_1xg", "batch"] == "b1" and "PROXY" in m.batch_source[0]
