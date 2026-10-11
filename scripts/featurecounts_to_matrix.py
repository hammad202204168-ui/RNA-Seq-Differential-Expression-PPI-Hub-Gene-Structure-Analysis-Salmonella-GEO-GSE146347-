#!/usr/bin/env python3
"""featureCounts table -> genes x samples raw-count CSV with TPM-table column names + metadata for batch_aware_de.py.

  python scripts/featurecounts_to_matrix.py data/rebuild/featureCounts_raw.txt data/sra/column_to_run_map.csv \
         [--batch_map data/sra/batch_assignment.csv] --out_counts data/rebuild/counts_matrix.csv --out_meta data/rebuild/sample_metadata_counts.csv

batch_assignment.csv (you write it AFTER reviewing batch_evidence.csv / hearing from the authors): columns tpm_column,batch.
Without it, batch defaults to the `_2`-suffix proxy and the output says so, because that is a guess, not a documented batch.
"""
import argparse, re, sys
import pandas as pd


def run(fc_path, map_path, batch_path, out_counts, out_meta):
    fc = pd.read_csv(fc_path, sep="\t", comment="#", index_col=0)
    coord_cols = [c for c in ("Chr", "Start", "End", "Strand", "Length") if c in fc.columns]
    fc[coord_cols].to_csv(out_counts.replace(".csv", "_gene_coords.csv"))   # join key to data/processed/gene_annotation.csv (replicon/start/end)
    fc = fc.drop(columns=coord_cols)
    fc.columns = [re.search(r"(SRR\d+)", c).group(1) for c in fc.columns]
    mp = pd.read_csv(map_path).dropna(subset=["tpm_column"])
    if set(mp.srr) - set(fc.columns):
        sys.exit(f"runs missing from featureCounts output: {sorted(set(mp.srr) - set(fc.columns))}")
    fc = fc[mp.srr.tolist()]
    fc.columns = mp.tpm_column.tolist()
    fc.to_csv(out_counts)
    meta = pd.DataFrame({"sample": mp.tpm_column.values,
                         "genotype": mp.genotype.replace({"wild type": "WT"}).values,
                         "culture": mp.orientation.values})
    if batch_path:
        b = pd.read_csv(batch_path).set_index("tpm_column")["batch"].astype(str)
        meta["batch"] = meta["sample"].map(b)
        if meta["batch"].isna().any():
            sys.exit("batch_assignment.csv does not cover all samples")
        meta["batch_source"] = "documented (" + batch_path + ")"
    else:
        meta["batch"] = meta["sample"].str.endswith("_2").map({True: "b2", False: "b1"})
        meta["batch_source"] = "PROXY: _2 suffix; meaning undocumented"
    meta.to_csv(out_meta, index=False)
    return fc, meta


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("featurecounts"); ap.add_argument("column_map")
    ap.add_argument("--batch_map")
    ap.add_argument("--out_counts", default="data/rebuild/counts_matrix.csv")
    ap.add_argument("--out_meta", default="data/rebuild/sample_metadata_counts.csv")
    a = ap.parse_args()
    run(a.featurecounts, a.column_map, a.batch_map, a.out_counts, a.out_meta)
