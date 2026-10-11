#!/usr/bin/env python3
"""
Parse GEO GSE146347 `fpkm_TPM_salonly` (Salmonella-only) into tidy files.

Input : data/raw/GSE146347_fpkm_TPM_salonly.txt.gz   (supplementary file from the GEO series page)
Output: data/processed/tpm_matrix.csv          genes x 12 samples, TPM
        data/processed/gene_annotation.csv     id, label, symbol, description, type, replicon, start, end
        data/processed/sample_metadata.csv     sample, genotype, culture, group

Notes
* The file holds FPKM (12 columns) and TPM (12 columns). Only TPM is used; neither is raw counts.
* Column names carry strain and culture ("WT2_1xg", "hfq1_LS"): 1xg = normal-gravity control orientation,
  LS = low-fluid-shear (LFS) orientation, matching `bacterial_growth_orientation` in the SRA run table.
* The file has no GSM/SRR identifiers, so each column is NOT linked to an individual GSM/SRR run. Only the
  group structure (WT/hfq x control/LFS, n = 3 each) is verified, against SraRunTable.csv.
"""
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "GSE146347_fpkm_TPM_salonly.txt.gz"
OUT = ROOT / "data" / "processed"


def sample_id(col):
    s = re.sub(r"[._](re)?merged?_STAR.*$", "", col)           # drop STAR/merge suffix
    s = re.sub(r"_R1(?=$|_)", "", s)
    return s


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(RAW, sep="\t")
    tpm_cols = [c for c in df.columns if c.endswith("_TPM")]
    assert len(tpm_cols) == 12, f"expected 12 TPM columns, found {len(tpm_cols)}"
    tpm = df.set_index("Gene")[tpm_cols].copy()
    tpm.columns = [sample_id(c[:-4]) for c in tpm_cols]
    assert tpm.columns.is_unique

    meta = pd.DataFrame({"sample": tpm.columns})
    meta["genotype"] = ["delta_hfq" if s.startswith("hfq") else "WT" for s in meta["sample"]]
    meta["culture"] = ["LFS" if "_LS" in s else "control" for s in meta["sample"]]
    meta["group"] = meta.genotype + "_" + meta.culture
    assert meta.group.value_counts().eq(3).all() and meta.group.nunique() == 4, meta.group.value_counts()
    assert (abs(tpm.sum() - 1e6) < 5).all(), "TPM columns should sum to ~1e6"

    ann = df[["Gene", "gene", "locus", "symbol", "description", "type"]].rename(columns={"Gene": "GeneID", "gene": "name"})
    loc = ann.locus.str.extract(r"^(?P<replicon>[^:]+):(?P<start>\d+)-(?P<end>\d+)$")
    ann = pd.concat([ann, loc], axis=1)
    dup = ann["name"].duplicated(keep=False)
    ann["label"] = ann["name"].where(~dup, ann["name"] + "|" + ann["GeneID"])      # readable and unique
    ann["type"] = ann["type"].where(ann["type"].isin(["protein-coding", "tRNA", "rRNA", "ncRNA", "pseudo", "other"]), "other/unspecified")

    tpm.index.name = "GeneID"
    tpm.to_csv(OUT / "tpm_matrix.csv")
    ann.set_index("GeneID").to_csv(OUT / "gene_annotation.csv")
    meta.to_csv(OUT / "sample_metadata.csv", index=False)
    print(f"{tpm.shape[0]} genes x {tpm.shape[1]} samples; groups: {meta.group.value_counts().to_dict()}")
    print("hfq TPM (WT vs delta_hfq mean):",
          tpm.loc[ann.loc[ann.name == 'hfq', 'GeneID']].T.groupby(meta.set_index('sample').genotype).mean().round(1).to_dict())


if __name__ == "__main__":
    main()
