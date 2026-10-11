# data/

| Folder / file | Content | Origin |
|---|---|---|
| `raw/` | `GSE146347_fpkm_TPM_salonly.txt.gz`, Salmonella-only FPKM/TPM table | GEO GSE146347 (supplementary file) |
| `sra/` | `SraRunTable_GSE146347.csv`, `SRR_Acc_List_bacterial.txt` (12 bacterial runs SRR11235557-68) | NCBI SRA, PRJNA610228 |
| `processed/` | `tpm_matrix.csv`, `gene_annotation.csv`, `sample_metadata.csv` | **generated** by `scripts/prepare_geo_tpm.py` |
| `string/` | STRING web exports (135 nodes, score >= 0.7) and UniProt ID-mapping of hub genes | STRING v12 / UniProt |
| `alphafold/` | 14 AlphaFold DB v6 models (`AF-<UniProt>-F1-model_v6.pdb`) | AlphaFold DB |
| `source_accessions.csv` | Provenance and verification status of every resource | this repository |

Raw data belong to the original authors (Barrila et al. 2022); please cite the paper and GEO accession when reusing them.
