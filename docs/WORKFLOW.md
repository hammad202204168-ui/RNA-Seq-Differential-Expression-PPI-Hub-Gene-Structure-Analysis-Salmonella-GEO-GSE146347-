# Workflow: step by step

Steps marked **manual** were done by hand with web tools and the exported files are stored in `data/`; steps marked **script** are Python scripts in `scripts/`.

| # | Step | Tool | How | Input -> output |
|---|---|---|---|---|
| 1 | Find the dataset and design | NCBI GEO, SRA Run Selector | manual | GSE146347, PRJNA610228 -> `data/sra/SraRunTable_GSE146347.csv`, `SRR_Acc_List_bacterial.txt` |
| 2 | Download the Salmonella TPM table | NCBI GEO (supplementary file) | manual | -> `data/raw/GSE146347_fpkm_TPM_salonly.txt.gz` |
| 3 | Parse table, sample metadata, annotation | Python (pandas) | script `prepare_geo_tpm.py` | -> `data/processed/` |
| 4 | QC, PCA, rRNA / depth / batch checks | Python (SciPy, scikit-learn, matplotlib) | script `real_data_analysis.py`, `robustness.py` | -> `results/figures/Fig2-3`, `robustness_summary.csv` |
| 5 | Differential expression (t-test on log2(TPM+1), BH-FDR) and high-confidence DEG set | Python (SciPy) | script `real_data_analysis.py` | -> `results/tables/de_*.csv`, `high_confidence_DEGs_*.csv` |
| 6 | Gene-set enrichment (rule-based) and co-expression network | Python (SciPy, NetworkX) | script `real_data_analysis.py` | -> `functional_enrichment*.csv`, Fig10-11 |
| 7 | PPI network for the 141 high-confidence DEGs | STRING web interface (*S.* Typhimurium LT2, score >= 0.7) | manual export | -> `data/string/*.tsv` |
| 8 | Hub-gene ranking (degree, weighted degree, betweenness; with/without text-mining) | Python (NetworkX) | script `ppi_hubs_structures.py` | -> `results/tables/ppi_*.csv`, Fig14 |
| 9 | Map hub genes to UniProt | UniProt ID mapping web tool | manual export | -> `data/string/uniprot_idmapping_hub_genes.tsv` |
| 10 | Download predicted structures for 14 targets | AlphaFold Protein Structure Database (v6) | manual download | -> `data/alphafold/*.pdb` |
| 11 | pLDDT summary and structure figures | Python (matplotlib) | script `ppi_hubs_structures.py` | -> `hub_structures.csv`, Fig15-16 |
| 12 | Cartoon renders (optional) | PyMOL / ChimeraX | scripts provided, **not tested** | `render_alphafold_cartoons.pml/.cxc` |

Not done: read trimming, alignment and raw-count generation (no raw reads were used; see `BATCH_AND_COUNTS_PLAN.md`), DESeq2/edgeR, GO/KEGG enrichment, Cytoscape/cytoHubba.

Re-running the scripted steps (3-6, 8, 11) is wrapped in `make run` and `make ppi`; the manual steps (1, 2, 7, 9, 10) need no re-running because their outputs are included in the repository.
