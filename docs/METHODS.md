# Methods

## Data
Input: `data/raw/GSE146347_fpkm_TPM_salonly.txt.gz`, the Salmonella-only supplementary table of GEO GSE146347 (Barrila *et al.* 2022). It contains
12 FPKM and 12 TPM columns for 4,743 genes annotated on NC_003197.1 and NC_003277.1. `scripts/prepare_geo_tpm.py` keeps the TPM columns, assigns
genotype (WT / delta_hfq) and culture (control = "1xg", LFS = "LS") from the column names, and checks 3 samples per group and column sums of 10^6.
Sample columns are not linked to individual GSM/SRR runs (the file has no run identifiers); group sizes agree with the SRA run table.

## Normalisation and filtering
TPM values are used as supplied, except that the 114 rRNA, tRNA and ncRNA genes are removed and each sample is rescaled to 10^6 over the remaining 4,629 genes,
because rRNA took 0–24 % of TPM depending on the library. Genes with TPM ≥ 1 in at least 3 of the 12 samples are tested (4,420 genes). Values are transformed as log2(TPM + 1).

## Statistics
A four-cell (genotype × culture) model on log2(TPM + 1): group means and a pooled within-group variance (8 residual df). Contrasts are tested with
Student's *t* statistics: Δ*hfq* vs WT in control and in LFS; LFS vs control in WT and in Δ*hfq*; and the genotype × culture interaction
([Δ*hfq*,LFS − Δ*hfq*,control] − [WT,LFS − WT,control]). log2FC is the difference of group means on the log2(TPM + 1) scale. p-values are adjusted per contrast
with Benjamini–Hochberg. Differentially expressed genes (DEGs) have FDR < 0.05 and |log2FC| ≥ 1. No variance moderation or count model is applied because raw counts are unavailable.

## Robustness checks (`scripts/robustness.py`)
Four data problems cannot be removed by analysis (unequal rRNA depletion, shallow and unequal depth, a possible library batch confounded with genotype, genes expressed in a single group); see README. Each contrast is therefore re-tested:
* **no_exclusive** – refit after removing genes with mean TPM ≥ 5 in one group and TPM < 1 in all samples of the other three;
* **detect_adj / rrna_adj / batch2_adj** – the four-group model plus one z-scored sample covariate: number of genes with TPM ≥ 1, rRNA share of the original TPM, or an indicator for a sample name ending in "_2" (7 residual df; condition numbers 4.5–7.8 are in `model_conditioning.csv`);
* **consistent** – the higher group must have TPM ≥ 1 in all three replicates and mean TPM ≥ 5.
A **high-confidence DEG** is a baseline DEG that is not group-exclusive, passes *consistent*, and keeps its direction with FDR < 0.05 and |log2FC| ≥ 1 in *detect_adj* and in *no_exclusive*. The rRNA and `_2` models are reported as sensitivity analyses only, because both covariates are nearly collinear with genotype and absorb genuine genotype effects (adjustment is then over-conservative: the deleted gene *hfq* fails). PC variance explained by genotype, culture and the `_2` suffix is in `pc_variance_by_factor.csv`. The covariate model with no covariate reproduces the baseline tables exactly (tested).

## Gene-set enrichment
Sets are defined from the supplied annotation (gene-name prefixes and description keywords; see `build_sets` in `scripts/real_data_analysis.py` and
`results/tables/gene_set_membership.csv`). Up- and down-regulated DEGs of each contrast (and, separately, the high-confidence DEGs) are tested separately against the 4,420 tested genes with a one-sided Fisher exact test,
sets of ≥ 5 genes, BH-adjusted within each list. This is not GO/KEGG enrichment.

## Co-expression network
The 60 most significant high-confidence DEGs of the LFS-in-Δ*hfq* contrast are connected when |Spearman ρ| ≥ 0.9 across the 12 samples; isolated genes are dropped. Hubs are the 5 genes with the best mean rank of
degree, weighted degree, betweenness and closeness (distance = 1 − |ρ|). Because correlations across 12 samples reflect the 2 × 2 group structure, this is descriptive and not an interaction network.

## PPI network and structures
Input: the 141 high-confidence LFS-in-Δ*hfq* DEGs. The network is a STRING web export (*S.* Typhimurium LT2, combined score >= 0.7; `data/string/`). Hubs are ranked by degree, weighted degree and betweenness; each ranking is repeated after removing text-mining evidence to flag `stable_hub` genes. UniProt accessions come from a UniProt ID-mapping export; predicted structures are AlphaFold DB v6 models (`data/alphafold/`) summarised by mean pLDDT. STRING associations are not proof of binding and AlphaFold monomers do not model assembled complexes. Run with `make ppi`.

## Figures
PCA and sample distances use log2(TPM + 1) of tested genes (PCA centred, unscaled). Heatmaps show row z-scores of the 40 most significant DEGs. Layouts follow Fig. 5–6 of Batut *et al.* (2021);
the annotated volcano chart follows the "Volcano Plot Interpretation" style.

## Reference
Batut B, van den Beek M, Doyle MA, Soranzo N (2021) RNA-Seq data analysis in Galaxy. Methods Mol Biol 2284:367–392.
