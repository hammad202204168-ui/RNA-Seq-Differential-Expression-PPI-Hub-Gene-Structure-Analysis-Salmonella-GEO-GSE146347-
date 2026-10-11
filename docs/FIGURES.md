# Figure legends

Fig0-Fig13 are written to `results/figures/` (600 dpi PNG + vector PDF) by `python scripts/real_data_analysis.py`; Fig14-Fig16 by `python scripts/ppi_hubs_structures.py` (`make ppi`). Data figures carry a footer stating that the analysis is exploratory and TPM-based.
Groups: WT/Δ*hfq* × control/LFS (n = 3 each). "DEG": FDR < 0.05 and |log2FC| ≥ 1.

| Figure | Legend |
|---|---|
| **Fig0** | Graphical summary. **a**, volcano of Δ*hfq* vs WT (control). **b**, DEG counts for the four contrasts (pale = all DEGs, dark = high-confidence). **c**, gene sets enriched at FDR < 0.05. |
| **Fig1** | Workflow. Dashed boxes = steps not performed (need raw reads). |
| **Fig2** | QC. **a**, log2(TPM + 1) per sample. **b**, genes with TPM ≥ 1. **c**, mean expression and tested genes. **d**, share of the original TPM taken by rRNA and tRNA/ncRNA (differs strongly between groups). **e**, sequenced bases per group from the SRA run table (group level). |
| **Fig3** | **a**, PCA. **b**, clustered Euclidean sample-distance matrix. Triangles = LFS, circles = control. |
| **Fig4** | log2(TPM + 1) of *hfq* (positive control) and the top three down- and two up-regulated **high-confidence** DEGs (mean TPM ≥ 10 in at least one group) for Δ*hfq* vs WT (control). |
| **Fig5** | **a**, row z-scores of the 40 most significant DEGs (* high-confidence, † expressed in one group only). **b**, volcano (top 12 labelled; black rings = high-confidence DEGs). |
| **Fig6** | Annotated "Volcano Plot Interpretation" chart for Δ*hfq* vs WT (control): negative/positive change relative to control, the zero point, the FDR thresholds and the meaning of each region. |
| **Fig7** | Volcano plots of all four contrasts, with up/down counts. |
| **Fig8** | **a**, p-value histograms. **b**, MA plots (mean log2(TPM + 1) vs log2FC) for each contrast. |
| **Fig9** | LFS effect in WT vs in Δ*hfq*; purple = genotype × culture FDR < 0.05. Nine of the 20 hits are genes expressed only in WT-control (see README, Data problems). |
| **Fig10** | Rule-based gene-set enrichment for all contrasts and directions. Dot size = −log10 FDR; the label under a dot is DEGs in set / set size; colour marks FDR < 0.05; a thick outline means the set is still FDR < 0.05 using high-confidence DEGs only. |
| **Fig11** | **a**, Spearman co-expression network (|ρ| ≥ 0.9) of the 60 most significant high-confidence LFS-in-Δ*hfq* DEGs; node colour = log2FC; ring = top-5 hub. **b**, top 8 genes by consensus rank. |
| **Fig13** | Robustness checks. **a**, all vs high-confidence DEGs per contrast. **b**, % of baseline DEGs retained under each alternative analysis (numbers in cells). **c**, variance of PC1/PC2 explained (η²) by genotype, culture and a sample-name ending "_2": PC1 tracks the suffix at least as well as genotype. |
| **Fig12** | Schematic volcano reading guide with **simulated points (not study data)**. |
| **Fig14** | STRING protein-protein association network (combined score >= 0.7, *S.* Typhimurium LT2) of the 141 high-confidence LFS-in-Δ*hfq* DEGs. Edges coloured by whether they survive removal of text-mining evidence. Also saved as SVG. |
| **Fig15** | AlphaFold DB v6 C-alpha traces of 14 hub/target proteins (predicted monomers, schematic). |
| **Fig16** | Per-residue pLDDT confidence bands for the same 14 models (very high > 90, confident 70-90, low 50-70, very low < 50). Also saved as SVG. |
## How to read a volcano plot
* **x-axis:** log2 fold change; 0 = no change, negative = lower in the first-named group.
* **y-axis:** −log10 of the adjusted p-value; higher = stronger evidence.
* Coloured points pass both the FDR and fold-change cut-offs. Statistical significance is not biological importance, and with n = 3 per group effect sizes are uncertain.
