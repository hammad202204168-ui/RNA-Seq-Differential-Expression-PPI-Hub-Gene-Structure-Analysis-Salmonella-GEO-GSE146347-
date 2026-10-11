# Figure legends

All figures are written to `results/figures/` as 600 dpi PNG and vector PDF by `python scripts/rna_seq_pipeline.py`.
Every data figure carries the footer *"Demonstration data – provenance unverified; not a biological result"* because the supplied
count matrix has no verified source (see `data/source_accessions.csv`). Set `DEMO_NOTE = None` in the script only after the matrix is replaced by verified data.

| Figure | Legend |
|---|---|
| **Fig0** | Graphical summary. **a**, volcano plot of the mutant vs control contrast. **b**, log2 fold change (± 1.96 SE) of five invasion-associated genes in the three contrasts (stress groups n = 2, exploratory). **c**, fraction of each curated gene set that is differentially expressed (30-gene background; not genome-wide GO). |
| **Fig1** | Workflow. Dashed boxes are example upstream steps that were not established or run. |
| **Fig2** | Quality control. **a**, total counts per sample. **b**, log2(count + 1) distributions. **c**, variance against mean per gene with the Poisson line (variance = mean). |
| **Fig3** | **a**, PCA of log2(normalised count + 1). **b**, clustered Euclidean sample-distance matrix (layout after Fig. 5 of Batut et al., 2021). |
| **Fig4** | PyDESeq2 diagnostics. **a**, gene-wise, fitted and final dispersions (values below 1e-4 fall off the axis). **b**, p-value histogram. **c**, MA plot. |
| **Fig5** | **a**, row z-scores of significant genes (padj < 0.05, \|log2FC\| ≥ 1) across all samples; colour strip = condition. **b**, volcano plot (layout after Fig. 6 of Batut et al., 2021). |
| **Fig6** | Normalised counts of the six significant genes; bars mark group means. |
| **Fig7** | log2 fold change of all 30 genes in the three contrasts (clipped at ±6). \* padj < 0.05, \*\* padj < 0.001. |
| **Fig8** | Over-representation of curated gene sets (one-sided Fisher exact test, Benjamini–Hochberg adjusted) against the 30-gene panel. |
| **Fig9** | **a**, STRING network (LT2, score ≥ 0.40, plus up to 15 first-shell interactors). Node colour = log2FC where the gene is in the count panel; ring = top-5 hub by mean rank of degree, weighted degree, betweenness and closeness. **b**, top 8 genes by consensus rank. |
| **Fig10** | Cα backbone of PDB 1HK9 (*E. coli* Hfq homolog; not a Salmonella-specific model), top and side views. Drawn only if the PDB header identifies Hfq. |
| **Fig11** | Annotated volcano plot of the actual mutant vs control contrast. The y-axis is symmetric-log so the significance bands stay readable. Dashed lines: adj. p = 0.05 and 0.001, and \|log2FC\| = 1. |
| **Fig12** | Schematic reading guide with **simulated points (not study data)**, showing how to interpret regions of a volcano plot. |
| **Fig13** | IGV-style view of the SPI-1 gene cluster on the **synthetic** reference (after Fig. 4a of Batut et al., 2021). Coverage (top), packed reads (blue = forward, red = reverse; first 45 rows) and gene models for wild type (**a**) and Δhfq (**b**), same y-scales. Reads were simulated from the count matrix. |
| **Fig14** | Sashimi-style plot of a **synthetic** three-exon toy locus (after Fig. 4b). Arcs show split reads across introns, with read counts, recovered from the BAM CIGAR strings. Illustrates the display only; Salmonella mRNAs are not spliced. |

## How to read a volcano plot
* **x-axis**: log2 fold change. 0 is "no change"; negative = lower in mutant than control, positive = higher.
* **y-axis**: −log10 of the (adjusted) p-value. Higher = stronger statistical evidence.
* Genes are coloured only if they pass both the significance threshold and the fold-change cut-off used by the pipeline (padj < 0.05 and \|log2FC\| ≥ 1).
* Statistical significance is not biological importance, and with n = 2–3 replicates per group fold-change estimates are uncertain.
