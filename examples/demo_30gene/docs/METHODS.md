# Methods

## Data
`data/count_matrix.csv`: 30 genes x 10 samples with declared labels for WT (n=3), delta-hfq (n=3), acid stress (n=2), and heat shock (n=2). The labels, sample-level biological provenance, sequencing type/instrument, and raw-read processing are **not established**. The previously listed GSE114853 was verified to be a human enteroendocrine-cell series, and SRR7219922 was verified to be a mouse FoxP3+ regulatory T-cell run; neither corresponds to this dataset. They are explicitly marked rejected in `data/source_accessions.csv`. NC_016810.1 is a verified SL1344 reference genome, but does not establish count provenance. Therefore all outputs from the shipped matrix are software demonstrations, not interpretable experimental results. A real candidate study, GEO GSE146347 (Barrila et al. 2022), was identified and its existence verified, but its data were not retrieved or analysed here; see `docs/REAL_DATA_GUIDE.md`.

## Differential expression
Median-of-ratios size-factor normalisation, gene-wise dispersion estimation with trend fitting and MAP dispersion shrinkage, and a Wald test, using **PyDESeq2** (a Python reimplementation of DESeq2; Love et al., 2014). Per-contrast p-values are Benjamini-Hochberg adjusted across genes. The pipeline's demonstration threshold is padj < 0.05 and |log2FC| >= 1. Each stress condition is compared with wild type in its own two-group model (acid and heat have n = 2). A selected 30-gene panel does not provide a representative genome-wide background; regardless of the software output, these contrasts are not biologically validated.

PyDESeq2 prints two warnings for this 30-gene panel: the parametric dispersion trend did not converge and the software falls back to a mean-based trend, and residual degrees of freedom are below 3 for some contrasts, which makes the dispersion prior poorly estimated. Both follow from the tiny gene panel and small replicate numbers; they are a further reason to treat p-values here as demonstrations.

## Visualisation
PCA and Euclidean sample distances use log2(normalised count + 1). Heatmaps show row z-scores. Layouts follow Fig. 5-6 of Batut et al. (2021).

## Functional enrichment
Gene Ontology annotation cannot be fetched offline, so `GENE_SETS` in `scripts/rna_seq_pipeline.py` holds hand-curated functional groups. Membership was checked against the Salmonella virulence/regulation review by Fàbrega and Vila and the type III secretion system review by Hueck:

| Curated group | Genes included | Membership rationale |
|---|---|---|
| SPI-1-associated invasion genes | `hilA`, `hilC`, `hilD`, `invA`, `invF`, `sipA`, `sipB`, `sipC`, `sipD` | Pathway-associated group, not a claim that every member is physically inside SPI-1. HilC/HilD regulate the SPI-1 program; Inv/Sip genes participate in SPI-1 invasion and secretion. |
| SPI-1 transcriptional regulators | `hilA`, `hilC`, `hilD`, `invF` | Regulatory members of the SPI-1 invasion program. |
| SPI-2 (intracellular survival) genes | `ssrA`, `ssrB`, `ssaG`, `ssaH`, `ssaE`, `sseB`, `sseC`, `sseD` | SsrA/SsrB regulatory pair and selected SPI-2 secretion/intracellular-survival components. |
| T3SS effectors and translocon proteins | `sipA`, `sipB`, `sipC`, `sipD`, `sopB`, `sopE`, `sseB`, `sseC`, `sseD` | SPI-1 Sip/Sop effectors or translocation components, together with SPI-2 Sse translocon components; this intentionally spans both secretion systems. |
| Two-component systems | `phoP`, `phoQ`, `ssrA`, `ssrB`, `ompR`, `envZ` | PhoP/PhoQ, SsrA/SsrB, and EnvZ/OmpR sensor-response-regulator pairs. |
| DNA maintenance (gyrase and repair) | `gyrA`, `gyrB`, `recA` | Gyrase subunits and the recombination/repair protein RecA. |
| Global regulators | `hfq`, `fnr`, `slyA`, `rpoS` | Curated global/transcriptional regulators. |

The groups are functional teaching sets, can overlap by design, and are intersected with the 30-gene panel before testing. Enrichment is a one-sided Fisher exact test against that panel, followed by Benjamini-Hochberg adjustment. The panel is a very small, targeted background: results are descriptive and must not be presented as genome-wide GO enrichment.

References:
- Fàbrega A, Vila J. *Salmonella enterica* serovar Typhimurium skills to succeed in the host: virulence and regulation. Clin Microbiol Rev. 2013;26(2):308-341. [doi:10.1128/CMR.00066-12](https://doi.org/10.1128/CMR.00066-12).
- Hueck CJ. Type III protein secretion systems in bacterial pathogens of animals and plants. Microbiol Mol Biol Rev. 1998;62(2):379-433. [doi:10.1128/MMBR.62.2.379-433.1998](https://doi.org/10.1128/MMBR.62.2.379-433.1998).
- Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol. 2014;15:550. [doi:10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8).
- Muzellec B, Teleńczuk M, Cabeli V, Andreux M. PyDESeq2: a python package for bulk RNA-seq differential expression analysis. Bioinformatics. 2023;39(9):btad547. [doi:10.1093/bioinformatics/btad547](https://doi.org/10.1093/bioinformatics/btad547).

## Synthetic alignments
`scripts/simulate_alignments.py` places the 30 genes head-to-tail on a random reference (seed 42) and simulates, for every sample, as many single-end 50-nt reads per gene as the count matrix specifies (start positions: a mixture of uniform and 3'-skewed; both strands). Reads are written with `pysam` as coordinate-sorted, indexed BAMs whose headers state that they are synthetic. Gene counts re-derived from the BAMs equal the matrix by construction, which serves as a software round-trip check only. A separate toy locus provides spliced reads (CIGAR `N`) for the junction display. `scripts/plot_alignments.py` draws Fig13-14.

## Network and hubs
`data/string_ppi_network_edges.csv` is analysed with NetworkX: degree, weighted degree, betweenness and closeness (distance = 1/score); hubs are the five genes with the best mean rank. STRING currently rejects the project strain's SL1344 taxon ID (216597), so the refresh script resolves these gene symbols against the supported *S. Typhimurium* LT2 reference (taxon 99287) and adds up to 15 first-shell interactors. STRING associations are network evidence, not proof of direct physical interaction.

## Structure
The bundled PDB 1HK9 is the *E. coli* Hfq crystal structure (Sauter et al., 2003), used as a homolog reference. Its header is checked before it is saved or plotted. It is not a Salmonella-specific structural model.

Reference: Sauter C, Basquin J, Suck D. SM-like proteins in Eubacteria: the crystal structure of the Hfq protein from *Escherichia coli*. Nucleic Acids Res. 2003;31(14):4091-4098. [doi:10.1093/nar/gkg480](https://doi.org/10.1093/nar/gkg480).

## Reference
Batut B, van den Beek M, Doyle MA, Soranzo N (2021) RNA-Seq data analysis in Galaxy. Methods Mol Biol 2284:367-392.
