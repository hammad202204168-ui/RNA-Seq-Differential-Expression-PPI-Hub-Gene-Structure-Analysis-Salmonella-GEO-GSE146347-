# Synthetic alignments

**These files are simulated. They are not sequencing data, and the reference is random sequence, not Salmonella.**

| File | What it is |
|---|---|
| `synthetic_reference.fa` | Random 52 kb sequence with the 30 panel genes laid out head-to-tail |
| `synthetic_genes.gtf`, `synthetic_gene_coordinates.csv` | Toy gene models (gene + exon records) |
| `<sample>.bam` / `.bai` | One coordinate-sorted BAM per sample: 50-nt single-end reads, number of reads per gene equal to `data/count_matrix.csv` |
| `spliced_locus_demo.*` | Toy 3-exon locus with spliced reads, used only to demonstrate a Sashimi-style display |

Regenerate with `python scripts/simulate_alignments.py` (seed 42, deterministic). The script re-counts reads from the BAMs and aborts if the
result differs from the count matrix. Each BAM header carries an `@CO` line stating that the reads are synthetic.

Because the reads are generated *from* the counts, the BAMs add no independent evidence; they exist so that the alignment-inspection
and counting steps of the reference workflow can be exercised and illustrated. Real data would replace them with STAR/HISAT2 output.
