#!/usr/bin/env python3
"""
Generate SYNTHETIC alignments (BAM) from the demonstration count matrix.

Why: the project has no verified raw reads, so real FASTQ/BAM files cannot be provided. To exercise the
alignment-inspection and counting steps of the reference workflow (Batut et al., 2021, Fig. 4), this script
simulates single-end 50-nt reads from `data/count_matrix.csv` onto a SYNTHETIC reference (random sequence, toy
gene models). Nothing here is Salmonella sequence or real sequencing output.

Outputs (data/synthetic_alignments/):
    synthetic_reference.fa / .fai        random 'genome' with the 30 panel genes laid out in matrix order
    synthetic_genes.gtf                  gene + exon records for those genes
    <sample>.bam / .bai                  one coordinate-sorted BAM per sample (read counts == matrix counts)
    spliced_locus_demo.{fa,gtf,bam,bai}  3-exon toy locus with spliced reads (display demo only)
results/tables/bam_derived_counts.csv    gene counts re-derived from the BAMs (round-trip check)

Deterministic (seed 42); re-running reproduces identical read positions.
"""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "synthetic_alignments"
TAB = ROOT / "results" / "tables"
SEED = 42
READ_LEN = 50
CONTIG = "synthetic_ref"
SPLICE_CONTIG = "spliced_locus_demo"
COMMENT = ("SYNTHETIC reads simulated from data/count_matrix.csv onto a random reference; "
           "not sequencing data and not Salmonella sequence")


def random_seq(rng, n):
    return "".join(rng.choice(list("ACGT"), size=n))


def build_reference(rng, genes):
    """Lay genes out head-to-tail with 400-700 bp intergenic gaps. Returns sequence and gene table."""
    rows, pos = [], 600
    for g in genes:
        length = int(rng.integers(700, 1600))
        strand = "+" if rng.random() < 0.5 else "-"
        rows.append(dict(gene=g, start=pos, end=pos + length, strand=strand))   # 0-based half-open
        pos += length + int(rng.integers(400, 700))
    seq = random_seq(rng, pos + 600)
    return seq, pd.DataFrame(rows).set_index("gene")


def write_fasta(path, name, seq):
    path.write_text(f">{name} synthetic random sequence (not a real genome)\n" +
                    "\n".join(seq[i:i + 70] for i in range(0, len(seq), 70)) + "\n")
    pysam.faidx(str(path))


def write_gtf(path, contig, table, exons=None):
    with open(path, "w") as fh:
        for g, r in table.iterrows():
            attr = f'gene_id "{g}"; gene_name "{g}";'
            fh.write(f"{contig}\tsynthetic\tgene\t{r.start + 1}\t{r.end}\t.\t{r.strand}\t.\t{attr}\n")
            for (s, e) in (exons or {}).get(g, [(r.start, r.end)]):
                fh.write(f"{contig}\tsynthetic\texon\t{s + 1}\t{e}\t.\t{r.strand}\t.\t{attr}\n")


def make_header(contig, length, sample):
    return {"HD": {"VN": "1.6", "SO": "unsorted"}, "SQ": [{"SN": contig, "LN": length}],
            "RG": [{"ID": sample, "SM": sample, "PL": "SIMULATED"}], "CO": [COMMENT]}


def new_read(name, seq, start, flag, cigar, sample):
    a = pysam.AlignedSegment()
    a.query_name, a.query_sequence, a.flag = name, seq, flag
    a.reference_id, a.reference_start, a.mapping_quality = 0, int(start), 60
    a.cigartuples = cigar
    a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
    a.set_tag("NH", 1); a.set_tag("NM", 0); a.set_tag("RG", sample)
    return a


def sort_and_index(tmp, final):
    pysam.sort("-o", str(final), str(tmp))
    pysam.index(str(final))
    Path(tmp).unlink()


def simulate_sample(sample, counts, table, seq, rng):
    tmp, final = OUT / f"{sample}.unsorted.bam", OUT / f"{sample}.bam"
    k = 0
    with pysam.AlignmentFile(str(tmp), "wb", header=make_header(CONTIG, len(seq), sample)) as bam:
        for gene, n in counts.items():
            g = table.loc[gene]; span = g.end - g.start - READ_LEN
            # start positions: mixture of uniform and 3'-skewed coverage (looks less artificial than pure uniform)
            offs = np.where(rng.random(int(n)) < 0.7, rng.random(int(n)), rng.beta(2.5, 1.0, int(n)))
            starts = g.start + (offs * span).astype(int)
            flags = np.where(rng.random(int(n)) < 0.5, 0, 16)
            for s, f in zip(starts, flags):
                k += 1
                bam.write(new_read(f"sim_{sample}_{k}", seq[s:s + READ_LEN], s, int(f), [(0, READ_LEN)], sample))
    sort_and_index(tmp, final)
    return final


def simulate_spliced_demo(rng):
    """Toy 3-exon locus with spliced reads: only to demonstrate junction / Sashimi-style displays."""
    seq = random_seq(rng, 4400)
    exons = [(500, 900), (1700, 2000), (3000, 3500)]
    table = pd.DataFrame([dict(gene="toyGene", start=500, end=3500, strand="+")]).set_index("gene")
    write_fasta(OUT / "spliced_locus_demo.fa", SPLICE_CONTIG, seq)
    write_gtf(OUT / "spliced_locus_demo.gtf", SPLICE_CONTIG, table, {"toyGene": exons})
    sample = "SIM_SPLICE_DEMO"
    tmp, final = OUT / "spliced_locus_demo.unsorted.bam", OUT / "spliced_locus_demo.bam"
    reads, k = [], 0
    for ex_start, ex_end, n in [(e[0], e[1], m) for e, m in zip(exons, (3200, 2400, 4200))]:     # unspliced exonic reads
        for s in rng.integers(ex_start, ex_end - READ_LEN, n):
            reads.append((int(s), [(0, READ_LEN)], seq[s:s + READ_LEN]))
    for (d, a_end, n) in [(900, 1700, 230), (2000, 3000, 200), (900, 3000, 50)]:             # junction reads
        for anchor in rng.integers(10, READ_LEN - 10, n):
            s = d - int(anchor)
            reads.append((s, [(0, int(anchor)), (3, a_end - d), (0, READ_LEN - int(anchor))],
                          seq[s:d] + seq[a_end:a_end + READ_LEN - int(anchor)]))
    with pysam.AlignmentFile(str(tmp), "wb", header=make_header(SPLICE_CONTIG, len(seq), sample)) as bam:
        for s, cig, sq in reads:
            k += 1
            bam.write(new_read(f"sim_splice_{k}", sq, s, 0, cig, sample))
    sort_and_index(tmp, final)
    return len(reads)


def bam_gene_counts(bam_path, table):
    """featureCounts-like: reads lying inside a gene's span are assigned to it (no overlaps in this layout)."""
    out = {}
    with pysam.AlignmentFile(str(bam_path)) as bam:
        for gene, r in table.iterrows():
            out[gene] = sum(1 for a in bam.fetch(CONTIG, int(r.start), int(r.end))
                            if a.reference_start >= r.start and a.reference_end <= r.end)
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True); TAB.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    counts = pd.read_csv(ROOT / "data" / "count_matrix.csv", index_col="GeneID")
    seq, table = build_reference(rng, counts.index.tolist())
    write_fasta(OUT / "synthetic_reference.fa", CONTIG, seq)
    write_gtf(OUT / "synthetic_genes.gtf", CONTIG, table)
    table.to_csv(OUT / "synthetic_gene_coordinates.csv")
    print(f"reference: {len(seq):,} bp, {len(table)} genes")

    derived = {}
    for sample in counts.columns:
        bam = simulate_sample(sample, counts[sample], table, seq, rng)
        derived[sample] = bam_gene_counts(bam, table)
        print(f"  {bam.name}: {int(counts[sample].sum()):,} reads")
    derived = pd.DataFrame(derived).loc[counts.index]
    derived.index.name = "GeneID"
    derived.to_csv(TAB / "bam_derived_counts.csv")
    ok = derived.equals(counts.astype(derived.dtypes.iloc[0]))
    print(f"round trip (BAM-derived counts == count_matrix.csv): {ok}")
    if not ok:
        sys.exit("ERROR: BAM-derived counts do not match the count matrix")
    n = simulate_spliced_demo(rng)
    print(f"  spliced_locus_demo.bam: {n:,} reads")
    for f in OUT.glob("*.unsorted.bam"):
        os.remove(f)
    print("Done.")


if __name__ == "__main__":
    main()
