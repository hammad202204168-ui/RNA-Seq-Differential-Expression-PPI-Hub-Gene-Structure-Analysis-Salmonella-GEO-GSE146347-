#!/usr/bin/env bash
# =====================================================================================================
# Rebuild RAW per-gene counts for the 12 bacterial samples of GSE146347 (BioProject PRJNA610228).
#
# STATUS: bash syntax-checked only (bash -n). NOT executed end-to-end: the environment that wrote it had no
# access to NCBI/ENA and no aligner. The Python steps it calls (featurecounts_to_matrix.py) are unit-tested.
# Run it on a machine with internet, ~10 GB disk and: STAR, samtools, subread (featureCounts), fastp, curl;
# optional sra-tools (prefetch, fasterq-dump) if you prefer SRA over ENA FASTQ.
#
#   bash scripts/rebuild_counts_from_sra.sh                  # LT2 reference (what the authors used: GCF_000006945.2)
#   REF=SL1344 GENOME_FA=... GFF=... bash scripts/rebuild_counts_from_sra.sh   # second pass, strain-matched
#
# Why two references: the sequenced strain is chi3339 (a mouse-passaged SL1344 derivative) but the published TPMs
# were mapped to LT2. Genes missing/divergent in LT2 look like "dropout" or "WT-only expression" whatever the
# biology. Compare the STAR unique-mapping rates printed at the end (data/rebuild/mapping_summary_<REF>.csv):
# a clear gain on SL1344 for specific samples means part of the dropout was a reference artefact. Verify the
# SL1344 accession set yourself before use; none is hard-coded here.
#
# STRAND: the library prep in the paper is NuGEN Ovation RNA-Seq System V2 + KAPA ligation, which is not strand-
# specific by design, so STRAND=0 is the default. Confirm with RSeQC infer_experiment.py before trusting it.
# =====================================================================================================
set -euo pipefail
REF="${REF:-LT2}"; THREADS="${THREADS:-8}"; STRAND="${STRAND:-0}"; FEATURE="${FEATURE:-gene}"; ATTR="${ATTR:-locus_tag}"
OUT="data/rebuild/${REF}"; mkdir -p "$OUT"/{fastq,trim,bam,index} data/rebuild

if [[ "$REF" == "LT2" && -z "${GENOME_FA:-}" ]]; then
  base=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2   # check this path if curl 404s
  curl -fL "$base/GCF_000006945.2_ASM694v2_genomic.fna.gz" | gunzip > "$OUT/genome.fa"
  curl -fL "$base/GCF_000006945.2_ASM694v2_genomic.gff.gz" | gunzip > "$OUT/genome.gff"
  GENOME_FA="$OUT/genome.fa"; GFF="$OUT/genome.gff"
fi
: "${GENOME_FA:?set GENOME_FA (and GFF) for REF=$REF}"; : "${GFF:?set GFF}"

STAR --runMode genomeGenerate --runThreadN "$THREADS" --genomeDir "$OUT/index" \
     --genomeFastaFiles "$GENOME_FA" --genomeSAindexNbases 10             # ~5 Mb genome: min(14, log2(5e6)/2-1) = 10

echo -e "srr\tinput_reads\tuniquely_mapped_pct\tmulti_mapped_pct\tunmapped_too_short_pct" > "data/rebuild/mapping_summary_${REF}.tsv"
while read -r SRR; do
  if command -v fasterq-dump >/dev/null; then
    prefetch "$SRR" -O "$OUT/fastq" && fasterq-dump "$OUT/fastq/$SRR" -O "$OUT/fastq" -e "$THREADS"
  else   # ENA mirror: ask ENA for the exact URL instead of guessing its directory scheme
    url=$(curl -fsL "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${SRR}&result=read_run&fields=fastq_ftp&format=tsv" | tail -n1)
    curl -fL "https://${url%%;*}" | gunzip > "$OUT/fastq/${SRR}.fastq"
  fi
  fastp -i "$OUT/fastq/${SRR}.fastq" -o "$OUT/trim/${SRR}.fq" -w "$THREADS" --json "$OUT/trim/${SRR}.json" --html /dev/null \
        --adapter_sequence AGATCGGAAGAGC --length_required 25
  STAR --runThreadN "$THREADS" --genomeDir "$OUT/index" --readFilesIn "$OUT/trim/${SRR}.fq" \
       --outFileNamePrefix "$OUT/bam/${SRR}_" --outSAMtype BAM SortedByCoordinate --alignIntronMax 1
  samtools index "$OUT/bam/${SRR}_Aligned.sortedByCoord.out.bam"
  log="$OUT/bam/${SRR}_Log.final.out"
  printf '%s\t%s\t%s\t%s\t%s\n' "$SRR" \
    "$(awk -F'|\t' '/Number of input reads/{gsub(/ /,"",$2);print $2}' "$log")" \
    "$(awk -F'|\t' '/Uniquely mapped reads %/{gsub(/[ %]/,"",$2);print $2}' "$log")" \
    "$(awk -F'|\t' '/% of reads mapped to multiple loci/{gsub(/[ %]/,"",$2);print $2}' "$log")" \
    "$(awk -F'|\t' '/% of reads unmapped: too short/{gsub(/[ %]/,"",$2);print $2}' "$log")" >> "data/rebuild/mapping_summary_${REF}.tsv"
  rm -f "$OUT/fastq/${SRR}.fastq" "$OUT/fastq/${SRR}.sra"; rm -rf "$OUT/fastq/${SRR}"
done < data/sra/SRR_Acc_List_bacterial.txt

featureCounts -T "$THREADS" -s "$STRAND" -t "$FEATURE" -g "$ATTR" -a "$GFF" -o "$OUT/featureCounts_raw.txt" "$OUT"/bam/*_Aligned.sortedByCoord.out.bam
echo "Done. Next:"
echo "  python scripts/collect_batch_evidence.py            # needs internet; writes data/sra/column_to_run_map.csv"
echo "  python scripts/featurecounts_to_matrix.py $OUT/featureCounts_raw.txt data/sra/column_to_run_map.csv [--batch_map data/sra/batch_assignment.csv] --out_counts data/rebuild/counts_${REF}.csv --out_meta data/rebuild/meta_${REF}.csv"
echo "  python scripts/batch_aware_de.py --counts data/rebuild/counts_${REF}.csv --meta data/rebuild/meta_${REF}.csv --outdir results/count_based_${REF}"
