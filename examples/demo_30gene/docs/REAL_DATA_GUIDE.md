# Moving from the demonstration matrix to real data

## Verified candidate source
**GEO GSE146347** (NASA GeneLab GLDS-277) is a real *S.* Typhimurium wild type vs Δ*hfq* RNA-seq study:
Barrila J *et al.* (2022) *Front Cell Infect Microbiol*, [doi:10.3389/fcimb.2022.705647](https://doi.org/10.3389/fcimb.2022.705647).
Its bacterial arm compares stationary-phase wild type and Δ*hfq* cells grown under low-shear modeled microgravity (LSMMG) or control culture; the series
also contains **human host** RNA-seq from 3-D co-cultures, which must be excluded.

What was and was not checked: the series' existence, organisms, title and abstract were confirmed through the OmicsDI mirror of the GEO record and the paper's
data-availability statement. **The GEO page, sample sheet and count files were not opened, and no counts were retrieved.**
This study is **not** the source of `data/count_matrix.csv`, and its design (LSMMG vs control) has no acid-stress or heat-shock arms.

## Steps (on a machine with internet)
1. `python scripts/geo_helper.py files GSE146347` to list supplementary files, and `python scripts/geo_helper.py samples GSE146347` to list samples
   (saved to `data/geo_samples_GSE146347.csv`). Identify the bacterial samples, their strain (wild type / Δ*hfq*) and culture (control / LSMMG).
2. Obtain **raw integer counts** per gene for the bacterial samples (a supplied counts file, or re-count from the SRA runs against the reference genome and
   annotation named in the paper's Methods; do not use TPM/FPKM/normalised values: DESeq2 needs raw counts).
3. Write `data/count_matrix.csv` (genes x samples, first column `GeneID`) and `data/sample_metadata.csv` (`SampleID`, `Condition`, ...).
4. Adapt the analysis to the new design. These are hard-coded for the demonstration and **must be changed**:
   * `COND_COL` / `COND_LAB` and the allowed `Condition` values in `scripts/rna_seq_pipeline.py` (currently Control, Mutant, AcidStress, HeatShock);
   * the contrast list in `main()` (currently Mutant, AcidStress, HeatShock vs Control); for a 2 x 2 design consider `~ strain + culture + strain:culture`;
   * gene names used in `fig_summary`, `fig_alignment_view` and `GENE_SETS` (30-gene panel), plus the STRING gene list (`scripts/fetch_string_network.py`);
   * replace the Galaxy-style panel-only enrichment with genome-wide GO/KEGG enrichment appropriate for the reference annotation.
5. Delete `data/synthetic_alignments/` (or regenerate from real BAMs) and set `DEMO_NOTE = None` **only after** provenance is documented in
   `data/source_accessions.csv` (add the run-to-sample mapping) and the README caution is replaced by a data description.
6. Re-run `python scripts/rna_seq_pipeline.py` and `pytest`; update `tests/` expectations (they currently assert demonstration-matrix values).

Reference genome: this repository does not assume a strain. Use the strain and assembly stated in the paper's Methods; `NC_016810.1` (SL1344) in
`source_accessions.csv` is only a reference-genome record, not a claim about GSE146347.
