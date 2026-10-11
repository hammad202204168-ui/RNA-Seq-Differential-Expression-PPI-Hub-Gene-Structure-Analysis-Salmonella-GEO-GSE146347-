# Separating genotype from batch, and handling dropout: what is done, what is blocked

## Status
| Item | Status | Where |
|---|---|---|
| Quantify how much genotype information survives a `_2` batch term | **Done** (real design) | `results/tables/design_identifiability.csv` |
| Batch-aware count model (PyDESeq2 `~ batch + group`), concordance with the unadjusted model | **Written, validated on simulated counts only** | `scripts/batch_aware_de.py`, `scripts/simulate_validation.py`, `results/tables/sim_validation_*.csv` |
| Dropout diagnosis and detection limits from counts | **Written, validated on simulated counts only** | same |
| Meaning of `_2` (authors / GEO / read headers) | **Blocked: not obtainable here.** The paper's Methods describe one library-prep protocol for all bacterial samples and never mention `_2`; GEO sample pages were not reachable | `scripts/collect_batch_evidence.py` (run with internet) |
| Raw counts from SRA | **Blocked: no NCBI/ENA access, no aligner in the build environment** | `scripts/rebuild_counts_from_sra.sh` (syntax-checked only) |

Nothing in `results/` other than the three tables above is changed by this work, and no result in the README was re-derived from counts.

## What the real design allows (batch = `_2` proxy)
Only `hfq1_1xg` is a Δhfq sample without `_2`, so once a batch term is fitted the Δhfq-vs-WT contrasts are informed by that one sample.

| Contrast | SE inflation | Effective n per group (nominal 3) |
|---|---:|---:|
| Δhfq vs WT, control | 1.41 | 1.5 |
| Δhfq vs WT, LFS | 1.80 | 0.9 |
| LFS vs control, WT | 1.00 | 3.0 |
| LFS vs control, Δhfq | 1.12 | 2.4 |

This is a property of the design. Raw counts and a correct batch label do not change it; only new samples (WT prepared in the `_2` batch, or Δhfq prepared outside it) would.

## What raw counts can and cannot do for "dropped-out" genes
* **Cannot:** turn an unsequenced transcript into a measurement. A zero at 0.5-1 M reads is information about a bound, not a value.
* **Can:** (1) replace log2(TPM+1) *t*-tests with a depth-aware negative-binomial model; (2) tell, per gene and group, whether all-zero observations are explainable by depth alone (`depth_limited_uninformative`) or not (`zero_not_explained_by_depth`, i.e. genuinely lower or absent); (3) give a 95 % upper detection limit of expression for an all-zero group; (4) test whether part of the dropout is a reference artefact by mapping to an SL1344-type reference as well as LT2 and comparing unique-mapping rates.

## Simulation results (3 seeds, `scripts/simulate_validation.py`)
Simulated counts, 12 samples, the real design and the real per-group read depths, genotype effect on 4 % of genes, batch effect on 12 %, 3 % of genes truly off in one WT group. Means over seeds, genotype contrast in control culture:

| | no batch term | batch term |
|---|---:|---:|
| Genes with a batch effect only, falsely called genotype-DE | 9.7 (7.7 %) | 0.3 (0.3 %) |
| Sensitivity for true genotype genes (FDR < 0.05, \|log2FC\| ≥ 1) | 31 % | 14.5 % |
| Median SE of log2FC | 0.49 | 0.68 |

Dropout diagnosis: 98.6 % of truly absent group-level zeros were flagged `zero_not_explained_by_depth`; 89 % of all flagged cells were truly below half of the gene's best-group expression (the rest are depth-borderline); the 95 % upper limit covered the true expression in 91 % of expressed all-zero cells (slightly liberal; treat as approximate). These are properties of the *method under the simulation's assumptions*, not of the real data.

## Decision table once the evidence is in
| What `batch_evidence.csv` / the authors say | Do this |
|---|---|
| `_2` = different flowcell/run or library-prep date for the Δhfq samples | Keep the batch term; report Δhfq-vs-WT as essentially unresolved; rely on LFS-within-genotype contrasts. Only new data fixes the rest. |
| `_2` = a second independent culture experiment | Same statistics (it is still nested in genotype); state it as a limitation. |
| `_2` is only a naming convention (no protocol, run or date differences) | Drop the batch term (`--meta` with a single batch level); keep detection/depth covariates; the rRNA and depth differences still need the count-based model. |
| Evidence silent | Ask the authors (corresponding authors: J. Barrila, C. A. Nickerson, ASU). Do not assume. |

## Run order (needs internet)
```bash
python scripts/collect_batch_evidence.py                 # GEO + ENA + read headers -> data/sra/batch_evidence.csv, column_to_run_map.csv
# review batch_evidence.csv; write data/sra/batch_assignment.csv (tpm_column,batch) if a batch is documented
bash scripts/rebuild_counts_from_sra.sh                  # LT2; then REF=SL1344 with GENOME_FA/GFF
python scripts/featurecounts_to_matrix.py data/rebuild/LT2/featureCounts_raw.txt data/sra/column_to_run_map.csv \
       --batch_map data/sra/batch_assignment.csv --out_counts data/rebuild/counts_LT2.csv --out_meta data/rebuild/meta_LT2.csv
python scripts/batch_aware_de.py --counts data/rebuild/counts_LT2.csv --meta data/rebuild/meta_LT2.csv --outdir results/count_based_LT2
```
