.PHONY: run bam test online
run:    ; python scripts/rna_seq_pipeline.py
bam:    ; python scripts/simulate_alignments.py && python scripts/plot_alignments.py
test:   ; python scripts/rna_seq_pipeline.py && python scripts/simulate_alignments.py && python scripts/plot_alignments.py && pytest -q
online: ; bash scripts/update_with_internet.sh
