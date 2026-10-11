.PHONY: install run ppi test sim demo clean

install:
	pip install -r requirements.txt

run:
	python scripts/prepare_geo_tpm.py && python scripts/real_data_analysis.py

ppi:
	python scripts/ppi_hubs_structures.py \
	  --string_edges data/string/string_interactions_short_135_nodes.tsv \
	  --string_annotations data/string/string_protein_annotations_135_nodes.tsv \
	  --uniprot_map data/string/uniprot_idmapping_hub_genes.tsv \
	  --structures_dir data/alphafold

test:
	python scripts/prepare_geo_tpm.py && python scripts/real_data_analysis.py && pytest -q

sim:
	python scripts/simulate_validation.py

demo:
	cd examples/demo_30gene && make test

clean:
	find . -name __pycache__ -type d -prune -exec rm -rf {} +
	rm -rf .pytest_cache
