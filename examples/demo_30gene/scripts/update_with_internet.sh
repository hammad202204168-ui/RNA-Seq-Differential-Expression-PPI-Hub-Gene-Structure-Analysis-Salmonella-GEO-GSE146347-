#!/usr/bin/env bash
# Run once on a machine with internet access: refreshes STRING network + Hfq structure, then re-runs everything.
set -euo pipefail
cd "$(dirname "$0")/.."
PYTHON="${PYTHON:-python}"
"$PYTHON" scripts/rna_seq_pipeline.py            # (1) DEGs -> results/tables/significant_degs.csv
"$PYTHON" scripts/fetch_string_network.py        # (2) STRING network for the current DEGs
"$PYTHON" scripts/fetch_hfq_structure.py 1HK9    # (3) Hfq structure (check the printed header!)
"$PYTHON" scripts/rna_seq_pipeline.py            # (4) regenerate hub analysis + Fig10
