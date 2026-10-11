#!/usr/bin/env python3
"""Download a validated Hfq structure into structure/ (run locally; needs internet)."""
import sys, urllib.request
from pathlib import Path

from structure_utils import is_hfq_pdb, pdb_header_text

pdb_id = sys.argv[1] if len(sys.argv) > 1 else "1HK9"
out = Path(__file__).resolve().parents[1] / "structure" / f"{pdb_id}.pdb"
out.parent.mkdir(exist_ok=True)
with urllib.request.urlopen(f"https://files.rcsb.org/download/{pdb_id}.pdb") as response:
    pdb_text = response.read().decode("ascii")

header = pdb_header_text(pdb_text)
print(header)
if not is_hfq_pdb(pdb_text):
    raise ValueError(f"PDB {pdb_id} header does not identify Hfq; refusing to save it.")

out.write_text(pdb_text, encoding="ascii")
print(f"saved validated Hfq structure: {out}")
