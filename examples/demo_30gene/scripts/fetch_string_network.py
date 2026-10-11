#!/usr/bin/env python3
"""Refresh the STRING network for current DEGs (run locally; needs internet).

STRING does not currently recognize the project strain SL1344 (taxon 216597), so the default
uses the supported Typhimurium LT2 reference (taxon 99287) to resolve the shared gene symbols.
"""
import io, json, sys, urllib.error, urllib.parse, urllib.request
from datetime import datetime, timezone
from pathlib import Path
import pandas as pd

root = Path(__file__).resolve().parents[1]
genes = pd.read_csv(root / "results/tables/significant_degs.csv")["GeneID"].tolist()
if not genes:
    raise ValueError("No significant DEGs found; run scripts/rna_seq_pipeline.py first.")
taxon = sys.argv[1] if len(sys.argv) > 1 else "99287"
if not taxon.isdigit():
    raise ValueError("STRING species taxon ID must be a numeric NCBI taxon identifier.")
q = urllib.parse.urlencode({"identifiers": "\r".join(genes), "species": taxon, "required_score": 400, "add_nodes": 15})
url = f"https://string-db.org/api/tsv/network?{q}"
try:
    with urllib.request.urlopen(url) as response:
        raw = response.read().decode()
except urllib.error.HTTPError as error:
    detail = error.read().decode(errors="replace").strip()
    raise RuntimeError(f"STRING request failed for taxon {taxon}: {detail}") from error

if not raw.strip():
    raise RuntimeError(f"STRING returned no network edges for taxon {taxon}.")
df = pd.read_csv(io.StringIO(raw), sep="\t")
required_columns = {"preferredName_A", "preferredName_B", "score"}
if not required_columns.issubset(df.columns):
    raise RuntimeError(f"STRING response is missing required columns: {sorted(required_columns - set(df.columns))}")
out = df[["preferredName_A", "preferredName_B", "score"]].drop_duplicates()
out.columns = ["source", "target", "score"]
if out.empty:
    raise RuntimeError(f"STRING returned no usable network edges for taxon {taxon}.")
if not (set(genes) & set(out.source) or set(genes) & set(out.target)):
    raise RuntimeError(f"STRING response for taxon {taxon} contains none of the requested DEG symbols.")
out.to_csv(root / "data/string_ppi_network_edges.csv", index=False)
metadata = {
    "source": "STRING API",
    "endpoint": "https://string-db.org/api/tsv/network",
    "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    "species_taxon_id": taxon,
    "species": ("Salmonella enterica serovar Typhimurium LT2"
                if taxon == "99287" else f"NCBI taxon {taxon}"),
    "requested_gene_symbols": genes,
    "required_score": 400,
    "first_shell_interactors_requested": 15,
    "edge_count": len(out),
}
(root / "data/string_network_metadata.json").write_text(
    json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
)
print(f"Saved {len(out)} STRING edges for {len(genes)} current DEGs plus first-shell interactors "
      f"(taxon {taxon}).")
print(out)
