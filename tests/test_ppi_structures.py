"""Offline tests for scripts/ppi_hubs_structures.py: parsing, hubs, pLDDT, plotting. Network calls are NOT tested."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import numpy as np, pandas as pd
import ppi_hubs_structures as P

IDS = "queryItem\tqueryNumber\tstringId\tncbiTaxonId\ttaxonName\tpreferredName\tannotation\nflgG\t0\t99287.STM1\t99287\tSalmonella\tflgG\tx\nflgL\t1\t99287.STM2\t99287\tS\tflgL\tx\nflgG\t0\t99287.STM1\t99287\tS\tflgG\tdup\n"
EDGES = ("stringId_A\tstringId_B\tpreferredName_A\tpreferredName_B\tncbiTaxonId\tscore\n"
         "a\tb\tflgG\tflgL\t99287\t0.99\nb\ta\tflgL\tflgG\t99287\t0.99\na\tc\tflgG\tcheA\t99287\t0.9\na\td\tflgG\tfliC\t99287\t0.8\nc\td\tcheA\tfliC\t99287\t0.7\ne\tf\tx1\tx2\t99287\t0.9\n")
UNI = "Entry\tReviewed\tGene Names\tProtein names\nQ1TREMBL\tunreviewed\tflgG\tFlagellar rod\nP0REV1\treviewed\tflgG STM1197\tFlagellar basal-body rod protein FlgG\n"

def pdb(n=6, b=95.0):
    L = []
    for i in range(n):
        L.append("ATOM  %5d  CA  ALA A%4d    %8.3f%8.3f%8.3f  1.00%6.2f           C" % (i + 1, i + 1, i * 3.8, 0.0, 0.0, b - 10 * i))
    return "\n".join(L)

def test_string_parsing_and_unresolved():
    ids = P.parse_string_ids(IDS, ["flgG", "flgL", "nope"])
    assert ids.stringId.notna().sum() == 2 and ids.loc[ids.gene == "nope", "stringId"].isna().all() and len(ids) == 3

def test_edges_deduplicated_and_hubs():
    e = P.parse_string_edges(EDGES)
    assert len(e) == 5                                   # reverse duplicate dropped
    h, G = P.hubs_from_edges(e, top=2)
    assert h.gene[0] == "flgG" and int(h.degree[0]) == 3 and h.hub.sum() == 2
    assert G.number_of_nodes() == 6

def test_uniprot_prefers_reviewed():
    assert P.pick_uniprot(UNI)["uniprot"] == "P0REV1" and P.pick_uniprot("") is None

def test_plddt():
    xyz, pl = P.plddt_from_pdb(pdb())
    assert xyz.shape == (6, 3) and pl[0] == 95.0
    s = P.plddt_summary(pl)
    assert s["n_residues"] == 6 and s["pct_residues_below_70"] == 50.0   # 95,85,75,65,55,45 -> 65,55,45 are below 70


def test_figures_render(tmp_path):
    e = P.parse_string_edges(EDGES)
    h, G = P.hubs_from_edges(e, top=2)
    P.fig_network(G, h, str(tmp_path / "n.png"), 2)
    xyz, pl = P.plddt_from_pdb(pdb())
    P.fig_structures([("flgG", xyz, pl), ("cheA", xyz, pl)], str(tmp_path / "s.png"))
    assert (tmp_path / "n.png").stat().st_size > 5000 and (tmp_path / "s.png").stat().st_size > 5000


ANN = pd.DataFrame({"node": ["invJ", "flgG", "mopB"], "other_names_and_aliases": ["spaN,STM2896", "STM1197", "groES,STM4332"]})

def test_alias_mapping():
    m = P.map_inputs_to_nodes(["flgG", "spaN", "groES", "pefI"], ANN).set_index("input_gene")
    assert m.loc["flgG", "match"] == "exact" and m.loc["spaN", "string_node"] == "invJ" and m.loc["groES", "match"] == "alias"
    assert m.loc["pefI", "match"] == "not_found_in_STRING" and pd.isna(m.loc["pefI", "string_node"])

def _edge(a, b, **ch):
    d = {c: 0.0 for c in P.STRING_CHANNELS}; d.update(ch)
    d.update(node1=a, node2=b); return d

def test_combined_without_textmining_drops_textmining_only_edges():
    tm = pd.Series(_edge("a", "b", automated_textmining=0.95)); ex = pd.Series(_edge("a", "c", experimentally_determined_interaction=0.95))
    assert P.combined_without(tm, []) > 0.9 and P.combined_without(tm, ["automated_textmining"]) < 0.1
    assert P.combined_without(ex, ["automated_textmining"]) > 0.9

def test_hub_robustness_changes_hub_when_textmining_removed():
    rows = [_edge("tm_hub", f"x{i}", automated_textmining=0.9) for i in range(5)] + [_edge("ex_hub", f"y{i}", experimentally_determined_interaction=0.9) for i in range(3)]
    e = pd.DataFrame(rows); e["combined_score"] = [P.combined_without(r, []) for _, r in e.iterrows()]
    rob = P.hub_robustness(e, top=1)
    assert rob["all_channels"].gene[0] == "tm_hub" and rob["no_textmining"].gene[0] == "ex_hub"


def test_find_structure_file_and_offline_summary(tmp_path):
    (tmp_path / "AF-P06179-F1-model_v4.pdb").write_text(pdb(10, 98.0))
    (tmp_path / "AF-P0A1I3-F1-model_v4.cif").write_text("not used")
    assert P.find_structure_file(str(tmp_path), "P06179").endswith("model_v4.pdb")
    assert P.find_structure_file(str(tmp_path), "P0A1I3") is None          # only PDB is read, CIF is ignored
    xyz, pl = P.plddt_from_pdb(open(P.find_structure_file(str(tmp_path), "P06179")).read())
    assert len(pl) == 10 and P.plddt_summary(pl)["n_residues"] == 10
