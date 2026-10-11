"""Helpers for identifying Hfq structures from PDB header records."""
import re


PDB_HEADER_RECORDS = ("HEADER", "TITLE", "COMPND")
HFQ_PATTERN = re.compile(r"\bHFQ\b|\bHOST FACTOR(?:-I| I)\b|\bHF-[I1]\b")


def pdb_header_text(pdb_text):
    return "\n".join(
        line[10:].strip()
        for line in pdb_text.splitlines()
        if line.startswith(PDB_HEADER_RECORDS)
    )


def is_hfq_pdb(pdb_text):
    return bool(HFQ_PATTERN.search(pdb_header_text(pdb_text).upper()))
