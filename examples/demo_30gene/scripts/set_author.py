#!/usr/bin/env python3
"""python scripts/set_author.py "Given Family" -> fills LICENSE and CITATION.cff"""
import sys, re
from pathlib import Path
root = Path(__file__).resolve().parents[1]
name = sys.argv[1]; given, _, family = name.rpartition(" ")
lic = root / "LICENSE"; lic.write_text(re.sub(r"Copyright \(c\) \d{4} .*", f"Copyright (c) 2026 {name}", lic.read_text()))
cff = root / "CITATION.cff"
cff.write_text(cff.read_text().replace("YOUR-LASTNAME", family).replace("YOUR-FIRSTNAME", given or family))
print("updated LICENSE and CITATION.cff for", name)
