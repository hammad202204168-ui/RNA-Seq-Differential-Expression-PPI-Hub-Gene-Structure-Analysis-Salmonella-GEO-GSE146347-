#!/usr/bin/env python3
"""Inspect a GEO series before building a count matrix (run locally; needs internet).

    python scripts/geo_helper.py files   GSE146347    # list supplementary files on the NCBI GEO server
    python scripts/geo_helper.py samples GSE146347    # sample titles/characteristics from the family SOFT file
                                                      # (also saved to data/geo_samples_<acc>.csv)

It only lists and parses metadata; it never builds a matrix or assigns conditions for you. Decide yourself which
samples are bacterial, which files hold raw (not normalised) counts, and map them to conditions.
"""
import gzip
import re
import sys
import urllib.request
from pathlib import Path

import pandas as pd

BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{acc}/"


def series_url(acc):
    return BASE.format(stub=acc[:-3] + "nnn", acc=acc)


def fetch(url):
    with urllib.request.urlopen(url, timeout=90) as r:
        return r.read()


def list_supplementary(acc):
    html = fetch(series_url(acc) + "suppl/").decode("utf-8", "replace")
    return [h for h in re.findall(r'href="([^"?/][^"]*)"', html) if not h.startswith(("http", "mailto"))]


def parse_soft_samples(text):
    """Parse '^SAMPLE' blocks of a GEO SOFT file into one dict per sample."""
    samples, cur = [], None
    for line in text.splitlines():
        if line.startswith("^SAMPLE"):
            cur = {"accession": line.split("=", 1)[1].strip()}
            samples.append(cur)
        elif cur is not None and line.startswith("!Sample_"):
            key, _, val = line[len("!Sample_"):].partition("=")
            key, val = key.strip(), val.strip()
            if key.startswith("characteristics"):
                ck, _, cv = val.partition(":")
                cur[f"char_{ck.strip()}" if cv else "char_other"] = (cv or ck).strip()
            elif key in ("title", "organism_ch1", "source_name_ch1", "library_strategy"):
                cur[key] = val
    return pd.DataFrame(samples)


def main():
    if len(sys.argv) != 3 or sys.argv[1] not in ("files", "samples"):
        sys.exit(__doc__)
    cmd, acc = sys.argv[1], sys.argv[2].upper()
    if cmd == "files":
        for f in list_supplementary(acc):
            print(series_url(acc) + "suppl/" + f)
    else:
        raw = fetch(series_url(acc) + f"soft/{acc}_family.soft.gz")
        df = parse_soft_samples(gzip.decompress(raw).decode("utf-8", "replace"))
        out = Path(__file__).resolve().parents[1] / "data" / f"geo_samples_{acc}.csv"
        df.to_csv(out, index=False)
        print(df.to_string(index=False)); print(f"\nsaved {out}")


if __name__ == "__main__":
    main()
