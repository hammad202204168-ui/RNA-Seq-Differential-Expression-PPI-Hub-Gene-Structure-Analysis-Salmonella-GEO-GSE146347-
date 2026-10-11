#!/usr/bin/env python3
"""Collect first-hand evidence for what the `_2` suffix means, and link TPM columns to GSM/SRR runs.

REQUIRES INTERNET (NCBI GEO + EBI ENA). NOT executed in the sandbox that produced this file: the parsers are
unit-tested offline (tests/test_collect_parsers.py); the network calls are not.

Evidence gathered per run (SRR11235557-68):
  1. GEO sample record (GSM4378545-56): title (-> TPM column name), protocols, library fields, SRA relation.
  2. ENA run record: instrument, first_created, library_name, read_count.
  3. FASTQ read names (first N reads): instrument / run / flowcell / lane. Different flowcell or run number
     = different sequencing batch. SRA-reformatted names (@SRR....1) carry no flowcell, which is reported
     as such rather than guessed; then use `fastq-dump -X N --origfmt` (needs sra-tools) via --use_sra_tools.
  4. Read-content QC that tracks library preparation: length distribution, adapter fraction, GC, duplication,
     5'-end composition.

Outputs (data/sra/):  column_to_run_map.csv   batch_evidence.csv   geo_samples_raw.json
IMPORTANT: this script reports evidence. It does not decide that `_2` == a batch. If the evidence is silent
(identical flowcell/protocol for all runs), ask the authors (Barrila / Nickerson lab, ASU) what `_2` denotes.
"""
import argparse, gzip, io, json, re, subprocess, sys, urllib.request
import collections
import pandas as pd

GSM_RANGE = [f"GSM{n}" for n in range(4378545, 4378557)]
ADAPTER = "AGATCGGAAGAGC"
CASAVA = re.compile(r"^@?(?P<instrument>[A-Za-z0-9_-]+):(?P<run>\d+):(?P<flowcell>[A-Za-z0-9-]+):(?P<lane>\d+):(?P<tile>\d+):\d+:\d+")


# ---------------------------------------------------------------------------------------------- parsers
def parse_soft(text):
    """Parse one GEO SOFT sample record -> dict of lists keyed by field name (without '!Sample_')."""
    d = collections.defaultdict(list)
    for line in text.splitlines():
        if line.startswith("!Sample_") and " = " in line:
            k, v = line[len("!Sample_"):].split(" = ", 1)
            d[k].append(v.strip())
    return dict(d)


def srr_from_soft(rec):
    for v in rec.get("relation", []):
        m = re.search(r"SRX\d+", v)
        if m:
            return m.group(0)
    return None


def parse_read_name(name):
    """Return dict(instrument, run, flowcell, lane) from an Illumina CASAVA-1.8+ read name, else None."""
    m = CASAVA.match(name.split()[0])
    return {k: m.group(k) for k in ("instrument", "run", "flowcell", "lane")} if m else None


def fastq_qc(lines, n_reads=50000):
    """Content QC from FASTQ text lines (iterable). Returns dict."""
    names, seqs = [], []
    it = iter(lines)
    for header in it:
        seq = next(it, "").strip(); next(it, ""); next(it, "")
        names.append(header.strip()); seqs.append(seq)
        if len(seqs) >= n_reads:
            break
    if not seqs:
        return {}
    lens = [len(s) for s in seqs]
    gc = sum(s.count("G") + s.count("C") for s in seqs) / max(sum(lens), 1)
    first = collections.Counter(s[:12] for s in seqs)
    parsed = [parse_read_name(n) for n in names[:1000]]
    parsed = [p for p in parsed if p]
    ids = collections.Counter((p["instrument"], p["run"], p["flowcell"], p["lane"]) for p in parsed)
    return dict(n_reads_sampled=len(seqs), mean_len=sum(lens) / len(lens), frac_len_lt_50=sum(l < 50 for l in lens) / len(lens),
                adapter_frac=sum(ADAPTER in s for s in seqs) / len(seqs), gc=gc,
                dup_frac=1 - len(set(seqs)) / len(seqs), top_12mer_frac=first.most_common(1)[0][1] / len(seqs),
                header_parsed_frac=len(parsed) / min(len(names), 1000),
                flowcells=";".join(f"{i}:{r}:{f}:L{l}" for (i, r, f, l), _ in ids.most_common(5)) or "NA (no Illumina read names stored)")


def normalise(x):
    return re.sub(r"[^a-z0-9]", "", x.lower())


def link_columns(tpm_columns, titles_by_gsm):
    """Map TPM column -> GSM by normalised title. Returns (dict, unmatched_columns)."""
    tnorm = {normalise(t): g for g, t in titles_by_gsm.items()}
    out, miss = {}, []
    for c in tpm_columns:
        hit = tnorm.get(normalise(c))
        if hit is None:  # titles may omit the trailing _2 or add a prefix; accept unique containment
            cand = [g for t, g in tnorm.items() if normalise(c) in t or t in normalise(c)]
            hit = cand[0] if len(cand) == 1 else None
        (out.__setitem__(c, hit) if hit else miss.append(c))
    return out, miss


# ---------------------------------------------------------------------------------------------- network
def get(url, timeout=60):
    with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "hfq-analysis/1.0"}), timeout=timeout) as r:
        return r.read()


def geo_record(gsm):
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}&targ=self&form=text&view=full"
    return parse_soft(get(url).decode("utf8", "replace"))


def ena_run(srr):
    f = "run_accession,experiment_accession,sample_alias,library_name,instrument_model,read_count,first_created,fastq_ftp"
    url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={srr}&result=read_run&fields={f}&format=tsv"
    rows = [l.split("\t") for l in get(url).decode().strip().splitlines()]
    return dict(zip(rows[0], rows[1])) if len(rows) > 1 else {}


def sample_fastq(ftp, n_reads):
    """Stream the head of a gzipped FASTQ over HTTP (no full download)."""
    url = "https://" + ftp.split(";")[0] if not ftp.startswith("http") else ftp
    req = urllib.request.Request(url, headers={"User-Agent": "hfq-analysis/1.0", "Range": f"bytes=0-{n_reads * 120}"})
    raw = urllib.request.urlopen(req, timeout=120).read()
    try:
        text = gzip.GzipFile(fileobj=io.BytesIO(raw)).read().decode("ascii", "replace")
    except (EOFError, OSError):  # truncated gzip is expected with Range requests
        d = gzip.decompressobj(31)
        text = d.decompress(raw).decode("ascii", "replace")
    return text.splitlines()[: n_reads * 4]


def sample_fastq_sra_tools(srr, n_reads):
    out = subprocess.run(["fastq-dump", "-X", str(n_reads), "--origfmt", "-Z", srr], capture_output=True, text=True, check=True).stdout
    return out.splitlines()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run_table", default="data/sra/SraRunTable_GSE146347.csv")
    ap.add_argument("--tpm", default="data/processed/tpm_matrix.csv")
    ap.add_argument("--outdir", default="data/sra")
    ap.add_argument("--n_reads", type=int, default=50000)
    ap.add_argument("--use_sra_tools", action="store_true", help="read original read names with fastq-dump --origfmt")
    a = ap.parse_args()

    rt = pd.read_csv(a.run_table)
    rt = rt[rt.Run.between("SRR11235557", "SRR11235568")].set_index("GEO_Accession (exp)")
    recs = {}
    for gsm in GSM_RANGE:
        recs[gsm] = geo_record(gsm)
        print(gsm, recs[gsm].get("title"))
    json.dump(recs, open(f"{a.outdir}/geo_samples_raw.json", "w"), indent=1)
    titles = {g: r.get("title", [""])[0] for g, r in recs.items()}

    cols = pd.read_csv(a.tpm, nrows=1, index_col=0).columns.tolist()
    mapping, missing = link_columns(cols, titles)
    rows = []
    for gsm in GSM_RANGE:
        srr = rt.loc[gsm, "Run"]
        ena = ena_run(srr)
        lines = sample_fastq_sra_tools(srr, a.n_reads) if a.use_sra_tools else sample_fastq(ena["fastq_ftp"], a.n_reads)
        qc = fastq_qc(lines, a.n_reads)
        rec = recs[gsm]
        rows.append(dict(gsm=gsm, srr=srr, geo_title=titles[gsm], tpm_column=next((c for c, g in mapping.items() if g == gsm), None),
                         genotype=rt.loc[gsm, "bacterial_genotype"], orientation=rt.loc[gsm, "bacterial_growth_orientation"],
                         bases=rt.loc[gsm, "Bases"], instrument=ena.get("instrument_model"), ena_first_created=ena.get("first_created"),
                         library_name=ena.get("library_name"), extract_protocol=" | ".join(rec.get("extract_protocol_ch1", []))[:300],
                         data_processing=" | ".join(rec.get("data_processing", []))[:300],
                         characteristics=" | ".join(rec.get("characteristics_ch1", [])), **qc))
    ev = pd.DataFrame(rows)
    ev.to_csv(f"{a.outdir}/batch_evidence.csv", index=False)
    ev[["tpm_column", "gsm", "srr", "genotype", "orientation", "bases"]].to_csv(f"{a.outdir}/column_to_run_map.csv", index=False)
    if missing:
        print("UNMATCHED TPM columns (link manually):", missing, file=sys.stderr)

    # does any recorded field co-vary with the `_2` suffix?
    ev["suffix_2"] = ev["tpm_column"].fillna("").str.endswith("_2")
    for c in ["instrument", "ena_first_created", "library_name", "extract_protocol", "data_processing", "flowcells", "adapter_frac", "mean_len"]:
        if c in ev:
            print(f"\n{c}\n", ev.groupby("suffix_2")[c].agg(lambda s: sorted(set(map(str, s)))[:4]).to_string())


if __name__ == "__main__":
    main()
