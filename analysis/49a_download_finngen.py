"""Download FinnGen summary statistics and keep only the variants needed for the replication (PRAEREG_finngen.md).

For each trait the public file is streamed, the SHA-256 of the original (gzip) file is logged, and only variants with
MAF < 1 % are kept that are pLoF (gnomAD v2 LOFTEE HC, GRCh38 position) or lie in a coding exon +-2 bp (refGene hg38,
as in script 25). Genes are not assigned here (script 49 does that). Chromosome 23 is renamed X.
Sources: FinnGen release R13 (height) and R12 (laboratory values), public Google Cloud buckets. The raw files are not
redistributed; please acknowledge "the participants and investigators of the FinnGen study" and cite Kurki et al. 2023.

Outputs: config.REPLICATION_DIR/finngen/<trait>.tsv.gz (columns ch, pos, ref, alt, beta, sebeta, af_alt); the log line
         per trait gives URL, SHA-256, size, number of variants read and kept.
Usage:   python 49a_download_finngen.py [trait ...]     (default: all 23 traits; existing files are skipped)
"""
import gzip
import hashlib
import sys
import time
import urllib.request
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

import config
from replication_tests import VariantAnnotation

OUT = config.REPLICATION_DIR / "finngen"
R12 = "https://storage.googleapis.com/finngen-public-data-r12/lab_values/summary_stats/finngen_R12_{}.gz"
SOURCES = {"height": "https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_HEIGHT_IRN.gz"}
SOURCES.update({k: R12.format(v) for k, v in {
    "hb": "3000963", "mcv": "3024731", "rbc": "3026361", "plt": "3007461", "wbc": "3010813", "lymph": "3019198", "creat": "3020564",
    "cysc": "3030366", "urate": "3026493", "urea": "3020779", "phosphate": "3003458", "alt": "3006923", "ast": "3013721", "ggt": "3026910",
    "bili": "3006140", "albumin": "3024561", "ldl": "3001308", "hdl": "3023602", "tg": "3048773", "chol": "3019900", "glucose": "3018251",
    "hba1c": "3004410"}.items()})
COLUMNS = {"ch": ["#chrom", "chrom", "CHR"], "pos": ["pos", "POS"], "ref": ["ref", "REF"], "alt": ["alt", "ALT"],
           "beta": ["beta", "BETA"], "sebeta": ["sebeta", "SE", "se"], "af_alt": ["af_alt", "A1FREQ", "af"]}


def variant_mask(build="hg38"):
    """pLoF positions and merged coding intervals (same annotation as script 25)."""
    ann = VariantAnnotation(gene_symbols=False)
    lof_pos = {}
    for ch, pos, _, _ in ann.lof_dict(build):
        lof_pos.setdefault(ch, set()).add(pos)
    cds = {}
    for ch, (v, _) in ann.cds[build].items():
        iv = sorted((s, e) for s, e, _ in v)
        st, en = [], []
        for s, e in iv:
            if st and s <= en[-1]:
                en[-1] = max(en[-1], e)
            else:
                st.append(s)
                en.append(e)
        cds[ch] = (np.array(st), np.array(en))
    return {ch: np.array(sorted(p)) for ch, p in lof_pos.items()}, cds


def keep_rows(chunk, LOF, CDS):
    """Boolean mask of variants in a coding interval or at a pLoF position."""
    keep = np.zeros(len(chunk), bool)
    pos = chunk.pos.to_numpy(np.int64)
    for c, ix in chunk.groupby("ch").indices.items():
        if c in CDS:
            st, en = CDS[c]
            j = np.searchsorted(st, pos[ix], side="left") - 1
            keep[ix] |= (j >= 0) & (pos[ix] <= en[np.clip(j, 0, None)])
        if c in LOF:
            keep[ix] |= np.isin(pos[ix], LOF[c])
    return keep


class HashingReader:
    """File-like wrapper that computes the SHA-256 and size of everything read."""

    def __init__(self, f):
        self.f, self.h, self.n = f, hashlib.sha256(), 0

    def read(self, n=-1):
        b = self.f.read(n)
        self.h.update(b)
        self.n += len(b)
        return b


def download(k):
    target = OUT / f"{k}.tsv.gz"
    if target.exists():
        return f"{k}: present, skipped"
    LOF, CDS = variant_mask()
    for attempt in range(4):
        try:
            t0 = time.time()
            raw = HashingReader(urllib.request.urlopen(urllib.request.Request(SOURCES[k], headers={"User-Agent": "research-script"}), timeout=600))
            gz = gzip.GzipFile(fileobj=raw)
            header = gz.readline().decode().rstrip("\n").split("\t")
            use = {}
            for col, candidates in COLUMNS.items():
                found = [c for c in candidates if c in header]
                if not found:
                    raise SystemExit(f"{k}: column {col} missing; header {header}")
                use[col] = found[0]
            parts, n_total = [], 0
            for chunk in pd.read_csv(gz, sep="\t", header=None, names=header, usecols=list(use.values()), chunksize=2_000_000,
                                     dtype={use["ch"]: str}):
                chunk = chunk.rename(columns={v: c for c, v in use.items()})
                n_total += len(chunk)
                chunk["ch"] = chunk.ch.str.replace("chr", "", regex=False).replace({"23": "X"})
                af = pd.to_numeric(chunk.af_alt, errors="coerce")
                chunk = chunk[np.minimum(af, 1 - af) < 0.01]
                parts.append(chunk[keep_rows(chunk, LOF, CDS)])
            while raw.read(1 << 20):   # read the rest for the complete checksum
                pass
            d = pd.concat(parts)
            tmp = target.with_suffix(".tmp.gz")
            d[list(COLUMNS)].to_csv(tmp, sep="\t", index=False)
            tmp.replace(target)
            return (f"{k}\t{SOURCES[k]}\tsha256 {raw.h.hexdigest()}\tbytes {raw.n}\tvariants {n_total}\tkept {len(d)}\t"
                    f"columns {use}\t{time.time() - t0:.0f} s")
        except SystemExit:
            raise
        except Exception as e:  # noqa: BLE001
            print(f"  {k}: attempt {attempt + 1} failed: {e}", flush=True)
            time.sleep(30 * (attempt + 1))
    return f"{k}: FAILED"


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    selection = sys.argv[1:] or list(SOURCES)
    print("Start", time.strftime("%Y-%m-%d %H:%M:%S"), "traits:", selection, flush=True)
    with ProcessPoolExecutor(4) as ex:
        for line in ex.map(download, selection):
            print(line, flush=True)
    print("done", time.strftime("%Y-%m-%d %H:%M:%S"))
