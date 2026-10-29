"""Download BioBank Japan summary statistics and keep only the variants needed for the replication (PRAEREG_bbj.md).

For each trait the archive is downloaded from the NBDC Human Database (accession hum0197, unrestricted access), its
SHA-256 is logged, and only variants with MAF < 1 % are kept that are pLoF (gnomAD v2 LOFTEE HC, GRCh37) or lie in a
coding exon +-2 bp (refGene hg19, as in script 25); the archive is then deleted. Chromosome 23 is renamed X.
Small description files (README) of each archive are kept in <output folder>/readme/<trait>/.
Please cite Sakaue et al. 2021 and the accession hum0197.

Outputs: config.REPLICATION_DIR/bbj/<trait>.tsv.gz (columns ch, pos, ref, alt, beta, sebeta, af_alt); the log line per
         trait gives URL, SHA-256, number of variants read and kept, and the columns used.
Usage:   python 49b_download_biobank_japan.py [trait ...]     (default: all 23 traits; existing files are skipped; needs curl)
"""
import gzip
import hashlib
import io
import subprocess
import sys
import time
import urllib.request
import zipfile
from concurrent.futures import ProcessPoolExecutor

from importlib import import_module

import numpy as np
import pandas as pd

import config

finngen = import_module("49a_download_finngen")   # same variant mask and filter

OUT = config.REPLICATION_DIR / "bbj"
TMP = config.REPLICATION_DIR / "bbj_download_tmp"
BASE = "https://humandbs.dbcls.jp/files/hum0197/"
SOURCES = {k: BASE + f"hum0197.v3.BBJ.{v}.v1.zip" for k, v in {
    "height": "Hei", "hb": "Hb", "mcv": "MCV", "rbc": "RBC", "plt": "PLT", "wbc": "WBC", "lymph": "LYM", "creat": "sCr", "urate": "UA",
    "urea": "BUN", "alt": "ALT", "ast": "AST", "ggt": "GGT", "bili": "TBil", "albumin": "Alb", "ldl": "LDLC", "hdl": "HDLC", "tg": "TG",
    "chol": "TC", "glucose": "Glucose", "hba1c": "HbA1c", "sbp": "SBP"}.items()}
SOURCES["phosphate"] = BASE + "hum0197.v5.gwas.P.v1.txt.gz"
COLUMNS = {"ch": ["CHR", "chr", "CHROM", "#CHROM", "chromosome"], "pos": ["POS", "pos", "BP", "position", "base_pair_location"],
           "ref": ["REF", "ref", "ALLELE0", "Allele1", "other_allele", "A2"],
           "alt": ["ALT", "alt", "ALLELE1", "Allele2", "effect_allele", "A1"],  # BOLT-LMM: ALLELE1 = effect allele (A1FREQ, BETA)
           "af_alt": ["Frq", "AF_Allele2", "ALT_Frq", "AF", "af_alt", "effect_allele_frequency", "A1FREQ", "FRQ"],
           "beta": ["BETA", "beta", "Beta"], "sebeta": ["SE", "se", "SEBETA", "sebeta", "standard_error"]}


def filter_file(fobj, LOF, CDS, k):
    header = fobj.readline().decode().rstrip("\r\n")
    sep = "\t" if "\t" in header else r"\s+"
    header = header.split("\t") if sep == "\t" else header.split()
    use = {}
    for col, candidates in COLUMNS.items():
        found = [c for c in candidates if c in header]
        if not found:
            raise SystemExit(f"{k}: column {col} missing; header {header}")
        use[col] = found[0]
    parts, n = [], 0
    for chunk in pd.read_csv(fobj, sep=sep, header=None, names=header, usecols=list(use.values()), chunksize=2_000_000,
                             dtype={use["ch"]: str}):
        chunk = chunk.rename(columns={v: c for c, v in use.items()})
        n += len(chunk)
        chunk["ch"] = chunk.ch.str.replace("chr", "", regex=False).replace({"23": "X"})
        af = pd.to_numeric(chunk.af_alt, errors="coerce")
        chunk = chunk[np.minimum(af, 1 - af) < 0.01]
        parts.append(chunk[finngen.keep_rows(chunk, LOF, CDS)])
    return pd.concat(parts), n, use


def download(k):
    target = OUT / f"{k}.tsv.gz"
    if target.exists():
        return f"{k}: present, skipped"
    LOF, CDS = finngen.variant_mask("hg19")
    url = SOURCES[k]
    archive = TMP / url.rsplit("/", 1)[1]
    for attempt in range(4):
        try:
            t0 = time.time()
            h = hashlib.sha256()
            # curl with resume (the server drops long connections), then check the size against Content-Length
            expected = int(urllib.request.urlopen(urllib.request.Request(url, method="HEAD", headers={"User-Agent": "research-script"}),
                                                  timeout=120).headers["Content-Length"])
            subprocess.run(["curl", "-sSL", "--retry", "20", "--retry-delay", "15", "--retry-all-errors", "--speed-limit", "10000",
                            "--speed-time", "120", "-C", "-", "-o", str(archive), url], check=False)
            if archive.stat().st_size != expected:
                raise IOError(f"incomplete: {archive.stat().st_size} of {expected} bytes")
            with open(archive, "rb") as f:
                while b := f.read(1 << 22):
                    h.update(b)
            parts, n, info = [], 0, []
            if archive.suffix == ".zip":
                with zipfile.ZipFile(archive) as z:
                    data_files = [i.filename for i in z.infolist() if i.file_size > 10_000_000]   # README etc. are small
                    info.append(f"archive: {z.namelist()}")
                    for i in z.infolist():   # keep description files (README)
                        if 0 < i.file_size <= 10_000_000:
                            z.extract(i, OUT / "readme" / k)
                    for m in data_files:
                        with z.open(m) as raw:
                            fo = gzip.GzipFile(fileobj=raw) if m.endswith(".gz") else io.BufferedReader(raw)
                            d, nn, use = filter_file(fo, LOF, CDS, k)
                            parts.append(d)
                            n += nn
                            info.append(f"{m}: {nn} variants, columns {use}")
            else:
                with gzip.open(archive, "rb") as fo:
                    d, nn, use = filter_file(fo, LOF, CDS, k)
                    parts.append(d)
                    n += nn
                    info.append(f"columns {use}")
            d = pd.concat(parts)
            tmp = target.with_suffix(".tmp.gz")
            d[list(COLUMNS)].to_csv(tmp, sep="\t", index=False)
            tmp.replace(target)
            archive.unlink()
            return f"{k}\t{url}\tsha256 {h.hexdigest()}\tvariants {n}\tkept {len(d)}\t{' | '.join(info)}\t{time.time() - t0:.0f} s"
        except SystemExit:
            raise
        except Exception as e:  # noqa: BLE001
            print(f"  {k}: attempt {attempt + 1} failed: {e}", flush=True)
            time.sleep(30 * (attempt + 1))
    return f"{k}: FAILED"


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    TMP.mkdir(parents=True, exist_ok=True)
    selection = sys.argv[1:] or list(SOURCES)
    print("Start", time.strftime("%Y-%m-%d %H:%M:%S"), "traits:", selection, flush=True)
    with ProcessPoolExecutor(2) as ex:
        for line in ex.map(download, selection):
            print(line, flush=True)
    print("done", time.strftime("%Y-%m-%d %H:%M:%S"))
