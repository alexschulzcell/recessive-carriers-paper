"""Download the Genebass gene-burden results for the pre-specified ICD-10 first-occurrence phenotypes.

The 122 ICD-10 codes and their Genebass analysis IDs were selected in advance (PRAEREG_diagnosen.md) and are listed in
config.ICD10_SELECTION. For every code and burden set (pLoF, synonymous, missense|LC) the gene-level results are fetched
from the Genebass API and written to config.GENEBASS_ICD_DIR/<ICD-10 code>_<set>.tsv.gz (set = plof, syn, mislc).
Files that already exist are skipped. Input for script 39.
Call: python 39a_download_genebass_icd10.py
"""
import gzip
import json
import sys
import time
import urllib.parse
import urllib.request

import pandas as pd

import config

API = "https://main.genebass.org/api"
OUT = config.GENEBASS_ICD_DIR
OUT.mkdir(parents=True, exist_ok=True)
SETS = {"plof": "pLoF", "syn": "synonymous", "mislc": "missense|LC"}


def get(url):
    req = urllib.request.Request(url, headers={"Accept-Encoding": "gzip", "User-Agent": "research-script"})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(req, timeout=300) as r:
                b = r.read()
            try:
                b = gzip.decompress(b)
            except OSError:
                pass
            return json.loads(b)
        except Exception as e:  # noqa: BLE001
            print(" retry", attempt, e, file=sys.stderr, flush=True)
            time.sleep(15 * (attempt + 1))
    raise RuntimeError(url)


S = pd.read_csv(config.ICD10_SELECTION, sep="\t")
for r in S.itertuples():
    for set_code, burden_set in SETS.items():
        f = OUT / f"{r.icd}_{set_code}.tsv.gz"
        if f.exists():
            continue
        d = pd.DataFrame(get(f"{API}/analysis/{r.analysis_id}/gene-manhattan?burdenSet={urllib.parse.quote(burden_set, safe='')}"))
        d.to_csv(f, sep="\t", index=False)
        print(r.icd, set_code, len(d), flush=True)
        time.sleep(2)
print("done")
