"""Download Genebass gene-burden results for organ-matched quantitative traits (generalisation across organ systems).

Same API as script 01 (https://main.genebass.org/api, UK Biobank 394,841 exomes); burden sets pLoF, synonymous (negative
control) and missense|LC; all traits are continuous UK Biobank fields, inverse-rank normalised, both sexes.

Output: config.GENEBASS_DIR / genebass_<trait>_<set>.tsv.gz (existing files are skipped).
Usage:  python 13_download_genebass_organ_traits.py
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
OUT = config.GENEBASS_DIR
OUT.mkdir(parents=True, exist_ok=True)
# trait code -> UK Biobank field id
TRAITS = {"wbc": "30000", "rbc": "30010", "hb": "30020", "mcv": "30040", "plt": "30080", "lymph": "30120", "albumin": "30600",
          "fev1": "3063", "fvc": "3062", "alt": "30620", "ast": "30650", "creat": "30700", "cysc": "30720", "ggt": "30730",
          "glucose": "30740", "hba1c": "30750", "hdl": "30760", "igf1": "30770", "phosphate": "30810", "bili": "30840",
          "tg": "30870", "urate": "30880", "sbp": "4080", "grip": "47", "fluidint": "20016", "reaction": "20023", "chol": "30690",
          "urea": "30670"}
SETS = {"plof": "pLoF", "syn": "synonymous", "mislc": "missense|LC"}


def get(url):
    """GET a JSON document (gzip-compressed or plain), with up to four attempts."""
    req = urllib.request.Request(url, headers={"Accept-Encoding": "gzip", "User-Agent": "research-script"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=180) as r:
                b = r.read()
            try:
                b = gzip.decompress(b)
            except OSError:
                pass
            return json.loads(b)
        except Exception as e:  # noqa: BLE001
            print(" retry", attempt, e, file=sys.stderr)
            time.sleep(10 * (attempt + 1))
    raise RuntimeError(url)


for trait, field in TRAITS.items():
    for set_key, set_name in SETS.items():
        f = OUT / f"genebass_{trait}_{set_key}.tsv.gz"
        if f.exists():
            continue
        d = pd.DataFrame(get(f"{API}/analysis/continuous-{field}-both_sexes--irnt/gene-manhattan?burdenSet={urllib.parse.quote(set_name, safe='')}"))
        d.to_csv(f, sep="\t", index=False)
        print(trait, set_key, len(d), flush=True)
        time.sleep(2)
print("done")
