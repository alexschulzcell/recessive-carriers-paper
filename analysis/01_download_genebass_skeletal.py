"""Download Genebass gene-burden results for skeletal phenotypes and the Genebass gene QC.

Genebass (Karczewski et al., Cell Genomics 2022) reports SAIGE-GENE burden tests for 394,841 UK Biobank exomes. This
script uses the public API of the Genebass browser (https://main.genebass.org/api). For each phenotype x burden set
(pLoF, missense|LC, synonymous) it stores all genes with BETA_Burden, Pvalue (SKAT-O), Pvalue_Burden and Pvalue_SKAT.
For each burden set it also stores the gene QC (CAF = cumulative allele frequency, carriers ~ 2 * CAF * N; keep flags).
Synonymous variants serve as a variant-level negative control, BMI and LDL as phenotype controls.

Output: config.GENEBASS_DIR / genebass_<trait>_<set>.tsv.gz and genebass_qc_<set>.tsv.gz (existing files are skipped).
Usage: python 01_download_genebass_skeletal.py
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

PHENOTYPES = {
    "height": "continuous-50-both_sexes--irnt",
    "sitheight": "continuous-20015-both_sexes--irnt",
    "ebmd": "continuous-bone_mineral_density_custom-both_sexes--custom",
    "oa_hip": "icd_first_occurrence-131870-both_sexes--",
    "oa_knee": "icd_first_occurrence-131872-both_sexes--",
    "oa_poly": "icd_first_occurrence-131868-both_sexes--",
    "osteopor": "icd_first_occurrence-131964-both_sexes--",
    "fract5y": "categorical-2463-both_sexes--",
    "alp": "continuous-30610-both_sexes--irnt",
    "calcium": "continuous-30680-both_sexes--irnt",
    "bmi": "continuous-21001-both_sexes--irnt",
    "ldl": "continuous-30780-both_sexes--irnt",
}
SETS = {"plof": "pLoF", "mislc": "missense|LC", "syn": "synonymous"}


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
            print("  retry", attempt, e, file=sys.stderr)
            time.sleep(10 * (attempt + 1))
    raise RuntimeError(url)


for set_key, set_name in SETS.items():
    f = OUT / f"genebass_qc_{set_key}.tsv.gz"
    if not f.exists():
        d = pd.DataFrame(get(f"{API}/gene-qc/burden-set/{urllib.parse.quote(set_name, safe='')}"))
        d.to_csv(f, sep="\t", index=False)
        print("qc", set_key, len(d))
        time.sleep(3)

for trait, analysis_id in PHENOTYPES.items():
    for set_key, set_name in SETS.items():
        f = OUT / f"genebass_{trait}_{set_key}.tsv.gz"
        if f.exists():
            continue
        d = pd.DataFrame(get(f"{API}/analysis/{analysis_id}/gene-manhattan?burdenSet={urllib.parse.quote(set_name, safe='')}"))
        d.to_csv(f, sep="\t", index=False)
        print(trait, set_key, len(d), flush=True)
        time.sleep(3)
print("done")
