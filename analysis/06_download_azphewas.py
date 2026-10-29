"""Download gene-level results per phenotype from the AstraZeneca PheWAS portal (UK Biobank 500k WGS v2 public release,
484,111 genomes).

API: https://azphewas.com/api/glr/v1/{version}/phenotype/{id}/associations?collapsing_model_ids=..&max_pvalue=1
Note: the portal stores only associations with p < 0.1, so it does not give a complete distribution. It is used for effect
sizes per carrier (beta, SE, yes_qv = number of carriers) and for the recessive model (rec = biallelic carriers).
All ancestries are stored. Collapsing models: ptv, flexdmg, rec, syn.

Output: config.AZPHEWAS_DIR / az_<trait>_<model>.tsv.gz (existing files are skipped).
Usage:  python 06_download_azphewas.py
"""
import json
import sys
import time
import urllib.request

import pandas as pd

import config

VERSION = "6319c068-fd59-46d8-85ee-82d82482eb14"
OUT = config.AZPHEWAS_DIR
OUT.mkdir(parents=True, exist_ok=True)
PHENOTYPES = {"height": "f92b34c1-d47c-44ed-ac75-dfaf636f320e", "sitheight": "422bc660-0895-423c-a6eb-6fd1c1613966",
              "ebmd": "c4174922-c168-4671-bd2d-e56a1b755f97", "alp": "fb454c13-d86a-4c2f-904e-6bcb8e0fd368",
              "oa_hip": "caa6489d-0d3b-4403-a1fd-5f9dec1197f5", "oa_knee": "2dc3f50b-a07a-4108-a7c4-ef67db3dabbd",
              "bmi": "6214eea9-05c9-4041-a1aa-6051e28f78e1", "ldl": "e9507546-8649-4c09-8387-5b732deb5718"}
MODELS = {"ptv": "70dd6e49-eed6-41b2-8911-4f85c8bfeaf9", "flexdmg": "84cc217e-4775-4b7d-941c-e226a7641ee5",
          "rec": "2fff89c7-db05-4b0a-95a1-bc9648e12050", "syn": "1eb77050-d98c-4cdf-a1a7-9dd7b6cf53f8"}

for trait, pheno_id in PHENOTYPES.items():
    for model, model_id in MODELS.items():
        f = OUT / f"az_{trait}_{model}.tsv.gz"
        if f.exists():
            continue
        url = (f"https://azphewas.com/api/glr/v1/{VERSION}/phenotype/{pheno_id}/associations"
               f"?collapsing_model_ids={model_id}&max_pvalue=1")
        for attempt in range(4):
            try:
                with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "research-script"}), timeout=300) as r:
                    d = json.loads(r.read())
                break
            except Exception as e:  # noqa: BLE001
                print(" retry", attempt, e, file=sys.stderr)
                time.sleep(15 * (attempt + 1))
        recs = d["data"]["associations"].get(model_id, [])
        df = pd.json_normalize(recs)
        df.columns = [c.replace("stats.", "") for c in df.columns]
        df.to_csv(f, sep="\t", index=False)
        print(trait, model, len(df), d["data"].get("traits_type"), flush=True)
        time.sleep(3)
print("done")
