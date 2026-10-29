"""Direction x organ specificity in Genes & Health 44k (pre-specified replication, PRAEREG_gnomad.md, addendum).

For the cells of the carrier atlas (script 36, without the cells excluded for clonal haematopoiesis) with |z| > 2 in
Genes & Health, the fraction pointing in the expected disease direction is compared between the top and the bottom
tissue-specificity tertile (T3 - T1). Null distribution: permutation of the tertile labels within organs (B = 5000),
one-sided p.

Inputs: 36_carrier_atlas.tsv, Genes & Health REGENIE gene-burden files (config.GH_BURDEN_DIR). Output: log only.
Call: python 37_direction_specificity_genes_and_health.py
"""
import re

import numpy as np
import pandas as pd

import common
import config

rng = np.random.default_rng(37)
B = 5000
GH_TRAITS = {"Height": "height", "Hb": "hb", "MCV": "mcv", "RBC": "rbc", "Platelets": "plt", "WBC": "wbc",
             "Lymphocytes": "lymph", "creatinine": "creat", "Urate": "urate", "Urea": "urea", "Phosphate": "phosphate",
             "ALT": "alt", "AST": "ast", "GGT": "ggt", "Bilirubin": "bili", "Albumin": "albumin", "LDL-C": "ldl",
             "HDL-C": "hdl", "Triglycerides": "tg", "Total_cholesterol": "chol", "Glucose_non_fasting": "glucose",
             "HbA1c": "hba1c", "FEV1": "fev1", "Systolic_BP": "sbp"}
MASK = {"C_plof_mis": "MASK_C.0.01", "A_plof": "MASK_A.0.01", "D_syn": "MASK_D.0.01"}
rows = []
for f in sorted(config.GH_BURDEN_DIR.glob("*_genetests_*.regenie.gz")):
    mt = re.search(r"2024_05_08_(.+?)\.residual_GNH", f.name)
    if not mt or mt.group(1) not in GH_TRAITS:
        continue
    d = pd.read_csv(f, sep=" ", comment="#", usecols=["ID", "A1FREQ", "N", "TEST", "BETA", "SE"])
    d = d[d.TEST == "ADD"]
    d["gene_id"] = d.ID.str.extract(r"\((ENSG\d+)\)")[0]
    d["mask"] = d.ID.str.extract(r"\)\.(MASK_.+)$")[0]
    d = d[(d.A1FREQ * 2 * d.N >= 5) & d.SE.gt(0)].dropna(subset=["gene_id"])
    for m, mk in MASK.items():
        x = d[d["mask"] == mk].drop_duplicates("gene_id")
        rows.append(pd.DataFrame({"gene_id": x.gene_id, "trait": GH_TRAITS[mt.group(1)], "mask": m, "z": x.BETA / x.SE}))
Z = pd.concat(rows)
X = pd.read_csv(config.results("36_carrier_atlas.tsv"), sep="\t")
X = (X[(X.status != "excluded_clonal_haematopoiesis") & X.specificity_tertile.notna()]
     [["gene_id", "symbol", "organ", "trait", "expected_direction", "specificity_tertile"]])
for m in MASK:
    k = X.merge(Z[Z["mask"] == m], on=["gene_id", "trait"])
    k = k[np.abs(k.z) > 2].copy()
    k["ok"] = (k.expected_direction * k.z) > 0
    t = k.specificity_tertile.values
    o = k.organ.values
    ok = k.ok.values
    obs = ok[t == 3].mean() - ok[t == 1].mean()
    # organs in the fixed order of the original analysis (keeps the permutation stream reproducible)
    idx = [np.where(o == g)[0] for g in common.sort_organs(np.unique(o))]
    dd = np.empty(B)
    for i in range(B):
        tp = t.copy()
        for ix in idx:
            tp[ix] = rng.permutation(tp[ix])
        dd[i] = ok[tp == 3].mean() - ok[tp == 1].mean()
    p = (1 + (dd >= obs).sum()) / (B + 1)
    per_tertile = " ".join(f"T{int(x)} {ok[t == x].mean():.3f}/{(t == x).sum()}" for x in (1, 2, 3))
    print(f"{m:11s} |z|>2 n {len(k):4d}: {per_tertile} | T3-T1 {obs:+.3f} (perm {dd.mean():+.3f}), one-sided p {p:.4f}")
