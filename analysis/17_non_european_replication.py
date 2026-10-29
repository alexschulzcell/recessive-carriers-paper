"""Replication in independent people: non-European ancestries in the AstraZeneca PheWAS portal (UK Biobank 500k WGS;
Genebass contains Europeans only). Exploratory.

Problem: the portal stores only p < 0.1. Therefore a direction test: for genes with a Genebass height (or eBMD) effect
(discovery, Europeans), count in how many non-European ancestries the gene appears with p < 0.1 and with which sign. Under
the null both signs are equally likely and the appearance rate is ~10 %. Comparison with control genes without a discovery
signal. Models: ptv (matched to Genebass pLoF) and flexdmg (matched to missense|LC).

Inputs: shet_merged.tsv.gz (script 11); config.AZPHEWAS_DIR / az_<trait>_<model>.tsv.gz (script 06)
Output: printed only.
Usage:  python 17_non_european_replication.py
"""
import numpy as np
import pandas as pd
from scipy.stats import norm, binomtest

import config

A = pd.read_csv(config.results("shet_merged.tsv.gz"), sep="\t")
for ph, zcol in [("height", "z_height"), ("ebmd", "z_ebmd")]:
    for m, s in [("ptv", "plof"), ("flexdmg", "mislc")]:
        az = pd.read_csv(config.AZPHEWAS_DIR / f"az_{ph}_{m}.tsv.gz", sep="\t")
        ne = az[az.ancestry != "European"]
        x = A[["symbol", "set", "category", f"{zcol}_{s}"]].dropna().rename(columns={f"{zcol}_{s}": "zd"})  # zd = discovery z
        x["pd"] = 2 * norm.sf(x.zd.abs())
        y = ne.merge(x, left_on="gene_name", right_on="symbol")
        y["same"] = np.sign(y.beta) == np.sign(y.zd)
        print(f"\n== {ph} {m}: non-European cells with p<0.1: {len(ne)}")
        for label, sel in [("AR_nos, discovery p<1e-3", (y.set == "AR_nos") & (y.pd < 1e-3)),
                           ("AR_nos, discovery p<0.05", (y.set == "AR_nos") & (y.pd < 0.05)),
                           ("all genes, discovery p<1e-4", y.pd < 1e-4), ("all genes, discovery p>0.5", y.pd > 0.5)]:
            z = y[sel]
            if len(z) == 0:
                continue
            bt = binomtest(int(z.same.sum()), len(z), 0.5, alternative="greater").pvalue
            print(f"  {label:30s}: {len(z):4d} cells, same direction {z.same.mean():.2f} (binomial p {bt:.1e}); "
                  f"ancestries: {z.ancestry.value_counts().to_dict()}")
        z = y[(y.set == "AR_nos") & (y.pd < 0.05)]
        if len(z):
            print(z.sort_values("pd")[["gene_name", "ancestry", "beta", "pvalue", "yes_qv", "zd", "same"]].round(3).head(20).to_string(index=False))
