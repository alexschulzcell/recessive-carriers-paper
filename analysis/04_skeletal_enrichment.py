"""Are recessive skeletal dysplasia genes enriched for carrier effects on skeletal phenotypes in UK Biobank?

Per gene: z = sign(BETA_Burden) * Phi^-1(1 - P_Burden / 2) from Genebass (UK Biobank, 394,841 exomes). Only genes that pass
the Genebass keep flags and have at least MIN carriers in the burden set are used.
Statistics per gene set: mean(z) (direction), mean(z^2) (inflation, independent of direction), fraction p < 0.01.
Comparison: matched draws (for every test gene one control gene from the same stratum, LOEUF decile x CDS-length quintile),
B = 2000 repetitions -> empirical p and z of the deviation (dz).
Controls: (a) pool = no Mendelian gene in HPO/OMIM and not in the nosology; (b) AR_nonskel = recessive Mendelian genes
without skeletal or growth HPO terms; AD_nos is tested for direct comparison. Negative controls: burden set synonymous;
phenotypes BMI and LDL.
Pre-specified reading (no stopping rule): a finding counts if AR_nos is raised against BOTH controls (a, b) for pLoF or
missense|LC, and not for synonymous.

Inputs:  gene_annotation.tsv.gz (script 03); config.GENEBASS_DIR / genebass_<trait>_<set>.tsv.gz (script 01)
Outputs: 04_skeletal_enrichment_min<MIN>.tsv (gene-set level), 04_gene_x_phenotype_min<MIN>.tsv (gene level for the
         nosology genes), both in RC_RESULTS_DIR.
Usage:   python 04_skeletal_enrichment.py [MIN]   (MIN = minimum number of carriers, default 10; the pipeline uses 10)
"""
import sys

import numpy as np
import pandas as pd
from scipy.stats import norm

import config

MIN_CARRIERS = int(sys.argv[1]) if len(sys.argv) > 1 else 10
B = 2000
rng = np.random.default_rng(1)
PHENOTYPES = ["height", "sitheight", "ebmd", "oa_hip", "oa_knee", "oa_poly", "osteopor", "fract5y", "alp", "calcium", "bmi", "ldl"]
SETS = ["plof", "mislc", "syn"]
TESTS = ["AR_nos", "AD_nos", "AR_skelsyndr", "AR_nonskel"]
CONTROLS = ["pool", "AR_nonskel"]

A = pd.read_csv(config.results("gene_annotation.tsv.gz"), sep="\t")
A = A[A.protein_coding & A.oe_lof_upper.notna() & A.cds_length.notna()].copy()
A["lb"] = pd.qcut(A.oe_lof_upper, 10, labels=False, duplicates="drop")
A["cb"] = pd.qcut(A.cds_length, 5, labels=False, duplicates="drop")
A["stratum"] = A.lb.astype(int) * 10 + A.cb.astype(int)

res, genes = [], []
for s in SETS:
    for ph in PHENOTYPES:
        f = config.GENEBASS_DIR / f"genebass_{ph}_{s}.tsv.gz"
        if not f.exists():
            continue
        g = pd.read_csv(f, sep="\t")[["gene_id", "BETA_Burden", "Pvalue_Burden", "Pvalue"]]
        g = g.merge(A, on="gene_id")
        g = g[(g[f"keep_{s}"] == True) & (g[f"carriers_{s}"] >= MIN_CARRIERS)  # noqa: E712
              & g.Pvalue_Burden.between(0, 1, inclusive="right")]
        g["z"] = np.sign(g.BETA_Burden) * norm.isf(g.Pvalue_Burden.clip(1e-300) / 2)
        g["z2"] = g.z ** 2
        g["sig"] = (g.Pvalue_Burden < 0.01).astype(float)
        nos = g[g.set.isin(["AR_nos", "AD_nos", "nos_other"])]
        genes.append(nos.assign(pheno=ph, bset=s)[["gene_id", "symbol", "set", "nos_groups", "pheno", "bset", f"carriers_{s}",
                                                   "BETA_Burden", "Pvalue_Burden", "Pvalue", "z"]]
                     .rename(columns={f"carriers_{s}": "carriers"}))
        for t in TESTS:
            T = g[g.set == t]
            if len(T) < 5:
                continue
            obs = T[["z", "z2", "sig"]].mean()
            need = T.stratum.value_counts()
            for c in CONTROLS:
                if c == t:
                    continue
                C = g[g.set == c]
                by_stratum = {k: v[["z", "z2", "sig"]].to_numpy() for k, v in C.groupby("stratum")}
                draws = np.zeros((B, 3))
                used = 0
                for k, n in need.items():
                    if k not in by_stratum:  # stratum missing among the controls -> nearest stratum
                        k = min(by_stratum, key=lambda x: abs(x - k))
                    arr = by_stratum[k]
                    idx = rng.integers(0, len(arr), size=(B, n))
                    draws += arr[idx].sum(axis=1)
                    used += n
                draws /= used
                for j, st in enumerate(["z", "z2", "sig"]):
                    mu, sd = draws[:, j].mean(), draws[:, j].std()
                    p_hi = (1 + (draws[:, j] >= obs.iloc[j]).sum()) / (B + 1)
                    p_lo = (1 + (draws[:, j] <= obs.iloc[j]).sum()) / (B + 1)
                    res.append(dict(bset=s, pheno=ph, test=t, ctrl=c, n=len(T), stat=st, obs=obs.iloc[j], ctrl_mean=mu,
                                    dz=(obs.iloc[j] - mu) / sd if sd > 0 else np.nan, p_hi=p_hi, p_lo=p_lo,
                                    median_carriers=T[f"carriers_{s}"].median()))
R = pd.DataFrame(res)
R.to_csv(config.results(f"04_skeletal_enrichment_min{MIN_CARRIERS}.tsv"), sep="\t", index=False)
pd.concat(genes).to_csv(config.results(f"04_gene_x_phenotype_min{MIN_CARRIERS}.tsv"), sep="\t", index=False)
pd.set_option("display.width", 250)
pd.set_option("display.max_rows", 500)
for st in ["z2", "z"]:
    P = R[(R.test == "AR_nos") & (R.stat == st)].pivot_table(index=["pheno"], columns=["bset", "ctrl"], values="dz").round(1)
    print(f"\nAR_nos against control, deviation in SD of the draws, statistic {st}:\n", P.to_string())
P = R[(R.stat == "z2") & (R.ctrl == "pool")].pivot_table(index="pheno", columns=["bset", "test"], values="dz").round(1)
print("\nAll test sets against pool, mean z^2:\n", P.to_string())
