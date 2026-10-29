"""pLoF replication with the official gene-burden results of Genes & Health 44k.

Genes & Health public gene-burden tables (gs://genesandhealth_publicdatasets, CC BY-SA 4.0; REGENIE; masks
A = LOFTEE high-confidence pLoF, C = pLoF + missense, D = synonymous; per gene AAF < 1 %, test ADD (burden)).
About 33k British South Asians without overlap with UK Biobank; 24 traits (height, blood count, kidney, liver, lipids,
glucose, FEV1, systolic blood pressure). Allele count per gene >= 5 (AAF * 2N). Burden z^2 = (BETA/SE)^2, capped at 50.
Tests as in script 25: T1 height (AR_nos versus matched controls), T2 organ diagonal (permutation of the organ labels),
T3 main-product rule W1/W2b (classes from script 22, without skeletal, permutation within organ; pre-specified in
PRAEREG_hauptprodukt.md).

Inputs: config.GH_BURDEN_DIR (*_genetests_*.regenie.gz), gene_annotation.tsv.gz, organ_fractions_recessive.tsv.gz,
22_main_product_units_plof.tsv (script 22).
Outputs: cache("replication_gene_stats")/GHburden_{plof,mis,syn}.tsv.gz (genes x traits z^2),
26_genes_and_health_burden_replication.tsv.
Call: python 26_genes_and_health_burden_replication.py
"""
import re
from collections import defaultdict

import numpy as np
import pandas as pd

import common
import config

OUT = config.cache("replication_gene_stats")
rng = np.random.default_rng(29)
# Genes & Health phenotype name -> trait code
TRAIT_CODES = {"Height": "height", "Hb": "hb", "MCV": "mcv", "RBC": "rbc", "Platelets": "plt", "WBC": "wbc",
               "Lymphocytes": "lymph", "creatinine": "creat", "Urate": "urate", "Urea": "urea", "Phosphate": "phosphate",
               "ALT": "alt", "AST": "ast", "GGT": "ggt", "Bilirubin": "bili", "Albumin": "albumin", "LDL-C": "ldl",
               "HDL-C": "hdl", "Triglycerides": "tg", "Total_cholesterol": "chol", "Glucose_non_fasting": "glucose",
               "HbA1c": "hba1c", "FEV1": "fev1", "Systolic_BP": "sbp"}
TRAIT_TO_ORGAN = {"height": "skeletal", "hb": "blood", "mcv": "blood", "rbc": "blood", "plt": "blood", "wbc": "immune",
                  "lymph": "immune", "creat": "kidney", "urate": "kidney", "urea": "kidney", "phosphate": "kidney",
                  "alt": "liver", "ast": "liver", "ggt": "liver", "bili": "liver", "albumin": "liver", "ldl": "lipids",
                  "hdl": "lipids", "tg": "lipids", "chol": "lipids", "glucose": "glucose", "hba1c": "glucose",
                  "fev1": "lung", "sbp": "circulatory"}
MASKS = {"plof": "MASK_A.0.01", "mis": "MASK_C.0.01", "syn": "MASK_D.0.01"}

Z = {m: {} for m in MASKS}
for f in sorted(config.GH_BURDEN_DIR.glob("*_genetests_*.regenie.gz")):
    trait = re.search(r"2024_05_08_(.+?)\.residual_GNH", f.name).group(1)
    if trait not in TRAIT_CODES:
        continue
    d = pd.read_csv(f, sep=" ", comment="#", usecols=["ID", "A1FREQ", "N", "TEST", "BETA", "SE"])
    d = d[d.TEST == "ADD"]
    d["gene_id"] = d.ID.str.extract(r"\((ENSG\d+)\)")[0]
    d["mask"] = d.ID.str.extract(r"\)\.(MASK_.+)$")[0]
    d["ac"] = d.A1FREQ * 2 * d.N
    for m, mask_name in MASKS.items():
        x = d[(d["mask"] == mask_name) & (d.ac >= 5) & d.SE.gt(0)].drop_duplicates("gene_id")
        Z[m][TRAIT_CODES[trait]] = pd.Series(((x.BETA / x.SE) ** 2).clip(upper=50).values, index=x.gene_id.values)
    print(trait, {m: len(Z[m][TRAIT_CODES[trait]]) for m in MASKS}, flush=True)
for m in MASKS:
    pd.DataFrame(Z[m]).to_csv(OUT / f"GHburden_{m}.tsv.gz", sep="\t")

A = pd.read_csv(config.results("gene_annotation.tsv.gz"), sep="\t")
A2 = A[A.protein_coding & A.oe_lof_upper.notna() & A.cds_length.notna()].copy()
A2["stratum"] = pd.qcut(A2.oe_lof_upper, 10, labels=False) * 10 + pd.qcut(A2.cds_length, 5, labels=False)


def matched(T, C, col, n_boot=5000):
    """Mean of col in T versus control genes C drawn with replacement within the same LOEUF x length stratum."""
    by_stratum = {k: v[col].to_numpy() for k, v in C.groupby("stratum")}
    tot = np.zeros(n_boot)
    for k, n in T.stratum.value_counts().items():
        arr = by_stratum[k] if k in by_stratum else by_stratum[min(by_stratum, key=lambda x: abs(x - k))]
        tot += arr[rng.integers(0, len(arr), size=(n_boot, n))].sum(axis=1)
    tot /= len(T)
    o = T[col].mean()
    return o, (o - tot.mean()) / tot.std(), (1 + (tot >= o).sum()) / (n_boot + 1)


rows = []
print("\nT1 height (G&H burden)")
for m in MASKS:
    H = A2.merge(Z[m]["height"].rename("z2"), left_on="gene_id", right_index=True)
    T = H[H.set == "AR_nos"]
    for c in ["pool", "AR_nonskel"]:
        o, dz, p = matched(T, H[H.set == c], "z2")
        print(f"  {m:4s} AR_nos n {len(T):3d} vs {c:10s}: z2 {o:.2f}, dz {dz:+.1f}, p {p:.4f}")
        rows.append(dict(test="T1", set=m, comparison=c, n=len(T), dz=dz, p=p))
    o, dz, p = matched(H[H.set == "AR_nonskel"], H[H.set == "pool"], "z2")
    print(f"  {m:4s} control AR_nonskel vs pool: dz {dz:+.1f}, p {p:.3f}")

OF = pd.read_csv(config.results("organ_fractions_recessive.tsv.gz"), sep="\t")
members = OF[OF.frac >= 0.2]


def centred(m, genes):
    """Doubly centred z^2 (traits with >= 50 genes), first per trait, then per gene."""
    M = pd.DataFrame(Z[m]).reindex(genes)
    M = M.loc[:, M.notna().sum() >= 50]
    Mc = M - M.mean(axis=0)
    return Mc.sub(Mc.mean(axis=1), axis=0)


print("\nT2 organ diagonal (G&H burden, permutation of organ labels B=2000)")
for m in MASKS:
    Zc = centred(m, sorted(set(members.gene_id)))
    Zc = Zc[Zc.notna().sum(axis=1) >= 3]
    organs = common.sort_organs({TRAIT_TO_ORGAN[t] for t in Zc.columns})   # fixed organ order of the original analysis
    genes = list(Zc.index)
    Lm = np.array([[1.0 if x in set(members[members.organ == o].gene_id) else 0.0 for o in organs] for x in genes])
    Cm = np.array([[1.0 if TRAIT_TO_ORGAN[t] == o else 0.0 for t in Zc.columns] for o in organs])
    V = Zc.values
    valid = ~np.isnan(V)
    V0 = np.nan_to_num(V)

    def stat(Lx, C=Cm):
        mk = ((Lx @ C) > 0) & valid
        return (V0 * mk).sum() / max(mk.sum(), 1)
    obs = stat(Lm)
    perm = np.array([stat(Lm[rng.permutation(len(Lm))]) for _ in range(2000)])
    dz, p = (obs - perm.mean()) / perm.std(), (1 + (perm >= obs).sum()) / 2001
    per = []
    for i, o in enumerate(organs):
        Lo, Co = Lm[:, [i]], Cm[[i], :]
        ob = stat(Lo, Co)
        pe = np.array([stat(Lo[rng.permutation(len(Lo))], Co) for _ in range(1000)])
        per.append(f"{o} {(ob - pe.mean()) / pe.std():+.1f}")
    print(f"  {m:4s}: {len(genes)} genes, {Zc.shape[1]} traits -> dz {dz:+.1f}, p {p:.4f} | {', '.join(per)}")
    rows.append(dict(test="T2", set=m, comparison="organ_permutation", n=len(genes), dz=dz, p=p))

print("\nT3 rule (classes W1/W2b from script 22, without skeletal, permutation within organ B=5000)")
U = pd.read_csv(config.results("22_main_product_units_plof.tsv"), sep="\t")[["gene_id", "organ", "W1", "W2b"]]
U = U[U.organ != "skeletal"]
for m in MASKS:
    Zc = centred(m, sorted(set(U.gene_id)))
    organ_traits = defaultdict(list)
    for t in Zc.columns:
        organ_traits[TRAIT_TO_ORGAN[t]].append(t)
    vals = [Zc.loc[g, organ_traits[o]].mean() if o in organ_traits and g in Zc.index else np.nan for g, o in zip(U.gene_id, U.organ)]
    X = U.assign(value=vals).dropna(subset=["value"])
    for route in ["W1", "W2b"]:
        lab = X[route].values.astype(bool)
        w = X.value.values
        # organs in order of first appearance in the unit table (as in the original analysis)
        idx = [np.where(X.organ.values == o)[0] for o in X.organ.unique()]
        obs = w[lab].mean() - w[~lab].mean()
        dd = np.empty(5000)
        for i in range(5000):
            Lx = lab.copy()
            for ix in idx:
                Lx[ix] = rng.permutation(Lx[ix])
            dd[i] = w[Lx].mean() - w[~Lx].mean()
        p = (1 + (dd >= obs).sum()) / 5001
        print(f"  {m:4s} {route}: diff {obs:+.3f} (n yes {int(lab.sum())}/{len(X)}; organs {common.sort_organs(X.organ.unique())}), p {p:.4f}")
        rows.append(dict(test=f"T3_{route}", set=m, comparison="within_organ", n=len(X), dz=obs, p=p))
pd.DataFrame(rows).to_csv(config.results("26_genes_and_health_burden_replication.tsv"), sep="\t", index=False)
