"""Continuous version of the mouse organ diagonal of script 38 (exploratory, not pre-specified).

Instead of the IMPC significance call (p < 1e-4), each gene x MP top-level system cell holds -log10 of the smallest
genotype p-value over the parameters of the system, corrected for the number of parameters k (Sidak: 1 - (1 - p)^k).
Everything else is identical to script 38 (same organ -> MP mapping, double centring, permutation of the organ labels
between genes). For heterozygous mice, the test is repeated without the nervous-system cells.

Inputs: as script 38. Output: log only.
Call: python 38b_mouse_organ_diagonal_continuous.py
"""
import importlib

import numpy as np
import pandas as pd

import impc_mouse

ORGAN_TO_MP = importlib.import_module("38_mouse_organ_diagonal").ORGAN_TO_MP
rng = np.random.default_rng(381)
B = 10000

cols = ["marker_symbol", "zygosity", "top_level_mp_term_name", "p_value", "status"]
parts = []
for ch in pd.read_csv(impc_mouse.STATISTICAL_RESULTS, usecols=cols, chunksize=500_000, dtype=str):
    ch = ch[ch.zygosity.isin(["heterozygote", "homozygote"]) & ch.top_level_mp_term_name.notna() & (ch.status == "Successful")]
    ch = ch.assign(sys=ch.top_level_mp_term_name.str.split(","),
                   p=pd.to_numeric(ch.p_value, errors="coerce")).explode("sys").dropna(subset=["p"])
    parts.append(ch.groupby(["marker_symbol", "zygosity", "sys"]).p.agg(["min", "size"]).reset_index())
R = pd.concat(parts).groupby(["marker_symbol", "zygosity", "sys"]).agg(pmin=("min", "min"), k=("size", "sum")).reset_index()
R["v"] = -np.log10(np.clip(1 - (1 - R.pmin.clip(1e-300, 1)) ** R.k, 1e-300, 1))   # Sidak-corrected minimum p
H = impc_mouse.one_to_one_orthologs()
members = impc_mouse.organ_membership()
for zygosity in ["heterozygote", "homozygote"]:
    x = R[R.zygosity == zygosity].merge(H[["mouse", "gene_id"]], left_on="marker_symbol", right_on="mouse")
    x = x[x.gene_id.isin(members.index)]
    Z = x.pivot_table(index="gene_id", columns="sys", values="v", aggfunc="max")
    Z = Z[[c for c in Z.columns if c != "mammalian phenotype"]]
    Zc = Z - Z.mean(axis=0)
    Zc = Zc.sub(Zc.mean(axis=1), axis=0)
    genes, systems, V = list(Zc.index), list(Zc.columns), Zc.values
    labels = [members[g] for g in genes]

    def mask(labs):
        return impc_mouse.organ_cell_mask(labs, systems, V, ORGAN_TO_MP)
    m = mask(labels)
    obs = V[m].mean()
    perm = np.array([V[mask([labels[i] for i in rng.permutation(len(labels))])].mean() for _ in range(B)])
    p = (1 + (perm >= obs).sum()) / (B + 1)
    print(f"{zygosity:12s} genes {len(genes)}, cells {m.sum()}: centred {obs:+.4f} (perm {perm.mean():+.4f} +- {perm.std():.4f}; "
          f"dz {(obs - perm.mean()) / perm.std():+.1f}) p {p:.4f}")
    per = []
    for o in ORGAN_TO_MP:
        mo = mask([{o} & lab for lab in labels])
        if mo.sum() >= 10:
            per.append(f"{o} {V[mo].mean():+.2f}/{mo.sum()}")
    print("   per organ (centred / n):", ", ".join(per))
    if zygosity == "heterozygote":
        mo = mask([lab - {"nervous"} for lab in labels])
        obs2 = V[mo].mean()
        perm2 = np.array([V[mask([(labels[i] - {"nervous"}) for i in rng.permutation(len(labels))])].mean() for _ in range(2000)])
        print(f"   without nervous-system cells: {obs2:+.4f} (perm {perm2.mean():+.4f}), p {(1 + (perm2 >= obs2).sum()) / 2001:.4f}")
