"""Directional organ diagonal on ICD-10 diagnoses (Genebass first-occurrence phenotypes).

Pre-specified test (PRAEREG_diagnosen.md). For 122 pre-selected ICD-10 first-occurrence codes, each assigned to one organ
system, the signed Genebass burden z of every recessive gene is averaged per organ after double centring. The primary test
asks whether carriers are diagnosed more often with diseases of the organ affected by the recessive disease of their gene
(directional, positive z = higher risk), with a permutation of the organ labels between genes. Secondary: the same test
on z^2 (direction-free), leave-one-organ-out, per-organ means and the strongest matching gene-diagnosis pairs.

Inputs: config.ICD10_SELECTION, config.GENEBASS_ICD_DIR/<code>_<set>.tsv.gz (script 39a), gene_annotation.tsv.gz,
organ_fractions_recessive.tsv.gz.
Output: 39_icd10_diagnoses_diagonal.tsv.
Call: python 39_icd10_diagnoses_diagonal.py
"""
import numpy as np
import pandas as pd
from scipy.stats import norm

import common
import config

rng = np.random.default_rng(39)
B = 10000
S = pd.read_csv(config.ICD10_SELECTION, sep="\t")
A = pd.read_csv(config.results("gene_annotation.tsv.gz"), sep="\t")
OF = pd.read_csv(config.results("organ_fractions_recessive.tsv.gz"), sep="\t")
OF = OF[OF.frac >= 0.2]
members = OF.groupby("gene_id").organ.apply(set)
GENES = sorted(members.index)
ORGANS = common.sort_organs(S.organ.unique())   # fixed organ order of the original analysis


def z_matrix(s):
    """Signed burden z (|z| from p, capped at sqrt(50); carriers >= 10), genes x ICD-10 codes, for burden set s."""
    out = {}
    for r in S.itertuples():
        f = config.GENEBASS_ICD_DIR / f"{r.icd}_{s}.tsv.gz"
        if not f.exists():
            continue
        g = pd.read_csv(f, sep="\t")[["gene_id", "BETA_Burden", "Pvalue_Burden"]].dropna().drop_duplicates("gene_id")
        g = g.merge(A[["gene_id", f"keep_{s}", f"carriers_{s}"]], on="gene_id")
        g = g[(g[f"keep_{s}"] == True) & (g[f"carriers_{s}"] >= 10)].set_index("gene_id")  # noqa: E712
        z = np.clip(norm.isf(g.Pvalue_Burden.clip(1e-300).values / 2), None, np.sqrt(50)) * np.sign(g.BETA_Burden.values)
        out[r.icd] = pd.Series(z, index=g.index)
    return pd.DataFrame(out)


def organ_values(Z, signed=True):
    """Doubly centred z (or z^2) averaged over the codes of each organ: genes x organs."""
    Z = Z.reindex(GENES)
    if not signed:
        Z = Z ** 2
    Zc = Z - Z.mean(axis=0)
    Zc = Zc.sub(Zc.mean(axis=1), axis=0)
    code_to_organ = dict(zip(S.icd, S.organ))
    return pd.DataFrame({o: Zc[[c for c in Zc.columns if code_to_organ[c] == o]].mean(axis=1) for o in ORGANS})


def test(G, name, drop_organ=None):
    organs = [o for o in ORGANS if o != drop_organ]
    V = G[organs].values
    labels = [members[g] & set(organs) for g in G.index]

    def mask(labs):
        x = np.zeros(V.shape, bool)
        for i, lab in enumerate(labs):
            for o in lab:
                x[i, organs.index(o)] = True
        return x & np.isfinite(V)
    mk = mask(labels)
    obs = V[mk].mean()
    perm = np.array([V[mask([labels[i] for i in rng.permutation(len(labels))])].mean()
                     for _ in range(B if drop_organ is None else 2000)])
    p = (1 + (perm >= obs).sum()) / (len(perm) + 1)
    print(f"  {name:34s} cells {mk.sum():4d}: {obs:+.4f} (perm {perm.mean():+.4f} +- {perm.std():.4f}; "
          f"dz {(obs - perm.mean()) / perm.std():+.1f}) p {p:.4f}", flush=True)
    return dict(test=name, n=int(mk.sum()), stat=obs, dz=(obs - perm.mean()) / perm.std(), p=p)


rows = []
Zs = {s: z_matrix(s) for s in ["plof", "mislc", "syn"]}
print("Codes loaded:", {s: z.shape[1] for s, z in Zs.items()})
print("\n=== Primary: directional (carriers diagnosed more often in their own organ) ===")
for s in ["plof", "mislc", "syn"]:
    rows.append(test(organ_values(Zs[s]), f"directional {s}"))
print("\n=== Secondary: direction-free (z2) ===")
for s in ["plof", "mislc", "syn"]:
    rows.append(test(organ_values(Zs[s], signed=False), f"z2 {s}"))
print("\n=== Leave-one-organ-out, directional pLoF ===")
G = organ_values(Zs["plof"])
for o in ORGANS:
    rows.append(test(G, f"without {o}", drop_organ=o))
print("\n=== Per organ (directional, mean of matching genes) ===")
for s in ["plof", "mislc"]:
    G = organ_values(Zs[s])
    print(f"  {s}: " + ", ".join(f"{o} {G.loc[[g for g in G.index if o in members[g]], o].mean():+.3f}/"
                                 f"{sum(o in members[g] for g in G.index)}" for o in ORGANS))
print("\n=== Strongest matching gene-diagnosis pairs (pLoF, z > 4) ===")
Z = Zs["plof"].reindex(GENES)
code_to_organ = dict(zip(S.icd, S.organ))
code_to_description = dict(zip(S.icd, S.description))
hits = [(A.set_index("gene_id").symbol.get(g), c, Z.loc[g, c], code_to_description[c][5:70]) for g in GENES for c in Z.columns
        if code_to_organ[c] in members[g] and np.isfinite(Z.loc[g, c]) and Z.loc[g, c] > 4]
for h in sorted(hits, key=lambda x: -x[2])[:25]:
    print(f"  {h[0]:10s} {h[1]} z {h[2]:+.1f}  {h[3]}")
pd.DataFrame(rows).to_csv(config.results("39_icd10_diagnoses_diagonal.tsv"), sep="\t", index=False, float_format="%.4f")
