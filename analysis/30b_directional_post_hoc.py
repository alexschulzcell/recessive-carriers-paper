"""Post hoc analysis of the directional organ-diagonal test of script 30 (exploratory, not pre-specified).

In script 30 the mean directional value d is clearly positive, but the fraction of d > 0 is only about 0.5: the signal
sits in the tail, not in a broad shift of all genes. Question here: when a carrier effect is detectable at all, does it
point in the disease direction? Statistic: fraction of matching cells with the expected sign, stratified by |z|
(raw signed z, not centred, so that the sign is the raw one).
Null distribution: permutation of organ labels between genes (B = 1000), as in script 30.

Inputs: gene_annotation.tsv.gz, organ_fractions_recessive.tsv.gz (results folder), Genebass burden files.
Output: 30b_directional_post_hoc.tsv.
Call: python 30b_directional_post_hoc.py
"""
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.stats import norm, binomtest

import common
import config

rng = np.random.default_rng(5)
B = 1000
UNAMBIGUOUS = {("blood", "hb"): -1, ("blood", "rbc"): -1, ("blood", "plt"): -1,
               ("immune", "wbc"): -1, ("immune", "lymph"): -1,
               ("kidney", "creat"): +1, ("kidney", "cysc"): +1, ("kidney", "urea"): +1,
               ("liver", "alt"): +1, ("liver", "ast"): +1, ("liver", "ggt"): +1, ("liver", "bili"): +1,
               ("liver", "albumin"): -1, ("lipids", "hdl"): -1, ("lipids", "tg"): +1,
               ("lung", "fev1"): -1, ("lung", "fvc"): -1, ("muscle", "grip"): -1,
               ("nervous", "fluidint"): -1, ("nervous", "reaction"): +1,
               ("skeletal", "height"): -1, ("skeletal", "sitheight"): -1}
# fixed organ order of the original analysis (keeps permutation streams reproducible)
ORGANS = common.sort_organs({o for o, _ in UNAMBIGUOUS})
TRAITS = sorted({t for _, t in UNAMBIGUOUS})
annotation = pd.read_csv(config.results("gene_annotation.tsv.gz"), sep="\t")
organ_fractions = pd.read_csv(config.results("organ_fractions_recessive.tsv.gz"), sep="\t")
members = organ_fractions[organ_fractions.frac >= 0.2]
M = {o: set(members[members.organ == o].gene_id) for o in ORGANS}
genes = sorted(set(members.gene_id))
gene_array = np.array(genes)
labels = {g: [o for o in ORGANS if g in M[o]] for g in genes}


def signed_z(trait, mask):
    f = config.GENEBASS_DIR / f"genebass_{trait}_{mask}.tsv.gz"
    if not f.exists():
        return None
    g = pd.read_csv(f, sep="\t")[["gene_id", "BETA_Burden", "Pvalue_Burden"]].dropna()
    g = g.merge(annotation[["gene_id", f"keep_{mask}", f"carriers_{mask}"]], on="gene_id")
    g = g[(g[f"keep_{mask}"] == True) & (g[f"carriers_{mask}"] >= 10)]  # noqa: E712
    g = g.drop_duplicates("gene_id").set_index("gene_id")
    z = np.clip(norm.isf(g.Pvalue_Burden.clip(1e-300).values / 2), None, np.sqrt(50))
    return pd.Series(np.sign(g.BETA_Burden.values) * z, index=g.index)


Z = {mask: pd.DataFrame({t: z for t in TRAITS if (z := signed_z(t, mask)) is not None}).reindex(genes)
     for mask in ["plof", "mislc", "syn"]}


def permuted_sets():
    shuffled = rng.permutation(gene_array)
    Mp = defaultdict(set)
    for g_new, g_old in zip(gene_array, shuffled):
        for o in labels[g_old]:
            Mp[o].add(g_new)
    return Mp


def fraction_expected(Zr, Mx, lo, hi):
    """Number of matching cells with the expected sign and number of cells, |z| in [lo, hi)."""
    hits = n = 0
    for (o, t), sign in UNAMBIGUOUS.items():
        if t not in Zr.columns:
            continue
        v = Zr.loc[sorted(Mx[o] & set(Zr.index)), t].dropna()
        v = v[(v.abs() >= lo) & (v.abs() < hi)]
        hits += (np.sign(v) == sign).sum()
        n += len(v)
    return hits, n


print("Fraction of matching gene x trait cells whose carrier effect points in the disease direction.")
print("Permutation of organ labels (B = 1000); one-sided p.\n")
rows = []
for mask in ["plof", "mislc", "syn"]:
    print(f"[{mask}]")
    for lo, hi, label in [(0, np.inf, "all"), (0, 2, "|z| < 2 (no signal)"),
                          (2, 3, "|z| 2-3"), (3, np.inf, "|z| > 3 (clear signal)"),
                          (2, np.inf, "|z| > 2")]:
        hits, n = fraction_expected(Z[mask], M, lo, hi)
        if n == 0:
            continue
        pm = np.array([(lambda r: r[0] / r[1] if r[1] else np.nan)(fraction_expected(Z[mask], permuted_sets(), lo, hi))
                       for _ in range(B)])
        obs = hits / n
        p = (1 + (pm >= obs).sum()) / (B + 1)
        bt = binomtest(hits, n, 0.5, alternative="greater").pvalue
        print(f"  {label:24s} {hits:5d}/{n:5d} = {obs:.3f}  (perm {np.nanmean(pm):.3f} +- {np.nanstd(pm):.3f}, "
              f"dz {(obs - np.nanmean(pm)) / np.nanstd(pm):+4.1f})  p_perm {p:.4f}  p_binom {bt:.2e}")
        rows.append(dict(set=mask, stratum=label, hits=hits, n=n, fraction=obs,
                         perm=np.nanmean(pm), p_perm=p, p_binom=bt))
    print()

# per organ, |z| > 2, pLoF
print("[plof] |z| > 2, per organ:")
for o in ORGANS:
    hits = n = 0
    for (oo, t), sign in UNAMBIGUOUS.items():
        if oo != o or t not in Z["plof"].columns:
            continue
        v = Z["plof"].loc[sorted(M[o] & set(Z["plof"].index)), t].dropna()
        v = v[v.abs() >= 2]
        hits += (np.sign(v) == sign).sum()
        n += len(v)
    if n:
        print(f"  {o:10s} {hits:3d}/{n:3d} = {hits / n:.2f}")
pd.DataFrame(rows).to_csv(config.results("30b_directional_post_hoc.tsv"), sep="\t", index=False, float_format="%.4f")
