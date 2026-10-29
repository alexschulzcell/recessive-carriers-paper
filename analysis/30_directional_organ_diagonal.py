"""Directional organ-diagonal test (pre-specified, PRAEREG_vorzeichen.md, parts A/B; tests T10-1 to T10-6).

The organ-diagonal tests of scripts 14/20 use z^2 and are therefore direction-free. Here every matching organ x trait
cell gets a pre-specified expected sign from disease biology. Directional value d = expected_sign * centred z;
d > 0 means "carriers shift in the direction of the disease".

Centring as in 14/20: first the trait mean, then the gene mean, on the raw signed z; only then multiplied by the
expected sign. Null distribution: permutation of organ labels between genes (B = 1000).

T10-1 main test pLoF | T10-2 missense | T10-3 synonymous (negative control) | T10-4 fraction with expected sign
T10-5 robustness (locus-adjusted, without top 10/25 genes, fixed primary trait, leave-one-organ-out)
T10-6 mixed cells separately

Inputs: gene_annotation.tsv.gz, organ_fractions_recessive.tsv.gz (results folder), Genebass burden files.
Output: 30_directional_organ_diagonal.tsv (one row per mask x cell).
Call: python 30_directional_organ_diagonal.py
"""
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.stats import norm, binomtest

import common
import config

rng = np.random.default_rng(5)
B = 1000

# --- pre-specified (PRAEREG_vorzeichen.md part A) ---
UNAMBIGUOUS = {  # (organ, trait): expected sign
    ("blood", "hb"): -1, ("blood", "rbc"): -1, ("blood", "plt"): -1,
    ("immune", "wbc"): -1, ("immune", "lymph"): -1,
    ("kidney", "creat"): +1, ("kidney", "cysc"): +1, ("kidney", "urea"): +1,
    ("liver", "alt"): +1, ("liver", "ast"): +1, ("liver", "ggt"): +1, ("liver", "bili"): +1,
    ("liver", "albumin"): -1,
    ("lipids", "hdl"): -1, ("lipids", "tg"): +1,
    ("lung", "fev1"): -1, ("lung", "fvc"): -1,
    ("muscle", "grip"): -1,
    ("nervous", "fluidint"): -1, ("nervous", "reaction"): +1,
    ("skeletal", "height"): -1, ("skeletal", "sitheight"): -1,
}
MIXED = {  # declared in advance as direction-mixed: tendency only, not part of the main test
    ("blood", "mcv"): -1, ("kidney", "urate"): +1, ("kidney", "phosphate"): -1,
    ("lipids", "ldl"): -1, ("lipids", "chol"): -1,
    ("glucose", "glucose"): +1, ("glucose", "hba1c"): +1,
    ("skeletal", "ebmd"): -1, ("circulatory", "sbp"): -1,
}
PRIMARY = {"blood": "hb", "immune": "wbc", "kidney": "creat", "liver": "bili", "lipids": "hdl",
           "lung": "fev1", "muscle": "grip", "nervous": "fluidint", "skeletal": "height"}
# fixed organ order of the original analysis (keeps permutation streams reproducible)
ORGANS = common.sort_organs({o for o, _ in UNAMBIGUOUS} | {o for o, _ in MIXED})
TRAITS = sorted({t for _, t in UNAMBIGUOUS} | {t for _, t in MIXED})

annotation = pd.read_csv(config.results("gene_annotation.tsv.gz"), sep="\t")
organ_fractions = pd.read_csv(config.results("organ_fractions_recessive.tsv.gz"), sep="\t")
members = organ_fractions[organ_fractions.frac >= 0.2]
M = {o: set(members[members.organ == o].gene_id) for o in ORGANS}
genes = sorted(set(members.gene_id))
print(f"Recessive genes with an organ profile: {organ_fractions.gene_id.nunique()}; in the test: {len(genes)}")
print("Assignment:", {o: len(M[o]) for o in ORGANS})


def signed_z(trait, mask):
    """Signed z = sign(BETA) * |z|, filters and cap as in 14/20."""
    f = config.GENEBASS_DIR / f"genebass_{trait}_{mask}.tsv.gz"
    if not f.exists():
        return None
    g = pd.read_csv(f, sep="\t")[["gene_id", "BETA_Burden", "Pvalue_Burden"]].dropna()
    g = g.merge(annotation[["gene_id", f"keep_{mask}", f"carriers_{mask}"]], on="gene_id")
    g = g[(g[f"keep_{mask}"] == True) & (g[f"carriers_{mask}"] >= 10)]  # noqa: E712
    g = g.drop_duplicates("gene_id").set_index("gene_id")
    z = norm.isf(g.Pvalue_Burden.clip(1e-300).values / 2)
    z = np.clip(z, None, np.sqrt(50))              # z^2 capped at 50, as in 14/20
    return pd.Series(np.sign(g.BETA_Burden.values) * z, index=g.index)


def centred(mask):
    Z = pd.DataFrame({t: z for t in TRAITS if (z := signed_z(t, mask)) is not None}).reindex(genes)
    Zc = Z - Z.mean(axis=0)                        # trait effect
    Zc = Zc.sub(Zc.mean(axis=1), axis=0)           # gene effect
    return Z, Zc


ZS = {mask: centred(mask) for mask in ["plof", "mislc", "syn"]}
Zsyn = ZS["syn"][1]
labels = {g: [o for o in ORGANS if g in M[o]] for g in genes}
gene_array = np.array(genes)


def permuted_sets():
    """One permuted organ -> genes mapping, keeping multiple assignments (as in 14/20)."""
    shuffled = rng.permutation(gene_array)
    Mp = defaultdict(set)
    for g_new, g_old in zip(gene_array, shuffled):
        for o in labels[g_old]:
            Mp[o].add(g_new)
    return Mp


def directional_values(Zc, Mx, cells, drop=frozenset(), locus=False, organs=None):
    """Directional values d of all matching cells."""
    out = []
    for (o, t), sign in cells.items():
        if t not in Zc.columns or (organs is not None and o not in organs):
            continue
        gs = (Mx[o] & set(Zc.index)) - drop
        v = Zc.loc[sorted(gs), t]
        if locus and t in Zsyn.columns:
            v = v[Zsyn.loc[v.index, t].abs() <= 1.96]
        out.append(sign * v.dropna().values)
    return np.concatenate([x for x in out if len(x)]) if out else np.array([])


def test(mask, cells=UNAMBIGUOUS, drop=frozenset(), locus=False, organs=None, name=""):
    Zc = ZS[mask][1]
    obs = directional_values(Zc, M, cells, drop, locus, organs)
    if len(obs) == 0:
        return None
    o = obs.mean()
    pm = np.array([directional_values(Zc, permuted_sets(), cells, drop, locus, organs).mean() for _ in range(B)])
    dz = (o - pm.mean()) / pm.std()
    p = (1 + (pm >= o).sum()) / (B + 1)
    print(f"  {name or mask:34s} d {o:+.3f} (perm {pm.mean():+.3f} +- {pm.std():.3f}) dz {dz:+5.1f}  p {p:.4f}  n {len(obs)}")
    return o, dz, p, obs


print("\n=== T10-1/2/3: main test, replication, negative control (22 unambiguous cells) ===")
res = {}
for mask, label in [("plof", "T10-1 pLoF"), ("mislc", "T10-2 missense"), ("syn", "T10-3 synonymous")]:
    res[mask] = test(mask, name=label)

print("\n=== T10-4: fraction with expected sign ===")
for mask in ["plof", "mislc", "syn"]:
    Zc = ZS[mask][1]
    d = res[mask][3]
    # non-matching cells: expected sign of the trait, but genes of other organs
    non_matching = []
    for (o, t), sign in UNAMBIGUOUS.items():
        if t not in Zc.columns:
            continue
        other = (set(Zc.index) - M[o])
        non_matching.append(sign * Zc.loc[sorted(other), t].dropna().values)
    non_matching = np.concatenate(non_matching)
    bt = binomtest((d > 0).sum(), len(d), 0.5)
    bu = binomtest((non_matching > 0).sum(), len(non_matching), 0.5)
    print(f"  {mask:6s} matching {(d > 0).mean():.3f} (n {len(d)}, p {bt.pvalue:.2e}) | "
          f"non-matching {(non_matching > 0).mean():.3f} (n {len(non_matching)}, p {bu.pvalue:.2f})")

print("\n=== T10-5: robustness (pLoF / missense) ===")
# top genes by |z| in the matching trait (pLoF)
Zp = ZS["plof"][1]
rank = []
for (o, t), _ in UNAMBIGUOUS.items():
    if t in Zp.columns:
        for g in sorted(M[o] & set(Zp.index)):
            if pd.notna(Zp.loc[g, t]):
                rank.append((abs(Zp.loc[g, t]), g))
top = [g for _, g in sorted(rank, reverse=True)]
top = list(dict.fromkeys(top))
primary_cells = {(o, t): UNAMBIGUOUS[(o, t)] for o, t in UNAMBIGUOUS if PRIMARY.get(o) == t}
for mask in ["plof", "mislc"]:
    print(f" [{mask}]")
    test(mask, locus=True, name="locus-adjusted")
    for k in (10, 25):
        test(mask, drop=frozenset(top[:k]), name=f"without top {k} genes")
    test(mask, cells=primary_cells, name="fixed primary trait")
    test(mask, locus=True, cells=primary_cells, name="primary trait + locus-adjusted")

print("\n  leave-one-organ-out (pLoF, 10 organs with unambiguous cells):")
organs_unambiguous = common.sort_organs({o for o, _ in UNAMBIGUOUS})
for o in organs_unambiguous:
    test("plof", organs=[x for x in organs_unambiguous if x != o], name=f"without {o}")

print("\n=== each organ separately (unambiguous cells) ===")
for o in organs_unambiguous:
    for mask in ["plof", "mislc"]:
        test(mask, organs=[o], name=f"{o} [{mask}]")

print("\n=== T10-6: mixed cells, separately, descriptive ===")
for mask in ["plof", "mislc", "syn"]:
    test(mask, cells=MIXED, name=f"mixed [{mask}]")

# table per cell
rows = []
for mask in ["plof", "mislc", "syn"]:
    Zc = ZS[mask][1]
    for cells, kind in [(UNAMBIGUOUS, "unambiguous"), (MIXED, "mixed")]:
        for (o, t), sign in cells.items():
            if t not in Zc.columns:
                continue
            v = (sign * Zc.loc[sorted(M[o] & set(Zc.index)), t]).dropna()
            rows.append(dict(set=mask, kind=kind, organ=o, trait=t, expected_sign=sign, n=len(v),
                             d_mean=v.mean(), fraction_positive=(v > 0).mean()))
T = pd.DataFrame(rows)
T.to_csv(config.results("30_directional_organ_diagonal.tsv"), sep="\t", index=False, float_format="%.4f")
pd.set_option("display.width", 200)
print("\n=== single cells (pLoF) ===")
print(T[T.set == "plof"].sort_values("d_mean", ascending=False)
      [["kind", "organ", "trait", "expected_sign", "n", "d_mean", "fraction_positive"]].round(3).to_string(index=False))
