"""Gene-specific expected direction of the carrier effect in UK Biobank (test T10-7, PRAEREG_vorzeichen.md part C).

Instead of one coarse direction per organ (script 30), every recessive gene gets its own expected direction per trait,
derived from the direction-bearing HPO terms of its OMIM diseases (direction_map.py; clinical descriptions, independent
of UK Biobank). Within one trait, different genes expect different directions, so the confounder found in script 30b
(the disease direction often coincides with the generic direction of a gene loss) no longer applies.

Matrix formulation: S (genes x traits, +1/-1/0) is the expectation map, Zc the centred signed Genebass burden z.
Statistic = mean of S*Zc over the cells with an expectation. Permutation = shuffling the rows of S between genes
(B = 1000). Also: robustness (locus-adjusted, without the top genes), sign fractions, and a table per trait x direction.

Inputs:  direction_map.py inputs (HPO, OMIM via HPO, Genebass), gene_annotation.tsv.gz, organ_fractions_recessive.tsv.gz
Outputs: 31_gene_specific_direction.tsv (per set x trait x expected sign: n, mean directed z, fraction positive)
Usage:   python 31_gene_specific_direction.py
"""
import numpy as np
import pandas as pd
from scipy.stats import binomtest

import config
from direction_map import build_direction_map

dm = build_direction_map()
GENES, S, TRAITS, RAW, CEN = dm.genes, dm.S, dm.traits, dm.raw, dm.centred
rng = np.random.default_rng(5)
B = 1000
Zsyn_c = CEN["syn"]


def stat(Z, Sm, mask):
    v = Sm * Z
    m = mask & (Sm != 0) & np.isfinite(Z)
    return np.nansum(np.where(m, v, 0)) / m.sum(), int(m.sum())


def test(s, locus=False, drop=None, name="", raw=False):
    Z = RAW[s] if raw else CEN[s]
    mask = np.isfinite(Z)
    if locus:   # locus adjustment: drop cells whose synonymous z is itself significant
        mask = mask & (np.abs(np.nan_to_num(Zsyn_c, nan=0.0)) <= 1.96)
    if drop is not None:
        mask = mask & ~drop[:, None]
    obs, n = stat(Z, S, mask)
    pm = np.empty(B)
    for b in range(B):
        pm[b] = stat(Z, S[rng.permutation(len(GENES))], mask)[0]
    dz = (obs - pm.mean()) / pm.std()
    p = (1 + (pm >= obs).sum()) / (B + 1)
    print(f"  {name or s:34s} d {obs:+.3f} (perm {pm.mean():+.3f} +- {pm.std():.3f}) dz {dz:+5.1f}  p {p:.4f}  n {n}", flush=True)
    return obs, dz, p


print("\n=== T10-7: gene-specific direction (main test), replication, negative control ===", flush=True)
for s, label in [("plof", "T10-7 pLoF"), ("mislc", "missense"), ("syn", "synonymous (control)")]:
    test(s, name=label)

print("\n=== Robustness ===", flush=True)
mx = np.nanmax(np.where(S != 0, np.abs(np.nan_to_num(CEN["plof"], nan=0.0)), np.nan), axis=1)
order = np.argsort(-np.nan_to_num(mx, nan=-1))
for s in ["plof", "mislc"]:
    print(f" [{s}]", flush=True)
    test(s, locus=True, name="locus-adjusted")
    for k in (10, 25):
        d = np.zeros(len(GENES), bool)
        d[order[:k]] = True
        test(s, drop=d, name=f"without top {k} genes")

print("\n=== Sign fraction (raw z; permutation of the direction rows) ===", flush=True)
for s in ["plof", "mislc", "syn"]:
    Z = RAW[s]
    for lo, hi, label in [(0, np.inf, "all"), (2, np.inf, "|z| > 2"), (3, np.inf, "|z| > 3")]:
        base = np.isfinite(Z) & (np.abs(Z) >= lo) & (np.abs(Z) < hi)

        def frac(Sm):
            m = base & (Sm != 0)
            if not m.sum():
                return np.nan, 0
            return ((Sm * Z)[m] > 0).mean(), int(m.sum())
        obs, n = frac(S)
        if not n:
            continue
        pm = np.array([frac(S[rng.permutation(len(GENES))])[0] for _ in range(B)])
        p = (1 + (pm >= obs).sum()) / (B + 1)
        bt = binomtest(int(round(obs * n)), n, 0.5, alternative="greater").pvalue
        print(f"  {s:6s} {label:10s} {obs:.3f} (n {n}; perm {pm.mean():.3f} +- {pm.std():.3f}, "
              f"dz {(obs - pm.mean()) / pm.std():+4.1f})  p_perm {p:.4f}  p_binom {bt:.2e}", flush=True)

# per trait x expected direction
rows = []
for s in ["plof", "mislc", "syn"]:
    Z = CEN[s]
    for j, t in enumerate(TRAITS):
        for sign in (+1, -1):
            m = (S[:, j] == sign) & np.isfinite(Z[:, j])
            if m.sum():
                v = sign * Z[m, j]
                rows.append(dict(set=s, trait=t, expected_sign=sign, n=int(m.sum()), d_mean=v.mean(), fraction_positive=(v > 0).mean()))
T = pd.DataFrame(rows)
T.to_csv(config.results("31_gene_specific_direction.tsv"), sep="\t", index=False, float_format="%.4f")
pd.set_option("display.width", 200)
print("\n=== Per trait x expected direction (pLoF), by contribution ===", flush=True)
print(T[T.set == "plof"].sort_values("d_mean", ascending=False).round(3).to_string(index=False), flush=True)
print("\n=== Traits with both directions, pLoF ===", flush=True)
print(T[(T.set == "plof") & T.trait.isin(dm.bidirectional_traits)].sort_values(["trait", "expected_sign"]).round(3).to_string(index=False), flush=True)
