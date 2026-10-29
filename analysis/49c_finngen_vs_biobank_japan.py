"""Cross-ancestry reading FinnGen against BioBank Japan (X1, X2; PRAEREG_bbj.md).

X1  Spearman correlation of gene scores (mean of expected sign x z over the direction cells of the gene-specific
    direction map, pLoF; Stouffer gene z from script 49) between FinnGen and BioBank Japan; permutation of gene labels
    (B = 2000, one-sided). Coding variants as a second line.
X2  Fisher combination of the directional tests (script 49, part 4) across FinnGen, BioBank Japan and Genes & Health
    (script 35, variant based).

Inputs:  cache/replication_signed_gene_stats_fg/, cache/replication_signed_gene_stats_bbj/ (written by script 49),
         49_finngen_4.tsv, 49_biobank_japan_4.tsv, 35_directional_variant_replication.tsv, direction_map.py inputs
Outputs: log only
Usage:   python 49c_finngen_vs_biobank_japan.py    (after script 49 for FG and BBJ)
"""
import numpy as np
import pandas as pd
from scipy.stats import combine_pvalues, spearmanr

import config
from direction_map import build_direction_map

dm = build_direction_map()
rng = np.random.default_rng(49)


def load_signed(cohort, s):
    out = {}
    for f in config.cache(f"replication_signed_gene_stats_{cohort.lower()}").glob(f"{cohort}_*.tsv.gz"):
        t = f.name[len(cohort) + 1: -len(".tsv.gz")]
        out[t] = pd.read_csv(f, sep="\t", index_col=0)[s].dropna()
    return out


def gene_scores(Z):   # as P2 in replication_tests.py
    out = {}
    for g, r in dm.directions.items():
        v = [sign * Z[t].get(g, np.nan) for t, sign in r.items() if t in Z]
        v = [a for a in v if np.isfinite(a)]
        if v:
            out[g] = np.mean(v)
    return pd.Series(out)


print("=== X1 gene scores FinnGen against BioBank Japan ===")
for s in ("plof", "coding"):
    a, b = gene_scores(load_signed("FG", s)), gene_scores(load_signed("BBJ", s))
    j = pd.concat([a, b], axis=1, keys=["a", "b"]).dropna()
    r = spearmanr(j.a, j.b).correlation
    null = np.array([spearmanr(j.a.values, rng.permutation(j.b.values)).correlation for _ in range(2000)])
    p = (1 + (null >= r).sum()) / 2001
    print(f"  {s}: rho {r:+.3f} (n {len(j)}), p {p:.4f}")

print("\n=== X2 Fisher across FinnGen, BioBank Japan, Genes & Health (directional test, primary) ===")
fg, bbj, gh = (pd.read_csv(config.results(f), sep="\t")
               for f in ("49_finngen_4.tsv", "49_biobank_japan_4.tsv", "35_directional_variant_replication.tsv"))
for s in ("plof", "coding"):
    ps = [x[(x.analysis == "primary") & (x.set == s) & (x.cohort == k)].p.iloc[0] for x, k in ((fg, "FG"), (bbj, "BBJ"), (gh, "GH"))]
    print(f"  {s}: FinnGen p {ps[0]:.4f}, BioBank Japan p {ps[1]:.4f}, G&H p {ps[2]:.4f} -> Fisher p {combine_pvalues(ps, method='fisher')[1]:.2g}")
