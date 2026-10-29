"""Generalised conservation: does an organ-matching carrier effect predict selection against heterozygotes (s_het) across all
organs, and better than an arbitrary (non-matching) carrier effect?

Exploratory (test register no. 12). For each recessive gene with an organ assignment from script 14 (frac >= 0.2):
matching = minimum burden p over the organ-matching traits (Bonferroni over their number); non-matching = the same over all
non-matching traits (Bonferroni). Groups: P (matching p_bonf < 0.01), U (only non-matching), N (neither). The groups are
compared for s_het (GeneBayes posterior mean) and LoF obs/exp, and by OLS of log10 s_het on both indicators (adjusted for
log10 CDS length). Missense|LC (variants other than those underlying s_het) is the main analysis, pLoF is shown in addition.

Inputs: organ_fractions_recessive.tsv.gz (14), shet_merged.tsv.gz (11), Genebass burden files. Output: log only.
Call: python 18_organ_matched_selection.py
"""
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import mannwhitneyu

import common
import config

ORGAN_TRAITS = common.ORGAN_TRAITS
ALL_TRAITS = [t for v in ORGAN_TRAITS.values() for t in v]
organ_fractions = pd.read_csv(config.results("organ_fractions_recessive.tsv.gz"), sep="\t")
members = organ_fractions[organ_fractions.frac >= 0.2]
A = pd.read_csv(config.results("shet_merged.tsv.gz"), sep="\t")[
    ["gene_id", "symbol", "set", "post_mean", "oe", "carriers_plof", "carriers_mislc", "oe_lof_upper", "cds_length"]]
for mask in ["mislc", "plof"]:
    P = {}
    for t in ALL_TRAITS:
        g = pd.read_csv(config.GENEBASS_DIR / f"genebass_{t}_{mask}.tsv.gz", sep="\t")[
            ["gene_id", "Pvalue_Burden"]].dropna().drop_duplicates("gene_id")
        P[t] = g.set_index("gene_id").Pvalue_Burden
    P = pd.DataFrame(P)
    rows = []
    for gid, sub in members.groupby("gene_id"):
        if gid not in P.index:
            continue
        matching_traits = sorted({t for o in sub.organ for t in ORGAN_TRAITS[o]})
        other_traits = [t for t in ALL_TRAITS if t not in matching_traits]
        pm = P.loc[gid, matching_traits].dropna()
        pu = P.loc[gid, other_traits].dropna()
        if len(pm) == 0 or len(pu) == 0:
            continue
        rows.append((gid, min(1, pm.min() * len(pm)), min(1, pu.min() * len(pu))))
    R = pd.DataFrame(rows, columns=["gene_id", "pm", "pu"]).merge(A, on="gene_id")
    R = R[R[f"carriers_{mask}"] >= 10].dropna(subset=["post_mean"])
    R["group"] = np.select([R.pm < 0.01, R.pu < 0.01], ["P matching", "U non-matching only"], "N none")
    print(f"\n[{mask}] recessive genes with an organ profile: {len(R)}")
    print(R.groupby("group").agg(n=("gene_id", "size"), s_het=("post_mean", "median"), obs_exp=("oe", "median"),
                                 carriers=(f"carriers_{mask}", "median")).round(4).to_string())
    a = R[R.group == "P matching"].post_mean
    u = R[R.group == "U non-matching only"].post_mean
    n = R[R.group == "N none"].post_mean
    print(f"  P vs U: MW p {mannwhitneyu(a, u, alternative='greater').pvalue:.1e} | "
          f"P vs N: {mannwhitneyu(a, n, alternative='greater').pvalue:.1e} | "
          f"U vs N: {mannwhitneyu(u, n, alternative='greater').pvalue:.1e}")
    R["lshet"] = np.log10(R.post_mean)
    R["ltr"] = np.log10(R[f"carriers_{mask}"])
    R["lcds"] = np.log10(R.cds_length)
    R["matching"] = (R.pm < 0.01).astype(int)
    R["non_matching"] = (R.pu < 0.01).astype(int)
    f = smf.ols("lshet ~ matching + non_matching + lcds", data=R).fit()
    print(f"  OLS log10 s_het: matching {f.params['matching']:+.3f} (p {f.pvalues['matching']:.1e}), "
          f"non-matching {f.params['non_matching']:+.3f} (p {f.pvalues['non_matching']:.1e})")
    # per organ
    R2 = R.merge(members[["gene_id", "organ"]], on="gene_id")
    per_organ = []
    for o in common.sort_organs(R2.organ.unique()):
        v = R2[R2.organ == o]
        if v.matching.sum() >= 3:
            per_organ.append(f"{o}: {v[v.matching == 1].post_mean.median():.3f}/"
                             f"{v[(v.matching == 0) & (v.non_matching == 0)].post_mean.median():.3f} ({int(v.matching.sum())})")
    print("  per organ (median s_het P / N, n):", " | ".join(per_organ))
