"""Organ diagonal without genes of clonal haematopoiesis (CHIP).

Exploratory sensitivity analysis. Re-implementation of an exploratory analysis whose original output was only a log.
CHIP genes: DNMT3A, ASXL1, JAK2, PPM1D, TET2. Of these, only TET2 is in the recessive organ table
(organ_fractions_recessive); its 11 gene-organ rows (all organs) are removed, so the gene drops out of centring and
permutation entirely. Otherwise unchanged from script 20:
  V0 all traits | V3 locus-adjusted | V1+V3 fixed primary trait + locus; each for pLoF and missense|LC, B = 1000,
  one random stream (seed 7) in this order.
A broader version (blood/immune rows only, extended gene list, gene-specific direction) is script 46.

Inputs: organ_fractions_recessive.tsv.gz (14), gene_annotation.tsv.gz (03), Genebass burden files. Output: log only.
Call: python 20b_organ_diagonal_without_chip_genes.py
"""
import numpy as np
import pandas as pd
from scipy.stats import norm

import common
import config

rng = np.random.default_rng(7)
B = 1000
CHIP = ["DNMT3A", "ASXL1", "JAK2", "PPM1D", "TET2"]
ORGAN_TRAITS = common.ORGAN_TRAITS
PRIMARY_TRAIT = {"blood": ["hb"], "immune": ["wbc"], "kidney": ["cysc"], "liver": ["alt"], "lipids": ["ldl"], "glucose": ["hba1c"],
                 "lung": ["fev1"], "muscle": ["grip"], "nervous": ["fluidint"], "skeletal": ["height"], "circulatory": ["sbp"]}
ALL_TRAITS = [t for v in ORGAN_TRAITS.values() for t in v]
organ_fractions = pd.read_csv(config.results("organ_fractions_recessive.tsv.gz"), sep="\t")
annotation = pd.read_csv(config.results("gene_annotation.tsv.gz"), sep="\t")
n0 = len(organ_fractions)
organ_fractions = organ_fractions[~organ_fractions.gene_id.isin(set(annotation[annotation.symbol.isin(CHIP)].gene_id))]
print(f"removed CHIP gene-organ rows: {n0 - len(organ_fractions)}")


def z_matrix(mask):
    """Signed burden z per gene x trait (sign of BETA, |z| from the burden p), carriers >= 10 and gene passes QC."""
    out = {}
    for t in ALL_TRAITS:
        g = pd.read_csv(config.GENEBASS_DIR / f"genebass_{t}_{mask}.tsv.gz", sep="\t")[
            ["gene_id", "BETA_Burden", "Pvalue_Burden"]].dropna().drop_duplicates("gene_id")
        g = g.merge(annotation[["gene_id", f"keep_{mask}", f"carriers_{mask}"]], on="gene_id")
        g = g[(g[f"keep_{mask}"] == True) & (g[f"carriers_{mask}"] >= 10)].set_index("gene_id")  # noqa: E712
        out[t] = pd.Series(np.sign(g.BETA_Burden.values) * norm.isf(g.Pvalue_Burden.clip(1e-300).values / 2), index=g.index)
    return pd.DataFrame(out)


ZS = {mask: z_matrix(mask) for mask in ["plof", "mislc", "syn"]}


def center(Z2):
    Z2 = Z2.clip(upper=50)
    Zc = Z2 - Z2.mean(axis=0)
    return Zc.sub(Zc.mean(axis=1), axis=0)


def run(mask, frac=0.2, traits=ORGAN_TRAITS, locus=False):
    """As run() in script 20 (without drop_top/drop_organ)."""
    members = organ_fractions[organ_fractions.frac >= frac]
    genes = sorted(set(members.gene_id) & set(ZS[mask].index))
    Z2 = ZS[mask].reindex(genes) ** 2
    if locus and mask != "syn":
        S = ZS["syn"].reindex(genes)
        for t in Z2.columns:
            ok = Z2[t].notna() & S[t].notna()
            bb = np.polyfit(S.loc[ok, t] ** 2, Z2.loc[ok, t], 1)[0] if ok.sum() > 50 else 0
            Z2[t] = Z2[t] - max(bb, 0) * S[t] ** 2
            Z2.loc[S[t].abs() > 1.96, t] = np.nan
    Zc = center(Z2)
    organs = list(traits)
    L = np.array([[1.0 if g in set(members[members.organ == o].gene_id) else 0.0 for o in organs] for g in genes])
    C = np.array([[1.0 if t in traits[o] else 0.0 for t in Zc.columns] for o in organs])
    V = Zc.values
    valid = ~np.isnan(V)
    V0 = np.nan_to_num(V)

    def stat(Lx):
        m = ((Lx @ C) > 0) & valid
        return (V0 * m).sum() / m.sum()

    obs = stat(L)
    perm = np.array([stat(L[rng.permutation(len(L))]) for _ in range(B)])
    return obs, (obs - perm.mean()) / perm.std(), (1 + (perm >= obs).sum()) / (B + 1), len(genes)


for name, kw in [("V0 without CHIP", {}), ("V3 locus without CHIP", {"locus": True}),
                 ("V1+V3 primary locus without CHIP", {"traits": PRIMARY_TRAIT, "locus": True})]:
    out = []
    for mask in ["plof", "mislc"]:
        obs, dz, p, n = run(mask, **kw)
        out.append(f"{mask} dz {dz:.1f} (p {p:.3f})")
    print(f"{name} | " + " | ".join(out), flush=True)
