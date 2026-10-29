"""Robustness of the organ diagonal of script 14.

Exploratory (test register no. 19). Each analysis is run for pLoF, missense|LC and synonymous variants; the organ labels are
permuted between genes (B = 1000):
  V0  as script 14 (all traits per organ, frac >= 0.2)
  V1  one FIXED primary trait per organ, chosen a priori: blood Hb, immune WBC, kidney cystatin C, liver ALT, lipids LDL,
      glucose HbA1c, lung FEV1, muscle grip strength, nervous fluid intelligence, skeletal height, circulatory systolic BP
  V2  without top genes: the k genes with the largest matching contribution (k = 10, 25, 50) are removed (also from the
      permutation)
  V3  locus adjustment: z2_adj = z2 - b * z2_syn of the same gene and trait (b per trait from a regression over all genes),
      and cells with |z_syn| > 1.96 are discarded (both; V3 does not apply to the synonymous set itself)
  V4  leave one organ out
  V5  organ threshold frac >= 0.1 / 0.35 / 0.5
Statistic: mean of the doubly centred z^2 (gene and trait effect removed) in matching cells; z^2 capped at 50.
Finally, per organ (primary trait, locus-adjusted): the single diagonal cell against the organ-label permutation.

Inputs: organ_fractions_recessive.tsv.gz (14), gene_annotation.tsv.gz (03), Genebass burden files.
Output: 20_organ_diagonal_robustness.tsv
Call: python 20_organ_diagonal_robustness.py
"""
import numpy as np
import pandas as pd
from scipy.stats import norm

import common
import config

rng = np.random.default_rng(7)
B = 1000
ORGAN_TRAITS = common.ORGAN_TRAITS
PRIMARY_TRAIT = {"blood": ["hb"], "immune": ["wbc"], "kidney": ["cysc"], "liver": ["alt"], "lipids": ["ldl"], "glucose": ["hba1c"],
                 "lung": ["fev1"], "muscle": ["grip"], "nervous": ["fluidint"], "skeletal": ["height"], "circulatory": ["sbp"]}
ALL_TRAITS = [t for v in ORGAN_TRAITS.values() for t in v]
organ_fractions = pd.read_csv(config.results("organ_fractions_recessive.tsv.gz"), sep="\t")
annotation = pd.read_csv(config.results("gene_annotation.tsv.gz"), sep="\t")


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


def run(mask, frac=0.2, traits=ORGAN_TRAITS, drop_top=0, drop_organ=None, locus=False):
    members = organ_fractions[organ_fractions.frac >= frac]
    if drop_organ:
        members = members[members.organ != drop_organ]
    traits = {o: ts for o, ts in traits.items() if o != drop_organ}
    genes = sorted(set(members.gene_id) & set(ZS[mask].index))
    Z = ZS[mask].reindex(genes)
    Z2 = Z ** 2
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
    if drop_top:
        mask_ = (L @ C) > 0
        contrib = np.where(mask_ & valid, V0, 0).sum(1) / np.maximum((mask_ & valid).sum(1), 1)
        contrib[(mask_ & valid).sum(1) == 0] = -np.inf
        keep = np.argsort(-contrib)[drop_top:]
        keep.sort()
        L, V0, valid = L[keep], V0[keep], valid[keep]
        genes = [genes[i] for i in keep]

    def stat(Lx):
        m = ((Lx @ C) > 0) & valid
        return (V0 * m).sum() / m.sum()

    obs = stat(L)
    perm = np.array([stat(L[rng.permutation(len(L))]) for _ in range(B)])
    return obs, (obs - perm.mean()) / perm.std(), (1 + (perm >= obs).sum()) / (B + 1), len(genes)


rows = []
analyses = [("V0 all traits", {}), ("V1 fixed primary trait", {"traits": PRIMARY_TRAIT}),
            ("V2 without top 10", {"drop_top": 10}), ("V2 without top 25", {"drop_top": 25}),
            ("V2 without top 50", {"drop_top": 50}),
            ("V1+V2 primary, without top 25", {"traits": PRIMARY_TRAIT, "drop_top": 25}),
            ("V3 locus-adjusted", {"locus": True}), ("V3 locus + without top 25", {"locus": True, "drop_top": 25}),
            ("V1+V3 primary, locus", {"traits": PRIMARY_TRAIT, "locus": True}),
            ("V5 frac>=0.1", {"frac": 0.1}), ("V5 frac>=0.35", {"frac": 0.35}), ("V5 frac>=0.5", {"frac": 0.5})]
analyses += [(f"V4 without {o}", {"drop_organ": o}) for o in ["lipids", "liver", "kidney", "blood", "immune", "skeletal"]]
for name, kw in analyses:
    r = {"analysis": name}
    for mask in ["plof", "mislc", "syn"]:
        if kw.get("locus") and mask == "syn":
            r["syn"] = "-"
            continue
        obs, dz, p, n = run(mask, **kw)
        r[mask] = f"{dz:5.1f} (p {p:.3f})"
        r["n_genes"] = n
    rows.append(r)
    print(r, flush=True)
R = pd.DataFrame(rows)
R.to_csv(config.results("20_organ_diagonal_robustness.tsv"), sep="\t", index=False)
print("\n" + R.to_string(index=False))

# Per organ with primary trait + locus adjustment: diagonal cell against the organ-label permutation (single p)
print("\nPer organ (primary trait, locus-adjusted, pLoF/missense): dz of the diagonal cell against permutation")
for mask in ["plof", "mislc"]:
    out = []
    for o, ts in PRIMARY_TRAIT.items():
        obs, dz, p, n = run(mask, traits={o: ts}, locus=True)
        out.append(f"{o} {dz:+.1f} (p {p:.3f})")
    print(f"  {mask}: " + " | ".join(out))
