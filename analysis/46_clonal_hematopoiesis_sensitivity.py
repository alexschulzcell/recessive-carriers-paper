"""Sensitivity to clonal haematopoiesis (CHIP): do the organ diagonal and the gene-specific direction hold without CHIP genes?

Exploratory, not pre-specified. CHIP genes: DNMT3A, ASXL1, JAK2, PPM1D, TET2 (as script 20b). Of these, only TET2 is in the
recessive gene set (organ_fractions_recessive). Main version (as 20b): the gene is removed entirely (all organ rows; it drops
out of centring and permutation). Further versions: (i) only the blood and immune rows/cells of the CHIP genes removed;
(ii) extended list: in addition ATM and MYD88 (the only further known clonal-haematopoiesis driver genes with an organ
fraction > 0 in the recessive table; for orientation only).
(a) Organ diagonal z^2 as in script 20: V0 all traits, V3 locus, V1+V3 primary trait + locus; B = 1000, seed 7, order as 20b.
    'before' = same random stream with the full gene set (V0 equals the log of 20; V3 and V1+V3 use a different random stream
    than 20 and can therefore differ by about +-0.1).
(b) T10-7 gene-specific direction as script 31: statistic mean of S * Zc, permutation of the rows of S (direction maps between
    genes), B = 1000, seed 5, order pLoF, missense, synonymous ('before' equals the log of 31). The direction map is rebuilt
    here for each gene set (same construction as direction_map.build_direction_map).

Inputs: HPO (hp.obo, phenotype.hpoa, genes_to_disease.txt), organ_fractions_recessive.tsv.gz (14), gene_annotation.tsv.gz (03),
Genebass burden files. Output: 46_clonal_hematopoiesis_sensitivity.tsv
Call: python 46_clonal_hematopoiesis_sensitivity.py
"""
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.stats import norm

import common
import config
import direction_map

REF = config.REF_DIR
B = 1000
CHIP = ["DNMT3A", "ASXL1", "JAK2", "PPM1D", "TET2"]
CHIP_EXTENDED = CHIP + ["ATM", "MYD88"]
BLOOD_IMMUNE = ["blood", "immune"]
ORGAN_TRAITS = common.ORGAN_TRAITS
PRIMARY_TRAIT = {"blood": ["hb"], "immune": ["wbc"], "kidney": ["cysc"], "liver": ["alt"], "lipids": ["ldl"], "glucose": ["hba1c"],
                 "lung": ["fev1"], "muscle": ["grip"], "nervous": ["fluidint"], "skeletal": ["height"], "circulatory": ["sbp"]}
BLOOD_IMMUNE_TRAITS = [t for o in BLOOD_IMMUNE for t in ORGAN_TRAITS[o]]
ALL_TRAITS = [t for v in ORGAN_TRAITS.values() for t in v]
OF0 = pd.read_csv(config.results("organ_fractions_recessive.tsv.gz"), sep="\t")
annotation = pd.read_csv(config.results("gene_annotation.tsv.gz"), sep="\t")
id2sym = dict(zip(annotation.gene_id, annotation.symbol))


def ids(symbols):
    return set(annotation[annotation.symbol.isin(symbols)].gene_id)


assigned = {g: common.sort_organs(OF0[(OF0.gene_id == g) & (OF0.frac >= 0.2)].organ) for g in ids(CHIP_EXTENDED) & set(OF0.gene_id)}
print("CHIP genes in the recessive organ table (organs with frac >= 0.2):",
      {id2sym[g]: o for g, o in sorted(assigned.items(), key=lambda x: id2sym[x[0]])})
print("not in the recessive gene set (enter no analysis):", sorted(set(CHIP_EXTENDED) - {id2sym[g] for g in assigned}))
BEFORE = "before (all genes)"
WITHOUT = "without CHIP, entirely (as 20b)"
WITHOUT_BLOOD_IMMUNE = "without CHIP, blood/immune only"
WITHOUT_EXTENDED = "without CHIP extended (+ATM, MYD88), entirely"
TABLES = {BEFORE: OF0,
          WITHOUT: OF0[~OF0.gene_id.isin(ids(CHIP))],
          WITHOUT_BLOOD_IMMUNE: OF0[~(OF0.gene_id.isin(ids(CHIP)) & OF0.organ.isin(BLOOD_IMMUNE))],
          WITHOUT_EXTENDED: OF0[~OF0.gene_id.isin(ids(CHIP_EXTENDED))]}
for k, t in TABLES.items():
    print(f"  {k}: removed gene-organ rows {len(OF0) - len(t)}")
rows = []

# ======================= (a) organ diagonal (as script 20) =======================


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


def diagonal(OF, rng, mask, frac=0.2, traits=ORGAN_TRAITS, locus=False):
    """As run() in script 20 (without drop_top/drop_organ)."""
    members = OF[OF.frac >= frac]
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


print("\n=== (a) organ diagonal z2 (as 20; dz and p of the organ-label permutation) ===", flush=True)
for k, OF in TABLES.items():
    rng = np.random.default_rng(7)
    for name, kw in [("V0 all traits", {}), ("V3 locus", {"locus": True}),
                     ("V1+V3 primary, locus", {"traits": PRIMARY_TRAIT, "locus": True})]:
        out = []
        for mask in ["plof", "mislc"]:
            obs, dz, p, n = diagonal(OF, rng, mask, **kw)
            out.append(f"{mask} dz {dz:5.1f} (p {p:.3f}; n {n})")
            rows.append(dict(part="a diagonal", analysis=name, gene_set=k, set=mask, stat=obs, dz=dz, p=p, n=n))
        print(f"  {k:46s} {name:22s} " + " | ".join(out), flush=True)

# ======================= (b) T10-7 gene-specific direction (as script 31) =======================
TRAITS = direction_map.TRAITS
children = direction_map.hpo_children(REF)
expected = defaultdict(lambda: defaultdict(list))
for term, trait_signs in direction_map.HPO_DIRECTION_MAP.items():
    for d in direction_map.descendants(term, children):
        for t, s in trait_signs.items():
            expected[d][t].append(s)
hpoa = pd.read_csv(REF / "phenotype.hpoa", sep="\t", comment="#", dtype=str)
hpoa = hpoa[hpoa.database_id.str.startswith("OMIM") & (hpoa.aspect == "P") & (hpoa.qualifier != "NOT")]
disease_terms = hpoa.groupby("database_id").hpo_id.apply(set)
g2d = pd.read_csv(REF / "genes_to_disease.txt", sep="\t", dtype=str)
g2d = g2d[g2d.association_type == "MENDELIAN"]
sym2id = dict(zip(annotation.symbol, annotation.gene_id))
gene_terms = g2d.groupby("gene_symbol").disease_id.apply(lambda ds: set().union(*[disease_terms.get(d, set()) for d in ds]))


def directions(recessive_genes):
    """{gene_id: {trait: +1/-1}} for the given recessive gene set (contradicting cells discarded)."""
    R = {}
    for sym, terms in gene_terms.items():
        gid = sym2id.get(sym)
        if gid is None or gid not in recessive_genes:
            continue
        acc = defaultdict(list)
        for t in terms:
            for trait, signs in expected.get(t, {}).items():
                acc[trait].extend(signs)
        r = {trait: signs[0] for trait, signs in acc.items() if len(set(signs)) == 1}
        if r:
            R[gid] = r
    return R


ALL_RECESSIVE = sorted(set(OF0.gene_id))
ZRAW = {mask: pd.DataFrame({t: z for t in TRAITS
                            if (z := direction_map.signed_ukb_z(annotation, t, mask, ALL_RECESSIVE)) is not None})
        for mask in ["plof", "mislc", "syn"]}


def t10_7(recessive_genes, rng, mask, drop_blood_immune=()):
    """T10-7 as script 31; drop_blood_immune: genes whose blood/immune cells (hb, mcv, plt, wbc, lymph) leave the direction map."""
    R = directions(recessive_genes)
    genes = sorted(R)
    S = pd.DataFrame(0.0, index=genes, columns=TRAITS)
    for g, r in R.items():
        for t, v in r.items():
            if not (g in drop_blood_immune and t in BLOOD_IMMUNE_TRAITS):
                S.loc[g, t] = v
    all_genes = sorted(recessive_genes)  # centring over all recessive genes (as 31)
    Z = ZRAW[mask].reindex(all_genes)
    Zc = Z - Z.mean(axis=0)
    Zc = Zc.sub(Zc.mean(axis=1), axis=0)
    Zc = Zc.reindex(columns=TRAITS).loc[genes].values
    Sm = S.values

    def stat(Sx):
        m = np.isfinite(Zc) & (Sx != 0)
        return np.nansum(np.where(m, Sx * Zc, 0)) / m.sum(), int(m.sum())

    obs, n = stat(Sm)
    pm = np.array([stat(Sm[rng.permutation(len(genes))])[0] for _ in range(B)])
    return obs, (obs - pm.mean()) / pm.std(), (1 + (pm >= obs).sum()) / (B + 1), n, len(genes)


print("\n=== (b) T10-7 gene-specific direction (as 31; permutation of the direction maps between genes) ===", flush=True)
VERSIONS_B = [(BEFORE, set(OF0.gene_id), set()),
              (WITHOUT, set(TABLES[WITHOUT].gene_id), set()),
              (WITHOUT_BLOOD_IMMUNE, set(OF0.gene_id), ids(CHIP)),
              (WITHOUT_EXTENDED, set(TABLES[WITHOUT_EXTENDED].gene_id), set())]
R0 = directions(set(OF0.gene_id))
print("  CHIP genes with a gene-specific direction:",
      {id2sym[g]: R0[g] for g in sorted(ids(CHIP_EXTENDED) & set(R0), key=lambda x: id2sym[x])})
for k, recessive_genes, dropped in VERSIONS_B:
    rng = np.random.default_rng(5)
    for mask in ["plof", "mislc", "syn"]:
        obs, dz, p, n, ng = t10_7(recessive_genes, rng, mask, dropped)
        rows.append(dict(part="b T10-7", analysis="main", gene_set=k, set=mask, stat=obs, dz=dz, p=p, n=n))
        print(f"  {k:46s} {mask:6s} d {obs:+.3f} dz {dz:+5.1f}  p {p:.4f}  n {n} (genes {ng})", flush=True)

R = pd.DataFrame(rows)
R.to_csv(config.results("46_clonal_hematopoiesis_sensitivity.tsv"), sep="\t", index=False, float_format="%.4f")
print("\n=== before -> after (main: CHIP genes removed entirely, as 20b) ===")
for (part, analysis, mask), x in R.groupby(["part", "analysis", "set"], sort=False):
    a = x[x.gene_set == BEFORE].iloc[0]
    b = x[x.gene_set == WITHOUT].iloc[0]
    print(f"  {part:12s} {analysis:22s} {mask:6s} dz {a.dz:+5.1f} -> {b.dz:+5.1f} | p {a.p:.4f} -> {b.p:.4f}")
