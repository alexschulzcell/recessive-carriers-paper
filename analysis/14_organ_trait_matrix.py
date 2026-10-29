"""Organ generalisation: do carriers of recessive genes show effects specifically in the traits of the disease organ?

Exploratory (test register no. 10, "organ diagonal"). Genes: purely recessive Mendelian genes (HPO: autosomal recessive, no
autosomal dominant or X-linked inheritance) plus the recessive skeletal-dysplasia genes of the nosology (AR_nos). Organ
assignment: fraction of the HPO phenotype terms (over all OMIM diseases of the gene) that lie in the HPO subtree of the organ;
a gene belongs to every organ with fraction >= FRAC (multiple assignment allowed). Traits: Genebass (UK Biobank, 394k exomes),
a priori organ-trait pairs (common.ORGAN_TRAITS): blood -> Hb, MCV, RBC, PLT; immune -> WBC, lymphocytes; kidney -> creatinine,
cystatin C, urate, urea, phosphate; liver -> ALT, AST, GGT, bilirubin, albumin; lipids -> LDL, HDL, TG, cholesterol;
glucose -> glucose, HbA1c; lung -> FEV1, FVC; muscle -> grip strength; nervous -> fluid intelligence, reaction time;
skeletal -> height, sitting height, eBMD; circulatory -> systolic blood pressure.
Test: mean of the doubly centred burden z^2 (first per trait, then per gene) in the matching organ-trait cells; the organ labels
are permuted between genes (B = 1000) for an empirical p. Also writes the organ x trait matrix (mean centred z^2 of the organ's
genes). Negative control: synonymous variants.

Inputs: HPO (hp.obo, phenotype.hpoa, genes_to_disease.txt), HGNC, gene_annotation.tsv.gz (script 03), Genebass burden files.
Outputs: organ_fractions_recessive.tsv.gz (gene x organ fraction of HPO terms; read by many later scripts),
14_organ_trait_matrix_{plof,mislc,syn}.tsv.
Call: python 14_organ_trait_matrix.py [FRAC]   (default FRAC = 0.2)
"""
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.stats import norm

import common
import config

REF = config.REF_DIR
FRAC = float(sys.argv[1]) if len(sys.argv) > 1 else 0.2
rng = np.random.default_rng(5)
ORGAN_HPO_ROOTS = common.ORGAN_HPO_ROOTS      # insertion order = organ order of the original analysis
ORGAN_TRAITS = common.ORGAN_TRAITS
TRAITS = [t for v in ORGAN_TRAITS.values() for t in v]

# HPO hierarchy: term -> direct children
children = defaultdict(set)
current = None
for line in open(REF / "hp.obo", encoding="utf-8"):
    if line.startswith("id: HP:"):
        current = line.split()[1]
    elif line.startswith("is_a:") and current:
        children[line.split()[1]].add(current)


def descendants(root):
    out, stack = {root}, [root]
    while stack:
        for c in children[stack.pop()]:
            if c not in out:
                out.add(c)
                stack.append(c)
    return out


SUBTREE = {o: set().union(*[descendants(r) for r in roots]) for o, roots in ORGAN_HPO_ROOTS.items()}
hpoa = pd.read_csv(REF / "phenotype.hpoa", sep="\t", comment="#", dtype=str)
hpoa = hpoa[hpoa.database_id.str.startswith("OMIM") & (hpoa.aspect == "P") & (hpoa.qualifier != "NOT")]
hgnc = pd.read_csv(REF / "hgnc.txt", sep="\t", dtype=str, usecols=["entrez_id", "ensembl_gene_id"]).dropna()
entrez2ensembl = dict(zip(hgnc.entrez_id, hgnc.ensembl_gene_id))
g2d = pd.read_csv(REF / "genes_to_disease.txt", sep="\t", dtype=str)
g2d = g2d[g2d.association_type == "MENDELIAN"]
g2d["gene_id"] = g2d.ncbi_gene_id.str.replace("NCBIGene:", "").map(entrez2ensembl)
disease_terms = hpoa.groupby("database_id").hpo_id.apply(set)
annotation = pd.read_csv(config.results("gene_annotation.tsv.gz"), sep="\t")
recessive = annotation[annotation.protein_coding
                       & (((annotation.hpo_ar == 1) & (annotation.hpo_ad != 1) & (annotation.hpo_xl != 1))
                          | (annotation.set == "AR_nos"))].copy()
gene_terms = (g2d.dropna(subset=["gene_id"]).groupby("gene_id").disease_id
              .apply(lambda ds: set().union(*[disease_terms.get(d, set()) for d in ds])))
rows = []
for gid in recessive.gene_id:
    t = gene_terms.get(gid, set())
    if len(t) < 3:
        continue
    for o, sub in SUBTREE.items():
        rows.append((gid, o, len(t & sub) / len(t)))
organ_fractions = pd.DataFrame(rows, columns=["gene_id", "organ", "frac"])
organ_fractions.to_csv(config.results("organ_fractions_recessive.tsv.gz"), sep="\t", index=False)
members = organ_fractions[organ_fractions.frac >= FRAC]
print(f"Recessive genes with an organ profile: {organ_fractions.gene_id.nunique()}; assignment (frac >= {FRAC}):",
      members.organ.value_counts().to_dict())


def z2(trait, mask):
    """Burden z^2 per gene (from the two-sided burden p), carriers >= 10 and gene passes QC."""
    g = pd.read_csv(config.GENEBASS_DIR / f"genebass_{trait}_{mask}.tsv.gz", sep="\t")[
        ["gene_id", "BETA_Burden", "Pvalue_Burden"]].dropna()
    g = g.merge(annotation[["gene_id", f"keep_{mask}", f"carriers_{mask}"]], on="gene_id")
    g = g[(g[f"keep_{mask}"] == True) & (g[f"carriers_{mask}"] >= 10)]  # noqa: E712
    g = g.drop_duplicates("gene_id").set_index("gene_id")
    return pd.Series(norm.isf(g.Pvalue_Burden.clip(1e-300).values / 2) ** 2, index=g.index)


genes = sorted(set(members.gene_id))
M = {o: set(members[members.organ == o].gene_id) for o in ORGAN_HPO_ROOTS}
for mask in ["plof", "mislc", "syn"]:
    Z = pd.DataFrame({t: z2(t, mask) for t in TRAITS
                      if (config.GENEBASS_DIR / f"genebass_{t}_{mask}.tsv.gz").exists()}).reindex(genes)
    Z = Z.clip(upper=50)                          # cap single huge z (e.g. Hb for HBB)
    Zc = Z - Z.mean(axis=0)                       # trait effect
    Zc = Zc.sub(Zc.mean(axis=1), axis=0)          # gene effect (like fixed effects in the balanced case)
    trait2organ = {t: o for o, ts in ORGAN_TRAITS.items() for t in ts if t in Z.columns}

    def stat(Mx):
        vals = [Zc.loc[list(Mx[o] & set(Z.index)), t].dropna().values for t, o in trait2organ.items()]
        v = np.concatenate([x for x in vals if len(x)])
        return v.mean()

    obs = stat(M)
    gene_array = np.array(genes)
    perm = []
    labels = {g: [o for o in ORGAN_HPO_ROOTS if g in M[o]] for g in genes}
    for b in range(1000):
        shuffled = rng.permutation(gene_array)
        Mp = defaultdict(set)
        for g_new, g_old in zip(gene_array, shuffled):
            for o in labels[g_old]:
                Mp[o].add(g_new)
        perm.append(stat(Mp))
    perm = np.array(perm)
    print(f"\n[{mask}] matching effect (mean centred z2 in matching cells) {obs:.3f}; permutation {perm.mean():.3f} +- {perm.std():.3f}; "
          f"dz {(obs - perm.mean()) / perm.std():.1f}; p {(1 + (perm >= obs).sum()) / 1001:.4f}")
    mat = pd.DataFrame({t: {o: Zc.loc[list(M[o] & set(Z.index)), t].mean() for o in ORGAN_HPO_ROOTS} for t in Z.columns})
    nrow = {o: len(M[o] & set(Z.dropna(how='all').index)) for o in ORGAN_HPO_ROOTS}
    mat.index = [f"{o} ({nrow[o]})" for o in mat.index]
    mat.to_csv(config.results(f"14_organ_trait_matrix_{mask}.tsv"), sep="\t", float_format="%.3f")
    if mask != "syn":
        pd.set_option("display.width", 260)
        print(mat.round(2).to_string())
    # per organ: diagonal cells against the rest of the organ's row
    out = []
    for o, ts in ORGAN_TRAITS.items():
        ts = [t for t in ts if t in Z.columns]
        if not ts:
            continue
        x = Zc.loc[list(M[o] & set(Z.index))]
        diag = x[ts].stack().mean()
        off = x.drop(columns=ts).stack().mean()
        out.append(f"{o}: {diag:+.2f}/{off:+.2f}")
    print("  matching/non-matching per organ: " + " | ".join(out))
