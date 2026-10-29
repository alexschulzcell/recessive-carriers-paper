"""Corrected skeletal confirmation for script 21 (PRAEREG_mechanismus.md, marked deviation).

Re-implementation of an exploratory analysis whose original output was only a log.
Deviation from the protocol (marked): in script 21 the skeletal confirmation was centred over the skeleton alone, which sets
every gene value to 0 (artefact). Here the gene values are computed as in 21, but centred over ALL organs (including the
skeleton); the genes whose main organ is skeletal (HPO-based, frac >= 0.2) are tested. In addition SECRETED - CATABOLIC over
all organs including the skeleton. Permutation as in 21 (within main organ, B = 10000, seed 11, one random stream).
GO classes, gene_values() and permtest() are unchanged from script 21.

Inputs: go-basic.obo, goa_human.gaf.gz, hgnc.txt, organ_fractions_recessive.tsv.gz (14), gene_annotation.tsv.gz (03),
Genebass burden files. Output: log only.
Call: python 21b_skeletal_mechanism_corrected.py
"""
import gzip
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.stats import norm

import common
import config

REF = config.REF_DIR
rng = np.random.default_rng(11)
ORGAN_TRAITS = common.ORGAN_TRAITS
# sorted order of the class labels in the original analysis (keeps the printed order identical)
CLASS_ORDER = ["CATABOLIC", "SECRETED", "CILIUM_ORGANS", "other"]

# --- GO classes (as in script 21)
go_children = defaultdict(set)
current = None
for line in open(REF / "go-basic.obo", encoding="utf-8"):
    if line.startswith("id: GO:"):
        current = line.split()[1]
    elif line.startswith(("is_a: GO:", "relationship: part_of GO:")) and current:
        go_children["GO:" + line.split("GO:")[1][:7]].add(current)


def descendants(root):
    out, stack = {root}, [root]
    while stack:
        for c in go_children[stack.pop()]:
            if c not in out:
                out.add(c)
                stack.append(c)
    return out


MECHANISM_CLASSES = {
    "SECRETED": set().union(*[descendants(x) for x in ["GO:0005576", "GO:0005198", "GO:0006888", "GO:0030198", "GO:0032964",
                                                       "GO:0006024", "GO:0030166"]]),
    "CATABOLIC": descendants("GO:0005764") | descendants("GO:0009056"),
    "CILIUM_ORGANS": descendants("GO:0005929"),
}
class_symbols = defaultdict(set)
with gzip.open(REF / "goa_human.gaf.gz", "rt", encoding="utf-8") as f:
    for line in f:
        if line.startswith("!"):
            continue
        c = line.split("\t")
        if "NOT" in c[3]:
            continue
        for k, terms in MECHANISM_CLASSES.items():
            if c[4] in terms:
                class_symbols[k].add(c[2])
hgnc = pd.read_csv(REF / "hgnc.txt", sep="\t", dtype=str, usecols=["symbol", "ensembl_gene_id", "enzyme_id"])
symbol2ensembl = dict(zip(hgnc.dropna(subset=["ensembl_gene_id"]).symbol, hgnc.dropna(subset=["ensembl_gene_id"]).ensembl_gene_id))
gene_class = {}
for k in ["SECRETED", "CATABOLIC", "CILIUM_ORGANS"]:  # precedence
    for sym in class_symbols[k]:
        e = symbol2ensembl.get(sym)
        if e and e not in gene_class:
            gene_class[e] = k

organ_fractions = pd.read_csv(config.results("organ_fractions_recessive.tsv.gz"), sep="\t")
annotation = pd.read_csv(config.results("gene_annotation.tsv.gz"), sep="\t")


def z_matrix(mask, traits):
    """|burden z| per gene x trait (from the two-sided burden p), carriers >= 10 and gene passes QC."""
    out = {}
    for t in traits:
        g = pd.read_csv(config.GENEBASS_DIR / f"genebass_{t}_{mask}.tsv.gz", sep="\t")[
            ["gene_id", "BETA_Burden", "Pvalue_Burden"]].dropna().drop_duplicates("gene_id")
        g = g.merge(annotation[["gene_id", f"keep_{mask}", f"carriers_{mask}"]], on="gene_id")
        g = g[(g[f"keep_{mask}"] == True) & (g[f"carriers_{mask}"] >= 10)].set_index("gene_id")  # noqa: E712
        out[t] = pd.Series(norm.isf(g.Pvalue_Burden.clip(1e-300).values / 2), index=g.index)
    return pd.DataFrame(out)


def gene_values(mask, organs):
    traits = [t for o in organs for t in ORGAN_TRAITS[o]]
    members = organ_fractions[(organ_fractions.frac >= 0.2) & organ_fractions.organ.isin(organs)]
    genes = sorted(set(members.gene_id))
    Z = z_matrix(mask, traits).reindex(genes)
    S = z_matrix("syn", traits).reindex(genes)
    Z2 = Z ** 2
    for t in traits:  # locus adjustment (V3 of script 20)
        ok = Z2[t].notna() & S[t].notna()
        b = np.polyfit(S.loc[ok, t] ** 2, Z2.loc[ok, t], 1)[0] if ok.sum() > 50 else 0
        Z2[t] = Z2[t] - max(b, 0) * S[t] ** 2
        Z2.loc[S[t].abs() > 1.96, t] = np.nan
    Z2 = Z2.clip(upper=50)
    Zc = Z2 - Z2.mean(axis=0)
    Zc = Zc.sub(Zc.mean(axis=1), axis=0)
    rows = []
    for g, sub in members.groupby("gene_id"):
        ts = [t for o in sub.organ for t in ORGAN_TRAITS[o]]
        v = Zc.loc[g, ts].dropna()
        if len(v):
            rows.append((g, v.mean(), sub.sort_values("frac").organ.iloc[-1]))
    G = pd.DataFrame(rows, columns=["gene_id", "value", "main_organ"]).merge(
        annotation[["gene_id", "symbol", f"carriers_{mask}", "oe_lof_upper"]], on="gene_id")
    G["gene_class"] = G.gene_id.map(gene_class).fillna("other")
    return G


def permtest(G, a, b, B=10000):
    """Difference of mean gene value between classes a and b; class labels permuted within the main organ."""
    x = G[G.gene_class.isin([a, b])].copy()
    obs = x[x.gene_class == a].value.mean() - x[x.gene_class == b].value.mean()
    lab = x.gene_class.values
    org = x.main_organ.values
    w = x.value.values
    idx = [np.where(org == o)[0] for o in common.sort_organs(np.unique(org))]  # fixed organ order (reproducible stream)
    d = np.empty(B)
    for i in range(B):
        L = lab.copy()
        for ix in idx:
            L[ix] = rng.permutation(L[ix])
        d[i] = w[L == a].mean() - w[L == b].mean()
    return obs, (1 + (d >= obs).sum()) / (B + 1), int((lab == a).sum()), int((lab == b).sum())


def in_class_order(df):
    return df.loc[[c for c in CLASS_ORDER if c in df.index]]


print("Deviation from the protocol: skeletal confirmation with centring over all organs (before: artefact = 0)")
for mask in ["plof", "mislc"]:
    G = gene_values(mask, list(ORGAN_TRAITS))
    sk = G[G.main_organ == "skeletal"]
    print(f"{mask} skeletal main organ: "
          f"{in_class_order(sk.groupby('gene_class').value.agg(['mean', 'median', 'size'])).round(3).to_dict('index')}")
    for a, b in [("SECRETED", "CATABOLIC"), ("SECRETED", "CILIUM_ORGANS")]:
        obs, p, na, nb = permtest(sk, a, b)
        print(f"  {a}-{b}: {obs:+.3f} (n {na}/{nb}) p {p:.4f}")
    obs, p, na, nb = permtest(G, "SECRETED", "CATABOLIC")
    print(f"  all organs including skeletal SECRETED-CATABOLIC: {obs:+.3f} p {p:.4f}", flush=True)
