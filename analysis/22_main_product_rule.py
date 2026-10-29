"""Main-product rule: are carrier effects larger for genes that make the organ's own main product? (pre-specified,
PRAEREG_hauptprodukt.md)

Pre-specified tests W1-W3 (test register nos. 21-24). Unit = gene x organ (recessive genes with organ fraction >= 0.2 from
script 14); unit value = mean doubly centred, locus-adjusted burden z^2 over the organ's traits. Four routes define whether
a gene makes the organ's main product:
  W1   gene annotated to the organ's main-product GO terms (MAIN_PRODUCT_GO, including descendants)
  W2a  absolute expression in the organ's tissues (HPA, MAIN_PRODUCT_TISSUES) in the upper tertile within the organ
  W2b  organ specificity of expression, log2((nTPM + 1) / (median over tissues + 1)), in the upper tertile within the organ
  W3   W1 and W2a
Test: difference of the mean unit value (yes - no), labels permuted within organ (B = 10000), skeletal units excluded (the
skeletal MATRIX classes of script 10 are the reference). Variants: primary, within-organ rank of the value, and without the
top 25 genes of script 20 (pLoF).

Inputs: go-basic.obo, goa_human.gaf.gz, hgnc.txt, hpa_tissue.zip, organ_fractions_recessive.tsv.gz (14),
gene_annotation.tsv.gz (03), Genebass burden files.
Outputs: 22_main_product_units_{plof,mislc}.tsv (read by 22b, 25, 26, 27, 29, 34, 49), 22_main_product_rule.tsv
Call: python 22_main_product_rule.py
"""
import gzip
import zipfile
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.stats import norm, rankdata

import common
import config

REF = config.REF_DIR
rng = np.random.default_rng(13)
B = 10000
ORGAN_TRAITS = common.ORGAN_TRAITS
MAIN_PRODUCT_GO = {"blood": ["GO:0042541", "GO:0006783", "GO:0055072", "GO:0030218", "GO:0030220"],
                   "immune": ["GO:0002521", "GO:0002377", "GO:0042582", "GO:0042581"],
                   "liver": ["GO:0072562", "GO:0006699", "GO:0015721", "GO:0042167"],
                   "lipids": ["GO:0042157", "GO:0034358", "GO:0030301"],
                   "kidney": ["GO:0070293", "GO:0097254", "GO:0015747", "GO:0003094", "GO:0035435"],
                   "muscle": ["GO:0030017", "GO:0006936"], "lung": ["GO:0043129", "GO:0007585"], "circulatory": ["GO:0008217"],
                   "nervous": ["GO:0099504"], "glucose": ["GO:0030073"],
                   "skeletal": ["GO:0032964", "GO:0006024", "GO:0030166", "GO:0006888", "GO:0030198"]}  # reference (MATRIX of script 10)
MAIN_PRODUCT_TISSUES = {"blood": ["bone marrow"], "immune": ["bone marrow", "spleen", "lymph node"], "liver": ["liver"],
                        "lipids": ["liver", "small intestine"], "kidney": ["kidney"], "muscle": ["skeletal muscle"],
                        "lung": ["lung"], "circulatory": ["heart muscle", "blood vessel"], "nervous": ["cerebral cortex"],
                        "glucose": ["pancreas"]}
# the 25 genes with the largest matching contribution in script 20 (pLoF)
TOP25 = {"KDM5B", "TMPRSS6", "SERPINH1", "CLDN10", "SLC13A1", "COL27A1", "CHSY1", "PAPPA2", "LCAT", "TFRC", "INPPL1", "TFR2",
         "LIPC", "CAD", "DDR2", "PRKG2", "SEC24D", "XYLT1", "TET2", "ABCC2", "LNPK", "SLC22A12", "POR", "ABCG5", "ADGRG6"}
# row order of the summary table printed at the end (sorted order of the original analysis labels)
ANALYSIS_PRINT_ORDER = ["rank", "without top 25", "primary"]

# GO
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


TERMS = {o: set().union(*[descendants(t) for t in ts]) for o, ts in MAIN_PRODUCT_GO.items()}
go_symbols = defaultdict(set)
with gzip.open(REF / "goa_human.gaf.gz", "rt", encoding="utf-8") as f:
    for line in f:
        if line.startswith("!"):
            continue
        c = line.split("\t")
        if "NOT" in c[3]:
            continue
        for o, terms in TERMS.items():
            if c[4] in terms:
                go_symbols[o].add(c[2])
hgnc = pd.read_csv(REF / "hgnc.txt", sep="\t", dtype=str, usecols=["symbol", "ensembl_gene_id"]).dropna()
symbol2ensembl = dict(zip(hgnc.symbol, hgnc.ensembl_gene_id))
GO_SET = {o: {symbol2ensembl[s] for s in v if s in symbol2ensembl} for o, v in go_symbols.items()}
# HPA tissue expression
z = zipfile.ZipFile(REF / "hpa_tissue.zip")
hpa = pd.read_csv(z.open(z.namelist()[0]), sep="\t").pivot_table(index="Gene", columns="Tissue", values="nTPM", aggfunc="max")
hpa_median = hpa.median(axis=1)

organ_fractions = pd.read_csv(config.results("organ_fractions_recessive.tsv.gz"), sep="\t")
annotation = pd.read_csv(config.results("gene_annotation.tsv.gz"), sep="\t")
ALL_TRAITS = [t for v in ORGAN_TRAITS.values() for t in v]


def z_matrix(mask):
    """|burden z| per gene x trait (from the two-sided burden p), carriers >= 10 and gene passes QC."""
    out = {}
    for t in ALL_TRAITS:
        g = pd.read_csv(config.GENEBASS_DIR / f"genebass_{t}_{mask}.tsv.gz", sep="\t")[
            ["gene_id", "Pvalue_Burden"]].dropna().drop_duplicates("gene_id")
        g = g.merge(annotation[["gene_id", f"keep_{mask}", f"carriers_{mask}"]], on="gene_id")
        g = g[(g[f"keep_{mask}"] == True) & (g[f"carriers_{mask}"] >= 10)].set_index("gene_id")  # noqa: E712
        out[t] = pd.Series(norm.isf(g.Pvalue_Burden.clip(1e-300).values / 2), index=g.index)
    return pd.DataFrame(out)


SYN = z_matrix("syn")


def units(mask):
    """Gene x organ units with value and the four route indicators."""
    members = organ_fractions[organ_fractions.frac >= 0.2]
    genes = sorted(set(members.gene_id))
    Z2 = z_matrix(mask).reindex(genes) ** 2
    S = SYN.reindex(genes)
    for t in ALL_TRAITS:
        ok = Z2[t].notna() & S[t].notna()
        b = np.polyfit(S.loc[ok, t] ** 2, Z2.loc[ok, t], 1)[0] if ok.sum() > 50 else 0
        Z2[t] = Z2[t] - max(b, 0) * S[t] ** 2
        Z2.loc[S[t].abs() > 1.96, t] = np.nan
    Z2 = Z2.clip(upper=50)
    Zc = Z2 - Z2.mean(axis=0)
    Zc = Zc.sub(Zc.mean(axis=1), axis=0)
    rows = []
    for r in members.itertuples():
        v = Zc.loc[r.gene_id, ORGAN_TRAITS[r.organ]].dropna()
        if len(v):
            rows.append((r.gene_id, r.organ, v.mean()))
    U = pd.DataFrame(rows, columns=["gene_id", "organ", "value"]).merge(annotation[["gene_id", "symbol"]], on="gene_id")
    U["W1"] = [g in GO_SET.get(o, set()) for g, o in zip(U.gene_id, U.organ)]
    expr, spec = [], []
    for g, o in zip(U.gene_id, U.organ):
        if o in MAIN_PRODUCT_TISSUES and g in hpa.index:
            e = hpa.loc[g, MAIN_PRODUCT_TISSUES[o]].max()
            expr.append(np.log10(e + 1))
            spec.append(np.log2((e + 1) / (hpa_median[g] + 1)))
        else:
            expr.append(np.nan)
            spec.append(np.nan)
    U["expr"] = expr
    U["specificity"] = spec
    U["W2a"] = U.groupby("organ").expr.transform(lambda x: x >= x.quantile(2 / 3)) & U.expr.notna()
    U["W2b"] = U.groupby("organ").specificity.transform(lambda x: x >= x.quantile(2 / 3)) & U.specificity.notna()
    U["W3"] = U.W1 & U.W2a
    return U


def perm(U, col, stat="mean", B=B):
    """Difference of mean value (yes - no); labels permuted within organ (organs in order of first appearance)."""
    x = U.copy()
    if stat == "rank":
        x["value"] = x.groupby("organ").value.transform(lambda v: rankdata(v) / len(v))
    lab = x[col].values.astype(bool)
    w = x.value.values
    idx = [np.where(x.organ.values == o)[0] for o in x.organ.unique()]
    obs = w[lab].mean() - w[~lab].mean()
    d = np.empty(B)
    for i in range(B):
        L = lab.copy()
        for ix in idx:
            L[ix] = rng.permutation(L[ix])
        d[i] = w[L].mean() - w[~L].mean()
    return obs, (1 + (d >= obs).sum()) / (B + 1), int(lab.sum())


out = []
for mask in ["plof", "mislc"]:
    U = units(mask)
    U.to_csv(config.results(f"22_main_product_units_{mask}.tsv"), sep="\t", index=False)
    X = U[U.organ != "skeletal"]
    for route in ["W1", "W2a", "W2b", "W3"]:
        Xw = X if route in ("W1",) else X[X.organ.isin(MAIN_PRODUCT_TISSUES)]
        for analysis, sub, st in [("primary", Xw, "mean"), ("rank", Xw, "rank"),
                                  ("without top 25", Xw[~Xw.symbol.isin(TOP25)], "mean")]:
            obs, p, n = perm(sub, route, st)
            out.append(dict(bset=mask, route=route, analysis=analysis, diff=round(obs, 3), p=p, n_yes=n, n=len(sub)))
            print(out[-1], flush=True)
    sk = U[U.organ == "skeletal"]
    obs, p, n = perm(sk, "W1")
    print(f"  reference skeletal MATRIX ({mask}): {obs:+.3f} p {p:.4f} n_yes {n}")
    # per organ (W1, W2a, W3)
    for route in ["W1", "W2a", "W3"]:
        txt = []
        for o in common.sort_organs(X.organ.unique()):
            x = X[X.organ == o]
            if x[route].sum() >= 3 and (~x[route]).sum() >= 3:
                txt.append(f"{o} {x[x[route]].value.mean() - x[~x[route]].value.mean():+.2f} ({int(x[route].sum())}/{len(x)})")
        print(f"  per organ {route} {mask}: " + " | ".join(txt))
    print(f"  top W1 units ({mask}):", ", ".join(X[X.W1].sort_values("value", ascending=False).head(15)
                                               .apply(lambda r: f"{r.symbol}({r.organ})", axis=1)))
R = pd.DataFrame(out)
R.to_csv(config.results("22_main_product_rule.tsv"), sep="\t", index=False)
summary = R.pivot_table(index=["route", "analysis"], columns="bset", values=["diff", "p"])
summary = summary.reindex(pd.MultiIndex.from_product([sorted(R.route.unique()), ANALYSIS_PRINT_ORDER], names=["route", "analysis"]))
print("\n" + summary.to_string())
