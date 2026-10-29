"""Matrix-flux hypothesis: GO classes of genes and the height carrier signal (exploratory; the hypothesis was derived from
scripts 05 and 09 and is therefore tested here only semi-independently).

Hypothesis: heterozygous effects on height occur in genes of anabolic matrix production and secretion (collagen
biosynthesis, glycosaminoglycan / proteoglycan synthesis, ER-to-Golgi export, extracellular matrix organisation), not in
catabolic (lysosomal), ciliary or osteoclast biology. Classes are defined a priori from GO (no manual gene selection):
  MATRIX = GO:0032964 collagen biosynthetic process | GO:0006024 glycosaminoglycan biosynthetic process | GO:0030166
           proteoglycan biosynthetic process | GO:0006888 ER to Golgi vesicle-mediated transport | GO:0030198 extracellular
           matrix organization (each with descendants, biological process)
  LYSO   = GO:0005764 lysosome (+ descendants, cellular component);  CILIUM = GO:0005929 cilium (+ descendants)
Test 1 (discovery, circular): within AR_nos. Test 2 (semi-independent, other genes): recessive Mendelian genes outside the
nosology (AR_skelsyndr + AR_nonskel) and pool -> does MATRIX show more height signal there as well? Synonymous variants
serve as negative control.

Inputs: gene_annotation.tsv.gz (script 03); config.REF_DIR / goa_human.gaf.gz, go-basic.obo, hgnc.txt;
        config.GENEBASS_DIR / genebass_<trait>_<set>.tsv.gz
Output: go_classes.tsv (gene_id, symbol, MATRIX, LYSO, CILIUM, category) in RC_RESULTS_DIR (used by later scripts,
        directly or through shet_merged.tsv.gz from script 11). The class definitions GO_CLASSES are read by script 47 (Table S3A).
Usage:  python 10_matrix_flux_go_classes.py
"""
import gzip
from collections import defaultdict

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

import config

# GO graph: term -> direct children (is_a and part_of)
children = defaultdict(set)
cur = None
for line in open(config.REF_DIR / "go-basic.obo", encoding="utf-8"):
    if line.startswith("id: GO:"):
        cur = line.split()[1]
    elif line.startswith(("is_a: GO:", "relationship: part_of GO:")) and cur:
        children["GO:" + line.split("GO:")[1][:7]].add(cur)


def descendants(root):
    """GO term and all its descendants."""
    out, stack = {root}, [root]
    while stack:
        for c in children[stack.pop()]:
            if c not in out:
                out.add(c)
                stack.append(c)
    return out


GO_CLASSES = {"MATRIX": set().union(*[descendants(x) for x in ["GO:0032964", "GO:0006024", "GO:0030166", "GO:0006888", "GO:0030198"]]),
              "LYSO": descendants("GO:0005764"), "CILIUM": descendants("GO:0005929")}
symbols = defaultdict(set)
with gzip.open(config.REF_DIR / "goa_human.gaf.gz", "rt", encoding="utf-8") as f:
    for line in f:
        if line.startswith("!"):
            continue
        c = line.split("\t")
        if "NOT" in c[3]:  # qualifier column
            continue
        for k, terms in GO_CLASSES.items():
            if c[4] in terms:
                symbols[k].add(c[2])
h = pd.read_csv(config.REF_DIR / "hgnc.txt", sep="\t", dtype=str, usecols=["symbol", "ensembl_gene_id"]).dropna()
s2e = dict(zip(h.symbol, h.ensembl_gene_id))
A = pd.read_csv(config.results("gene_annotation.tsv.gz"), sep="\t")
A = A[A.protein_coding & A.oe_lof_upper.notna()].copy()
for k in GO_CLASSES:
    A[k] = A.gene_id.isin({s2e[s] for s in symbols[k] if s in s2e}).astype(int)
# one category per gene, priority MATRIX > LYSO > CILIUM
A["category"] = np.select([A.MATRIX == 1, A.LYSO == 1, A.CILIUM == 1], ["MATRIX", "LYSO", "CILIUM"], "other")
A[["gene_id", "symbol", "MATRIX", "LYSO", "CILIUM", "category"]].to_csv(config.results("go_classes.tsv"), sep="\t", index=False)
print("Categories in AR_nos:", A[A.set == "AR_nos"].category.value_counts().to_dict())
print("MATRIX genes in AR_nos:", ", ".join(sorted(A[(A.set == "AR_nos") & (A.MATRIX == 1)].symbol)))

for ph in ["height", "sitheight", "ebmd"]:
    for s in ["plof", "mislc", "syn"]:
        g = pd.read_csv(config.GENEBASS_DIR / f"genebass_{ph}_{s}.tsv.gz", sep="\t")[["gene_id", "Pvalue_Burden"]].merge(A, on="gene_id")
        g = g[(g[f"keep_{s}"] == True) & (g[f"carriers_{s}"] >= 10) & g.Pvalue_Burden.notna()].copy()  # noqa: E712
        g["sig"] = (g.Pvalue_Burden < 0.05).astype(int)
        g["ltr"] = np.log10(g[f"carriers_{s}"])  # log10 carriers
        out = []
        for label, sel in [("AR_nos", g.set == "AR_nos"), ("AR_other", g.set.isin(["AR_skelsyndr", "AR_nonskel"])),
                           ("pool", g.set == "pool")]:
            x = g[sel]
            fr = x.groupby("category").sig.agg(["mean", "size"])
            txt = " ".join(f"{k} {r['mean']:.2f}({int(r['size'])})" for k, r in fr.iterrows())
            try:
                m = smf.logit("sig ~ MATRIX + ltr + oe_lof_upper", data=x).fit(disp=0)
                txt += f" | OR MATRIX adj {np.exp(m.params['MATRIX']):.2f} p {m.pvalues['MATRIX']:.1e}"
            except Exception as ex:  # noqa: BLE001
                txt += f" | {ex}"
            out.append(f"{label}: {txt}")
        print(f"\n{ph} {s}\n  " + "\n  ".join(out))
