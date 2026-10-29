"""Prediction of new skeletal dysplasia genes from the carrier effect on height (exploratory).

(A) Matrix rule outside the nosology, adjusted for linkage disequilibrium. (B) Time split: does the carrier effect in UK
Biobank data (Genebass, data freeze 2021/22) predict genes that were annotated as dysplasia genes only AFTER the 2023
nosology? (C) Candidate list of new dysplasia genes.

Protection against LD ('synthetic association' with common height variants): synonymous z of the same gene as covariate
and the filter |z_syn| < 1.96; candidates must also replicate in the missense set (same direction).
New dysplasia genes (default mode "date"): OMIM disease whose earliest HPO biocuration entry is on or after 2023-01-01,
with a term below HP:0002652 (skeletal dysplasia) or HP:0003498 (disproportionate short stature); gene not in the 2023
nosology and without an older dysplasia diagnosis.
Mode "omim": OMIM number >= 620000 (assigned from about the end of 2022) and a broader skeletal definition that also
includes short stature together with >= 5 terms below HP:0000924 (skeletal system).

Inputs: shet_merged.tsv.gz (script 11); config.GENEBASS_DIR / genebass_height_syn.tsv.gz; config.REF_DIR / hp.obo,
        phenotype.hpoa, hgnc.txt, genes_to_disease.txt
Output: 12_candidate_genes.tsv in RC_RESULTS_DIR (written in both modes; the run in "omim" mode overwrites it, as in the
        original pipeline, which runs "date" first and "omim" second; the two versions differ only in the column `new`)
Usage:  python 12_prediction_new_dysplasia_genes.py [omim]
"""
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import norm, fisher_exact
from sklearn.metrics import roc_auc_score
from statsmodels.stats.multitest import multipletests

import config

REF = config.REF_DIR
MODE = sys.argv[1] if len(sys.argv) > 1 else "date"

A = pd.read_csv(config.results("shet_merged.tsv.gz"), sep="\t")
g = pd.read_csv(config.GENEBASS_DIR / "genebass_height_syn.tsv.gz", sep="\t")[["gene_id", "BETA_Burden", "Pvalue_Burden"]].dropna()
g["z_height_syn"] = np.sign(g.BETA_Burden) * norm.isf(g.Pvalue_Burden.clip(1e-300) / 2)
A = A.merge(g[["gene_id", "z_height_syn"]], on="gene_id", how="left")
A = A.dropna(subset=["z_height_plof", "z_height_syn"]).copy()
A["sig"] = (2 * norm.sf(A.z_height_plof.abs()) < 0.05).astype(int)
A["z2syn"] = A.z_height_syn ** 2
A["ltr"] = np.log10(A.carriers_plof)  # log10 carriers
A["MATRIX"] = (A.category == "MATRIX").astype(int)
A["synclean"] = A.z_height_syn.abs() < 1.96

print("(A) Matrix rule outside the nosology (pLoF height, fraction p<0.05):")
for label, sel in [("AR outside the nosology", A.set.isin(["AR_skelsyndr", "AR_nonskel"])),
                   ("no Mendelian gene (pool)", A.set == "pool"),
                   ("all outside the nosology", ~A.set.isin(["AR_nos", "AD_nos", "nos_other"]))]:
    for clean in [False, True]:
        x = A[sel & (A.synclean if clean else True)]
        f = smf.logit("sig ~ MATRIX + ltr + oe_lof_upper + z2syn", data=x).fit(disp=0)
        print(f"  {label:26s} {'syn-clean ' if clean else 'all       '}: MATRIX {x[x.MATRIX == 1].sig.mean():.3f} (n {x.MATRIX.sum()}) vs rest "
              f"{x[x.MATRIX == 0].sig.mean():.3f} | adj. OR {np.exp(f.params['MATRIX']):.2f} p {f.pvalues['MATRIX']:.1e}")

# --- (B) new dysplasia genes
children = defaultdict(set)
cur = None
for line in open(REF / "hp.obo", encoding="utf-8"):
    if line.startswith("id: HP:"):
        cur = line.split()[1]
    elif line.startswith("is_a:") and cur:
        children[line.split()[1]].add(cur)


def descendants(root):
    """HPO term and all its descendants."""
    out, stack = {root}, [root]
    while stack:
        for c in children[stack.pop()]:
            if c not in out:
                out.add(c)
                stack.append(c)
    return out


SD = descendants("HP:0002652") | descendants("HP:0003498")
ha = pd.read_csv(REF / "phenotype.hpoa", sep="\t", comment="#", dtype=str)
ha = ha[ha.database_id.str.startswith("OMIM")]
ha["date"] = ha.biocuration.str.extract(r"\[(\d{4}-\d{2}-\d{2})\]")[0]
first = ha.groupby("database_id").date.min()
sd_dis = set(ha[(ha.aspect == "P") & ha.hpo_id.isin(SD) & (ha.qualifier != "NOT")].database_id)
inh = ha[ha.aspect == "I"].groupby("database_id").hpo_id.apply(set)
if MODE == "omim":
    # OMIM number >= 620000 (assigned from about the end of 2022); broader skeletal definition: short stature plus
    # >= 5 terms of skeletal morphology
    SK = descendants("HP:0000924")
    short_stature = set(ha[(ha.aspect == "P") & ha.hpo_id.isin(descendants("HP:0004322")) & (ha.qualifier != "NOT")].database_id)
    sk = ha[(ha.aspect == "P") & ha.hpo_id.isin(SK) & (ha.qualifier != "NOT")].groupby("database_id").size()
    sd_dis = sd_dis | (short_stature & set(sk[sk >= 5].index))
    new_dis = {d for d in sd_dis if int(d[5:]) >= 620000}
else:
    new_dis = {d for d in sd_dis if isinstance(first.get(d), str) and first[d] >= "2023-01-01"}
old_dis = sd_dis - new_dis
h = pd.read_csv(REF / "hgnc.txt", sep="\t", dtype=str, usecols=["entrez_id", "ensembl_gene_id"]).dropna()
e2e = dict(zip(h.entrez_id, h.ensembl_gene_id))
g2d = pd.read_csv(REF / "genes_to_disease.txt", sep="\t", dtype=str)
g2d = g2d[g2d.association_type == "MENDELIAN"]
g2d["gene_id"] = g2d.ncbi_gene_id.str.replace("NCBIGene:", "").map(e2e)
new_genes = set(g2d[g2d.disease_id.isin(new_dis)].gene_id.dropna())
old_sd_genes = set(g2d[g2d.disease_id.isin(old_dis)].gene_id.dropna())
nos = set(A[A.set.isin(["AR_nos", "AD_nos", "nos_other"])].gene_id)
cand_new = (new_genes - old_sd_genes) - nos
newAR = {gg for gg in cand_new
         if any("HP:0000007" in inh.get(d, set()) for d in g2d[(g2d.gene_id == gg) & g2d.disease_id.isin(new_dis)].disease_id)}
print(f"\n(B) New dysplasia genes after 2023-01-01 (not in the nosology, no older dysplasia diagnosis): {len(cand_new)}, of which AR {len(newAR)}")
U = A[~A.gene_id.isin(nos)].copy()
U["new"] = U.gene_id.isin(cand_new).astype(int)
U["newAR"] = U.gene_id.isin(newAR).astype(int)
print("   in the Genebass analysis:", U.new.sum(), "AR:", U.newAR.sum())
print("   ", ", ".join(sorted(U[U.new == 1].symbol)))
for label in ["new", "newAR"]:
    x = U
    a, b = x[x[label] == 1], x[x[label] == 0]
    orr, p = fisher_exact([[a.sig.sum(), len(a) - a.sig.sum()], [b.sig.sum(), len(b) - b.sig.sum()]])
    f = smf.logit(f"sig ~ {label} + ltr + oe_lof_upper + z2syn", data=x).fit(disp=0)
    auc = roc_auc_score(x[label], x.z_height_plof.abs())
    auc_mis = roc_auc_score(x[x.z_height_mislc.notna()][label], x[x.z_height_mislc.notna()].z_height_mislc.abs())
    auc_syn = roc_auc_score(x[label], x.z_height_syn.abs())
    print(f"   {label:5s}: p<0.05 {a.sig.mean():.3f} (n {len(a)}) vs {b.sig.mean():.3f} | OR {orr:.2f} p {p:.1e} | adj. OR {np.exp(f.params[label]):.2f} "
          f"p {f.pvalues[label]:.1e} | AUC |z_pLoF| {auc:.3f}, |z_mis| {auc_mis:.3f}, |z_syn| {auc_syn:.3f} (negative control)")
    # combination: carrier effect + matrix + s_het
    x = x.dropna(subset=["post_mean"]).copy()
    x["score"] = x.z_height_plof.abs().rank(pct=True) + x.MATRIX + x.post_mean.rank(pct=True)
    print(f"          AUC combined (|z| + MATRIX + s_het): {roc_auc_score(x[label], x.score):.3f} | s_het only {roc_auc_score(x[label], x.post_mean):.3f}")

# --- (C) candidates
U["q"] = multipletests(2 * norm.sf(U.z_height_plof.abs()), method="fdr_bh")[1]
C = U[(U.q < 0.05) & U.synclean & (np.sign(U.z_height_plof) == np.sign(U.z_height_mislc))].copy()
C = C.sort_values("q")[["symbol", "set", "category", "carriers_plof", "z_height_plof", "z_height_mislc", "z_height_syn", "post_mean",
                        "new", "hpo_n_dis"]]
C.to_csv(config.results("12_candidate_genes.tsv"), sep="\t", index=False)
print(f"\n(C) Candidates outside the nosology (pLoF FDR<0.05 genome-wide, syn clean, missense same direction): {len(C)}; "
      f"of which MATRIX {int((C.category == 'MATRIX').sum())}, already described as dysplasia gene after 2023: {int(C.new.sum())}, "
      f"without any OMIM disease: {int(C.hpo_n_dis.isna().sum())}")
pd.set_option("display.width", 220)
print(C.head(60).round(4).to_string(index=False))
