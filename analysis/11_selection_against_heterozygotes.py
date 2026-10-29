"""Selection against heterozygotes: does the carrier effect on height explain the selection against heterozygous loss of
function in recessive dysplasia genes? (exploratory)

s_het = selection coefficient against heterozygous LoF (GeneBayes, Zeng ... Pritchard, Nat Genet 2024; Zenodo 10403680).
Additional measure without gene features: LoF observed/expected from the same source.
Circularity check: s_het is based on LoF frequency, and LoF burden power falls with s_het (fewer carriers), so an
association 'effect <-> high s_het' in the pLoF set would run against the power direction (conservative). Cleaner: the
missense|LC effect (other variants).
Tests: (1) s_het of AR_nos genes with vs without a height effect (FDR < 0.05 within AR_nos, pLoF or missense);
(2) Spearman |z_missense| ~ s_het in AR_nos vs AR_nonskel vs pool (organ-specific?);
(3) s_het by GO category (MATRIX / CILIUM / LYSO) in AR_nos;
(4) regression log10 s_het ~ z2_missense_height + log10 cds_length, per gene set.

Inputs: gene_annotation.tsv.gz (script 03), go_classes.tsv (script 10), config.GENEBAYES_FILE,
        config.GENEBASS_DIR / genebass_<trait>_<set>.tsv.gz
Output: shet_merged.tsv.gz in RC_RESULTS_DIR (annotation + GO category + s_het + Genebass z for height, eBMD and BMI;
        used by later scripts)
Usage:  python 11_selection_against_heterozygotes.py
"""
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import norm, mannwhitneyu, spearmanr
from statsmodels.stats.multitest import multipletests

import config

A = pd.read_csv(config.results("gene_annotation.tsv.gz"), sep="\t")
A = A[A.protein_coding].merge(pd.read_csv(config.results("go_classes.tsv"), sep="\t")[["gene_id", "category"]], on="gene_id", how="left")
S = pd.read_csv(config.GENEBAYES_FILE, sep="\t").rename(columns={"ensg": "gene_id"})
S["oe"] = (S.obs_lof + 0.5) / (S.exp_lof + 0.5)
A = A.merge(S[["gene_id", "post_mean", "oe", "exp_lof"]], on="gene_id", how="inner")
A["lshet"] = np.log10(A.post_mean)
for ph in ["height", "ebmd", "bmi"]:
    for s in ["plof", "mislc"]:
        g = pd.read_csv(config.GENEBASS_DIR / f"genebass_{ph}_{s}.tsv.gz", sep="\t")[["gene_id", "BETA_Burden", "Pvalue_Burden"]]
        g = g.merge(A[["gene_id", f"keep_{s}", f"carriers_{s}"]], on="gene_id")
        g = g[(g[f"keep_{s}"] == True) & (g[f"carriers_{s}"] >= 10) & g.Pvalue_Burden.notna()]  # noqa: E712
        g[f"z_{ph}_{s}"] = np.sign(g.BETA_Burden) * norm.isf(g.Pvalue_Burden.clip(1e-300) / 2)
        A = A.merge(g[["gene_id", f"z_{ph}_{s}"]], on="gene_id", how="left")
A.to_csv(config.results("shet_merged.tsv.gz"), sep="\t", index=False)

print("(1) s_het (median) of AR_nos genes with/without a height effect (FDR<0.05 within AR_nos):")
for s in ["plof", "mislc"]:
    x = A[(A.set == "AR_nos") & A[f"z_height_{s}"].notna()].copy()
    x["q"] = multipletests(2 * norm.sf(x[f"z_height_{s}"].abs()), method="fdr_bh")[1]
    hit, no = x[x.q < 0.05], x[x.q >= 0.05]
    p = mannwhitneyu(hit.post_mean, no.post_mean, alternative="greater").pvalue
    p2 = mannwhitneyu(hit.oe, no.oe, alternative="less").pvalue
    print(f"  {s:5s}: hits n={len(hit)} s_het {hit.post_mean.median():.4f} (obs/exp {hit.oe.median():.2f}) | rest n={len(no)} s_het {no.post_mean.median():.4f} "
          f"(obs/exp {no.oe.median():.2f}) | MW p {p:.1e}, obs/exp p {p2:.1e}")

print("\n(2) Spearman |z| (missense|LC) against s_het, per gene set and trait:")
for ph in ["height", "ebmd", "bmi"]:
    row = []
    for st in ["AR_nos", "AD_nos", "AR_nonskel", "AR_skelsyndr", "pool"]:
        x = A[(A.set == st) & A[f"z_{ph}_mislc"].notna()]
        r, p = spearmanr(x[f"z_{ph}_mislc"].abs(), x.post_mean)
        row.append(f"{st} {r:+.2f} (p {p:.0e}, n {len(x)})")
    print(f"  {ph:7s}: " + " | ".join(row))

print("\n(3) s_het by GO category in AR_nos (median, n):")
x = A[A.set == "AR_nos"]
print("  " + " | ".join(f"{k} {v.post_mean.median():.4f} (obs/exp {v.oe.median():.2f}, n {len(v)})" for k, v in x.groupby("category")))
m, o = x[x.category == "MATRIX"].post_mean, x[x.category.isin(["CILIUM", "LYSO"])].post_mean
print(f"  MATRIX against CILIUM+LYSO: MW p {mannwhitneyu(m, o, alternative='greater').pvalue:.1e}")

print("\n(4) Regression log10 s_het ~ z2 height (missense) + log10 CDS:")
A["z2h"] = A.z_height_mislc ** 2
A["lcds"] = np.log10(A.cds_length)
for st in ["AR_nos", "AR_nonskel", "pool"]:
    x = A[(A.set == st)].dropna(subset=["z2h", "lcds", "lshet"])
    f = smf.ols("lshet ~ z2h + lcds", data=x).fit()
    print(f"  {st:10s}: beta z2 {f.params['z2h']:+.4f} (p {f.pvalues['z2h']:.1e}), n {len(x)}")
x = A[A.set == "AR_nos"].dropna(subset=["z2h", "lcds", "lshet", "category"])
f = smf.ols("lshet ~ z2h + lcds + C(category, Treatment('other'))", data=x).fit()
print(f.summary().tables[1])
