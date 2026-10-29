"""Conservation across species: are heterozygous knockout mice (IMPC) of recessive dysplasia genes shorter, and does the
matrix rule also hold in the mouse? (exploratory)

IMPC statistical results, zygosity heterozygote, 'Body length' (linear mixed model including body weight) and
'Bone Mineral Density (excluding skull)'. Note: IMPC mostly phenotypes heterozygotes only when homozygotes are not viable,
so the selection of genes is not random. Ortholog mapping is simplified via the upper-case symbol. z = estimate / SE,
averaged per gene over centres.
Tests: (1) AR_nos against all other genes with heterozygous mouse data (fraction p < 0.05, mean z); (2) GO category
within AR_nos and overall; (3) concordance human (Genebass pLoF z, height / eBMD) vs mouse z.

Inputs: shet_merged.tsv.gz (script 11); config.IMPC_API_DIR / impc_het_Body.tsv.gz, impc_het_Bone.tsv.gz
Output: printed only.
Usage:  python 15_mouse_heterozygous_knockouts.py
"""
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu, fisher_exact

import config

A = pd.read_csv(config.results("shet_merged.tsv.gz"), sep="\t")
A["SYM"] = A.symbol.str.upper()
for f, label, human_z in [("impc_het_Body.tsv.gz", "body length", "z_height_plof"),
                          ("impc_het_Bone.tsv.gz", "bone mineral density", "z_ebmd_plof")]:
    m = pd.read_csv(config.IMPC_API_DIR / f, sep="\t")
    m = m[(m.status == "Successful") & m.genotype_effect_stderr_estimate.gt(0)].copy()
    m["z"] = m.genotype_effect_parameter_estimate / m.genotype_effect_stderr_estimate
    m["SYM"] = m.marker_symbol.str.upper()
    g = m.groupby("SYM").agg(zm=("z", "mean"), pmin=("p_value", "min"), n=("z", "size")).reset_index()
    g["sig"] = (g.pmin < 0.05).astype(int)
    x = g.merge(A[["SYM", "set", "category", human_z, "post_mean"]], on="SYM", how="left")
    x["set"] = x.set.fillna("unknown")
    print(f"\n== Mouse heterozygous, {label}: {len(g)} genes; AR_nos {int((x.set == 'AR_nos').sum())}")
    for st in ["AR_nos", "AD_nos", "nos_other", "AR_nonskel", "AR_skelsyndr", "pool"]:
        y = x[x.set == st]
        if len(y):
            print(f"  {st:12s} n {len(y):4d}: p<0.05 {y.sig.mean():.3f}, mean z {y.zm.mean():+.2f}, fraction z<0 {(y.zm < 0).mean():.2f}")
    a, b = x[x.set == "AR_nos"], x[~x.set.isin(["AR_nos", "AD_nos", "nos_other"])]
    orr, p = fisher_exact([[a.sig.sum(), len(a) - a.sig.sum()], [b.sig.sum(), len(b) - b.sig.sum()]])
    print(f"  AR_nos vs non-nosology: OR {orr:.2f} p {p:.3f}; z smaller (MW) p {mannwhitneyu(a.zm, b.zm, alternative='less').pvalue:.3f}")
    for label2, sel in [("AR_nos", x.set == "AR_nos"), ("all", x.set.notna())]:
        y = x[sel]
        print(f"  [{label2}] by GO category: "
              + " | ".join(f"{k} p<0.05 {v.sig.mean():.2f} z {v.zm.mean():+.2f} (n {len(v)})" for k, v in y.groupby("category")))
    y = x.dropna(subset=[human_z])
    r, p = spearmanr(y[human_z], y.zm)
    yy = y[y.set == "AR_nos"]
    r2, p2 = spearmanr(yy[human_z], yy.zm) if len(yy) > 5 else (np.nan, np.nan)
    print(f"  Concordance human z vs mouse z: all rho {r:+.3f} (p {p:.1e}, n {len(y)}); AR_nos rho {r2:+.3f} (p {p2:.2f}, n {len(yy)})")
    if label == "body length":
        print("  AR_nos MATRIX genes:", yy[yy.category == "MATRIX"][["SYM", "zm", "pmin", human_z]].round(2).to_string(index=False))
