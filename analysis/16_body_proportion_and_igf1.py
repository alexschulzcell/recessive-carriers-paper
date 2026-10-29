"""Body proportion (leg vs trunk) and IGF-1 (local vs systemic) of the carrier effects (exploratory).

Proportion: Genebass BETA_Burden for height and sitting height (both inverse-rank normalised, same variant weights, so the
ratio is comparable). In cm: dH = b_h * SD_h, dS = b_s * SD_s (UK Biobank: SD height ~9.3 cm, SD sitting height ~4.8 cm,
both sexes pooled); leg fraction = (dH - dS) / dH. Reference: all genes with pLoF height p < 1e-3 (genome-wide) or pool.
IGF-1: carrier effect on IGF-1 (Genebass field 30770). Local growth-plate genes should not lower IGF-1 (if anything raise it
through feedback); genes of the GH/IGF axis (IGF1R, PAPPA, IRS1, GHSR, GHRH, ...) serve as the opposite pole.

Inputs: shet_merged.tsv.gz (script 11); config.GENEBASS_DIR / genebass_<trait>_<set>.tsv.gz (scripts 01, 13)
Output: printed only.
Usage:  python 16_body_proportion_and_igf1.py
"""
import numpy as np
import pandas as pd
from scipy.stats import norm, mannwhitneyu

import config

SDH, SDS = 9.3, 4.8  # SD of height and sitting height in cm
A = pd.read_csv(config.results("shet_merged.tsv.gz"), sep="\t")


def load(t, s="plof"):
    """Genebass burden results of one trait x burden set with signed z (index gene_id)."""
    g = pd.read_csv(config.GENEBASS_DIR / f"genebass_{t}_{s}.tsv.gz", sep="\t")[["gene_id", "BETA_Burden", "Pvalue_Burden"]].dropna()
    g["z"] = np.sign(g.BETA_Burden) * norm.isf(g.Pvalue_Burden.clip(1e-300) / 2)
    return g.set_index("gene_id")


for s in ["plof", "mislc"]:
    h, si = load("height", s), load("sitheight", s)
    x = A[["gene_id", "symbol", "set", "category", f"keep_{s}", f"carriers_{s}"]].merge(
        h[["BETA_Burden", "Pvalue_Burden", "z"]], left_on="gene_id", right_index=True).merge(
        si[["BETA_Burden"]].rename(columns={"BETA_Burden": "b_s"}), left_on="gene_id", right_index=True)
    x = x[(x[f"keep_{s}"] == True) & (x[f"carriers_{s}"] >= 10)]  # noqa: E712
    x["dH"] = x.BETA_Burden * SDH
    x["dS"] = x.b_s * SDS
    x["leg_fraction"] = (x.dH - x.dS) / x.dH
    hit = x[x.Pvalue_Burden < 1e-3]
    print(f"\n[{s}] leg fraction of the height effect (median), genes with height p < 1e-3:")
    for label, sel in [("all hits genome-wide", hit.index == hit.index), ("AR_nos", hit.set == "AR_nos"),
                       ("AR_nos MATRIX", (hit.set == "AR_nos") & (hit.category == "MATRIX")),
                       ("MATRIX outside the nosology", (hit.set != "AR_nos") & (hit.category == "MATRIX")),
                       ("AD_nos+nos_other", hit.set.isin(["AD_nos", "nos_other"]))]:
        y = hit[sel]
        print(f"  {label:28s} n {len(y):3d}: leg fraction {y.leg_fraction.median():.2f} "
              f"(IQR {y.leg_fraction.quantile(.25):.2f}-{y.leg_fraction.quantile(.75):.2f})")
    y = hit[hit.set == "AR_nos"].sort_values("Pvalue_Burden")
    print(y[["symbol", "category", "dH", "dS", "leg_fraction"]].round(2).head(20).to_string(index=False))
    a = hit[hit.set.isin(["AR_nos", "AD_nos", "nos_other"])].leg_fraction
    b = hit[~hit.set.isin(["AR_nos", "AD_nos", "nos_other"])].leg_fraction
    print(f"  nosology genes vs other hits: leg fraction {a.median():.2f} vs {b.median():.2f}, MW p {mannwhitneyu(a, b).pvalue:.3f}")

print("\nIGF-1: carrier effect (z) among height hits")
if (config.GENEBASS_DIR / "genebass_igf1_plof.tsv.gz").exists():
    ig = load("igf1")
    h = load("height")
    x = A[["gene_id", "symbol", "set", "category"]].merge(h[["z", "Pvalue_Burden"]].rename(columns={"z": "z_h"}), left_on="gene_id", right_index=True).merge(
        ig[["z"]].rename(columns={"z": "z_igf"}), left_on="gene_id", right_index=True)
    hit = x[x.Pvalue_Burden < 1e-3]
    axis = hit[hit.symbol.isin(["IGF1R", "PAPPA", "IRS1", "GHSR", "GHRH", "IGF2BP2", "IGFALS", "GHR", "STAT5B", "PAPPA2", "IGF1"])]
    print(axis[["symbol", "z_h", "z_igf"]].round(2).to_string(index=False))
    for label, sel in [("AR_nos", hit.set == "AR_nos"), ("MATRIX all", hit.category == "MATRIX"), ("all hits", hit.index == hit.index)]:
        y = hit[sel]
        conc = (np.sign(y.z_h) == np.sign(y.z_igf)).mean()
        print(f"  {label:12s} n {len(y):3d}: |z_IGF1| > 1.96 in {(y.z_igf.abs() > 1.96).mean():.2f}; same direction as height {conc:.2f}; "
              f"median z_IGF1 x sign(z_h) {(y.z_igf * np.sign(y.z_h)).median():+.2f}")
    print(hit[hit.set == "AR_nos"].sort_values("Pvalue_Burden")[["symbol", "category", "z_h", "z_igf"]].round(2).head(20).to_string(index=False))
