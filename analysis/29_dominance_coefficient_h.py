"""Dominance coefficient h = beta_het / beta_hom (pre-specified, PRAEREG_stufe9.md, part C).

beta_hom is the patient SDS from the literature table (literature_patient_values.tsv, scripts 32/32b; genes chosen by
data availability). beta_het (SD per carrier) is the AstraZeneca PheWAS PTV estimate (UK Biobank 500k, European
ancestry) where reported; otherwise the Genebass z_pLoF / sqrt(carriers), calibrated by a regression through the origin
on the AZ values (genes with both estimates, pooled over traits).
Rules: |beta_hom| >= 0.5; h clipped to [-1, 2]; negative h are reported.
H1 median h << 0.5; H2 h(main product gene: W1 or W2b from script 22) > h(other genes), one-sided Mann-Whitney;
H3 Spearman(h, s_het) > 0. Validation: heterozygotes reported in the literature -> h_lit = sds_het_lit / sds_patients.

Inputs: literature_patient_values.tsv, gene_annotation.tsv.gz, 22_main_product_units_plof.tsv (results folder),
Genebass pLoF burden files, AZ PheWAS PTV files, GeneBayes s_het.
Output: 29_dominance_h.tsv (read by 29b).
Call: python 29_dominance_coefficient_h.py
"""
import numpy as np
import pandas as pd
from scipy.stats import norm, mannwhitneyu, spearmanr, wilcoxon

import common
import config

L = pd.read_csv(config.results("literature_patient_values.tsv"), sep="\t")
L.columns = [c.strip() for c in L.columns]
for c in ["patient_mean", "control_or_reference_mean", "reference_sd", "sds_patients", "het_literature_mean", "patients_n"]:
    L[c] = pd.to_numeric(L[c], errors="coerce")
L["trait"] = L.trait.str.strip().str.lower().replace({"neut": "wbc", "height_sds": "height"})
# recalculate missing SDS where a reference is available
m = L.sds_patients.isna() & L.patient_mean.notna() & L.control_or_reference_mean.notna() & L.reference_sd.gt(0)
L.loc[m, "sds_patients"] = (L.loc[m, "patient_mean"] - L.loc[m, "control_or_reference_mean"]) / L.loc[m, "reference_sd"]
L = L.dropna(subset=["sds_patients"])
L = L.sort_values("patients_n", ascending=False).drop_duplicates(["symbol", "trait"])
print(f"Literature pairs with patient SDS: {len(L)}; per organ: {L.organ.value_counts().to_dict()}")

annotation = pd.read_csv(config.results("gene_annotation.tsv.gz"), sep="\t")
shet = pd.read_csv(config.GENEBAYES_FILE, sep="\t").rename(columns={"ensg": "gene_id"})[["gene_id", "post_mean"]]
GENEBASS_TRAIT = {"height": "height", "hdl": "hdl", "ldl": "ldl", "tg": "tg", "bili": "bili", "alt": "alt", "ggt": "ggt",
                  "hb": "hb", "mcv": "mcv", "plt": "plt", "urate": "urate", "creat": "creat", "wbc": "wbc", "lymph": "lymph"}
rows = []
for t in sorted(set(L.trait) & set(GENEBASS_TRAIT)):
    g = (pd.read_csv(config.GENEBASS_DIR / f"genebass_{GENEBASS_TRAIT[t]}_plof.tsv.gz", sep="\t")
         [["gene_id", "BETA_Burden", "Pvalue_Burden"]].dropna().drop_duplicates("gene_id"))
    g = g.merge(annotation[["gene_id", "symbol", "carriers_plof", "keep_plof"]], on="gene_id")
    g = g[(g.carriers_plof >= 5)]
    g["b_gb"] = np.sign(g.BETA_Burden) * norm.isf(g.Pvalue_Burden.clip(1e-300) / 2) / np.sqrt(g.carriers_plof)
    f = config.AZPHEWAS_DIR / f"az_{t}_ptv.tsv.gz"
    if f.exists():
        az = pd.read_csv(f, sep="\t")
        az = az[az.ancestry == "European"][["gene_name", "beta"]].rename(columns={"gene_name": "symbol", "beta": "b_az"})
        g = g.merge(az, on="symbol", how="left")
    else:
        g["b_az"] = np.nan
    g["trait"] = t
    rows.append(g)
G = pd.concat(rows)
both = G.dropna(subset=["b_az"])
k = (both.b_az * both.b_gb).sum() / (both.b_gb ** 2).sum()
print(f"Calibration Genebass -> AZ: factor {k:.3f} (n {len(both)} gene x trait; r {both.b_az.corr(both.b_gb):.2f})")
G["b_het"] = G.b_az.fillna(k * G.b_gb)
G["het_source"] = np.where(G.b_az.notna(), "AZ", "Genebass calibrated")
X = L.merge(G[["symbol", "trait", "gene_id", "b_het", "het_source", "carriers_plof"]], on=["symbol", "trait"], how="left")
print("without carrier estimate (no Genebass pLoF):", ", ".join(X[X.b_het.isna()].symbol + ":" + X[X.b_het.isna()].trait))
X = X.dropna(subset=["b_het"])
X = X[X.sds_patients.abs() >= 0.5].copy()
X["h"] = (X.b_het / X.sds_patients).clip(-1, 2)
units = pd.read_csv(config.results("22_main_product_units_plof.tsv"), sep="\t")[["gene_id", "organ", "W1", "W2b"]]
X = X.merge(units, on=["gene_id", "organ"], how="left").merge(shet, on="gene_id", how="left")
X["main"] = (X.W1.fillna(False).astype(bool) | X.W2b.fillna(False).astype(bool))
X.to_csv(config.results("29_dominance_h.tsv"), sep="\t", index=False)
pd.set_option("display.width", 220)
print(f"\nPairs with h: {len(X)} (organs {X.organ.value_counts().to_dict()})")
print(X.sort_values("h", ascending=False)[["symbol", "organ", "trait", "patients_n", "sds_patients", "b_het", "het_source",
                                           "h", "main", "post_mean"]].round(3).to_string(index=False))
print(f"\nH1: median h {X.h.median():.3f} (IQR {X.h.quantile(.25):.3f}-{X.h.quantile(.75):.3f}); fraction h > 0: {(X.h > 0).mean():.2f}; "
      f"Wilcoxon h < 0.5 p {wilcoxon(X.h - 0.5, alternative='less').pvalue:.1e}; "
      f"Wilcoxon h > 0 p {wilcoxon(X.h, alternative='greater').pvalue:.1e}")
organ_median = X.groupby("organ").h.median().round(3)
print("   median h per organ:", organ_median.reindex(common.sort_organs(organ_median.index)).to_dict())
a, b = X[X.main].h, X[~X.main].h
print(f"H2: main product h {a.median():.3f} (n {len(a)}) vs other {b.median():.3f} (n {len(b)}); "
      f"one-sided MW p {mannwhitneyu(a, b, alternative='greater').pvalue:.3f}")
y = X.dropna(subset=["post_mean"])
r, p = spearmanr(y.h, y.post_mean)
print(f"H3: Spearman h ~ s_het rho {r:+.3f}, one-sided p {p / 2 if r > 0 else 1 - p / 2:.3f} (n {len(y)})")
V = X.dropna(subset=["het_literature_mean"]).copy()
if len(V):
    V["sds_het_lit"] = (V.het_literature_mean - V.control_or_reference_mean) / V.reference_sd
    V["h_lit"] = V.sds_het_lit / V.sds_patients
    print(f"\nValidation with heterozygotes from the literature (n {len(V)}): Spearman h vs h_lit "
          f"{spearmanr(V.h, V.h_lit)[0]:+.2f}; median h_lit {V.h_lit.median():.3f}")
    print(V[["symbol", "trait", "h", "h_lit", "het_n"]].round(3).to_string(index=False))
