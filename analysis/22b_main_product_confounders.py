"""Confounders of the main-product rule (script 22, PRAEREG_hauptprodukt.md).

Exploratory follow-up. Re-implementation of an exploratory analysis whose original output was only a log.
Gene x organ units from script 22 (without skeletal units); covariates from shet_merged (script 11): number of carriers
(log10), LOEUF (oe_lof_upper), s_het (log10 GeneBayes posterior mean), CDS length (log10). Units with a missing covariate are
dropped.
  (1) median number of carriers and s_het for route yes/no (W1, W2b)
  (2) OLS value ~ W + ltr + oe_lof_upper + lshet + lcds + C(organ), robust SE (HC1)
  (3) W2b per organ (difference yes - no) and dose gradient: mean value per specificity tertile within the organ (same
      boundaries as W2b)

Inputs: 22_main_product_units_{plof,mislc}.tsv (22), shet_merged.tsv.gz (11). Output: log only.
Call: python 22b_main_product_confounders.py
"""
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

import common
import config

M = pd.read_csv(config.results("shet_merged.tsv.gz"), sep="\t")


def organ_categorical(x):
    """Organ as a categorical in the level order of the original analysis (identical regression design)."""
    return x.assign(organ=pd.Categorical(x.organ, categories=common.sort_organs(x.organ.unique())))


for mask in ["plof", "mislc"]:
    U = pd.read_csv(config.results(f"22_main_product_units_{mask}.tsv"), sep="\t")
    X = U[U.organ != "skeletal"].merge(M[["gene_id", "post_mean", "lshet", f"carriers_{mask}", "oe_lof_upper", "cds_length"]],
                                       on="gene_id", how="left")
    X["ltr"] = np.log10(X[f"carriers_{mask}"])
    X["lcds"] = np.log10(X.cds_length)
    X = X.dropna(subset=["lshet", "oe_lof_upper", "lcds"])
    for route in ["W1", "W2b"]:
        Y = X if route == "W1" else X[X.specificity.notna()]
        carriers = Y.groupby(route)[f"carriers_{mask}"].median()
        shet = Y.groupby(route).post_mean.median()
        print(f"{mask} {route} carriers median yes/no: {carriers[True]} {carriers[False]} | "
              f"s_het yes/no: {round(shet[True], 4)} {round(shet[False], 4)}")
        Y = Y.assign(**{route: Y[route].astype(int)})
        f = smf.ols(f"value ~ {route} + ltr + oe_lof_upper + lshet + lcds + C(organ)", data=organ_categorical(Y)).fit(cov_type="HC1")
        print(f"   adjusted: {route} {f.params[route]:+.3f} (p {f.pvalues[route]:.4f}); "
              f"ltr {f.params['ltr']:+.3f} (p {f.pvalues['ltr']:.3f})")
    Y = X[X.specificity.notna()].copy()
    per_organ = []
    for o in common.sort_organs(Y.organ.unique()):
        x = Y[Y.organ == o]
        per_organ.append(f"{o} {x[x.W2b].value.mean() - x[~x.W2b].value.mean():+.2f} ({int(x.W2b.sum())}/{len(x)})")
    print("   per organ W2b: " + " | ".join(per_organ))
    q = Y.groupby("organ").specificity
    Y["tertile"] = (1 + (Y.specificity >= q.transform(lambda v: v.quantile(1 / 3))).astype(int)
                    + (Y.specificity >= q.transform(lambda v: v.quantile(2 / 3))).astype(int))
    print(f"   dose gradient (specificity tertiles): {Y.groupby('tertile').value.mean().round(3).to_dict()}", flush=True)
