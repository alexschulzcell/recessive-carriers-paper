"""Cell-type resolution of the main-product rule: producer cell instead of tissue (pre-specified, PRAEREG_stufe9.md, part A, W4).

Pre-specified test W4 (test register no. 31). Uses the gene x organ units of script 22 (locus-adjusted value, tissue
specificity). For each unit, the cell-type specificity is log2((max nCPM over the organ's producer cells + 1) / (median nCPM
over all cell types + 1)) from the HPA single-cell data (PRODUCER_CELLS); W4 = upper tertile within the organ. Test: difference
of mean value W4 yes - no, labels permuted within organ (B = 10000); the same for W2b on the same units for comparison; an OLS
of the value on both standardised specificities with organ fixed effects (robust SE), and per-organ differences.

Inputs: hpa_singlecell.zip, 22_main_product_units_{plof,mislc}.tsv (22).
Outputs: 27_producer_cell_units_{plof,mislc}.tsv
Call: python 27_producer_cell_type_w4.py
"""
import zipfile

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

import common
import config

rng = np.random.default_rng(31)
PRODUCER_CELLS = {"blood": ["erythrocyte progenitors", "megakaryocyte-erythroid progenitors", "megakaryocytes"],
                  "immune": ["neutrophils", "neutrophil progenitors", "t-cells", "b-cells", "monocytes", "nk-cells"],
                  "liver": ["hepatocytes", "cholangiocytes"], "lipids": ["hepatocytes", "enterocytes"],
                  "kidney": ["proximal tubule cells", "loop of henle epithelial cells", "distal convoluted tubule cells", "podocytes"],
                  "lung": ["alveolar cells type 1", "alveolar cells type 2"], "muscle": ["myonuclei"],
                  "circulatory": ["cardiomyocytes", "vascular smooth muscle cells", "vascular endothelial cells"],
                  "nervous": ["brain excitatory neurons", "brain inhibitory neurons"], "glucose": ["pancreatic islet cells"]}
z = zipfile.ZipFile(config.REF_DIR / "hpa_singlecell.zip")
sc = pd.read_csv(z.open(z.namelist()[0]), sep="\t").pivot_table(index="Gene", columns="Cell type", values="nCPM", aggfunc="max")
sc_median = sc.median(axis=1)


def perm(X, col, B=10000):
    """Difference of mean value (yes - no); labels permuted within organ (organs in order of first appearance)."""
    lab = X[col].values.astype(bool)
    w = X.value.values
    idx = [np.where(X.organ.values == o)[0] for o in X.organ.unique()]
    obs = w[lab].mean() - w[~lab].mean()
    d = np.empty(B)
    for i in range(B):
        L = lab.copy()
        for ix in idx:
            L[ix] = rng.permutation(L[ix])
        d[i] = w[L].mean() - w[~L].mean()
    return obs, (1 + (d >= obs).sum()) / (B + 1), int(lab.sum())


for mask in ["plof", "mislc"]:
    U = pd.read_csv(config.results(f"22_main_product_units_{mask}.tsv"), sep="\t")
    U = U[U.organ.isin(PRODUCER_CELLS)].copy()
    U["cell_specificity"] = [np.log2((sc.loc[g, PRODUCER_CELLS[o]].max() + 1) / (sc_median[g] + 1)) if g in sc.index else np.nan
                             for g, o in zip(U.gene_id, U.organ)]
    U = U.dropna(subset=["cell_specificity", "specificity"])
    U["W4"] = U.groupby("organ").cell_specificity.transform(lambda x: x >= x.quantile(2 / 3))
    obs, p, n = perm(U, "W4")
    obs2, p2, n2 = perm(U, "W2b")
    print(f"[{mask}] W4 cell type: diff {obs:+.3f} (n yes {n}/{len(U)}), p {p:.4f} | "
          f"W2b tissue (same units): diff {obs2:+.3f}, p {p2:.4f}")
    for c in ["specificity", "cell_specificity"]:
        U[c + "_z"] = U.groupby("organ")[c].transform(lambda x: (x - x.mean()) / x.std())
    reg = U.assign(organ=pd.Categorical(U.organ, categories=common.sort_organs(U.organ.unique())))  # level order of the original
    f = smf.ols("value ~ specificity_z + cell_specificity_z + C(organ)", data=reg).fit(cov_type="HC1")
    print(f"   regression: tissue {f.params['specificity_z']:+.3f} (p {f.pvalues['specificity_z']:.4f}) | "
          f"cell type {f.params['cell_specificity_z']:+.3f} (p {f.pvalues['cell_specificity_z']:.4f}); "
          f"correlation tissue/cell specificity r {U.specificity.corr(U.cell_specificity):.2f}")
    per_organ = []
    for o in common.sort_organs(U.organ.unique()):
        x = U[U.organ == o]
        per_organ.append(f"{o} {x[x.W4].value.mean() - x[~x.W4].value.mean():+.2f} ({int(x.W4.sum())}/{len(x)})")
    print("   per organ W4:", " | ".join(per_organ))
    U.to_csv(config.results(f"27_producer_cell_units_{mask}.tsv"), sep="\t", index=False)
