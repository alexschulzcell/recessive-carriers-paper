"""Does organ specificity in the disease organ predict selection against heterozygous LoF (gnomAD LOEUF), without any UK
Biobank phenotype? (pre-specified, PRAEREG_gnomad.md)

Pre-specified test. Genes: purely recessive genes with a main organ (largest HPO organ fraction from script 14, >= 0.2;
skeletal excluded because it has no tissue definition), gnomAD expected LoF >= 10. Primary test P1: OLS of the rank of
-LOEUF on the standardised tissue specificity of the gene in its main organ (HPA, MAIN_PRODUCT_TISSUES of script 22),
adjusted for the tissue specificity index tau, median expression, CDS length and organ fixed effects; the main-organ labels
are permuted between genes (B = 10000), one-sided p. Secondary: cell-type specificity (producer cells of script 27),
s_het (circular, since s_het uses LoF constraint), obs/exp; per-organ Spearman correlations; unadjusted correlations; and a
bridge to the phenotype: within-organ rank of the carrier effect (units of script 22) against within-organ rank of -LOEUF.

Inputs: hpa_tissue.zip, hpa_singlecell.zip, organ_fractions_recessive.tsv.gz (14), shet_merged.tsv.gz (11),
22_main_product_units_{plof,mislc}.tsv (22).
Outputs: 34_gnomad_constraint_prediction.tsv, 34_genes.tsv
Call: python 34_gnomad_constraint_prediction.py
"""
import zipfile

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

import config

REF = config.REF_DIR
rng = np.random.default_rng(34)
B = 10000
TISSUES = {"blood": ["bone marrow"], "immune": ["bone marrow", "spleen", "lymph node"], "liver": ["liver"],
           "lipids": ["liver", "small intestine"], "kidney": ["kidney"], "muscle": ["skeletal muscle"], "lung": ["lung"],
           "circulatory": ["heart muscle", "blood vessel"], "nervous": ["cerebral cortex"], "glucose": ["pancreas"]}  # as in 22
PRODUCER_CELLS = {"blood": ["erythrocyte progenitors", "megakaryocyte-erythroid progenitors", "megakaryocytes"],
                  "immune": ["neutrophils", "neutrophil progenitors", "t-cells", "b-cells", "monocytes", "nk-cells"],
                  "liver": ["hepatocytes", "cholangiocytes"], "lipids": ["hepatocytes", "enterocytes"],
                  "kidney": ["proximal tubule cells", "loop of henle epithelial cells", "distal convoluted tubule cells", "podocytes"],
                  "lung": ["alveolar cells type 1", "alveolar cells type 2"], "muscle": ["myonuclei"],
                  "circulatory": ["cardiomyocytes", "vascular smooth muscle cells", "vascular endothelial cells"],
                  "nervous": ["brain excitatory neurons", "brain inhibitory neurons"], "glucose": ["pancreatic islet cells"]}  # as in 27
ORGANS = list(TISSUES)


def load_pivot(zip_name, col, value_col):
    z = zipfile.ZipFile(REF / zip_name)
    return pd.read_csv(z.open(z.namelist()[0]), sep="\t").pivot_table(index="Gene", columns=col, values=value_col, aggfunc="max")


hpa = load_pivot("hpa_tissue.zip", "Tissue", "nTPM")
hpa_median = hpa.median(axis=1)
sc = load_pivot("hpa_singlecell.zip", "Cell type", "nCPM")
sc_median = sc.median(axis=1)
L = np.log2(hpa + 1)
tau = ((1 - L.div(L.max(axis=1), axis=0)).sum(axis=1) / (L.shape[1] - 1)).rename("tau")
SPEC = pd.DataFrame({o: np.log2((hpa[TISSUES[o]].max(axis=1) + 1) / (hpa_median + 1)) for o in ORGANS})
CELL_SPEC = pd.DataFrame({o: np.log2((sc[PRODUCER_CELLS[o]].max(axis=1) + 1) / (sc_median + 1)) for o in ORGANS})

organ_fractions = pd.read_csv(config.results("organ_fractions_recessive.tsv.gz"), sep="\t")
main = organ_fractions.sort_values("frac", ascending=False).drop_duplicates("gene_id")
main = main[(main.frac >= 0.2) & main.organ.isin(ORGANS)]
A = pd.read_csv(config.results("shet_merged.tsv.gz"), sep="\t")  # gene annotation + GeneBayes s_het (post_mean), obs/exp (oe), exp_lof
G = main.merge(A[["gene_id", "symbol", "oe_lof_upper", "oe", "exp_lof", "post_mean", "cds_length"]], on="gene_id")
G = G[(G.exp_lof >= 10) & G.oe_lof_upper.notna() & G.gene_id.isin(hpa.index)].copy()
G["tau"] = tau.reindex(G.gene_id).values
G["lmed"] = np.log10(hpa_median.reindex(G.gene_id).values + 1)
G["llen"] = np.log10(G.cds_length)
G = G.dropna(subset=["tau", "lmed", "llen"]).reset_index(drop=True)
print(f"Genes: {len(G)} purely recessive, main organ not skeletal, exp_lof >= 10 | per organ: {G.organ.value_counts().to_dict()}")


def test(name, y, M):
    """OLS y ~ specificity in own organ (z) + tau + lmed + llen + organ fixed effects; organ labels permuted between genes."""
    ok = ~np.isnan(y) & M.reindex(G.gene_id).notna().all(axis=1).values
    gid = G.gene_id.values[ok]
    org = G.organ.values[ok]
    yy = y[ok]
    Mg = M.reindex(gid).values
    col = {o: i for i, o in enumerate(ORGANS)}
    sub = G[ok]

    def coef(organ_labels):
        s = Mg[np.arange(len(gid)), [col[o] for o in organ_labels]]
        sz = (s - s.mean()) / s.std()
        X = np.column_stack([sz, sub.tau, sub.lmed, sub.llen] + [(organ_labels == o).astype(float) for o in ORGANS])
        return np.linalg.lstsq(X, yy, rcond=None)[0][0]

    obs = coef(org)
    d = np.array([coef(rng.permutation(org)) for _ in range(B)])
    p = (1 + (d >= obs).sum()) / (B + 1)
    print(f"  {name:42s} n {ok.sum():4d}: coef {obs:+.4f} (perm {d.mean():+.4f} +- {d.std():.4f}; "
          f"dz {(obs - d.mean()) / d.std():+.1f}), p one-sided {p:.4f}")
    return dict(test=name, n=int(ok.sum()), coef=obs, perm_mean=d.mean(), perm_sd=d.std(), p=p)


rows = []
yL = rankdata(-G.oe_lof_upper.values) / len(G)
yS = np.where(G.post_mean.notna(), G.post_mean.rank(pct=True).values, np.nan)
yO = rankdata(-G.oe.fillna(G.oe.median()).values) / len(G)
print("\n[Primary] tissue specificity in the disease organ -> LOEUF")
rows.append(test("P1 tissue specificity -> -LOEUF", yL, SPEC))
print("[Secondary]")
rows.append(test("S1 cell type specificity -> -LOEUF", yL, CELL_SPEC))
rows.append(test("S2 tissue specificity -> s_het (circular?)", yS, SPEC))
rows.append(test("S3 tissue specificity -> -obs/exp", yO, SPEC))
for o in ORGANS:
    m = G.organ.values == o
    if m.sum() >= 20:
        r = spearmanr(SPEC.reindex(G.gene_id[m]).values[:, ORGANS.index(o)], -G.oe_lof_upper[m])
        print(f"    {o:11s} n {m.sum():3d}: rho(specificity, -LOEUF) {r.correlation:+.2f} (p {r.pvalue:.3f})")
# unadjusted association (context: tissue-specific genes are in general less constrained)
own_spec = SPEC.reindex(G.gene_id).values[np.arange(len(G)), [ORGANS.index(o) for o in G.organ]]
print(f"  unadjusted: rho(specificity, -LOEUF) {spearmanr(own_spec, -G.oe_lof_upper).correlation:+.3f}; "
      f"rho(tau, -LOEUF) {spearmanr(G.tau, -G.oe_lof_upper).correlation:+.3f}")

# Bridge to the phenotype: carrier effect (locus-adjusted organ value, script 22) vs -LOEUF within organ
print("\n[Bridge] carrier effect (22) ~ -LOEUF within organ")
for mask in ["plof", "mislc"]:
    U = pd.read_csv(config.results(f"22_main_product_units_{mask}.tsv"), sep="\t").merge(
        A[["gene_id", "oe_lof_upper"]], on="gene_id").dropna(subset=["oe_lof_upper", "value"])
    U["rw"] = U.groupby("organ").value.rank(pct=True)
    U["rl"] = U.groupby("organ").oe_lof_upper.rank(pct=True, ascending=False)
    r = spearmanr(U.rw, U.rl)
    print(f"  {mask:5s} n {len(U)}: rho(within-organ rank carrier effect, within-organ rank -LOEUF) {r.correlation:+.3f} "
          f"(p {r.pvalue:.2g})")
    rows.append(dict(test=f"bridge_{mask}", n=len(U), coef=r.correlation, p=r.pvalue))
pd.DataFrame(rows).to_csv(config.results("34_gnomad_constraint_prediction.tsv"), sep="\t", index=False)
G.assign(own_specificity=own_spec).to_csv(config.results("34_genes.tsv"), sep="\t", index=False)
