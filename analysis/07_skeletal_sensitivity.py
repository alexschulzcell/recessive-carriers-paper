"""Sensitivity analyses of the height carrier signal of the recessive skeletal dysplasia genes.

S1 Strictly recessive: AR_nos genes without any dominant OMIM disease (HPO hpo_ad == 0); enrichment (fraction p < 0.05,
   mean z^2) against matched pool and AR_nonskel genes.
S2 Internal replication: pLoF and missense|LC carriers are largely different people. Correlation z_pLoF ~ z_missense per
   gene, AR_nos against control sets; and: do genes with pLoF p < 0.05 show the same direction in the missense set?
S3 Synonymous per gene: are pLoF z and synonymous z correlated (hint of linkage disequilibrium with common height variants)?

Inputs: gene_annotation.tsv.gz (script 03); config.GENEBASS_DIR / genebass_<trait>_<set>.tsv.gz (script 01)
Output: printed only.
Usage:  python 07_skeletal_sensitivity.py
"""
import numpy as np
import pandas as pd
from scipy.stats import norm, spearmanr, binomtest

import config

rng = np.random.default_rng(3)
A = pd.read_csv(config.results("gene_annotation.tsv.gz"), sep="\t")
A = A[A.protein_coding & A.oe_lof_upper.notna() & A.cds_length.notna()].copy()
A["stratum"] = pd.qcut(A.oe_lof_upper, 10, labels=False) * 10 + pd.qcut(A.cds_length, 5, labels=False)


def z(ph, s):
    """Signed burden z per gene (index gene_id) for genes passing the keep flags with >= 10 carriers."""
    g = pd.read_csv(config.GENEBASS_DIR / f"genebass_{ph}_{s}.tsv.gz", sep="\t")[["gene_id", "BETA_Burden", "Pvalue_Burden"]].merge(
        A[["gene_id", f"keep_{s}", f"carriers_{s}"]], on="gene_id")
    g = g[(g[f"keep_{s}"] == True) & (g[f"carriers_{s}"] >= 10)]  # noqa: E712
    return g.set_index("gene_id").pipe(lambda d: np.sign(d.BETA_Burden) * norm.isf(d.Pvalue_Burden.clip(1e-300) / 2))


def matched_p(T, C, col, B=2000):
    """Observed mean, matched control mean, dz and empirical one-sided p (LOEUF x CDS-length strata)."""
    by_stratum = {k: v[col].to_numpy() for k, v in C.groupby("stratum")}
    tot = np.zeros(B)
    for k, n in T.stratum.value_counts().items():
        arr = by_stratum[k] if k in by_stratum else by_stratum[min(by_stratum, key=lambda x: abs(x - k))]
        tot += arr[rng.integers(0, len(arr), size=(B, n))].sum(axis=1)
    tot /= len(T)
    o = T[col].mean()
    return o, tot.mean(), (o - tot.mean()) / tot.std(), (1 + (tot >= o).sum()) / (B + 1)


print("S1 strictly recessive (no dominant OMIM disease):")
for ph in ["height", "sitheight", "ebmd", "oa_hip", "alp"]:
    for s in ["plof", "mislc", "syn"]:
        zz = z(ph, s).rename("z")
        g = A.merge(zz, left_on="gene_id", right_index=True)
        g["sig"] = (2 * norm.sf(g.z.abs()) < 0.05).astype(float)
        g["z2"] = g.z ** 2
        T = g[(g.set == "AR_nos") & (g.hpo_ad != 1)]
        for c in ["pool", "AR_nonskel"]:
            o, m, dz, p = matched_p(T, g[g.set == c], "sig")
            o2, m2, dz2, p2 = matched_p(T, g[g.set == c], "z2")
            print(f"  {ph:9s} {s:5s} n={len(T):3d} vs {c:10s}: p<0.05 {o:.3f} against {m:.3f} (dz {dz:.1f}, p {p:.4f}); "
                  f"z^2 {o2:.2f} against {m2:.2f} (dz {dz2:.1f}, p {p2:.4f})")

print("\nS2 internal replication pLoF vs missense|LC (height):")
for ph in ["height", "sitheight", "ebmd"]:
    zp, zm = z(ph, "plof").rename("zp"), z(ph, "mislc").rename("zm")
    g = A.merge(zp, left_on="gene_id", right_index=True).merge(zm, left_on="gene_id", right_index=True)
    for st in ["AR_nos", "AR_nonskel", "pool", "AD_nos", "nos_other"]:
        x = g[g.set == st]
        r, p = spearmanr(x.zp, x.zm)
        hit = x[2 * norm.sf(x.zp.abs()) < 0.05]
        same = (np.sign(hit.zp) == np.sign(hit.zm)).sum()
        bt = binomtest(int(same), len(hit), 0.5, alternative="greater").pvalue if len(hit) else np.nan
        print(f"  {ph:9s} {st:10s} n={len(x):5d}: rho {r:+.3f} (p {p:.1e}); pLoF hits p<0.05: {len(hit)}, "
              f"same direction in missense {same} (binomial p {bt:.3f})")

print("\nS3 pLoF z against synonymous z per gene (height):")
zp, zs = z("height", "plof").rename("zp"), z("height", "syn").rename("zs")
g = A.merge(zp, left_on="gene_id", right_index=True).merge(zs, left_on="gene_id", right_index=True)
for st in ["AR_nos", "pool"]:
    x = g[g.set == st]
    r, p = spearmanr(x.zp, x.zs)
    print(f"  {st}: rho {r:+.3f} (p {p:.2f}); |z_syn|>1.96 among pLoF hits: "
          f"{((x.zs.abs() > 1.96) & (x.zp.abs() > 1.96)).sum()} of {(x.zp.abs() > 1.96).sum()}")
