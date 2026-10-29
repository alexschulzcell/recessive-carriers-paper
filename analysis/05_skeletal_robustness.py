"""Robustness of the height / skeletal carrier signal of the recessive skeletal dysplasia genes (AR_nos).

R1 Breadth: fraction p < 0.05 and mean z^2 in AR_nos against matched pool genes after removing the top k genes
   (k = 0, 5, 10, 20).
R2 Matching on the allele count in the burden set (log CAF decile) instead of LOEUF x CDS length.
R3 Power: minimal detectable effect per gene (SD units, 80 % power, alpha 0.05 / number of genes) from the pLoF carrier count.
R4 Nosology groups: mean z^2 per group (>= 5 recessive genes) for height, sitting height, eBMD and ALP, against matched pool.
R5 Fisher: AR_nos against AR_nonskel and AR_skelsyndr directly (fraction p < 0.05), unmatched.

Inputs:  gene_annotation.tsv.gz (script 03); config.GENEBASS_DIR / genebass_<trait>_<set>.tsv.gz (script 01)
Outputs: 05_robustness.tsv (R1/R2), 05_power_recessive.tsv (R3), 05_nosology_groups.tsv (R4), in RC_RESULTS_DIR.
Usage:   python 05_skeletal_robustness.py
"""
import numpy as np
import pandas as pd
from scipy.stats import norm, fisher_exact

import config

rng = np.random.default_rng(2)
B = 2000
A = pd.read_csv(config.results("gene_annotation.tsv.gz"), sep="\t")
A = A[A.protein_coding & A.oe_lof_upper.notna() & A.cds_length.notna()].copy()
A["lb"] = pd.qcut(A.oe_lof_upper, 10, labels=False, duplicates="drop")
A["cb"] = pd.qcut(A.cds_length, 5, labels=False, duplicates="drop")
A["stratum"] = A.lb * 10 + A.cb


def load(ph, s, min_carriers=10):
    """Genebass burden results of one trait x burden set, merged with the annotation and filtered."""
    g = pd.read_csv(config.GENEBASS_DIR / f"genebass_{ph}_{s}.tsv.gz", sep="\t")[["gene_id", "BETA_Burden", "Pvalue_Burden"]].merge(A, on="gene_id")
    g = g[(g[f"keep_{s}"] == True) & (g[f"carriers_{s}"] >= min_carriers)].copy()  # noqa: E712
    g["z"] = np.sign(g.BETA_Burden) * norm.isf(g.Pvalue_Burden.clip(1e-300) / 2)
    g["sig"] = (g.Pvalue_Burden < 0.05).astype(float)
    g["z2"] = g.z ** 2
    g["cafb"] = pd.qcut(np.log10(g[f"caf_{s}"]), 10, labels=False, duplicates="drop")
    return g


def matched(T, C, col, key, B=B):
    """B matched draws of the mean of `col`: for every test gene one control gene from the same stratum `key`."""
    by_stratum = {k: v[col].to_numpy() for k, v in C.groupby(key)}
    tot = np.zeros(B)
    for k, n in T[key].value_counts().items():
        arr = by_stratum[k] if k in by_stratum else by_stratum[min(by_stratum, key=lambda x: abs(x - k))]
        tot += arr[rng.integers(0, len(arr), size=(B, n))].sum(axis=1)
    return tot / len(T)


out = []
print("R1/R2: fraction p<0.05 and mean z^2, AR_nos against pool, after removing the top k genes")
for ph in ["height", "sitheight", "ebmd", "oa_hip", "oa_knee", "alp", "bmi"]:
    for s in ["plof", "mislc", "syn"]:
        g = load(ph, s)
        T0 = g[g.set == "AR_nos"].sort_values("Pvalue_Burden")
        C = g[g.set == "pool"]
        for k in [0, 5, 10, 20]:
            T = T0.iloc[k:]
            for key in ["stratum", "cafb"]:
                for col in ["sig", "z2"]:
                    d = matched(T, C, col, key)
                    o = T[col].mean()
                    # only the test genes lose their top genes; the controls keep theirs
                    out.append(dict(pheno=ph, bset=s, drop_top=k, match=key, stat=col, n=len(T), obs=o, ctrl=d.mean(),
                                    dz=(o - d.mean()) / d.std(), p=(1 + (d >= o).sum()) / (B + 1)))
R = pd.DataFrame(out)
R.to_csv(config.results("05_robustness.tsv"), sep="\t", index=False)
pd.set_option("display.width", 250)
print(R[R.stat == "sig"].pivot_table(index=["pheno", "bset"], columns=["match", "drop_top"], values="dz").round(1).to_string())
print("\nobs/ctrl fraction p<0.05 (matching LOEUF x length):")
print(R[(R.stat == "sig") & (R.match == "stratum")].pivot_table(index=["pheno", "bset"], columns="drop_top",
                                                               values=["obs", "ctrl"]).round(3).to_string())

# R3 power
T = A[A.set == "AR_nos"].copy()
ntest = (T.carriers_plof >= 10).sum()
zc = norm.isf(0.05 / ntest / 2) + norm.ppf(0.8)
T["mde_sd_plof"] = zc / np.sqrt(T.carriers_plof.clip(lower=1))
print(f"\nR3 power (pLoF, quantitative trait, alpha 0.05/{ntest}): carriers per AR gene median {T.carriers_plof.median():.0f}, "
      f"IQR {T.carriers_plof.quantile(.25):.0f}-{T.carriers_plof.quantile(.75):.0f}; genes with >=100 carriers: {(T.carriers_plof >= 100).sum()}, "
      f">=500: {(T.carriers_plof >= 500).sum()}, <10: {(T.carriers_plof < 10).sum()}")
print("minimal detectable effect (SD) median {:.2f}; genes with MDE<0.2 SD: {}; <0.1 SD: {}".format(
    T.mde_sd_plof.median(), (T.mde_sd_plof < 0.2).sum(), (T.mde_sd_plof < 0.1).sum()))
T[["gene_id", "symbol", "nos_groups", "oe_lof_upper", "carriers_plof", "mde_sd_plof"]].sort_values("carriers_plof").to_csv(
    config.results("05_power_recessive.tsv"), sep="\t", index=False)

# R4 nosology groups
rows = []
for ph in ["height", "sitheight", "ebmd", "alp"]:
    for s in ["plof", "mislc"]:
        g = load(ph, s)
        C = g[g.set == "pool"]
        T = g[g.set.isin(["AR_nos", "AD_nos", "nos_other"])].copy()
        T["grp"] = T.nos_groups.str.split("|")
        T = T.explode("grp")
        for (grp, st), sub in T.groupby(["grp", "set"]):
            if len(sub) < 5 or st != "AR_nos":
                continue
            d = matched(sub, C, "z2", "stratum")
            rows.append(dict(pheno=ph, bset=s, grp=grp[:60], n=len(sub), mean_z2=sub.z2.mean(), mean_z=sub.z.mean(),
                             dz=(sub.z2.mean() - d.mean()) / d.std(),
                             top=",".join(sub.sort_values("Pvalue_Burden").symbol.head(3))))
G = pd.DataFrame(rows)
G.to_csv(config.results("05_nosology_groups.tsv"), sep="\t", index=False)
print("\nR4 nosology groups (AR genes only, >=5), dz of mean z^2 against matched pool; height:")
print(G[G.pheno == "height"].pivot_table(index="grp", columns="bset", values=["dz", "n"]).round(1)
      .sort_values(("dz", "plof"), ascending=False).to_string())

# R5 Fisher, direct comparison
print("\nR5 Fisher (p<0.05), AR_nos against other recessive gene sets, unmatched:")
for ph in ["height", "sitheight", "ebmd", "oa_hip"]:
    for s in ["plof", "mislc", "syn"]:
        g = load(ph, s)
        a = g[g.set == "AR_nos"].sig
        for c in ["AR_nonskel", "AR_skelsyndr"]:
            b = g[g.set == c].sig
            orr, p = fisher_exact([[a.sum(), len(a) - a.sum()], [b.sum(), len(b) - b.sum()]])
            print(f"  {ph:9s} {s:5s} vs {c:12s}: {a.mean():.3f} against {b.mean():.3f}, OR {orr:.2f}, p {p:.1e}")
