"""Main gene x phenotype table of the skeletal nosology genes, fraction of non-silent genes, and dose continuum
(heterozygous vs biallelic carriers).

T1 Gene x phenotype (AR_nos, and AD_nos / nos_other genes as positive controls): Genebass z (pLoF, missense|LC) for height,
   sitting height, eBMD, hip and knee osteoarthritis; FDR (Benjamini-Hochberg) within the AR_nos set per phenotype x burden
   set (columns qAR_*); AstraZeneca PheWAS (UK Biobank 500k WGS, European ancestry) beta / SE / p / carriers for the ptv
   and rec (biallelic) models.
T2 Fraction of non-silent genes: Storey pi1 (lambda 0.5) of the burden p-values, AR_nos against control sets.
T3 Dose: for genes with an AZ ptv and an AZ rec result (European ancestry): beta heterozygous vs beta biallelic (height).

Inputs: gene_annotation.tsv.gz (script 03); config.GENEBASS_DIR / genebass_<trait>_<set>.tsv.gz (script 01);
        config.AZPHEWAS_DIR / az_<trait>_<model>.tsv.gz (script 06)
Output: 08_gene_x_phenotype_table.tsv in RC_RESULTS_DIR
Usage:  python 08_gene_x_phenotype_table_and_dose.py
"""
import numpy as np
import pandas as pd
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

import config

A = pd.read_csv(config.results("gene_annotation.tsv.gz"), sep="\t")
A = A[A.protein_coding].copy()
PHENOTYPES = ["height", "sitheight", "ebmd", "oa_hip", "oa_knee"]

wide = A[A.set.isin(["AR_nos", "nos_other", "AD_nos"])][["gene_id", "symbol", "set", "nos_inh", "nos_groups", "hpo_ad",
                                                         "oe_lof_upper", "carriers_plof"]].copy()
for ph in PHENOTYPES:
    for s in ["plof", "mislc"]:
        g = pd.read_csv(config.GENEBASS_DIR / f"genebass_{ph}_{s}.tsv.gz", sep="\t")[["gene_id", "BETA_Burden", "Pvalue_Burden"]]
        g = g.merge(A[["gene_id", f"keep_{s}", f"carriers_{s}", "set"]], on="gene_id")
        g = g[(g[f"keep_{s}"] == True) & (g[f"carriers_{s}"] >= 10) & g.Pvalue_Burden.notna()].copy()  # noqa: E712
        g["z"] = np.sign(g.BETA_Burden) * norm.isf(g.Pvalue_Burden.clip(1e-300) / 2)
        ar = g.set == "AR_nos"
        g.loc[ar, "q"] = multipletests(g.loc[ar, "Pvalue_Burden"], method="fdr_bh")[1]
        wide = wide.merge(g[["gene_id", "z", "q"]].rename(columns={"z": f"z_{ph}_{s}", "q": f"qAR_{ph}_{s}"}), on="gene_id", how="left")

az = {}
for ph in ["height", "sitheight", "ebmd"]:
    for m in ["ptv", "flexdmg", "rec"]:
        f = config.AZPHEWAS_DIR / f"az_{ph}_{m}.tsv.gz"
        if not f.exists():
            continue
        a = pd.read_csv(f, sep="\t")
        a = a[a.ancestry == "European"][["gene_name", "pvalue", "beta", "beta_std", "yes_qv"]]
        az[(ph, m)] = a
        if ph in ("height", "ebmd") and m in ("ptv", "rec"):
            wide = wide.merge(a[["gene_name", "beta", "beta_std", "pvalue", "yes_qv"]].rename(columns={
                "gene_name": "symbol", "beta": f"az_{ph}_{m}_beta", "beta_std": f"az_{ph}_{m}_se", "pvalue": f"az_{ph}_{m}_p",
                "yes_qv": f"az_{ph}_{m}_n"}), on="symbol", how="left")
wide = wide.sort_values(["set", "z_height_plof"])
wide.to_csv(config.results("08_gene_x_phenotype_table.tsv"), sep="\t", index=False, float_format="%.4g")

pd.set_option("display.width", 250)
pd.set_option("display.max_rows", 200)
hits = wide[(wide.set == "AR_nos") & ((wide.filter(like="qAR_") < 0.05).any(axis=1))]
cols = (["symbol", "nos_groups", "carriers_plof", "hpo_ad"] + [c for c in wide.columns if c.startswith("z_")]
        + ["az_height_ptv_beta", "az_height_ptv_n", "az_height_rec_beta", "az_height_rec_n"])
h = hits[cols].copy()
h["nos_groups"] = h.nos_groups.str.replace("Group ", "G").str[:28]
print(f"T1 AR genes with FDR<0.05 (within the AR set) in any phenotype: {len(hits)}")
print(h.round(2).to_string(index=False))

print("\nT2 Storey pi1 (fraction of non-null genes), burden p, carriers >= 10:")
for ph in ["height", "sitheight", "ebmd", "oa_hip", "oa_knee", "bmi"]:
    for s in ["plof", "mislc", "syn"]:
        g = pd.read_csv(config.GENEBASS_DIR / f"genebass_{ph}_{s}.tsv.gz", sep="\t")[["gene_id", "Pvalue_Burden"]].merge(
            A[["gene_id", "set", f"keep_{s}", f"carriers_{s}"]], on="gene_id")
        g = g[(g[f"keep_{s}"] == True) & (g[f"carriers_{s}"] >= 10) & g.Pvalue_Burden.notna()]  # noqa: E712
        out = []
        for st in ["AR_nos", "AR_nonskel", "AR_skelsyndr", "pool", "AD_nos", "nos_other"]:
            p = g[g.set == st].Pvalue_Burden
            pi1 = 1 - min(1, (p > 0.5).mean() / 0.5)
            out.append(f"{st} {pi1:.2f}")
        print(f"  {ph:9s} {s:5s}: " + " | ".join(out))

print("\nT3 dose heterozygous (ptv) vs biallelic (rec), height, AZ European ancestry, nosology genes:")
if ("height", "ptv") in az and ("height", "rec") in az:
    m = az[("height", "ptv")].merge(az[("height", "rec")], on="gene_name", suffixes=("_het", "_bi"))
    m = m.merge(A[["symbol", "set", "nos_groups"]].rename(columns={"symbol": "gene_name"}), on="gene_name")
    m = m[m.set.isin(["AR_nos", "nos_other", "AD_nos"])]
    print(m[["gene_name", "set", "beta_het", "beta_std_het", "pvalue_het", "yes_qv_het", "beta_bi", "beta_std_bi", "pvalue_bi",
             "yes_qv_bi"]].round(3).sort_values("pvalue_het").to_string(index=False))
