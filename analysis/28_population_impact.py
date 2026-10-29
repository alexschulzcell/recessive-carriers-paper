"""Population impact of carrier effects (pre-specified, descriptive; PRAEREG_stufe9.md, part B).

Pre-specified descriptive analysis (test register no. 32; deviation: TG, platelets, creatinine and WBC were added and were
not pre-specified). Per trait: recessive genes of the matching organ (height: recessive nosology genes, AR_nos) with a Genebass
pLoF burden FDR < 0.05 within the gene set; effect per carrier from AZ PheWAS (European ancestry, PTV, SD units); carrier
frequency f = carriers / N (Genebass, pLoF). Summary measures: fraction of people with >= 1 pLoF carrier genotype;
attributable fraction at the tails < -2 SD and > +2 SD; variance explained sum f * beta^2. Conservative version: only genes
with AZ p < 1e-4. In addition the carrier fraction over ALL genes of the set (silent or not).

Inputs: gene_annotation.tsv.gz (03), organ_fractions_recessive.tsv.gz (14), Genebass and AZ PheWAS files.
Output: 28_population_impact.tsv
Call: python 28_population_impact.py
"""
import numpy as np
import pandas as pd
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

import config

N = 394841
annotation = pd.read_csv(config.results("gene_annotation.tsv.gz"), sep="\t")
organ_fractions = pd.read_csv(config.results("organ_fractions_recessive.tsv.gz"), sep="\t")
members = organ_fractions[organ_fractions.frac >= 0.2]
# trait (AZ code) -> (organ, Genebass trait code)
TRAITS = {"height": ("skeletal", "height"), "hdl": ("lipids", "hdl"), "tg": ("lipids", "tg"), "hb": ("blood", "hb"),
          "mcv": ("blood", "mcv"), "plt": ("blood", "plt"), "urate": ("kidney", "urate"), "creat": ("kidney", "creat"),
          "bili": ("liver", "bili"), "wbc": ("immune", "wbc")}
rows = []
for t, (organ, gb) in TRAITS.items():
    genes = set(annotation[annotation.set == "AR_nos"].gene_id) if t == "height" else set(members[members.organ == organ].gene_id)
    g = pd.read_csv(config.GENEBASS_DIR / f"genebass_{gb}_plof.tsv.gz", sep="\t")[["gene_id", "Pvalue_Burden"]].dropna().merge(
        annotation[["gene_id", "symbol", "keep_plof", "carriers_plof"]], on="gene_id")
    g = g[g.gene_id.isin(genes) & (g.keep_plof == True) & (g.carriers_plof >= 10)].copy()  # noqa: E712
    g["q"] = multipletests(g.Pvalue_Burden, method="fdr_bh")[1]
    az = pd.read_csv(config.AZPHEWAS_DIR / f"az_{t}_ptv.tsv.gz", sep="\t")
    az = az[az.ancestry == "European"][["gene_name", "beta", "pvalue", "yes_qv"]]
    x = g[g.q < 0.05].merge(az, left_on="symbol", right_on="gene_name")
    all_f = (g.carriers_plof / N).clip(upper=0.5)
    fraction_all = 1 - np.prod(1 - all_f)
    for analysis, sel in [("FDR<0.05", x.index == x.index), ("conservative AZ p<1e-4", x.pvalue < 1e-4)]:
        y = x[sel]
        f = (y.carriers_plof / N).values
        b = y.beta.values
        carrier_fraction = 1 - np.prod(1 - f)
        lo = (f * (norm.cdf(-2 - b) - norm.cdf(-2)))[b < 0].sum() / norm.cdf(-2)
        hi = (f * (norm.sf(2 - b) - norm.sf(2)))[b > 0].sum() / norm.sf(2)
        rows.append(dict(trait=t, analysis=analysis, n_genes=len(y), carrier_fraction=carrier_fraction,
                         one_in=round(1 / carrier_fraction) if carrier_fraction > 0 else np.nan,
                         AF_lower=lo, AF_upper=hi, variance_explained=(f * b ** 2).sum(), carriers_all_AR_genes=fraction_all,
                         genes=", ".join(y.sort_values("beta").symbol.head(8))))
R = pd.DataFrame(rows)
pd.set_option("display.width", 250)
pd.set_option("display.max_colwidth", 80)
print(R.round(4).to_string(index=False))
R.to_csv(config.results("28_population_impact.tsv"), sep="\t", index=False)
