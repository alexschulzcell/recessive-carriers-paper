"""Pre-specified consistency checks P2 and P3 (PRAEREG_stufe19.md).

P2  Consistency across cohorts at the gene level: gene score = mean of expected sign x z over the gene's direction
    cells (gene-specific direction map, pLoF), UK Biobank (Genebass) against Genes & Health (REGENIE gene burden, mask A).
    Spearman correlation, one-sided permutation p (B = 2000).
P3  Yeast-human at the gene level: one-to-one best-hit orthologs (Alliance), yeast gene diagonal of Hoepfner et al.
    (script 40) and of Hillenmeyer et al. (computed as in script 41) against the UK Biobank gene score.
    Spearman correlation, one-sided permutation p (B = 5000).

Inputs:  direction_map.py inputs, config.GH_BURDEN_DIR, reference/alliance_orthology.tsv.gz, reference/SGD_features.tab,
         reference/hgnc.txt, 40_yeast_hoepfner_genes.tsv, config.YEAST_HILLENMEYER_DIR
Outputs: 45_cross_cohort_and_yeast_human_consistency.tsv
Usage:   python 45_cross_cohort_and_yeast_human_consistency.py
"""
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import config
import replication_tests as rt
from direction_map import build_direction_map
from yeast_data import read_hillenmeyer_z

dm = build_direction_map()
rng = np.random.default_rng(45)
rows = []

# ---------------- P2
print("\n=== P2 consistency across cohorts (gene score = mean of s*z over direction cells, pLoF) ===")
ukb = rt.ukb_gene_scores(dm)
Zg = rt.genes_and_health_burden_z()
ghv = rt.gene_scores(dm.directions, lambda g, t: Zg[t].get(g, np.nan) if t in Zg else np.nan)
rt.cohort_consistency([(ukb, ghv, "UKB", "G&H")], rng, rows)

# ---------------- P3
print("\n=== P3 yeast-human at the gene level (orthologs, best hits in both directions) ===")
o = pd.read_csv(config.REF_DIR / "alliance_orthology.tsv.gz", sep="\t", comment="#", dtype=str)
o = o[(o.Gene1SpeciesTaxonID == "NCBITaxon:9606") & (o.Gene2SpeciesTaxonID == "NCBITaxon:559292")
      & (o.IsBestScore == "Yes") & (o.IsBestRevScore == "Yes")]
sgd = pd.read_csv(config.REF_DIR / "SGD_features.tab", sep="\t", header=None, dtype=str)
sgd = sgd[sgd[1] == "ORF"]
sgd_to_orf = dict(zip("SGD:" + sgd[0], sgd[3]))
o["orf"] = o.Gene2ID.map(sgd_to_orf)
hgnc = pd.read_csv(config.REF_DIR / "hgnc.txt", sep="\t", dtype=str, usecols=["hgnc_id", "ensembl_gene_id"]).dropna()
o = o.merge(hgnc, left_on="Gene1ID", right_on="hgnc_id").dropna(subset=["orf"])
print(f"  ortholog pairs (1:1 best in both directions): {len(o)}, with a direction map: {o.ensembl_gene_id.isin(dm.directions).sum()}")
hoepfner = pd.read_csv(config.results("40_yeast_hoepfner_genes.tsv"), sep="\t").set_index("gene").diagonal

# Hillenmeyer gene diagonal as in script 41: mean centred heterozygous z over the conditions in which the homozygote is hit
het, hom = read_hillenmeyer_z("het"), read_hillenmeyer_z("hom")
by_compound = lambda d_: d_.T.groupby([c.split(":")[1] for c in d_.columns]).mean().T  # noqa: E731
Hs, Os = by_compound(het), by_compound(hom)
compounds = sorted(set(Hs.columns) & set(Os.columns))
genes = sorted(set(Hs.dropna(how="all").index) & set(Os.dropna(how="all").index))
Hs, Os = Hs.loc[genes, compounds], Os.loc[genes, compounds]
Mk = (Os >= 3).values & np.isfinite(Hs.values)
ok = Mk.sum(axis=1) >= 2
Hc = Hs - Hs.mean(axis=0)
Hc = Hc.sub(Hc.mean(axis=1), axis=0)
hillenmeyer = pd.Series([np.nanmean(Hc.values[i][Mk[i]]) for i in np.where(ok)[0]], index=np.array(genes)[ok])

for yeast_name, yeast_values in [("Hoepfner", hoepfner), ("Hillenmeyer", hillenmeyer)]:
    x = o.assign(y=o.orf.map(yeast_values), h=o.ensembl_gene_id.map(ukb)).dropna(subset=["y", "h"]).drop_duplicates("ensembl_gene_id")
    if len(x) < 15:
        print(f"  {yeast_name}-UKB: n {len(x)} too small")
        continue
    r = spearmanr(x.y, x.h).correlation
    null = np.array([spearmanr(x.y.values, rng.permutation(x.h.values)).correlation for _ in range(5000)])
    p = (1 + (null >= r).sum()) / 5001
    print(f"  {yeast_name:11s} ~ human UKB: rho {r:+.3f} (n {len(x)}), p {p:.4f}")
    rows.append(dict(test=f"P3 {yeast_name}-UKB", n=len(x), stat=r, p=p))

pd.DataFrame(rows).to_csv(config.results("45_cross_cohort_and_yeast_human_consistency.tsv"), sep="\t", index=False, float_format="%.4f")
