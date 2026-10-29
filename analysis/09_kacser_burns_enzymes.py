"""Kacser-Burns test: are carriers of recessive dysplasia genes that encode enzymes silent, while structural, matrix,
secretion and signalling genes show a heterozygous effect?

Hypothesis (Kacser & Burns 1981, Genetics 97:639): enzyme defects are mostly recessive because metabolic flux is buffered
against a loss of 50 % activity. Prediction: among the AR dysplasia genes, the carrier effect on height is stronger for
non-enzymes than for enzymes. Enzyme = HGNC field enzyme_id (EC number) present. Comparison within AR_nos (same class of
disease severity), with a check that the difference is not explained by carrier count or LOEUF
(logistic regression sig ~ enzyme + log10(carriers) + LOEUF). Exploratory.

Inputs:  gene_annotation.tsv.gz (script 03); config.REF_DIR / hgnc.txt; config.GENEBASS_DIR / genebass_<trait>_<set>.tsv.gz
Output:  enzyme_genes.tsv (gene_id, enzyme = 0/1) in RC_RESULTS_DIR; results are printed.
Usage:   python 09_kacser_burns_enzymes.py
"""
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import norm, fisher_exact, mannwhitneyu

import config

h = pd.read_csv(config.REF_DIR / "hgnc.txt", sep="\t", dtype=str, usecols=["ensembl_gene_id", "enzyme_id", "gene_group"])
enzymes = set(h[h.enzyme_id.notna()].ensembl_gene_id)
A = pd.read_csv(config.results("gene_annotation.tsv.gz"), sep="\t")
A = A[A.protein_coding].copy()
A["enzyme"] = A.gene_id.isin(enzymes).astype(int)
A[["gene_id", "enzyme"]].to_csv(config.results("enzyme_genes.tsv"), sep="\t", index=False)
for ph in ["height", "sitheight", "ebmd"]:
    for s in ["plof", "mislc", "syn"]:
        g = pd.read_csv(config.GENEBASS_DIR / f"genebass_{ph}_{s}.tsv.gz", sep="\t")[["gene_id", "BETA_Burden", "Pvalue_Burden"]].merge(A, on="gene_id")
        g = g[(g[f"keep_{s}"] == True) & (g[f"carriers_{s}"] >= 10) & g.Pvalue_Burden.notna()].copy()  # noqa: E712
        g["z"] = np.sign(g.BETA_Burden) * norm.isf(g.Pvalue_Burden.clip(1e-300) / 2)
        g["sig"] = (g.Pvalue_Burden < 0.05).astype(int)
        g["z2"] = g.z ** 2
        g["ltr"] = np.log10(g[f"carriers_{s}"])  # log10 carriers
        line = []
        for st in ["AR_nos", "AR_nonskel", "pool"]:
            x = g[(g.set == st) & g.oe_lof_upper.notna()]
            e, n = x[x.enzyme == 1], x[x.enzyme == 0]
            orr, p = fisher_exact([[n.sig.sum(), len(n) - n.sig.sum()], [e.sig.sum(), len(e) - e.sig.sum()]])
            mw = mannwhitneyu(n.z2, e.z2, alternative="greater").pvalue
            line.append(f"{st}: non-enzyme {n.sig.mean():.3f} (n {len(n)}) / enzyme {e.sig.mean():.3f} (n {len(e)}), "
                        f"OR {orr:.2f} p {p:.3f}, MW-z2 p {mw:.3f}")
            if st == "AR_nos":
                try:
                    m = smf.logit("sig ~ enzyme + ltr + oe_lof_upper", data=x).fit(disp=0)
                    line[-1] += f" | adj. OR enzyme {np.exp(m.params['enzyme']):.2f} (p {m.pvalues['enzyme']:.3f})"
                except Exception as ex:  # noqa: BLE001
                    line[-1] += f" | regression: {ex}"
        print(f"{ph:9s} {s:5s}\n   " + "\n   ".join(line))

x = A[A.set == "AR_nos"]
print("\nEnzymes among AR_nos:", x.enzyme.sum(), "of", len(x), "| among AD_nos:", A[A.set == "AD_nos"].enzyme.sum(),
      "of", (A.set == "AD_nos").sum())
