"""Shared inputs of the IMPC mouse organ-diagonal analyses (scripts 38 and 38b).

The organ -> MP top-level system mapping (ORGAN_TO_MP) is defined in 38_mouse_organ_diagonal.py itself (it is the
pre-specified mapping and is also read by script 47 for Table S3B); script 38b imports it from there.

Inputs in config.IMPC_RELEASE_DIR: statistical-results-ALL.csv.gz (IMPC release, one row per gene x zygosity x
parameter) and HMD_HumanPhenotype.rpt (MGI human-mouse orthologs).
"""
import numpy as np
import pandas as pd

import config

STATISTICAL_RESULTS = config.IMPC_RELEASE_DIR / "statistical-results-ALL.csv.gz"


def one_to_one_orthologs():
    """1:1 human-mouse orthologs (MGI HMD) with the human gene_id: columns human, mouse, gene_id, symbol."""
    H = pd.read_csv(config.IMPC_RELEASE_DIR / "HMD_HumanPhenotype.rpt", sep="\t", header=None, usecols=[0, 2],
                    names=["human", "mouse"], dtype=str).dropna()
    H = H[~H.human.duplicated(keep=False) & ~H.mouse.duplicated(keep=False)]
    A = pd.read_csv(config.results("gene_annotation.tsv.gz"), sep="\t", usecols=["gene_id", "symbol"]).drop_duplicates("symbol")
    return H.merge(A, left_on="human", right_on="symbol")


def organ_membership():
    """{gene_id: set of organs} of the recessive genes (organ fraction >= 0.2, as in scripts 14/20)."""
    OF = pd.read_csv(config.results("organ_fractions_recessive.tsv.gz"), sep="\t")
    OF = OF[OF.frac >= 0.2]
    return OF.groupby("gene_id").organ.apply(set)


def organ_cell_mask(labels, systems, V, organ_to_mp):
    """Boolean genes x systems mask of the cells that match the organ labels of each gene (finite cells of V only).

    labels: one set of organs per row of V; systems: column names of V (MP top-level systems).
    """
    m = np.zeros(V.shape, bool)
    for i, organs in enumerate(labels):
        for o in organs:
            for s in organ_to_mp[o]:
                if s in systems:
                    m[i, systems.index(s)] = True
    return m & np.isfinite(V)
