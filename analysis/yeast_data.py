"""Readers for the yeast chemogenomic screens (scripts 40, 40b, 40c, 41, 41b, 45).

Hoepfner et al. 2014 (config.YEAST_HOEPFNER_DIR): HIP_scores.txt / HOP_scores.txt, genes x experiments.
Hillenmeyer et al. 2008 (config.YEAST_HILLENMEYER_DIR): het.z_result_nm.pub / hom.z_result_nm.pub, genes x experiments.

Reading the large text files is slow, so a Parquet copy is read instead when it is present next to the text file
(HIP_scores.parquet, HOP_scores.parquet, het.z_result_nm.parquet, hom.z_result_nm.parquet). The Parquet copies contain
exactly the tables that the text fallback below produces (identical values, row and column order), so both routes give
identical results.
"""
import pandas as pd

import config


def read_hoepfner_scores(screen):
    """Hoepfner score table for screen "HIP" (heterozygous) or "HOP" (homozygous); index = systematic gene name."""
    parquet = config.YEAST_HOEPFNER_DIR / f"{screen}_scores.parquet"
    if parquet.exists():
        return pd.read_parquet(parquet)
    return pd.read_csv(config.YEAST_HOEPFNER_DIR / f"{screen}_scores.txt", sep="\t", index_col=0,
                       low_memory=False).apply(pd.to_numeric, errors="coerce")


def read_hillenmeyer_z(zygosity):
    """Hillenmeyer z-score table for "het" or "hom"; index = ORF (rows of the same ORF averaged)."""
    parquet = config.YEAST_HILLENMEYER_DIR / f"{zygosity}.z_result_nm.parquet"
    if parquet.exists():
        return pd.read_parquet(parquet)
    d = pd.read_csv(config.YEAST_HILLENMEYER_DIR / f"{zygosity}.z_result_nm.pub", sep="\t", index_col=0, low_memory=False)
    d.index = d.index.str.split(":").str[0]          # "YAL001C:chr1_1" -> "YAL001C"
    return d.groupby(level=0).mean()
