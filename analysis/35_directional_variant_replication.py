"""Directed replication with variant-level summary statistics (test T10-8b, PRAEREG_vorzeichen.md addendum B).

Sources as in script 25: Genes & Health variant level (GH), FinMetSeq (FMS), exome-chip consortia before UK Biobank (CHIP).
Pre-specified aggregation: z_gene = sum z_i / sqrt(k) (Stouffer), z_i = BETA/SE of the alternative allele (oriented to the
LoF allele for pLoF), |z_i| <= sqrt(50), MAF < 1 %. The gene-specific direction map is the unchanged map of script 31
(direction_map.py); statistic, centring and permutation as in scripts 31/33. Only traits of the direction map are used
(not eGFR; as a sensitivity analysis eGFR -> creat and eGFRcys -> cysc with reversed sign).

Inputs:  as script 25, plus the direction_map.py inputs
Outputs: 35_directional_variant_replication.tsv; signed gene statistics cached in cache/replication_signed_gene_stats/
Usage:   python 35_directional_variant_replication.py
"""
import numpy as np

import config
import replication_tests as rt
from direction_map import build_direction_map

ann = rt.VariantAnnotation()
sources = rt.independent_cohort_sources()
dm = build_direction_map()
rng = np.random.default_rng(35)

signed = rt.load_gene_stats(sources, ann, config.cache("replication_signed_gene_stats"), signed=True)
R = rt.directional_tests(rt.signed_by_cohort_and_set(signed), list(sources), dm, rng)
R.to_csv(config.results("35_directional_variant_replication.tsv"), sep="\t", index=False, float_format="%.4f")
