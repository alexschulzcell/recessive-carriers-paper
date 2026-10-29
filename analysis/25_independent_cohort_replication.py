"""Replication of the three core findings in cohorts without UK Biobank participants (undirected, z2).

Cohorts (freely downloadable summary statistics, see data/README.md):
    GH    Genes & Health 44k (exomes of British South Asians; GWAS Catalog, GRCh38, variant level, MAC >= 5): 12 traits
    FMS   FinMetSeq (exomes of 19k Finns; hg19): height, creatinine, eGFR, LDL, HDL, triglycerides
    CHIP  exome-chip consortia before UK Biobank (hg19): GIANT height, Blood Cell Consortium Hb/MCV/RBC/PLT/WBC/lymphocytes,
          CKDGen eGFRcrea/eGFRcys (European ancestry), GLGC LDL/HDL/TG
Gene statistics (pLoF and coding variants, MAF < 1 %): see replication_tests.py.

Tests
    T1  height: recessive skeletal-dysplasia genes (AR_nos) against matched controls (pool) and recessive non-skeletal genes.
    T2  organ diagonal: doubly centred z2, permutation of organ labels, per cohort x variant set.
    T3  pre-specified producer rule (PRAEREG_hauptprodukt.md): W1/W2b classes from script 22 (gene x organ), value from the
        replication traits, permutation within organ, skeletal organ excluded, one-sided.
Supports the independent replication reported in the Results (Genes & Health, FinMetSeq, exome-chip consortia).

Inputs:  config.REPLICATION_DIR (annot, gh, finmetseq, bcx, ckdgen), config.EXOME_CHIP_DIR, refGene (config.REF_DIR),
         gene_annotation.tsv.gz, organ_fractions_recessive.tsv.gz, 22_main_product_units_plof.tsv
Outputs: 25_independent_cohort_replication.tsv; gene statistics cached in cache/replication_gene_stats/
Usage:   python 25_independent_cohort_replication.py
"""
import numpy as np
import pandas as pd

import config
import replication_tests as rt

rng = np.random.default_rng(23)
annotation = pd.read_csv(config.results("gene_annotation.tsv.gz"), sep="\t")
ann = rt.VariantAnnotation()
sources = rt.independent_cohort_sources()
stats = rt.load_gene_stats(sources, ann, config.cache("replication_gene_stats"))

rows = []
rt.t1_height_enrichment(stats, sources, rt.matched_strata_annotation(annotation), rng, rows)
organ_fractions = pd.read_csv(config.results("organ_fractions_recessive.tsv.gz"), sep="\t")
rt.t2_organ_diagonal(stats, sources, organ_fractions[organ_fractions.frac >= 0.2], rt.TRAIT_TO_ORGAN, rng, rows)
rt.t3_producer_rule(stats, sources, rt.load_producer_units(), rt.TRAIT_TO_ORGAN, rng, rows)
pd.DataFrame(rows).to_csv(config.results("25_independent_cohort_replication.tsv"), sep="\t", index=False)
