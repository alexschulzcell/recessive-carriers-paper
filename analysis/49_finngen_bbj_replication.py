"""Pre-specified replication in FinnGen (F1-F8, PRAEREG_finngen.md) and BioBank Japan (B1-B8, PRAEREG_bbj.md).

The methods are those of the earlier replication analyses, applied unchanged to one cohort (replication_tests.py):
    1-3  (F1-F3 / B1-B3)  tests T1-T3 of script 25 (undirected z2: height enrichment, organ diagonal, producer rule W1/W2b)
    4    (F4 / B4)        test T10-8b of script 35 (gene-specific direction, Stouffer gene z); a second line combines the
                          result with Genes & Health (script 35, no sample overlap) by Fisher's method
    5    (F5 / B5)        P1: frozen atlas prediction (Table S5)
    6    (F6 / B6)        P2: gene-level consistency with UK Biobank and with Genes & Health
    7    (F7 / B7)        P4: skeletal mechanism for height (matrix genes vs cilium + lysosome genes)
    8    (F8 / B8)        P5: dominance coefficient h for height, with the carrier effect beta_het per gene taken as the
                          inverse-variance weighted mean of the pLoF variant betas (oriented to the LoF allele)
Only the reader is new: the pre-filtered files written by 49a (FinnGen) and 49b (BioBank Japan).
Trait -> organ assignment a priori as in scripts 14/19.

Inputs:  config.REPLICATION_DIR/finngen or /bbj (<trait>.tsv.gz), inputs of scripts 25 and 35, the frozen prediction,
         literature_patient_values.tsv, 29b_dominance_h_main_test.tsv, 35_directional_variant_replication.tsv, config.GH_BURDEN_DIR
Outputs: 49_<cohort>_1-3.tsv, 49_<cohort>_4.tsv, 49_<cohort>_5-8.tsv (cohort = finngen or biobank_japan);
         gene statistics cached in cache/replication_gene_stats_<fg|bbj>/ and cache/replication_signed_gene_stats_<fg|bbj>/
Usage:   python 49_finngen_bbj_replication.py FG      (FinnGen)
         python 49_finngen_bbj_replication.py BBJ     (BioBank Japan)
"""
import sys

import numpy as np
import pandas as pd
from scipy.stats import combine_pvalues

import common
import config
import replication_tests as rt
from direction_map import build_direction_map

COHORT = sys.argv[1]
BASE_TRAITS = ["height", "hb", "mcv", "rbc", "plt", "wbc", "lymph", "creat", "urate", "urea", "phosphate", "alt", "ast", "ggt",
               "bili", "albumin", "ldl", "hdl", "tg", "chol", "glucose", "hba1c"]
SETTINGS = {"FG": ("finngen", "finngen", "hg38", BASE_TRAITS + ["cysc"], "F"),
            "BBJ": ("bbj", "biobank_japan", "hg19", BASE_TRAITS + ["sbp"], "B")}
folder, out_name, BUILD, TRAITS, test_prefix = SETTINGS[COHORT]
cohort_dir = config.REPLICATION_DIR / folder
missing = [k for k in TRAITS if not (cohort_dir / f"{k}.tsv.gz").exists()]
assert not missing, f"missing files: {missing}"
print(f"Cohort {COHORT}: {folder}, {BUILD}, {len(TRAITS)} traits", flush=True)

# trait -> organ, a priori as in scripts 14/19 (PRAEREG_finngen.md, PRAEREG_bbj.md)
TRAIT_TO_ORGAN = {"height": "skeletal", "hb": "blood", "mcv": "blood", "rbc": "blood", "plt": "blood", "wbc": "immune",
                  "lymph": "immune", "creat": "kidney", "cysc": "kidney", "urate": "kidney", "urea": "kidney", "phosphate": "kidney",
                  "alt": "liver", "ast": "liver", "ggt": "liver", "bili": "liver", "albumin": "liver", "ldl": "lipids", "hdl": "lipids",
                  "tg": "lipids", "chol": "lipids", "glucose": "glucose", "hba1c": "glucose", "sbp": "circulatory"}
sources = {COHORT: [(k, BUILD, lambda k=k: rt.read_finngen_bbj(cohort_dir, k)) for k in TRAITS]}
cache_suffix = COHORT.lower()
ann = rt.VariantAnnotation()
annotation = pd.read_csv(config.results("gene_annotation.tsv.gz"), sep="\t")

# ======================= 1-3: tests T1-T3 of script 25 =======================
print(f"\n########## {COHORT} {test_prefix}1-{test_prefix}3 (method of script 25: height, organ diagonal z2, rule W1/W2b) ##########", flush=True)
rng = np.random.default_rng(23)
stats = rt.load_gene_stats(sources, ann, config.cache(f"replication_gene_stats_{cache_suffix}"))
rows = []
rt.t1_height_enrichment(stats, sources, rt.matched_strata_annotation(annotation), rng, rows)
organ_fractions = pd.read_csv(config.results("organ_fractions_recessive.tsv.gz"), sep="\t")
rt.t2_organ_diagonal(stats, sources, organ_fractions[organ_fractions.frac >= 0.2], TRAIT_TO_ORGAN, rng, rows)
rt.t3_producer_rule(stats, sources, rt.load_producer_units(), TRAIT_TO_ORGAN, rng, rows)
pd.DataFrame(rows).to_csv(config.results(f"49_{out_name}_1-3.tsv"), sep="\t", index=False)

# ======================= 4: test T10-8b of script 35 =======================
print(f"\n########## {COHORT} {test_prefix}4 (method of script 35: gene-specific direction, Stouffer) ##########", flush=True)
dm = build_direction_map()
rng = np.random.default_rng(35)
signed = rt.signed_by_cohort_and_set(
    rt.load_gene_stats(sources, ann, config.cache(f"replication_signed_gene_stats_{cache_suffix}"), signed=True))
R4 = rt.directional_tests(signed, [COHORT], dm, rng)
R4.to_csv(config.results(f"49_{out_name}_4.tsv"), sep="\t", index=False, float_format="%.4f")
gh35 = pd.read_csv(config.results("35_directional_variant_replication.tsv"), sep="\t")
print("\nFisher combination with Genes & Health only (variant based, script 35; no sample overlap):")
for s in ("plof", "coding"):
    a = R4[(R4.analysis == "primary") & (R4.set == s)].p
    b = gh35[(gh35.analysis == "primary") & (gh35.cohort == "GH") & (gh35.set == s)].p
    if len(a) and len(b):
        print(f"  {s}: {COHORT} p {a.iloc[0]:.4f}, G&H p {b.iloc[0]:.4f} -> Fisher p {combine_pvalues([a.iloc[0], b.iloc[0]], method='fisher')[1]:.2g}")

# ======================= 5-8: P1, P2, P4, P5 =======================
print(f"\n########## {COHORT} {test_prefix}5-{test_prefix}8 (P1, P2, P4, P5) ##########", flush=True)
cohort_z = {"plof": signed[(COHORT, "plof")], "plofmis": signed[(COHORT, "coding")]}   # plofmis = pLoF + coding variants

# beta_het per gene: pLoF variants with MAF < 1 %, oriented to the LoF allele (as in script 35), inverse-variance weighted
lofd = ann.lof_dict(BUILD)
d = rt.read_finngen_bbj(cohort_dir, "height")
d = d[np.isfinite(d.z) & (d.maf < 0.01)]
acc = {}
for ch, pos, ref, alt, b, se in zip(d.ch.astype(str), d.pos.astype(int), d.ref.astype(str), d.alt.astype(str), d.beta, d.se):
    g = lofd.get((ch, pos, ref, alt))
    sign = 1
    if g is None and lofd.get((ch, pos, alt, ref)):
        g = lofd[(ch, pos, alt, ref)]
        sign = -1
    if g:
        w = 1 / se ** 2
        a = acc.setdefault(g, [0.0, 0.0])
        a[0] += w * sign * b
        a[1] += w
beta_het = {g: s / w for g, (s, w) in acc.items()}
print(f"  beta_het for {len(beta_het)} genes with pLoF variants (height)")

rng = np.random.default_rng(45)
rows = []
print(f"=== {COHORT} {test_prefix}5 (P1) frozen prediction ===")
rt.frozen_prediction_test(common.load_frozen_prediction(), cohort_z, rng, rows)

print(f"\n=== {COHORT} {test_prefix}6 (P2) consistency across cohorts (gene score = mean of s*z over direction cells, pLoF) ===")
ukb = rt.ukb_gene_scores(dm)
cohort_scores = rt.gene_scores(dm.directions, lambda g, t: cohort_z["plof"][t].get(g, np.nan) if t in cohort_z["plof"] else np.nan)
Zg = rt.genes_and_health_burden_z()
gh_scores = rt.gene_scores(dm.directions, lambda g, t: Zg[t].get(g, np.nan) if t in Zg else np.nan)
rt.cohort_consistency([(ukb, cohort_scores, "UKB", COHORT), (cohort_scores, gh_scores, COHORT, "G&H")], rng, rows)

print(f"\n=== {COHORT} {test_prefix}7 (P4) skeletal mechanism (height) ===")
rt.skeletal_mechanism_test({"plof": signed[(COHORT, "plof")]["height"], "coding": signed[(COHORT, "coding")]["height"]}, rows)

print(f"\n=== {COHORT} {test_prefix}8 (P5) h (height, beta_het from pLoF betas) ===")
rt.dominance_h_height(beta_het, rows)
pd.DataFrame(rows).to_csv(config.results(f"49_{out_name}_5-8.tsv"), sep="\t", index=False, float_format="%.4f")
print("\ndone")
