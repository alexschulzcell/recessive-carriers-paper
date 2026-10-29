> English translation (21 October 2026) of `prespecification/PRAEREG_vorzeichen.md`. The German original is the pre-specification document and is unchanged; this translation was made after the analyses and is provided for readability only.

# Pre-specification 4 (stage 10): directed diagonal test

Specified on 19 October 2026, **before** any computation. Known: stages 1–9; that the sign convention of Genebass `BETA_Burden` is biologically correct (checked on ACAN, IHH, NPR2, FBN1, LDLR, PCSK9, APOB, LCAT, ABCA1, CETP, UGT1A1, TMPRSS6, HBB — signs of single genes only, no group values). Not known: any value of this test.

## Why

All main analyses so far (stages 1–9) use z², i.e. they are **directionless**. The claim "carriers show a mild form of the same disease", however, requires the **correct sign**: carriers of a dysplasia gene must be *shorter*, not merely different. Without sign, the diagonal is also consistent with "carriers of recessive genes are more variable in the organ trait" — a much weaker claim. This test closes the gap.

## Data

Genebass (UKB 394k, inverse-rank-normal transformed), signed z = sign(BETA_Burden) · |z| from `Pvalue_Burden`; filters as in scripts 14/20 (keep, carriers ≥ 10, z² capped at 50). Organ assignment from `daten/organ_anteile_AR.tsv.gz` (HPO fraction ≥ 0.2), unchanged from stage 2. Centring as in script 14: first trait mean, then gene mean, **in the raw signed z**, only then multiplied by the expected sign.

Directed value per cell: `d = s_expected · z_centred`. d > 0 means: the carrier shifts towards the disease.

## A) Expected sign per organ × trait (main test, unambiguous cells)

Derived from the disease biology of the organ's recessive diseases, not from data. s = −1 means "carriers expected lower".

| Organ | Trait | s | Rationale (class of recessive disease) |
|---|---|---|---|
| Blood | hb | −1 | anaemias (thalassaemia, Fanconi, DBA, CDA, spherocytosis) |
| Blood | rbc | −1 | the same |
| Blood | plt | −1 | congenital thrombocytopenias |
| Immune | wbc | −1 | SCID, agammaglobulinaemia, neutropenias, CGD |
| Immune | lymph | −1 | the same |
| Kidney | creat | +1 | nephrotic syndrome, nephronophthisis, cystic kidneys: GFR falls |
| Kidney | cysc | +1 | the same |
| Kidney | urea | +1 | the same |
| Liver | alt | +1 | cholestasis, metabolic liver disease, Wilson: hepatocyte damage |
| Liver | ast | +1 | the same |
| Liver | ggt | +1 | the same |
| Liver | bili | +1 | cholestasis and conjugation defects (UGT1A1, ABCC2, ABCB11) |
| Liver | albumin | −1 | synthesis disorder |
| Lipids | hdl | −1 | LCAT, ABCA1 (Tangier), APOA1 deficiency |
| Lipids | tg | +1 | chylomicronaemia (LPL, APOC2, GPIHBP1, LMF1) |
| Lung | fev1 | −1 | CF, ciliary dyskinesia, surfactant defects, ILD |
| Lung | fvc | −1 | the same |
| Muscle | grip | −1 | myopathies, muscular dystrophies |
| Nervous | fluidint | −1 | intellectual disability, neurodegeneration |
| Nervous | reaction | +1 | reaction *time*: worse = longer = higher value |
| Skeletal | height | −1 | skeletal dysplasias: short stature |
| Skeletal | sitheight | −1 | the same |

22 unambiguous cells across 10 organs.

### Declared in advance as mixed-direction (NOT in the main test, reported separately)

Here, recessive diseases exist in both directions; an a priori direction would be a guess.

| Organ | Trait | Tendency | Counterexample |
|---|---|---|---|
| Blood | mcv | −1 (microcytic) | megaloblastic (TCN2, MTR, folate transport) → +1 |
| Kidney | urate | +1 (GFR falls) | hypouricaemia SLC22A12, SLC2A9 → −1 |
| Kidney | phosphate | −1 (tubular loss) | CKD hyperphosphataemia → +1 |
| Lipids | ldl | −1 (hypobetalipoproteinaemia, MTTP, APOB) | LDLRAP1, recessive FH → +1 |
| Lipids | chol | −1 | the same split |
| Glucose | glucose | +1 (neonatal diabetes) | hyperinsulinism, glycogen storage diseases → −1 |
| Glucose | hba1c | +1 | the same |
| Skeletal | ebmd | −1 (osteogenesis imperfecta) | osteopetrosis, sclerosing dysplasias → +1 |
| Circulatory | sbp | −1 (cardiomyopathy) | weak proxy, both possible |

## B) Tests and criteria

**T10-1 main test (pLoF, unambiguous cells).** Statistic = mean of d over all matching gene × trait cells. Null distribution: permutation of organ labels between genes (B = 1000, multiple membership preserved, as in scripts 14/20). One-sided, criterion **P < 0.05**.

**T10-2 replication.** The same with missense (mislc), criterion P < 0.05.

**T10-3 negative control.** The same with synonymous. Expectation: null. P < 0.05 for synonymous is reported as a warning, and the main test is then read only locus-adjusted.

**T10-4 sign fraction (illustrative).** Fraction of matching cells with d > 0, binomial test against 0.5, two-sided. Prediction: > 0.5. Plus the same measure for non-matching cells (expectation 0.5).

**T10-5 hardening, pre-specified as in stage 3.**
- locus-adjusted: cells with |z_syn| > 1.96 discarded;
- without the top 10 and the top 25 genes (by |z| in the matching trait);
- one fixed primary trait per organ (hb, wbc, creat, bili, hdl, fev1, grip, fluidint, height);
- leave-one-organ-out over the 10 organs with unambiguous cells.

Criterion: the sign remains positive and P < 0.05 both locus-adjusted and with a fixed primary trait.

**T10-6 mixed cells (separately).** Main-test statistic over the 9 mixed cells only, with the tendency stated above. No prediction; purely descriptive.

## C) Gene-specific expected direction from HPO (secondary test)

The organ direction in A is coarse. More precise: each gene receives its **own** direction from the HPO terms of its OMIM diseases — described clinically and independent of UKB.

**Term → trait × direction, specified in advance** (HPO term present for the gene → this expectation):

- **Blood/immune:** HP:0001903 Anemia → hb −1; HP:0001935 Microcytic anemia → mcv −1, hb −1; HP:0001889 Macrocytic anemia → mcv +1, hb −1; HP:0001873 Thrombocytopenia → plt −1; HP:0001894 Thrombocytosis → plt +1; HP:0001882 Leukopenia → wbc −1; HP:0001875 Neutropenia → wbc −1; HP:0001974 Leukocytosis → wbc +1; HP:0001888 Lymphopenia → lymph −1
- **Kidney:** HP:0003259 Elevated circulating creatinine → creat +1, cysc +1; HP:0012622 Chronic kidney disease → creat +1, cysc +1, urea +1; HP:0002905 Hyperphosphatemia → phosphate +1; HP:0002148 Hypophosphatemia → phosphate −1; HP:0002149 Hyperuricemia → urate +1; HP:0003537 Hypouricemia → urate −1
- **Liver:** HP:0002904 Hyperbilirubinemia, HP:0002908 Conjugated hyperbilirubinemia, HP:0008282 Unconjugated hyperbilirubinemia, HP:0003573 Increased total bilirubin → bili +1; HP:0002910 Elevated hepatic transaminase → alt +1, ast +1; HP:0003234 Elevated GGT → ggt +1; HP:0003073 Hypoalbuminemia → albumin −1
- **Lipids:** HP:0003233 Decreased HDL cholesterol → hdl −1; HP:0002155 Hypertriglyceridemia → tg +1; HP:0003124 Hypercholesterolemia → chol +1, ldl +1; HP:0003563 Hypobetalipoproteinemia, HP:0008181 Abetalipoproteinemia → ldl −1, chol −1
- **Glucose:** HP:0003074 Hyperglycemia, HP:0000857 Neonatal insulin-dependent diabetes mellitus → glucose +1, hba1c +1; HP:0001943 Hypoglycemia → glucose −1
- **Lung:** HP:0002097 Emphysema, HP:0006536 Airway obstruction, HP:0002113 Pulmonary infiltrates, HP:0002206 Pulmonary fibrosis, HP:0002205 Recurrent respiratory infections → fev1 −1, fvc −1
- **Muscle:** HP:0003701 Proximal muscle weakness, HP:0003324 Generalized muscle weakness, HP:0003202 Skeletal muscle atrophy → grip −1
- **Nervous:** HP:0001249 Intellectual disability, HP:0001263 Global developmental delay, HP:0002376 Developmental regression → fluidint −1, reaction +1
- **Skeletal:** HP:0004322 Short stature, HP:0003510 Severe short stature, HP:0008873 Disproportionate short-limb short stature, HP:0003498 Disproportionate short stature → height −1, sitheight −1; HP:0000098 Tall stature → height +1, sitheight +1; HP:0004349 Reduced bone mineral density, HP:0000939 Osteoporosis → ebmd −1; HP:0011001 Increased bone mineral density, HP:0011002 Osteopetrosis → ebmd +1

Directionless parent terms (HP:0003119 Abnormal circulating lipid concentration, HP:0003116 Abnormal circulating cholesterol, HP:0012100 Abnormal circulating creatinine and similar) are discarded. Terms are counted **with their HPO descendants**. If a gene has terms of both directions for the same trait, this cell is dropped (conflicting).

**T10-7:** main-test statistic on these gene-specific cells, pLoF, permutation of gene labels. Criterion P < 0.05, missense replication P < 0.05, synonymous null.

**Prediction:** the effect (dz) is **larger** than in T10-1 because the direction is more precise. If it is smaller, this is reported as an argument for the coarser organ level, not concealed.

## D) What makes the test fail

A mean d ≤ 0 or P ≥ 0.05 in T10-1 means: the diagonal is directionless, i.e. only "more variance in the matching organ" and **not** "mild form of the disease". This would refute the claim of the title and would be reported as such; the z² results of stages 1–9 would remain valid, but their interpretation would have to be rewritten.

---

# Addendum (stage 10b): directed replication outside UKB — T10-8

Specified on 19 October 2026, **before** any computation. Known: the results of T10-1 to T10-7 in UKB; that the replication tables stored in stages 6–8 **contain only Σz² and therefore no signs**, so the signs have to be retrieved again from the source data. Not known: any value of this addendum.

**The direction map is not touched.** The same gene-specific HPO map as in part C applies, unchanged, including the rule that conflicting cells are dropped. No traits, genes or directions are chosen anew. This is therefore a true replication without degrees of freedom.

## A) Main source: Genes & Health 44k, REGENIE burden (gene level, signed)

Source as in stage 8 (`gs://genesandhealth_publicdatasets/results_44k_ExWAS`, CC BY-SA 4.0), 33k British South Asians, no UKB overlap. Masks A (LOFTEE high-confidence pLoF), C (pLoF + missense), D (synonymous), test ADD, AAF < 1%, allele count ≥ 5 — all as in stage 8. The only difference: instead of (BETA/SE)², **sign(BETA) · |BETA/SE|** is kept.

Centring and statistic exactly as in T10-7: trait mean, then gene mean in the raw signed z, then multiplied by the gene-specific direction; null distribution by permutation of direction maps between genes (B = 1000), one-sided.

**T10-8a:** criterion **P < 0.05** with mask A (pLoF). Mask C (pLoF + missense) as second row, mask D (synonymous) as negative control, expectation null.

**Expectation in advance:** stage 8 showed that pLoF in 33k reaches about 8% of the UKB strength and that the diagonal and rule replicated there only with missense. **It is therefore to be expected that mask A is too weak and only mask C reaches the criterion.** This case is defined in advance as "partially replicated", neither as success nor as failure.

## B) Secondary sources: variant level, with a pre-specified sign aggregation

Sources as in stage 7: FinMetSeq (19k Finns), pre-UKB exome-chip consortia (GIANT height, BCX blood count, CKDGen, GLGC lipids), Genes & Health variant level. All free, without login, to be downloaded again (`work/24_download_replikation.sh`).

**Aggregation per gene × trait (new, specified here):** for the qualifying variants of a gene (MAF < 1%, assignment as in stage 7; pLoF = gnomAD LOFTEE high confidence, coding = CDS), `z_gene = Σ z_i / √k` is formed (Stouffer), where z_i = BETA_i / SE_i refers to the **alternative allele** and k is the number of variants. Rationale: for LoF variants the alternative allele is the damaging one, so signs are oriented consistently across variants. Capping of z_i at ±√50 as everywhere.

**Known limitation, stated in advance:** Stouffer without LD adjustment is anticonservative in magnitude for correlated variants. This affects the **strength**, not the direction; the test is read against the permutation of direction maps, which carries the same correlation structure.

**T10-8b:** the same statistic, per source, criterion P < 0.05; additionally Fisher's combination across sources. Synonymous or non-coding variants as negative control where available.

## C) What counts as replication

- **Replicated:** T10-8a with mask A or C reaches P < 0.05, synonymous negative control remains null.
- **Partially replicated:** only mask C, or only the secondary sources.
- **Not replicated:** no mean > 0 with P < 0.05 in any source. The direction finding is then considered UKB-specific and reported as such; T10-7 would remain valid, but without external support.
