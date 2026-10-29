> English translation (21 October 2026) of `prespecification/PRAEREG_gnomad.md`. The German original is the pre-specification document and is unchanged; this translation was made after the analyses and is provided for readability only.

# Pre-specified: does the rule predict selection pressure – without a single phenotype? (stage 11, `work/34_*`)

Specified on 20 October 2026, **before** any computation. Internal specification, no public time stamp.

## Question
The rule from stages 5/9A states: half the dose shows in the carrier where the gene works **organ-specifically** in its disease organ. If this is true, exactly these recessive genes should also be **under stronger selection against heterozygous LoF** in the population. This can be tested from allele frequencies in gnomAD, without any UKB traits – a third, independent level of measurement (population genetics).

## Unit and data
- Genes: exclusively recessive genes from `organ_anteile_AR.tsv.gz` (as in stage 5), **main organ** = organ with the largest HPO fraction, fraction ≥ 0.2.
- Without skeleton (HPA has no cartilage/bone; as in scripts 22/27). Organ → tissue and organ → producer cell exactly as in `22_hauptprodukt.py` (TIS) and `27_zelltyp_W4.py` (PROD).
- Specificity in the disease organ: `spez` = log2((nTPM in the organ tissue + 1)/(median over all tissues + 1)) as in script 22; cell type `zspez` as in script 27.
- General tissue specificity as confounder: `tau` (Yanai) over all HPA tissues, plus log10 median expression.
- Only genes with exp_lof ≥ 10 (LOEUF otherwise unreliable).

## Outcomes
- **Primary: LOEUF** (gnomAD v2.1.1, `oe_lof_upper`). Based purely on LoF frequency, contains no expression features. gnomAD v2 contains no UKB.
- Secondary: s_het (GeneBayes, Zeng 2024). **Caution, circularity:** the GeneBayes prior uses gene features, including expression. Secondary outcome only.
- Secondary: point estimate obs/exp LoF.

## Statistics
- Regression: rank(−LOEUF) ~ spez + tau + log10 median expression + log10 CDS length + organ (fixed). Test statistic: coefficient of `spez`.
- Null: permute organ labels between genes (B = 10,000), recompute `spez`, re-estimate the coefficient. One-sided P (expectation positive).
- Primary criterion: tissue `spez`, LOEUF, P < 0.05.
- Secondary: cell-type `zspez` (same logic); s_het; and, as a bridge to the phenotype, the correlation of the carrier effect (locus-adjusted organ value from script 22) with −LOEUF within organs.

## Interpretation, specified in advance
- Positive: the rule predicts from expression alone which recessive genes are treated "as semi-dominant" by selection. Connects stage 5 with Judd 2026/Fridman 2021.
- Null: the rule describes the phenotype, but selection follows other pleiotropy (consistent with stage 2, test 18). This is also an answer.
- Negative (specific genes less constrained): argues for the known rule "tissue-specific = less essential", which we would then not have adjusted for sufficiently.

---

# Addendum 20 October 2026: direction × specificity (stage 12, `work/36`, `work/37`)

**Post hoc finding in UKB (from the atlas, script 36):** in cells with |z| > 2, carrier effects in tissue-specificity tertile 3 point towards disease in 0.677 (n 257), in tertile 1 in 0.513 (n 199); permutation of tertiles within organ P 0.0002; synonymous flat. This was **not** pre-specified.

**Pre-specified replication in Genes & Health 44k (before computing):** the same direction map (script 31), the same specificity tertiles (from `ATLAS_traeger.tsv`, unchanged), signed REGENIE z as in script 33. Statistic: fraction in the disease direction in tertile 3 minus tertile 1, cells with |z| > 2; null by permutation of tertiles within organ (B = 5,000), one-sided. Primary mask C (pLoF + missense, more cells), mask A second row, mask D (synonymous) negative control. Criterion: mask C P < 0.05 and synonymous not significant. Expectation: same direction, but clearly fewer cells because of 33k; "not replicated" will be reported as such.
