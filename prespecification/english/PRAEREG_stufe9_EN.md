> English translation (21 October 2026) of `prespecification/PRAEREG_stufe9.md`. The German original is the pre-specification document and is unchanged; this translation was made after the analyses and is provided for readability only.

# Pre-specification 3 (stage 9): dominance h, cell-type resolution, population level

Specified on 19 October 2026, **before** any computation. Known: stages 1–8 (in particular the W2b results at tissue level); the HPA cell-type list was checked for names only.

## A) Cell-type resolution of the rule (W4)
- **Hypothesis:** the carrier effect on the organ trait follows the specificity in the **producer cell** of the organ at least as well as tissue specificity (W2b).
- **Producer cells** (HPA single cell type, nCPM; maximum over the listed types):
  - Blood: erythrocyte progenitors, megakaryocyte-erythroid progenitors, megakaryocytes
  - Immune: neutrophils, neutrophil progenitors, t-cells, b-cells, monocytes, nk-cells
  - Liver: hepatocytes, cholangiocytes
  - Lipids: hepatocytes, enterocytes
  - Kidney: proximal tubule cells, loop of henle epithelial cells, distal convoluted tubule cells, podocytes
  - Lung: alveolar cells type 1 and 2
  - Muscle: myonuclei
  - Circulatory: cardiomyocytes, vascular smooth muscle cells, vascular endothelial cells
  - Nervous: brain excitatory/inhibitory neurons
  - Glucose: pancreatic islet cells
  - Skeletal: no chondrocytes in HPA → excluded
- **Class W4:** cell specificity = log2((producer + 1) / (median of all cell types + 1)) in the top tertile per organ.
- **Primary test:** pLoF, without skeleton, units and value as in script 22. Permutation within organ, one-sided, P < 0.05; replication missense P < 0.05.
- **Secondary:** regression value ~ z(tissue specificity) + z(cell specificity) + organ. Prediction: the coefficient of cell specificity is > 0 given tissue specificity.

## B) Population level
- **Gene sets per organ:** recessive genes with an organ-matched carrier effect, FDR < 0.05 in Genebass pLoF (within the recessive genes, per trait).
- **Effect** per carrier from the AstraZeneca PheWAS (UKB 500k, Europeans, PTV, s.d.); **carrier frequency** from Genebass (2 × CAF, pLoF).
- **Measures:**
  - fraction of the population with ≥ 1 such pLoF variant;
  - fraction of the lower or upper tail (< −2 s.d. or > +2 s.d.) of the trait attributable to carriers: attributable fraction Σ f·[Φ(−2−β) − Φ(−2)] / Φ(−2);
  - explained variance Σ f·β².
- **Traits:** height (main case), HDL, Hb, MCV, urate, bilirubin.
- **Interpretation:** descriptive, no test. Selection on significance → effects biased upwards (winner's curse). Therefore additionally a conservative variant with AstraZeneca effects at AstraZeneca P < 1 × 10⁻⁴ only.

## C) Dominance coefficient h
- **Patient values:** agent-based literature extraction. Genes selected by data availability (≥ 3 biallelic patients), not by carrier effect.
- **Carrier effect β_het** in s.d.: AstraZeneca PTV (Europeans) where reported (P < 0.1). Otherwise from Genebass: z_pLoF / √(carriers), calibrated over genes with both values (linear regression without intercept).
- **Calculation:** h = β_het / β_hom (both in s.d., same trait).
- **Hypotheses:**
  - **H1:** median h lies clearly below 0.5 (partially recessive). Expectation 0.02–0.2.
  - **H2:** h is higher for main-product genes (W1 OR W2b) than for others; Mann–Whitney, one-sided.
  - **H3:** h correlates positively with s_het (Spearman, one-sided).
- **Validation:** where the literature contains values of heterozygous parents, the literature h is compared with our h.
- **Quality rules:**
  - genes with |β_hom| < 0.5 s.d. are excluded (unstable denominator);
  - h is bounded to [−1, 2];
  - opposite signs (h < 0) are reported, not removed.
