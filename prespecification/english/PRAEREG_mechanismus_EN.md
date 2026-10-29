> English translation (21 October 2026) of `prespecification/PRAEREG_mechanismus.md`. The German original is the pre-specification document and is unchanged; this translation was made after the analyses and is provided for readability only.

# Pre-specification: mechanism rule across organs (StilleTraeger, stage 4)

Specified on 19 October 2026, **before** any class analysis outside the skeleton was computed. After this, nothing is changed in hypothesis, classes, traits or test; deviations are marked as such.

**What was known beforehand (honest):** for the skeleton, the rule is derived from the data (matrix/secretion dosage-sensitive, cilium/lysosome silent; scripts 05/10/11). I have seen the 25 strongest genes of the organ diagonal (among them LCAT, LIPC, ABCG5, TMPRSS6, TFR2, CLDN10, SLC22A12, ABCC2). I have not computed class membership or class comparisons outside the skeleton.

## Hypothesis
Recessivity "holds" (carriers silent) when the gene product works in a buffered flux with reserve capacity: catabolic/lysosomal enzymes, ciliary machinery. Recessivity "breaks" (carriers show an attenuated form of the disease) when the gene product is made, secreted or structurally incorporated in large amounts: secretory and structural high throughput.

## Classes (GO, descendants included, no manual selection; order = precedence for multiple membership)
1. **SECRETED** = extracellular region (GO:0005576, CC) OR structural molecule activity (GO:0005198, MF) OR ER-to-Golgi transport (GO:0006888) OR ECM organization (GO:0030198) OR collagen/GAG/proteoglycan biosynthesis (GO:0032964, GO:0006024, GO:0030166)
2. **CATABOLIC** = lysosome (GO:0005764, CC) OR catabolic process (GO:0009056, BP)
3. **CILIUM** = cilium (GO:0005929, CC)
4. **other**

## Genes, organs, traits
Exclusively recessive Mendelian genes with an organ profile (HPO fraction ≥ 0.2; as in scripts 14/20), traits per organ as in script 14 (all assigned traits). Genebass UKB, carriers ≥ 10, keep flags. Organ-matched gene value = mean of the double-centred, **locus-adjusted** z² (V3 from script 20) over the matching traits.

## Primary test
- **Data set:** pLoF, **without the skeleton**. Skeletal organ and skeletal traits are removed so that the test is independent of the derivation.
- **Statistic:** mean organ-matched gene value in SECRETED minus CATABOLIC.
- **Significance:** permutation of class labels within the main organ (organ with the largest HPO fraction), B = 10,000, one-sided (SECRETED > CATABOLIC).
- **Prediction:** difference > 0, P < 0.05.

## Secondary
- Missense|LC, same test.
- SECRETED against CILIUM and SECRETED against other.
- Per organ (blood, immune, kidney, liver, lipids, lung, muscle, nervous, circulatory), where ≥ 5 genes per class are available.
- Regression: gene value ~ class + log10 carriers + LOEUF + main organ.
- Skeleton alone as confirmation of the derivation (not independent).
- Kacser–Burns directly: enzyme (EC according to HGNC) against non-enzyme across all organs (expectation from script 09: no difference).

## Interpretation
- **Primary test positive:** the rule generalizes across organs.
- **Primary test negative:** the rule is skeleton-specific (growth plate as a high-throughput matrix organ). This is also a statement, only a smaller one.
