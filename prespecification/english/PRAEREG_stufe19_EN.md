> English translation (October 21, 2026) of `prespecification/PRAEREG_stufe19.md`, made after the analyses and provided for
> readability only. **Excerpt:** limited to the parts used in the paper (P2 and P3, for the cohorts analyzed there); omitted
> text is marked [...].

# Pre-specified: stage 19 – five hardening tests for Nature Genetics (`work/45_*`)

Specified on 20 October 2026, **before** the respective numbers were computed. [...]

## P1 [...]

## P2 – consistency across cohorts of different ancestry [...]
Gene value = mean of s_expected · z over the cells of the direction map (script 31), pLoF, for each cohort: UKB (Genebass, European), Genes & Health (REGENIE mask A, South Asian), [...]. Spearman of gene values [...] UKB–G&H [...] over genes with data in both; null: permute gene labels in one cohort (B = 2,000), one-sided. Criterion: all three ρ > 0, at least two with P < 0.05.

## P3 – yeast and human at gene level
Human–yeast orthologs from the Alliance (ORTHOLOGY-ALLIANCE_COMBINED, IsBestScore = Yes and IsBestRevScore = Yes), SGD ID → ORF via SGD_features. Yeast value: gene diagonal from `work/40_hefe_gene.tsv` (Hoepfner, primary); Hillenmeyer as second row (gene diagonal from script 41). Human value: gene value from P2 (UKB pLoF) [...]. Only exclusively recessive human genes with a direction map. Test: Spearman yeast value ~ human value, one-sided positive, permutation (B = 5,000). Criterion P < 0.05 with UKB [...]. Interpretation: positive = the same genes are "half-sensitive" in both species. A null result is reported as "analogy, but no gene concordance".

## P4 [...]

## P5 [...]
