# Pre-specified: independent replication in FinnGen (stage 21, `work/49_*`)

English translation of `PRAEREG_finngen.md`; the German original is authoritative.

Specified on October 28, 2026, **before** any FinnGen trait file was downloaded.
Known: file listing and manifest of the buckets `finngen-public-data-r12/lab_values/` and `finngen-public-data-r13/summary_stats/`,
the Kanta documentation (`FinnGen_kanta_lab_documentation_v2.pdf`, `Kanta_labs_GWAS_results_v2_summary.txt`: lab values, n, units,
outlier filters), and all previous results (UKB, All of Us, Genes & Health, FinMetSeq, exome chip).
Not known: any FinnGen association value (no beta, no P, no variant list).

Reason: the All of Us replication (stages 18/19) used queries of the public All by All browser, whose download the All of Us FAQ does not
permit. It is being re-derived through the Researcher Workbench. FinnGen is freely available; it does not replace All of Us but is a further
independent cohort.

## Data
FinnGen public GWAS summary statistics (registration October 28, 2026; terms: acknowledgment "We want to acknowledge the participants and
investigators of the FinnGen study", citation Kurki et al., Nature 2023). Finnish ancestry, GRCh38, imputed (SISU panel), REGENIE v3.3,
covariates age, sex, genotyping batch, PC1–10; quantitative traits inverse-rank normalized (SD units).
FinnGen applies its own filters: MAC > 30 (lab values) and imputation INFO > 0.6. These filters are kept unchanged.
**No UK Biobank overlap.**

Traits (FinnGen ID; selection rule: the quantitative, rank-normalized trait for the same analyte in serum/plasma or blood, the largest n if
several; chosen from the documentation only):

| Trait | Source | ID | n (documentation) |
|---|---|---|---|
| height | R13 endpoint | HEIGHT_IRN | 364,515 |
| hb | R12 Kanta | 3000963 | 447,271 |
| mcv | R12 Kanta | 3024731 | 447,040 |
| rbc | R12 Kanta | 3026361 | 440,562 |
| plt | R12 Kanta | 3007461 | 441,342 |
| wbc | R12 Kanta | 3010813 | 440,693 |
| lymph | R12 Kanta | 3019198 | 183,481 |
| creat | R12 Kanta | 3020564 | 426,104 |
| cysc | R12 Kanta | 3030366 | 10,125 |
| urate | R12 Kanta | 3026493 | 81,374 |
| urea | R12 Kanta | 3020779 | 31,957 |
| phosphate | R12 Kanta | 3003458 | 46,763 |
| alt | R12 Kanta | 3006923 | 413,996 |
| ast | R12 Kanta | 3013721 | 124,818 |
| ggt | R12 Kanta | 3026910 | 257,813 |
| bili | R12 Kanta | 3006140 (total) | 190,407 |
| albumin | R12 Kanta | 3024561 | 153,962 |
| ldl | R12 Kanta | 3001308 | 368,976 |
| hdl | R12 Kanta | 3023602 | 361,195 |
| tg | R12 Kanta | 3048773 (fasting) | 347,199 |
| chol | R12 Kanta | 3019900 | 365,484 |
| glucose | R12 Kanta | 3018251 (fasting) | 349,417 |
| hba1c | R12 Kanta | 3004410 | 332,429 |

Not available: sitheight, ebmd, sbp, fev1, fvc, grip, fluidint, reaction. Tests that would need these traits run without them.

Download (`work/49a_finngen_download.py`): each file is read as a stream, the SHA-256 of the original file is logged, and only variants with
MAF < 1% that are pLoF or lie in a coding exon (± 2 bp) are stored. Expected columns: `#chrom pos ref alt … beta sebeta af_alt`. If names
differ, columns are matched by meaning and this is logged; nothing else changes. Chromosome 23 is read as X.

## Variants and gene statistics (unchanged from stages 7 and 13)
- z_i = beta/sebeta for the alternative allele; MAF = min(af_alt, 1 − af_alt) < 1%; |z_i| ≤ √50.
- pLoF: gnomAD v2 LOFTEE high-confidence variants from the GeneBayes file (`gnomad_v2_lofs.tsv.gz`, SHA 310be244…), matched by GRCh38
  position and alleles in both orientations; where the file's alternative allele is the reference allele, z is reoriented to the LoF allele (as in 35).
- Coding: RefSeq NM_ CDS exon ± 2 bp (`refGene_hg38.txt.gz`, SHA 58fbc18d…), as in 25.
- Undirected: Σz_i² ~ χ²(k) → z² per gene (cap 50), as in 25. Directional: z_gene = Σz_i/√k (Stouffer), as in 35.
- For F8 only: β_het per gene = inverse-variance weighted mean of pLoF betas (oriented to the LoF allele), SD units.
- Organ assignment of traits a priori as in 14/19 (blood: hb, mcv, rbc, plt; immune: wbc, lymph; kidney: creat, cysc, urate, urea, phosphate;
  liver: alt, ast, ggt, bili, albumin; lipids: ldl, hdl, tg, chol; glucose: glucose, hba1c; skeletal: height). Script 25 listed only the subset
  of traits available there; for FinnGen this assignment is applied to all available traits.

## Tests (code unchanged; only the reader is new)
For every test: pLoF primary, coding as second line. One-sided; permutation and random seed as in the original script (25: 23, 35: 35, 45: 45).

- **F1 Height (25 T1):** AR_nos against matched pool and against AR_nonskel, z², B = 5,000. Criterion: both P < 0.05 (pLoF).
  Control: AR_nonskel against pool, expected null.
- **F2 Organ diagonal (25 T2):** double-centred z², permutation of organ labels, B = 2,000. Criterion P < 0.05 (pLoF).
- **F3 Rule W1/W2b (25 T3):** classes from `22_einheiten_plof.tsv` (SHA 585e264d…), skeleton excluded, permutation within organ, B = 5,000.
  Criterion P < 0.05 per route (pLoF).
- **F4 Gene-specific direction (35, T10-8b):** direction map from 31 unchanged (script SHA d4c61868…), Stouffer z, centring and statistic as in 35,
  permutation of direction maps between genes, B = 1,000. **Primary criterion of the stage: F4 with pLoF P < 0.05.**
  Secondary: fraction in the expected direction among |z| > 2.
- **F5 Frozen prediction (P1):** `ATLAS_VORHERSAGE_2026-10-20.tsv` (SHA-256 7694042a…, frozen October 20, 2026, 08:37).
  d = expected direction · z_FinnGen (Stouffer, pLoF); statistic mean d; null: directions permuted between genes within the same trait,
  B = 2,000. Criterion P < 0.05 (pLoF). Secondary: fraction d > 0 among |z| > 2 (B = 500), specificity tertiles (tertile 3 expected strongest).
  Unlike in All of Us, the overall result F4 is also unknown at the time of this specification.
- **F6 Cross-cohort consistency (P2):** gene value = mean s · z over the direction cells (pLoF). Spearman FinnGen–UKB (criterion ρ > 0,
  P < 0.05) and FinnGen–G&H (second line), permutation of gene labels, B = 2,000. The pair FinnGen–All of Us is added with the same statistic
  once the Workbench data are available.
- **F7 Skeletal mechanism (P4):** GO classes from `go_kategorien.tsv` (SHA 0ed94db9…), AR_nos, height: fraction |z| > 1.96 in MATRIX against
  CILIUM + LYSO, one-sided Fisher test. Criterion P < 0.05 (pLoF).
- **F8 h for height (P5):** β_hom from `h_literatur.tsv` (SHA 521f5a68…, corrected by 32b), h = β_het/β_hom per gene. Median and IQR,
  sign test h > 0 (binomial, one-sided), Wilcoxon h < 0.5; concordance with UKB h (Spearman, `29b_h_haupttest.tsv`).
  Criterion: sign test P < 0.05.

Minimum sizes: F1 ≥ 20 AR_nos genes, F2/F4 ≥ 20 cells (as in 25/35), F5 ≥ 100 cells, F7 ≥ 5 genes per class, F8 ≥ 10 genes.
Below these, the test is "not testable", not "failed".

## Overlap and combination
FinnGen includes FINRISK participants of the THL Biobank, as do FinMetSeq and parts of the exome-chip consortia. FinnGen is therefore reported
**on its own** and **not** combined with FMS or exome-chip data by Fisher's method. A Fisher combination with Genes & Health only (no overlap)
is allowed and reported as a second line.

## Reading
- **Replicated:** F4 with pLoF P < 0.05.
- **Partly replicated:** F4 with coding variants only P < 0.05.
- **Not replicated:** F4 P ≥ 0.05 with both variant sets.
F1–F3 and F5–F8 are reported with their criteria regardless of outcome. Failed tests enter the test register and the manuscript.

## Limitations stated in advance
- Variant level and imputation: very rare pLoF variants are missing (MAC > 30, INFO > 0.6). Fewer genes than in an exome cohort are expected,
  with Finnish-enriched variants.
- Stouffer without LD adjustment is anticonservative in magnitude; it is read against the permutation, which carries the same correlation (as in 35).
- No synonymous negative control at variant level (as in 25/35); the control is the permutation and, in F1, AR_nonskel against pool.
- No stratification by sex (summary statistics).

## Deviations
Every deviation from this document is recorded in FINDINGS with reason and time and reported in the manuscript. Traits, variant sets, criteria
and the direction map are not re-chosen after download.

## Checksums of the methods (at specification)
`25_replikation_unabhaengig.py` 3fc66f6a…; `31_vorzeichen_genspezifisch.py` d4c61868…; `35_t108b_gerichtet.py` 48465d43…;
`45_haertungen.py` 0516eaec…; `14_organmatrix.py` f14760aa…; `genannotation.tsv.gz` b052ec97…; `organ_anteile_AR.tsv.gz` 9b3c85b7….
