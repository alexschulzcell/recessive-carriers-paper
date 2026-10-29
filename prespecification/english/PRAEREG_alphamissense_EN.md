> English translation (21 October 2026) of `prespecification/PRAEREG_alphamissense.md`. The German original is the pre-specification document and is unchanged; this translation was made after the analyses and is provided for readability only.

# Pre-specified: dose curve within carriers – does the disease direction scale with predicted deleteriousness? (stage 17, `work/42_*`)

Specified on 20 October 2026, **before** AlphaMissense scores were linked to association data. Known: stages 1–16, in particular T10-8b (variant level, coding: G&H P 0.031, exome chip P 0.009, FMS P 0.14). Not known: any number of this test.

## Question
If carriers show a mild form of the disease because their gene dose falls, the effect must **scale with the loss of function of the variant**: harmless missense ≈ 0, damaging missense → towards the disease, close to the pLoF effect. This tests the mechanism within the same genes and traits.

## Data
- Per-variant associations (alternative allele, BETA/SE): Genes & Health variant level (hg38), FinMetSeq (hg19), pre-UKB exome-chip consortia (hg19); exactly the readers of script 25.
- AlphaMissense (Cheng et al. 2023, Zenodo 8208688), canonical transcripts, hg19/hg38 matched to the source; key chromosome–position–ref–alt. If the AlphaMissense alternative allele is the reference allele of the file, z is flipped (oriented to the missense allele).
- pLoF variants as in script 35 (gnomAD LOFTEE high confidence, oriented to the LoF allele) as reference point for the "full half dose".
- MAF < 1% only; |z| ≤ √50. Gene assignment via the AlphaMissense transcript → gene symbol (HGNC) or CDS as in script 25.
- Genes × traits: gene-specific HPO direction map from script 31, **unchanged**. d_v = s_expected · z_v.

## Statistics
- **Primary T17-1 (within cell):** regression d_v ~ am_pathogenicity with a fixed effect for each (source × gene × trait); slope β. Null: permute am_pathogenicity **within** each cell (B = 1,000), one-sided (expectation β > 0). All three sources pooled. Criterion P < 0.05.
- **T17-2 (classes, descriptive + test):** mean d in the AlphaMissense classes likely benign (< 0.34), ambiguous, likely pathogenic (> 0.564), plus pLoF; expectation monotonic benign < ambiguous < pathogenic ≤ pLoF; test pathogenic − benign within cell as in T17-1.
- **Negative control T17-3:** the same slope with direction maps permuted between genes (B = 200) → expectation 0. In addition, the slope without direction sign on |z| (descriptive: more damaging variants have larger effects, regardless of direction).
- Secondary: each source separately; only genes with ≥ 5 missense variants in the cell.
- Side row (descriptive, selected): Genebass UKB single variants with P < 1 × 10⁻⁴ (API; only these are unbinned): fraction in the disease direction by AlphaMissense class.

## Interpretation, specified in advance
- β > 0 and monotonic: dose–response within carriers; supports "dose" as mechanism and makes AlphaMissense a predictor of the carrier phenotype.
- Null: either too little power at variant level (33k/19k cohorts) or the direction does not arise through dose. Report power via the number of variants.

---

# Addendum 20 October 2026: T17 with UKB single variants from Backman et al. 2021 (Regeneron, 450k exomes; `work/43`)

Specified **before** loading the data. Known: result of T17 in G&H/FMS/exome chip (T17-1 P 0.032; T17-3 control almost as large; post hoc bidirectional traits P 0.08).

- Source: GWAS Catalog, PMID 34662886, single-variant studies (not "gene-based burden"), UKB_Freeze_450_EUR, REGENIE, traits rank-inverse-normal transformed, hg38, CC0. Matching via the UKB field number (identical to Genebass: height 50, sitheight 20015, blood count 30000–30120, biochemistry 306xx–308xx, fev1 3063, fvc 3062, grip 47, fluidint 20016, reaction 20023, sbp 4080, ebmd 3148). If a field is missing, the trait is dropped.
- **UKB overlaps with Genebass** (the discovery data of the direction map, however, are HPO, not UKB). For T17 this is unproblematic because T17 asks a new question (dose within carriers) and is not a replication.
- Statistics **identical to script 42** (T17-1 to T17-3, code unchanged, only source "BKM" added). Additionally pre-specified: **T17-4 (now the primary criterion for gene specificity):** only traits with genes in both directions; slope against permutation of directions between genes **within the trait** (B = 500), one-sided, criterion P < 0.05. This is the control found post hoc, now pre-specified on new data.
