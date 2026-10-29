# Pre-specified: independent replication in BioBank Japan (stage 21, addendum to `PRAEREG_finngen.md`, `work/49_*`)

English translation of `PRAEREG_bbj.md`; the German original is authoritative.

Specified on October 28, 2026, **before** any BBJ trait file was downloaded.
Known: the download list of PheWeb.jp (traits, n, file names), the NBDC pages hum0014/hum0197 and the methods of Sakaue et al.
(Nat Genet 2021; imputation with Minimac3 using BBJ WGS n = 1,037 + 1KGP3, Rsq ≥ 0.7, 13.5 million variants, GRCh37).
At the time of specification the FinnGen download (`49a`) is running; **no FinnGen result has been computed or looked at.**
Not known: any association value from BBJ.

Purpose: a second large, freely available cohort independent of UKB, of a different ancestry (East Asian, Japan). It also tests whether a
FinnGen result depends on Finnish-enriched founder variants: BBJ carries different rare variants in the same genes.

## Data
BioBank Japan, GWAS summary statistics from hum0197 (PheWeb.jp "Downloads", only the **BBJ column**, never the meta-analyses with UKB or
FinnGen), open without restriction (Sakaue 2021: "openly shared without any restrictions"; NBDC policy). Citation Sakaue et al. 2021.
Imputed, GRCh37, linear regression; no overlap with UKB, FinnGen or G&H.

| Trait | File | n |
|---|---|---|
| height | hum0197.v3.BBJ.Hei.v1.zip | 165,056 |
| hb | hum0197.v3.BBJ.Hb.v1.zip | 152,447 |
| mcv | hum0197.v3.BBJ.MCV.v1.zip | 129,832 |
| rbc | hum0197.v3.BBJ.RBC.v1.zip | 153,512 |
| plt | hum0197.v3.BBJ.PLT.v1.zip | 148,623 |
| wbc | hum0197.v3.BBJ.WBC.v1.zip | 154,355 |
| lymph | hum0197.v3.BBJ.LYM.v1.zip | 95,717 |
| creat | hum0197.v3.BBJ.sCr.v1.zip | 150,266 |
| urate | hum0197.v3.BBJ.UA.v1.zip | 129,405 |
| urea | hum0197.v3.BBJ.BUN.v1.zip | 148,767 |
| phosphate | hum0197.v5.gwas.P.v1.txt.gz | 46,848 |
| alt | hum0197.v3.BBJ.ALT.v1.zip | 150,545 |
| ast | hum0197.v3.BBJ.AST.v1.zip | 150,068 |
| ggt | hum0197.v3.BBJ.GGT.v1.zip | 133,471 |
| bili | hum0197.v3.BBJ.TBil.v1.zip | 124,341 |
| albumin | hum0197.v3.BBJ.Alb.v1.zip | 120,539 |
| ldl | hum0197.v3.BBJ.LDLC.v1.zip | 72,866 |
| hdl | hum0197.v3.BBJ.HDLC.v1.zip | 74,970 |
| tg | hum0197.v3.BBJ.TG.v1.zip | 111,667 |
| chol | hum0197.v3.BBJ.TC.v1.zip | 135,808 |
| glucose | hum0197.v3.BBJ.Glucose.v1.zip | 133,336 |
| hba1c | hum0197.v3.BBJ.HbA1c.v1.zip | 71,221 |
| sbp | hum0197.v3.BBJ.SBP.v1.zip | 145,505 |

Not available: cysc, sitheight, ebmd, fev1, fvc, grip, fluidint, reaction. Phosphate comes from a later release (v5); it is used as
delivered.

Download (`work/49b_bbj_download.py`): download the archive, log the SHA-256, assign columns by meaning (chromosome, position, REF, ALT,
frequency of the ALT allele, BETA, SE; name in the log); only variants with MAF < 1 % are stored that are pLoF (GRCh37 position) or lie in a
coding exon ± 2 bp (refGene hg19); afterwards the archive is deleted. Chromosome 23 is read as X. If a column is missing, the script
stops instead of guessing.

## Variants, gene statistics, tests
**Everything as in `PRAEREG_finngen.md`**, with GRCh37 matching (`lof19`, `refGene_hg19`, as FinMetSeq in 25/35). Tests B1–B8 correspond to F1–F8
with the same criteria and minimum sizes. Organ assignment additionally sbp → circulation (as 14/19).
**Primary criterion for BBJ: B4 (gene-specific direction) with pLoF P < 0.05.**
B8 (h) only if the source (README in the archive or Sakaue 2021) reports height as rank-normalised or in SD units; otherwise "not testable".

## Cross-ancestry reading ("founder bias")
1. F4 (FinnGen) and B4 (BBJ) are reported **separately**.
2. **X1:** Spearman correlation of the gene values (mean s · z over the direction cells, pLoF) FinnGen against BBJ, permutation of the gene labels
   (B = 2,000, one-sided), criterion ρ > 0, P < 0.05. In addition BBJ–UKB and BBJ–G&H as F6.
3. **X2:** Fisher combination of the F4/B4 P values with Genes & Health (35, variant based) across the three non-overlapping cohorts
   FinnGen, BBJ, G&H, per variant set. Second line.
4. Reading: if the direction replicates in FinnGen **and** BBJ (both pLoF or both coding P < 0.05), the result does not depend on one
   ancestry. If it replicates only in FinnGen, this is reported as a possible consequence of Finnish founder variants.

## Limitations named in advance
- Imputation with Rsq ≥ 0.7: very rare variants are largely missing; fewer pLoF genes are to be expected in BBJ than in FinnGen.
- Quantitative traits in Sakaue 2021 are not consistently documented as rank-normalised; z = BETA/SE does not depend on this, β_het (B8) does.
- Lipids and HbA1c, with n ≈ 71,000–75,000, are considerably smaller.
