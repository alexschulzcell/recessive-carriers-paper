> English translation (21 October 2026) of `prespecification/PRAEREG_maus.md`. The German original is the pre-specification document and is unchanged; this translation was made after the analyses and is provided for readability only.

# Pre-specified: organ diagonal in heterozygous knockout mice (stage 14, `work/38_*`)

Specified on 20 October 2026, **before** IMPC data were analysed. Known: stages 2/15 (IMPC heterozygous, body length/BMD only, group null, human–mouse direction ρ 0.28). Not known: any number of this test.

## Question
Do heterozygous knockout mice of recessive human disease genes show their phenotype preferentially in the **organ system of the human disease**? This would be the organ diagonal as a **causal experiment** in a second species: defined null alleles, controlled environment, no linkage disequilibrium, no stratification.

## Data
IMPC data release (latest, ftp.ebi.ac.uk/pub/databases/impc/all-data-releases/latest/results): `statistical-results-ALL` (all tested parameters with P value, zygosity, top-level MP) and `genotype-phenotype-assertions-ALL` (significant hits). Orthologues: MGI HMD_HumanPhenotype (1:1). Genes: exclusively recessive human genes with organ assignment (`organ_anteile_AR.tsv.gz`, fraction ≥ 0.2; as in stage 2).

## Organ → MP top level (a priori)
| Organ | MP top level |
|---|---|
| Blood | hematopoietic system phenotype |
| Immune | immune system phenotype |
| Kidney | renal/urinary system phenotype |
| Liver | liver/biliary system phenotype |
| Lipids, glucose | homeostasis/metabolism phenotype |
| Lung | respiratory system phenotype |
| Muscle | muscle phenotype |
| Nervous | nervous system phenotype; behavior/neurological phenotype |
| Skeletal | skeleton phenotype; limbs/digits/tail phenotype; growth/size/body region phenotype |
| Circulatory | cardiovascular system phenotype |

## Unit and statistics
- Cell = gene × MP top-level system tested for this gene in **heterozygous** animals (at least one parameter in statistical-results).
- Value = 1 if a significant assertion (genotype-phenotype-assertions, zygosity heterozygote) exists in the system, otherwise 0.
- Test statistic: mean value of matching cells minus mean value of non-matching cells, centred by gene and system (double centring as in stage 2). Null: permute organ labels between genes (B = 10,000), one-sided.
- **Primary criterion:** heterozygous, P < 0.05.
- **Positive control:** the same statistic for homozygous lines (expectation: strongly positive; shows that mapping and statistics work).
- **Negative control:** non-recessive genes without Mendelian disease have no organ labels; the permutation itself serves as null instead.
- Secondary: only genes with both het and hom data; ratio of hit rates het/hom in matching cells as a rough "mouse h".

## Interpretation, specified in advance
- Positive: the diagonal holds across species and causally. Half the dose also shows in the mouse in the disease organ.
- Null with positive hom control: mouse hets (mostly ~8 animals per line, systematically only for lethal lines) too weak; "no, with a number", report power.
- Caution, selection: IMPC phenotypes hets systematically mainly for homozygous lethal lines. Therefore analyse separately by viability as well.
