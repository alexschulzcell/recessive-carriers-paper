> English translation (21 October 2026) of `prespecification/PRAEREG_diagnosen.md`. The German original is the pre-specification document and is unchanged; this translation was made after the analyses and is provided for readability only.

# Pre-specified: organ diagonal on clinical diagnoses (stage 15, `work/39_*`)

Specified on 20 October 2026, **before** diagnosis data were loaded. Known: stages 1–13 (laboratory values/measurements), stage 1 (osteoarthritis/osteoporosis/fracture null).

## Question
Do carriers of recessive disease genes have an increased risk of **clinical diagnoses in the organ of their recessive disease**? So far all endpoints were measurements. Diagnoses turn "mild form of the disease" into a clinically visible statement (framed biologically, no clinical validation).

## Data
Genebass (UKB 394k exomes) ICD-10 first occurrence (UKB fields 130000 onwards), three-character codes with ≥ 1,000 cases, sets pLoF, missense|LC, synonymous. Organ assignment of ICD blocks a priori (before loading): blood D50–D77; immune D80–D89; kidney N00–N29; liver K70–K83; lipids E78; glucose E10–E16; lung J40–J99; muscle G70–G73, M60–M63; nervous G00–G69, G74–G99, F00–F03, F70–F79; skeletal M15–M19, M40–M43, M80–M94, Q65–Q79; circulatory I00–I99. This gives 122 codes (`work/39_icd_auswahl.tsv`). Genes and organ labels as in stage 2 (exclusively recessive, HPO fraction ≥ 0.2).

## Statistics
- **Directed** (diagnosis = disease direction, expectation: carriers more often diagnosed): signed z = sign(BETA_Burden)·|z|, genes with ≥ 10 carriers, |z| ≤ √50.
- Centring by code and gene (as in scripts 30/31). Gene value per organ = mean of the centred z over the organ's codes.
- Statistic: mean of the matching gene × organ cells; null: permute organ labels between genes (B = 10,000), one-sided.
- **Primary:** pLoF, P < 0.05. Missense second row, synonymous negative control (expectation null).
- Secondary: directionless (z²) as in stage 2; leave-one-organ-out; without circulatory (largest block, weakest link).
- Descriptive: odds ratio per carrier for the strongest gene–diagnosis pairs.

## Interpretation, specified in advance
- Positive: carrier status becomes visible down to the level of diagnoses in the disease organ.
- Null: effects remain subclinical (consistent with h ≈ 1%), an honest statement with power reported.
