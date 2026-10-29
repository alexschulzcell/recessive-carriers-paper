> English translation (21 October 2026) of `prespecification/PRAEREG_hauptprodukt.md`. The German original is the pre-specification document and is unchanged; this translation was made after the analyses and is provided for readability only.

# Pre-specification 2: the organ's own main product (StilleTraeger, stage 5)

Specified on 19 October 2026, **before** any computation. A new, separate test; no modification of stage 4 (where the primary test was negative and remains so).

**Prior knowledge (honest):** I know the 25 strongest genes of the organ diagonal (among them TMPRSS6, TFR2, TFRC → blood; LCAT, LIPC, ABCG5 → lipids; CLDN10, SLC22A12 → kidney; ABCC2 → liver) and the class means from stage 4 (broad GO classes). The classes here are defined exclusively by GO terms or HPA expression, not by gene names. I checked the HPA file only for format (columns, list of tissues).

## Hypothesis
Within an organ, carriers of recessive genes show an effect on the organ trait when the gene is involved in the **high-throughput main product of this organ** (erythrocyte: haemoglobin/haem/iron; liver: plasma proteins, bile acids, bilirubin; lipids: lipoprotein assembly and remodelling; kidney: tubular transport; muscle: sarcomere; lung: surfactant/gas exchange). Other recessive genes of the same organ are buffered.

## Unit and value
Unit = gene × organ (exclusively recessive Mendelian genes, HPO organ fraction ≥ 0.2; as in scripts 14/20/21). Value = mean of the double-centred, locus-adjusted z² (V3) over the traits of **this** organ (centring over all 31 traits). Genebass UKB, carriers ≥ 10.

## Routes (classes "in the main product" yes/no, per organ)
- **W1 GO main product** (descendants included):
  - Blood: GO:0042541 haemoglobin biosynthesis, GO:0006783 haem biosynthesis, GO:0055072 iron homeostasis, GO:0030218 erythrocyte differentiation, GO:0030220 platelet formation
  - Immune: GO:0002521 leukocyte differentiation, GO:0002377 immunoglobulin production, GO:0042582 azurophil granule, GO:0042581 specific granule
  - Liver: GO:0072562 blood microparticle (plasma proteins), GO:0006699 bile acid biosynthesis, GO:0015721 bile acid transport, GO:0042167 haem catabolism
  - Lipids: GO:0042157 lipoprotein metabolism, GO:0034358 plasma lipoprotein particle, GO:0030301 cholesterol transport
  - Kidney: GO:0070293 renal absorption, GO:0097254 renal tubular secretion, GO:0015747 urate transport, GO:0003094 glomerular filtration, GO:0035435 phosphate transport
  - Muscle: GO:0030017 sarcomere, GO:0006936 muscle contraction
  - Lung: GO:0043129 surfactant homeostasis, GO:0007585 gas exchange
  - Circulatory: GO:0008217 regulation of blood pressure
  - Nervous: GO:0099504 synaptic vesicle cycle
  - Glucose: GO:0030073 insulin secretion
- **W2a absolute expression:** HPA nTPM in the organ tissue in the top tertile among the recessive genes of this organ. Tissue assignment: blood → bone marrow; immune → max(bone marrow, spleen, lymph node); liver → liver; lipids → max(liver, small intestine); kidney → kidney; muscle → skeletal muscle; lung → lung; circulatory → max(heart muscle, blood vessel); nervous → cerebral cortex; glucose → pancreas.
- **W2b organ-specific expression:** ratio nTPM(organ tissue) / median nTPM over all tissues in the top tertile (per organ).
- **W3 combination:** W1 AND W2a.

## Test
- **Primary:** pLoF, **without skeleton**, all other organs pooled.
- **Statistic:** mean(value | class yes) − mean(value | class no).
- **Significance:** permutation of class labels within organ, B = 10,000, one-sided.
- **Multiple testing:** 4 routes → Bonferroni α = 0.0125 per route.
- **Replication:** missense|LC (other variants, largely other carriers), same direction with P < 0.05.
- **Secondary:**
  - per organ;
  - robust statistic (rank test within organ, mean of ranks);
  - without the 25 top genes of the diagonal (list from script 20);
  - skeleton with the MATRIX definition as reference (not independent).

## Interpretation
- **At least one route significant in the primary test and replicated with missense:** the main-product rule generalizes.
- **Robustness:** if the result holds without the top 25 genes, it does not rest only on textbook genes.
- **All routes negative:** the rule remains skeleton-specific (end of this line of analysis).
