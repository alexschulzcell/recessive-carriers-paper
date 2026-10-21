# Vorab festgelegt: Organ-Diagonale bei heterozygoten Knockout-Mäusen (Stufe 14, `work/38_*`)

Festgelegt am 20.10.2026, **bevor** IMPC-Daten ausgewertet wurden. Bekannt: Stufe 2/15 (IMPC heterozygot, nur Körperlänge/BMD, Gruppe null,
Richtung Mensch–Maus ρ 0,28). Nicht bekannt: jede Zahl dieses Tests.

## Frage
Zeigen heterozygote Knockout-Mäuse rezessiver menschlicher Krankheitsgene ihren Phänotyp bevorzugt im **Organsystem der menschlichen Krankheit**?
Das wäre die Organ-Diagonale als **kausales Experiment** in einer zweiten Art: definierte Nullallele, kontrollierte Umwelt, kein LD, keine Stratifikation.

## Daten
IMPC Data Release (latest, ftp.ebi.ac.uk/pub/databases/impc/all-data-releases/latest/results): `statistical-results-ALL` (alle getesteten
Parameter mit p-Wert, Zygotie, Top-Level-MP) und `genotype-phenotype-assertions-ALL` (signifikante Treffer). Orthologe: MGI HMD_HumanPhenotype (1:1).
Gene: rein rezessive menschliche Gene mit Organ-Zuordnung (`organ_anteile_AR.tsv.gz`, Anteil ≥ 0,2; wie Stufe 2).

## Organ → MP-Top-Level (a priori)
| Organ | MP-Top-Level |
|---|---|
| Blut | hematopoietic system phenotype |
| Immun | immune system phenotype |
| Niere | renal/urinary system phenotype |
| Leber | liver/biliary system phenotype |
| Lipide, Glukose | homeostasis/metabolism phenotype |
| Lunge | respiratory system phenotype |
| Muskel | muscle phenotype |
| Nerven | nervous system phenotype; behavior/neurological phenotype |
| Skelett | skeleton phenotype; limbs/digits/tail phenotype; growth/size/body region phenotype |
| Kreislauf | cardiovascular system phenotype |

## Einheit und Statistik
- Zelle = Gen × MP-Top-Level-System, das für dieses Gen in **heterozygoten** Tieren getestet wurde (mindestens ein Parameter in statistical-results).
- Wert = 1, wenn eine signifikante Assertion (genotype-phenotype-assertions, Zygotie heterozygote) im System vorliegt, sonst 0.
- Teststatistik: mittlerer Wert der passenden Zellen minus mittlerer Wert der nicht passenden Zellen, zentriert nach Gen und System
  (doppelte Zentrierung wie Stufe 2). Null: Organ-Etiketten zwischen Genen permutieren (B = 10.000), einseitig.
- **Primärkriterium:** heterozygot, p < 0,05.
- **Positivkontrolle:** dieselbe Statistik für homozygote Linien (Erwartung: stark positiv; zeigt, dass Mapping und Statistik funktionieren).
- **Negativkontrolle:** nicht-rezessive Gene ohne Mendel-Krankheit haben keine Organ-Etiketten; stattdessen Permutation selbst als Null.
- Sekundär: nur Gene mit sowohl het- als auch hom-Daten; Verhältnis der Trefferraten het/hom in passenden Zellen als grobes „Maus-h“.

## Lesart vorab
- Positiv: Die Diagonale gilt artübergreifend und kausal. Halbe Dosis zeigt sich auch in der Maus im Krankheitsorgan.
- Null bei positiver Hom-Kontrolle: Maus-Hets (meist ~8 Tiere je Linie, nur letale Linien systematisch) zu schwach; „nein mit Zahl“, Power angeben.
- Achtung Selektion: IMPC phänotypisiert Hets systematisch vor allem bei homozygot letalen Linien. Daher Auswertung auch getrennt nach Viabilität.
