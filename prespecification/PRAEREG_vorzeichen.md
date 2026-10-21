# Vorab-Festlegung 4 (Stufe 10): Gerichteter Diagonalen-Test

Festgelegt am 19.10.2026, **bevor** gerechnet wurde. Bekannt sind: Stufen 1–9; dass die Signkonvention von Genebass `BETA_Burden` biologisch stimmt (geprüft an ACAN, IHH, NPR2, FBN1, LDLR, PCSK9, APOB, LCAT, ABCA1, CETP, UGT1A1, TMPRSS6, HBB — nur Vorzeichen einzelner Gene, keine Gruppenwerte). Nicht bekannt: jeder Wert dieses Tests.

## Warum

Alle bisherigen Hauptanalysen (Stufen 1–9) laufen auf z², also **richtungsfrei**. Die Behauptung „Träger sind eine milde Form derselben Krankheit" verlangt aber das **richtige Vorzeichen**: Träger eines Dysplasie-Gens müssen *kleiner* sein, nicht bloß auffällig. Ohne Vorzeichen ist die Diagonale auch mit „Träger rezessiver Gene sind im Organmerkmal variabler" vereinbar — einer deutlich schwächeren Aussage. Dieser Test schließt die Lücke.

## Datengrundlage

Genebass (UKB 394k, irnt), signiertes z = sign(BETA_Burden) · |z| aus `Pvalue_Burden`; Filter wie in 14/20 (keep, carriers ≥ 10, z² bei 50 gekappt). Organ-Zuordnung aus `daten/organ_anteile_AR.tsv.gz` (HPO-Anteil ≥ 0,2), unverändert aus Stufe 2. Zentrierung wie in 14: zuerst Merkmals-Mittel, dann Gen-Mittel, **im rohen signierten z**, erst danach mit dem Erwartungsvorzeichen multipliziert.

Gerichteter Wert je Zelle: `d = s_erwartet · z_zentriert`. d > 0 heißt: der Träger verschiebt sich in Krankheitsrichtung.

## A) Erwartungsvorzeichen je Organ × Merkmal (Haupttest, eindeutige Zellen)

Abgeleitet aus der Krankheitsbiologie der rezessiven Erkrankungen des Organs, nicht aus Daten. s = −1 heißt „Träger erwartungsgemäß niedriger".

| Organ | Merkmal | s | Begründung (rezessive Krankheitsklasse) |
|---|---|---|---|
| Blut | hb | −1 | Anämien (Thalassämie, Fanconi, DBA, CDA, Sphärozytose) |
| Blut | rbc | −1 | dieselbe |
| Blut | plt | −1 | kongenitale Thrombozytopenien |
| Immun | wbc | −1 | SCID, Agammaglobulinämie, Neutropenien, CGD |
| Immun | lymph | −1 | dieselbe |
| Niere | creat | +1 | nephrotisches Syndrom, Nephronophthise, Zystennieren: GFR sinkt |
| Niere | cysc | +1 | dieselbe |
| Niere | urea | +1 | dieselbe |
| Leber | alt | +1 | Cholestase, Stoffwechselleber, Wilson: Hepatozytenschaden |
| Leber | ast | +1 | dieselbe |
| Leber | ggt | +1 | dieselbe |
| Leber | bili | +1 | Cholestase und Konjugationsdefekte (UGT1A1, ABCC2, ABCB11) |
| Leber | albumin | −1 | Synthesestörung |
| Lipide | hdl | −1 | LCAT-, ABCA1- (Tangier), APOA1-Defizienz |
| Lipide | tg | +1 | Chylomikronämie (LPL, APOC2, GPIHBP1, LMF1) |
| Lunge | fev1 | −1 | CF, Ziliendyskinesie, Surfactant-Defekte, ILD |
| Lunge | fvc | −1 | dieselbe |
| Muskel | grip | −1 | Myopathien, Muskeldystrophien |
| Nerven | fluidint | −1 | Intelligenzminderung, Neurodegeneration |
| Nerven | reaction | +1 | Reaktions*zeit*: schlechter = länger = höherer Wert |
| Skelett | height | −1 | Skelettdysplasien: Kleinwuchs |
| Skelett | sitheight | −1 | dieselbe |

22 eindeutige Zellen über 10 Organe.

### Vorab als richtungs-gemischt erklärt (NICHT im Haupttest, getrennt berichtet)

Hier gibt es rezessive Krankheiten in beide Richtungen; eine a-priori-Richtung wäre geraten.

| Organ | Merkmal | Tendenz | Gegenbeispiel |
|---|---|---|---|
| Blut | mcv | −1 (mikrozytär) | megaloblastisch (TCN2, MTR, Folattransport) → +1 |
| Niere | urate | +1 (GFR sinkt) | Hypourikämie SLC22A12, SLC2A9 → −1 |
| Niere | phosphate | −1 (tubulärer Verlust) | CKD-Hyperphosphatämie → +1 |
| Lipide | ldl | −1 (Hypobetalipoproteinämie, MTTP, APOB) | LDLRAP1, rezessive FH → +1 |
| Lipide | chol | −1 | dieselbe Spaltung |
| Glukose | glucose | +1 (neonataler Diabetes) | Hyperinsulinismus, Glykogenosen → −1 |
| Glukose | hba1c | +1 | dieselbe |
| Skelett | ebmd | −1 (Osteogenesis imperfecta) | Osteopetrose, sklerosierende Dysplasien → +1 |
| Kreislauf | sbp | −1 (Kardiomyopathie) | schwacher Stellvertreter, beides möglich |

## B) Tests und Kriterien

**T10-1 Haupttest (pLoF, eindeutige Zellen).** Statistik = Mittel von d über alle passenden Gen × Merkmal-Zellen. Nullverteilung: Permutation der Organ-Etiketten zwischen Genen (B = 1000, Mehrfachzuordnung erhalten, wie 14/20). Einseitig, Kriterium **p < 0,05**.

**T10-2 Replikation.** Dasselbe mit Missense (mislc), Kriterium p < 0,05.

**T10-3 Negativkontrolle.** Dasselbe mit synonym. Erwartung: null. Ein p < 0,05 bei synonym wird als Warnung berichtet und der Haupttest dann nur locus-bereinigt gelesen.

**T10-4 Vorzeichen-Anteil (anschaulich).** Anteil der passenden Zellen mit d > 0, Binomialtest gegen 0,5, zweiseitig. Vorhersage: > 0,5. Dazu dieselbe Kennzahl für die unpassenden Zellen (Erwartung 0,5).

**T10-5 Härtung, vorab festgelegt wie in Stufe 3.**
- locus-bereinigt: Zellen mit |z_syn| > 1,96 verworfen;
- ohne die Top-10 und die Top-25 Gene (nach |z| im passenden Merkmal);
- ein festes Primärmerkmal je Organ (hb, wbc, creat, bili, hdl, fev1, grip, fluidint, height);
- leave-one-organ-out über die 10 Organe mit eindeutigen Zellen.

Kriterium: Das Vorzeichen bleibt positiv und p < 0,05 sowohl locus-bereinigt als auch bei festem Primärmerkmal.

**T10-6 Gemischte Zellen (getrennt).** Haupttest-Statistik nur über die 9 gemischten Zellen, mit der oben genannten Tendenz. Keine Vorhersage; rein beschreibend.

## C) Gen-eigene Erwartungsrichtung aus HPO (Sekundärtest)

Die Organ-Richtung in A ist grob. Präziser: jedes Gen bekommt seine **eigene** Richtung aus den HPO-Termen seiner OMIM-Krankheiten — klinisch beschrieben und von der UKB unabhängig.

**Term → Merkmal × Richtung, vorab festgelegt** (HPO-Term beim Gen vorhanden → diese Erwartung):

- **Blut/Immun:** HP:0001903 Anemia → hb −1; HP:0001935 Microcytic anemia → mcv −1, hb −1; HP:0001889 Macrocytic anemia → mcv +1, hb −1; HP:0001873 Thrombocytopenia → plt −1; HP:0001894 Thrombocytosis → plt +1; HP:0001882 Leukopenia → wbc −1; HP:0001875 Neutropenia → wbc −1; HP:0001974 Leukocytosis → wbc +1; HP:0001888 Lymphopenia → lymph −1
- **Niere:** HP:0003259 Elevated circulating creatinine → creat +1, cysc +1; HP:0012622 Chronic kidney disease → creat +1, cysc +1, urea +1; HP:0002905 Hyperphosphatemia → phosphate +1; HP:0002148 Hypophosphatemia → phosphate −1; HP:0002149 Hyperuricemia → urate +1; HP:0003537 Hypouricemia → urate −1
- **Leber:** HP:0002904 Hyperbilirubinemia, HP:0002908 Conjugated hyperbilirubinemia, HP:0008282 Unconjugated hyperbilirubinemia, HP:0003573 Increased total bilirubin → bili +1; HP:0002910 Elevated hepatic transaminase → alt +1, ast +1; HP:0003234 Elevated GGT → ggt +1; HP:0003073 Hypoalbuminemia → albumin −1
- **Lipide:** HP:0003233 Decreased HDL cholesterol → hdl −1; HP:0002155 Hypertriglyceridemia → tg +1; HP:0003124 Hypercholesterolemia → chol +1, ldl +1; HP:0003563 Hypobetalipoproteinemia, HP:0008181 Abetalipoproteinemia → ldl −1, chol −1
- **Glukose:** HP:0003074 Hyperglycemia, HP:0000857 Neonatal insulin-dependent diabetes mellitus → glucose +1, hba1c +1; HP:0001943 Hypoglycemia → glucose −1
- **Lunge:** HP:0002097 Emphysema, HP:0006536 Airway obstruction, HP:0002113 Pulmonary infiltrates, HP:0002206 Pulmonary fibrosis, HP:0002205 Recurrent respiratory infections → fev1 −1, fvc −1
- **Muskel:** HP:0003701 Proximal muscle weakness, HP:0003324 Generalized muscle weakness, HP:0003202 Skeletal muscle atrophy → grip −1
- **Nerven:** HP:0001249 Intellectual disability, HP:0001263 Global developmental delay, HP:0002376 Developmental regression → fluidint −1, reaction +1
- **Skelett:** HP:0004322 Short stature, HP:0003510 Severe short stature, HP:0008873 Disproportionate short-limb short stature, HP:0003498 Disproportionate short stature → height −1, sitheight −1; HP:0000098 Tall stature → height +1, sitheight +1; HP:0004349 Reduced bone mineral density, HP:0000939 Osteoporosis → ebmd −1; HP:0011001 Increased bone mineral density, HP:0011002 Osteopetrosis → ebmd +1

Richtungslose Oberbegriffe (HP:0003119 Abnormal circulating lipid concentration, HP:0003116 Abnormal circulating cholesterol, HP:0012100 Abnormal circulating creatinine u. ä.) werden verworfen. Terme werden **mit ihren HPO-Nachkommen** gezählt. Hat ein Gen für dasselbe Merkmal Terme beider Richtungen, fällt diese Zelle weg (widersprüchlich).

**T10-7:** Haupttest-Statistik auf diesen gen-eigenen Zellen, pLoF, Permutation der Gen-Etiketten. Kriterium p < 0,05, Missense-Replikation p < 0,05, synonym null.

**Vorhersage:** Der Effekt (dz) ist **größer** als in T10-1, weil die Richtung präziser ist. Wird er kleiner, wird das als Argument für die gröbere Organ-Ebene berichtet, nicht verschwiegen.

## D) Was den Test scheitern lässt

Ein Mittel von d ≤ 0 oder p ≥ 0,05 in T10-1 heißt: die Diagonale ist richtungsfrei, also nur „mehr Varianz im passenden Organ" und **nicht** „milde Form der Krankheit". Das würde den Titelanspruch widerlegen und wäre so zu berichten; die z²-Ergebnisse der Stufen 1–9 blieben gültig, aber ihre Deutung müsste umgeschrieben werden.

---

# Nachtrag (Stufe 10b): gerichtete Replikation außerhalb der UKB — T10-8

Festgelegt am 19.10.2026, **bevor** gerechnet wurde. Bekannt sind: die Ergebnisse von T10-1 bis T10-7 in der UKB; dass die in Stufe 6–8 gespeicherten Replikationstabellen **nur Σz² enthalten, also keine Vorzeichen**, und die Vorzeichen deshalb neu aus den Quelldaten geholt werden müssen. Nicht bekannt: jeder Wert dieses Nachtrags.

**Die Richtungskarte wird nicht angetastet.** Es gilt dieselbe gen-eigene HPO-Karte wie in Teil C, unverändert, inklusive der Regel, dass widersprüchliche Zellen wegfallen. Es werden keine Merkmale, Gene oder Richtungen neu gewählt. Damit ist dies eine echte Replikation ohne Freiheitsgrade.

## A) Hauptquelle: Genes & Health 44k, REGENIE-Burden (gen-ebene, signiert)

Quelle wie Stufe 8 (`gs://genesandhealth_publicdatasets/results_44k_ExWAS`, CC BY-SA 4.0), 33k britische Südasiaten, keine UKB-Überlappung. Masken A (LOFTEE-HC-pLoF), C (pLoF + Missense), D (synonym), Test ADD, AAF < 1 %, Allelzahl ≥ 5 — alles wie Stufe 8. Einziger Unterschied: statt (BETA/SE)² wird **sign(BETA) · |BETA/SE|** behalten.

Zentrierung und Statistik genau wie T10-7: Merkmals-Mittel, dann Gen-Mittel im rohen signierten z, danach mit der gen-eigenen Richtung multiplizieren; Nullverteilung durch Permutation der Richtungskarten zwischen den Genen (B = 1000), einseitig.

**T10-8a:** Kriterium **p < 0,05** mit Maske A (pLoF). Maske C (pLoF + Missense) als zweite Zeile, Maske D (synonym) als Negativkontrolle, Erwartung null.

**Erwartung vorab:** Stufe 8 hat gezeigt, dass pLoF in 33k etwa 8 % der UKB-Stärke erreicht und dort Diagonale und Regel nur mit Missense repliziert haben. Es ist deshalb **damit zu rechnen, dass Maske A zu schwach ist und nur Maske C das Kriterium erreicht.** Dieser Fall ist vorab als „teilweise repliziert" definiert, nicht als Erfolg und nicht als Scheitern.

## B) Zweitquellen: variantenbasiert, mit vorab festgelegter Vorzeichen-Aggregation

Quellen wie Stufe 7: FinMetSeq (19k Finnen), Exom-Chip-Verbund vor UKB (GIANT-Größe, BCX-Blutbild, CKDGen, GLGC-Lipide), Genes & Health variant-level. Alle frei, ohne Login, neu zu laden (`work/24_download_replikation.sh`).

**Aggregation je Gen × Merkmal (neu, hier festgelegt):** Für die qualifizierenden Varianten eines Gens (MAF < 1 %, Zuordnung wie Stufe 7; pLoF = gnomAD-LOFTEE-HC, kodierend = CDS) wird
`z_gen = Σ z_i / √k` gebildet (Stouffer), wobei z_i = BETA_i / SE_i bezogen auf das **Alternativallel** ist und k die Variantenzahl. Begründung: bei LoF-Varianten ist das Alternativallel das schädigende, die Vorzeichen sind damit über Varianten hinweg gleich orientiert. Kappung von z_i auf ±√50 wie überall.

**Bekannte Grenze, vorab benannt:** Stouffer ohne LD-Bereinigung ist bei korrelierten Varianten antikonservativ im Betrag. Das betrifft die **Stärke**, nicht die Richtung; der Test wird gegen die Permutation der Richtungskarten gelesen, die dieselbe Korrelationsstruktur trägt.

**T10-8b:** dieselbe Statistik, je Quelle, Kriterium p < 0,05; zusätzlich Fisher-Kombination über die Quellen. Synonym bzw. nicht-kodierend als Negativkontrolle, wo verfügbar.

## C) Was zählt als Replikation

- **Repliziert:** T10-8a mit Maske A oder C erreicht p < 0,05, Negativkontrolle synonym bleibt null.
- **Teilweise repliziert:** nur Maske C, oder nur die Zweitquellen.
- **Nicht repliziert:** kein Mittel > 0 bei p < 0,05 in irgendeiner Quelle. Dann gilt der Richtungsbefund als UKB-spezifisch und ist so zu berichten; T10-7 bliebe gültig, aber ohne externe Stütze.
