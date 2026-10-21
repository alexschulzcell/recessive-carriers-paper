# Vorab-Festlegung 3 (Stufe 9): Dominanz h, Zelltyp-Auflösung, Bevölkerungsebene

Festgelegt am 19.10.2026, **bevor** gerechnet wurde. Bekannt sind: Stufen 1–8 (insbesondere die W2b-Ergebnisse auf Gewebeebene); die HPA-Zelltypliste wurde nur auf Namen geprüft.

## A) Zelltyp-Auflösung der Regel (W4)
- **Hypothese:** Der Trägereffekt im Organ-Merkmal folgt der Spezifität in der **Produzenten-Zelle** des Organs, mindestens so gut wie der Gewebe-Spezifität (W2b).
- **Produzenten-Zellen** (HPA single cell type, nCPM; Maximum über die genannten Typen):
  - Blut: erythrocyte progenitors, megakaryocyte-erythroid progenitors, megakaryocytes
  - Immun: neutrophils, neutrophil progenitors, t-cells, b-cells, monocytes, nk-cells
  - Leber: hepatocytes, cholangiocytes
  - Lipide: hepatocytes, enterocytes
  - Niere: proximal tubule cells, loop of henle epithelial cells, distal convoluted tubule cells, podocytes
  - Lunge: alveolar cells type 1 und 2
  - Muskel: myonuclei
  - Kreislauf: cardiomyocytes, vascular smooth muscle cells, vascular endothelial cells
  - Nerven: brain excitatory/inhibitory neurons
  - Glukose: pancreatic islet cells
  - Skelett: keine Chondrozyten in HPA → ausgeschlossen
- **Klasse W4:** Zellspezifität = log2((Produzent + 1) / (Median aller Zelltypen + 1)) im oberen Tertil je Organ.
- **Primärtest:** pLoF, ohne Skelett, Einheiten und Wert wie in 22. Permutation innerhalb Organ, einseitig, p < 0,05; Replikation Missense p < 0,05.
- **Sekundär:** Regression Wert ~ z(Gewebe-Spez.) + z(Zell-Spez.) + Organ. Vorhersage: Der Koeffizient der Zell-Spezifität ist > 0 bei gegebener Gewebe-Spezifität.

## B) Bevölkerungsebene
- **Gen-Sets je Organ:** AR-Gene mit organ-passendem Trägereffekt, FDR < 0,05 in Genebass pLoF (innerhalb der AR-Gene je Merkmal).
- **Effekt** je Träger aus AZ-PheWAS (UKB 500k, Europäer, PTV, SD); **Trägerfrequenz** aus Genebass (2 × CAF, pLoF).
- **Kennzahlen:**
  - Anteil der Bevölkerung mit ≥ 1 solchen pLoF;
  - Anteil am unteren bzw. oberen Rand (< −2 SD bzw. > +2 SD) des Merkmals, der auf Träger entfällt: attributabler Anteil Σ f·[Φ(−2−β) − Φ(−2)] / Φ(−2);
  - erklärte Varianz Σ f·β².
- **Merkmale:** Größe (Hauptfall), HDL, Hb, MCV, Harnsäure, Bilirubin.
- **Lesart:** beschreibend, kein Test. Auswahl über Signifikanz → Effekte nach oben verzerrt (Winner's curse). Daher zusätzlich eine konservative Variante nur mit AZ-Effekten in AZ-p < 1e-4.

## C) Dominanzkoeffizient h
- **Patientenwerte:** Agent, Literatur. Gene ausgewählt nach Datenverfügbarkeit (≥ 3 biallelische Patienten), nicht nach Trägereffekt.
- **Trägereffekt β_het** in SD: AZ-PTV (Europäer), wenn berichtet (p < 0,1). Sonst aus Genebass: z_pLoF / √(Träger), kalibriert über Gene mit beiden Angaben (lineare Regression ohne Achsenabschnitt).
- **Berechnung:** h = β_het / β_hom (beide in SD, gleiches Merkmal).
- **Hypothesen:**
  - **H1:** Median h liegt deutlich unter 0,5 (partiell rezessiv). Erwartung 0,02–0,2.
  - **H2:** h ist für Hauptprodukt-Gene (W1 ODER W2b) höher als für andere; Mann-Whitney einseitig.
  - **H3:** h korreliert positiv mit s_het (Spearman, einseitig).
- **Validierung:** Wo die Literatur Werte heterozygoter Eltern enthält, wird das Literatur-h mit unserem h verglichen.
- **Qualitätsregeln:**
  - Gene mit |β_hom| < 0,5 SD fallen heraus (Nenner instabil).
  - h wird auf [−1, 2] begrenzt.
  - Gegenläufige Vorzeichen (h < 0) werden berichtet, nicht entfernt.
