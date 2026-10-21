# Vorab festgelegt: Dosiskurve innerhalb der Träger – skaliert die Krankheitsrichtung mit der vorhergesagten Schädlichkeit? (Stufe 17, `work/42_*`)

Festgelegt am 20.10.2026, **bevor** AlphaMissense-Scores mit Assoziationsdaten verknüpft wurden. Bekannt: Stufen 1–16, insbesondere T10-8b
(variantenbasiert kodierend: GH p 0,031, CHIP p 0,009, FMS p 0,14). Nicht bekannt: jede Zahl dieses Tests.

## Frage
Wenn Träger eine milde Form der Krankheit sind, weil ihre Gen-Dosis sinkt, dann muss der Effekt **mit dem Funktionsverlust der Variante skalieren**:
harmlose Missense ≈ 0, schädliche Missense → in Krankheitsrichtung, nahe am pLoF-Effekt. Das prüft den Mechanismus innerhalb derselben Gene und Merkmale.

## Daten
- Assoziationen je Variante (Alternativallel, BETA/SE): Genes & Health variant-level (hg38), FinMetSeq (hg19), Exom-Chip-Verbund vor UKB (hg19); exakt die Reader aus `25`.
- AlphaMissense (Cheng et al. 2023, Zenodo 8208688), kanonische Transkripte, hg19/hg38 passend zur Quelle; Schlüssel Chrom-Pos-Ref-Alt.
  Steht das AlphaMissense-Alt als Datei-Ref, wird z umgedreht (auf das Missense-Allel orientiert).
- pLoF-Varianten wie 35 (gnomAD-LOFTEE-HC, auf LoF-Allel orientiert) als Referenzpunkt „volle halbe Dosis“.
- Nur MAF < 1 %; |z| ≤ √50. Gen-Zuordnung über das AlphaMissense-Transkript → Gen-Symbol (HGNC) bzw. CDS wie 25.
- Gene × Merkmale: gen-eigene HPO-Richtungskarte aus 31, **unverändert**. d_v = s_erwartet · z_v.

## Statistik
- **Primär T17-1 (innerhalb Zelle):** Regression d_v ~ am_pathogenicity mit festem Effekt je (Quelle × Gen × Merkmal); Steigung β.
  Null: am_pathogenicity **innerhalb** jeder Zelle permutieren (B = 1.000), einseitig (Erwartung β > 0). Alle drei Quellen gepoolt. Kriterium p < 0,05.
- **T17-2 (Klassen, beschreibend + Test):** Mittel d in AlphaMissense-Klassen likely benign (< 0,34), ambiguous, likely pathogenic (> 0,564), dazu pLoF;
  Erwartung monoton benign < ambiguous < pathogenic ≤ pLoF; Test pathogenic − benign innerhalb Zelle wie T17-1.
- **Negativkontrolle T17-3:** dieselbe Steigung mit zwischen Genen permutierten Richtungskarten (B = 200) → Erwartung 0. Zusätzlich die Steigung
  ohne Richtungsvorzeichen auf |z| (beschreibend: schädlichere Varianten haben größere Effekte, unabhängig von der Richtung).
- Sekundär: je Quelle getrennt; nur Gene mit ≥ 5 Missense-Varianten in der Zelle.
- Nebenzeile (beschreibend, selektiert): Genebass-UKB-Einzelvarianten mit p < 1e-4 (API, nur diese sind ungebinnt): Anteil Krankheitsrichtung je AlphaMissense-Klasse.

## Lesart vorab
- β > 0 und Monotonie: Dosis-Wirkung innerhalb der Träger; stützt „Dosis“ als Mechanismus und macht AlphaMissense zu einem Prädiktor des Trägerphänotyps.
- Null: Entweder zu wenig Power auf Variantenebene (33k/19k-Kohorten) oder die Richtung entsteht nicht über Dosis. Power über Zahl der Varianten berichten.

---

# Nachtrag 20.10.2026: T17 mit UKB-Einzelvarianten aus Backman et al. 2021 (Regeneron, 450k Exome; `work/43`)

Festgelegt **vor** dem Laden der Daten. Bekannt: Ergebnis von T17 in GH/FMS/CHIP (T17-1 p 0,032; T17-3-Kontrolle fast gleich groß;
post hoc bidirektionale Merkmale p 0,08).

- Quelle: GWAS Catalog, PMID 34662886, Einzelvarianten-Studien (nicht „Gene-based burden“), UKB_Freeze_450_EUR, REGENIE, Merkmale RINT, hg38,
  CC0. Zuordnung über die UKB-Feldnummer (identisch zu Genebass: height 50, sitheight 20015, Blutbild 30000–30120, Biochemie 306xx–308xx,
  fev1 3063, fvc 3062, grip 47, fluidint 20016, reaction 20023, sbp 4080, ebmd 3148). Fehlt ein Feld, fällt das Merkmal weg.
- **UKB überlappt mit Genebass** (Entdeckungsdaten der Richtungskarte sind aber HPO, nicht UKB). Für T17 ist das unproblematisch, weil T17 eine neue
  Frage ist (Dosis innerhalb der Träger), keine Replikation.
- Statistik **identisch zu 42** (T17-1 bis T17-3, Code unverändert, nur Quelle „BKM“ ergänzt). Zusätzlich vorab:
  **T17-4 (jetzt Primärkriterium für Gen-Spezifität):** nur Merkmale mit Genen in beiden Richtungen; Steigung gegen Permutation der Richtungen
  zwischen Genen **innerhalb des Merkmals** (B = 500), einseitig, Kriterium p < 0,05. Das ist die post hoc gefundene Kontrolle, jetzt an neuen Daten vorab.
