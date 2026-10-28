# Vorab festgelegt: Unabhängige Replikation in FinnGen (Stufe 21, `work/49_*`)

Festgelegt am 28.10.2026, **bevor** eine einzige FinnGen-Merkmalsdatei geladen wurde.
Bekannt: Dateiliste und Manifest der Buckets `finngen-public-data-r12/lab_values/` und `finngen-public-data-r13/summary_stats/`,
die Kanta-Dokumentation (`FinnGen_kanta_lab_documentation_v2.pdf`, `Kanta_labs_GWAS_results_v2_summary.txt`: Laborwerte, n, Einheiten,
Ausreißerfilter), alle bisherigen Ergebnisse (UKB, All of Us, Genes & Health, FinMetSeq, Exom-Chip).
Nicht bekannt: jeder Assoziationswert aus FinnGen (kein Beta, kein P, keine Variantenliste).

Anlass: Die All-of-Us-Replikation (Stufe 18/19) beruht auf Abfragen des öffentlichen All-by-All-Browsers, dessen Download die
All-of-Us-FAQ nicht erlaubt. Sie wird über die Researcher Workbench neu bezogen. FinnGen ist frei verfügbar und ersetzt sie nicht,
sondern ist eine weitere unabhängige Kohorte.

## Daten
FinnGen, öffentliche GWAS-Summary-Statistiken (Registrierung 28.10.2026; Bedingungen: Danksagung „We want to acknowledge the participants
and investigators of the FinnGen study“, Zitat Kurki et al., Nature 2023). Finnische Abstammung, GRCh38, imputiert (SISU-Panel),
REGENIE v3.3, Kovariaten Alter, Geschlecht, Genotypisierungs-Batch, PC1–10; quantitative Merkmale inverse-rang-normalisiert (SD-Einheit).
FinnGen filtert selbst: MAC > 30 (Laborwerte) und Imputations-INFO > 0,6. Diese Filter werden übernommen, nicht verändert.
**Keine UKB-Überlappung.**

Merkmale (FinnGen-ID; Auswahlregel: das quantitative, rangnormalisierte Merkmal mit demselben Analyten in Serum/Plasma bzw. Blut,
bei mehreren das mit dem größten n; festgelegt nur nach Dokumentation):

| Merkmal | Quelle | ID | n (Doku) |
|---|---|---|---|
| height | R13 Endpunkt | HEIGHT_IRN | 364.515 |
| hb | R12 Kanta | 3000963 | 447.271 |
| mcv | R12 Kanta | 3024731 | 447.040 |
| rbc | R12 Kanta | 3026361 | 440.562 |
| plt | R12 Kanta | 3007461 | 441.342 |
| wbc | R12 Kanta | 3010813 | 440.693 |
| lymph | R12 Kanta | 3019198 | 183.481 |
| creat | R12 Kanta | 3020564 | 426.104 |
| cysc | R12 Kanta | 3030366 | 10.125 |
| urate | R12 Kanta | 3026493 | 81.374 |
| urea | R12 Kanta | 3020779 | 31.957 |
| phosphate | R12 Kanta | 3003458 | 46.763 |
| alt | R12 Kanta | 3006923 | 413.996 |
| ast | R12 Kanta | 3013721 | 124.818 |
| ggt | R12 Kanta | 3026910 | 257.813 |
| bili | R12 Kanta | 3006140 (gesamt) | 190.407 |
| albumin | R12 Kanta | 3024561 | 153.962 |
| ldl | R12 Kanta | 3001308 | 368.976 |
| hdl | R12 Kanta | 3023602 | 361.195 |
| tg | R12 Kanta | 3048773 (nüchtern) | 347.199 |
| chol | R12 Kanta | 3019900 | 365.484 |
| glucose | R12 Kanta | 3018251 (nüchtern) | 349.417 |
| hba1c | R12 Kanta | 3004410 | 332.429 |

Nicht vorhanden: sitheight, ebmd, sbp, fev1, fvc, grip, fluidint, reaction. Tests, die diese Merkmale bräuchten, laufen ohne sie.

Laden (`work/49a_finngen_download.py`): Datei als Strom lesen, SHA-256 der Originaldatei protokollieren, gespeichert werden nur Varianten mit
MAF < 1 %, die pLoF sind oder in einem kodierenden Exon (± 2 bp) liegen. Erwartete Spalten: `#chrom pos ref alt … beta sebeta af_alt`.
Weichen die Namen ab, wird nach Bedeutung zugeordnet und das im Log vermerkt; sonst ändert sich nichts. Chromosom 23 wird als X gelesen.

## Varianten und Gen-Statistiken (unverändert aus Stufe 7 und 13)
- z_i = beta/sebeta bezogen auf das Alternativallel; MAF = min(af_alt, 1 − af_alt) < 1 %; |z_i| ≤ √50.
- pLoF: gnomAD-v2-LOFTEE-HC aus der GeneBayes-Datei (`gnomad_v2_lofs.tsv.gz`, SHA 310be244…), Abgleich über GRCh38-Position und Allele in
  beiden Orientierungen; ist das Datei-Alternativallel das Referenzallel, wird z auf das LoF-Allel umgedreht (wie 35).
- Kodierend: RefSeq-NM_-CDS-Exon ± 2 bp (`refGene_hg38.txt.gz`, SHA 58fbc18d…), wie 25.
- Ungerichtet: Σz_i² ~ χ²(k) → z² je Gen (Kappe 50), wie 25. Gerichtet: z_gen = Σz_i/√k (Stouffer), wie 35.
- Nur für F8: β_het je Gen = inversvarianzgewichtetes Mittel der pLoF-Betas (auf das LoF-Allel orientiert), SD-Einheit.
- Organzuordnung der Merkmale a priori wie 14/19 (Blut: hb, mcv, rbc, plt; Immun: wbc, lymph; Niere: creat, cysc, urate, urea, phosphate;
  Leber: alt, ast, ggt, bili, albumin; Lipide: ldl, hdl, tg, chol; Glukose: glucose, hba1c; Skelett: height). In 25 stand nur die Teilmenge
  der dort vorhandenen Merkmale; für FinnGen wird diese Zuordnung auf alle vorhandenen Merkmale angewendet.

## Tests (Code unverändert; neu ist nur der Einleser)
Für jeden Test: pLoF primär, kodierend als zweite Zeile. Einseitig, Permutation und Zufallsstartwert wie im Ursprungsskript
(25: 23, 35: 35, 45: 45).

- **F1 Größe (25 T1):** AR_nos gegen gematchten pool und gegen AR_nonskel, z², B = 5.000. Kriterium: beide p < 0,05 (pLoF).
  Kontrolle: AR_nonskel gegen pool, Erwartung null.
- **F2 Organ-Diagonale (25 T2):** doppelt zentriertes z², Permutation der Organ-Etiketten, B = 2.000. Kriterium p < 0,05 (pLoF).
- **F3 Regel W1/W2b (25 T3):** Klassen aus `22_einheiten_plof.tsv` (SHA 585e264d…), ohne Skelett, Permutation innerhalb Organ, B = 5.000.
  Kriterium je Route p < 0,05 (pLoF).
- **F4 Gen-eigene Richtung (35, T10-8b):** Richtungskarte aus 31 unverändert (Skript-SHA d4c61868…), Stouffer-z, Zentrierung und Statistik
  wie 35, Permutation der Richtungskarten zwischen Genen, B = 1.000. **Primärkriterium der Stufe: F4 mit pLoF p < 0,05.**
  Sekundär: Anteil richtig bei |z| > 2.
- **F5 Eingefrorene Vorhersage (P1):** `ATLAS_VORHERSAGE_2026-10-20.tsv` (SHA-256 7694042a…, eingefroren 20.10.2026 08:37).
  d = richtung_erwartet · z_FinnGen (Stouffer, pLoF), Statistik Mittel d, Null: Richtungen innerhalb desselben Merkmals zwischen Genen
  permutiert, B = 2.000. Kriterium p < 0,05 (pLoF). Sekundär: Anteil d > 0 bei |z| > 2 (B = 500), Spezifitäts-Tertile (Erwartung Tertil 3 am
  stärksten). Anders als in All of Us ist hier auch das Gesamtergebnis F4 bei der Festlegung unbekannt.
- **F6 Konsistenz über Kohorten (P2):** Gen-Wert = Mittel s · z über die Richtungszellen (pLoF). Spearman FinnGen–UKB (Kriterium ρ > 0,
  p < 0,05) und FinnGen–G&H (zweite Zeile), Permutation der Gen-Etiketten, B = 2.000. Das Paar FinnGen–All of Us kommt mit derselben
  Statistik hinzu, sobald die Workbench-Daten vorliegen.
- **F7 Skelett-Mechanismus (P4):** GO-Klassen aus `go_kategorien.tsv` (SHA 0ed94db9…), AR_nos, Größe: Anteil |z| > 1,96 MATRIX gegen
  CILIUM + LYSO, Fisher einseitig. Kriterium p < 0,05 (pLoF).
- **F8 h für Größe (P5):** β_hom aus `h_literatur.tsv` (SHA 521f5a68…, nach 32b korrigiert), h = β_het/β_hom je Gen. Median und IQR,
  Vorzeichentest h > 0 (binomial, einseitig), Wilcoxon h < 0,5; Konkordanz mit UKB-h (Spearman, `29b_h_haupttest.tsv`).
  Kriterium: Vorzeichentest p < 0,05.

Mindestgrößen: F1 ≥ 20 AR_nos-Gene, F2/F4 ≥ 20 Zellen (wie 25/35), F5 ≥ 100 Zellen, F7 ≥ 5 Gene je Klasse, F8 ≥ 10 Gene.
Darunter heißt der Test „nicht prüfbar“, nicht „gescheitert“.

## Überlappung und Kombination
FinnGen enthält FINRISK-Teilnehmer der THL-Biobank; FinMetSeq und Teile der Exom-Chip-Konsortien ebenfalls. FinnGen wird deshalb **für sich**
berichtet und **nicht** per Fisher mit FMS oder Exom-Chip kombiniert. Eine Fisher-Kombination nur mit Genes & Health (keine Überlappung) ist
erlaubt und wird als zweite Zeile berichtet.

## Lesart
- **Repliziert:** F4 mit pLoF p < 0,05.
- **Teilweise repliziert:** F4 nur kodierend p < 0,05.
- **Nicht repliziert:** F4 mit beiden Varianten-Sets p ≥ 0,05.
F1–F3 und F5–F8 werden mit ihren Kriterien berichtet, unabhängig vom Ausgang. Gescheiterte Tests kommen ins Testregister und ins Manuskript.

## Vorab benannte Grenzen
- Variantenebene und Imputation: Sehr seltene pLoF-Varianten fehlen (MAC > 30, INFO > 0,6). Erwartet werden weniger Gene als in einer
  Exom-Kohorte, mit finnisch angereicherten Varianten.
- Stouffer ohne LD-Bereinigung ist im Betrag antikonservativ; gelesen wird gegen die Permutation, die dieselbe Korrelation trägt (wie 35).
- Keine synonyme Negativkontrolle auf Variantenebene (wie 25/35); die Kontrolle ist die Permutation und in F1 der Vergleich AR_nonskel gegen pool.
- Keine Aufteilung nach Geschlecht (Summary-Statistiken).

## Abweichungen
Jede Abweichung von diesem Dokument wird in FINDINGS mit Grund und Zeitpunkt vermerkt und im Manuskript genannt. Merkmale, Varianten-Sets,
Kriterien und Richtungskarte werden nach dem Laden nicht neu gewählt.

## Prüfsummen der Methoden (Stand der Festlegung)
`25_replikation_unabhaengig.py` 3fc66f6a…; `31_vorzeichen_genspezifisch.py` d4c61868…; `35_t108b_gerichtet.py` 48465d43…;
`45_haertungen.py` 0516eaec…; `14_organmatrix.py` f14760aa…; `genannotation.tsv.gz` b052ec97…; `organ_anteile_AR.tsv.gz` 9b3c85b7….
