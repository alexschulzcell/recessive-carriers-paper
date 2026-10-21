# Vorab-Festlegung 2: Organ-eigenes Hauptprodukt (StilleTraeger, Stufe 5)

Festgelegt am 19.10.2026, **bevor** gerechnet wurde. Neuer, eigener Test; kein Umbau von Stufe 4 (dort Primärtest negativ, bleibt so stehen).

**Vorwissen (ehrlich):** Ich kenne die 25 stärksten Gene der Organ-Diagonale (u. a. TMPRSS6, TFR2, TFRC → Blut; LCAT, LIPC, ABCG5 → Lipide;
CLDN10, SLC22A12 → Niere; ABCC2 → Leber) und die Klassenmittel aus Stufe 4 (breite GO-Klassen). Die Klassen hier werden ausschließlich über GO-Begriffe bzw.
HPA-Expression definiert, nicht über Gennamen. Die HPA-Datei habe ich nur auf Format geprüft (Spalten, Gewebeliste).

## Hypothese
Innerhalb eines Organs zeigen Träger rezessiver Gene dann einen Effekt im Organ-Merkmal, wenn das Gen am **Hauptprodukt dieses Organs mit hohem Durchsatz** beteiligt ist
(Erythrozyt: Hämoglobin/Häm/Eisen; Leber: Plasmaproteine, Gallensäuren, Bilirubin; Lipide: Lipoprotein-Aufbau und -Umbau; Niere: tubulärer Transport;
Muskel: Sarkomer; Lunge: Surfactant/Gasaustausch). Andere rezessive Gene desselben Organs sind gepuffert.

## Einheit und Wert
Einheit = Gen × Organ (rein rezessive Mendel-Gene, HPO-Organanteil ≥ 0,2; wie 14/20/21). Wert = Mittel des doppelt zentrierten, locus-bereinigten z² (V3)
über die Merkmale **dieses** Organs (Zentrierung über alle 31 Merkmale). Genebass UKB, Träger ≥ 10.

## Wege (Klassen „im Hauptprodukt“ ja/nein, je Organ)
- **W1 GO-Hauptprodukt** (Nachfahren eingeschlossen):
  - Blut: GO:0042541 Hämoglobin-Biosynthese, GO:0006783 Häm-Biosynthese, GO:0055072 Eisen-Homöostase, GO:0030218 Erythrozyten-Differenzierung, GO:0030220 Plättchenbildung
  - Immun: GO:0002521 Leukozyten-Differenzierung, GO:0002377 Immunglobulin-Produktion, GO:0042582 azurophile Granula, GO:0042581 spezifische Granula
  - Leber: GO:0072562 Blut-Mikropartikel (Plasmaproteine), GO:0006699 Gallensäure-Biosynthese, GO:0015721 Gallensäure-Transport, GO:0042167 Häm-Abbau
  - Lipide: GO:0042157 Lipoprotein-Stoffwechsel, GO:0034358 Plasma-Lipoproteinpartikel, GO:0030301 Cholesterin-Transport
  - Niere: GO:0070293 renale Absorption, GO:0097254 tubuläre Sekretion, GO:0015747 Urat-Transport, GO:0003094 glomeruläre Filtration, GO:0035435 Phosphat-Transport
  - Muskel: GO:0030017 Sarkomer, GO:0006936 Muskelkontraktion
  - Lunge: GO:0043129 Surfactant-Homöostase, GO:0007585 Gasaustausch
  - Kreislauf: GO:0008217 Blutdruck-Regulation
  - Nerven: GO:0099504 synaptischer Vesikelzyklus
  - Glukose: GO:0030073 Insulinsekretion
- **W2a Expression absolut:** HPA-nTPM im Organgewebe im oberen Tertil unter den AR-Genen dieses Organs.
  Gewebezuordnung: Blut → Knochenmark; Immun → max(Knochenmark, Milz, Lymphknoten); Leber → Leber; Lipide → max(Leber, Dünndarm); Niere → Niere;
  Muskel → Skelettmuskel; Lunge → Lunge; Kreislauf → max(Herzmuskel, Blutgefäß); Nerven → Großhirnrinde; Glukose → Pankreas.
- **W2b Expression organ-spezifisch:** Verhältnis nTPM(Organgewebe) / Median nTPM über alle Gewebe im oberen Tertil (je Organ).
- **W3 Kombination:** W1 UND W2a.

## Test
- **Primär:** pLoF, **ohne Skelett**, alle anderen Organe gepoolt.
- **Statistik:** Mittel(Wert | Klasse ja) − Mittel(Wert | Klasse nein).
- **Signifikanz:** Permutation der Klassen-Etiketten innerhalb des Organs, B = 10.000, einseitig.
- **Mehrfachtest:** 4 Wege → Bonferroni α = 0,0125 je Weg.
- **Replikation:** Missense|LC (andere Varianten, weitgehend andere Träger), gleiche Richtung mit p < 0,05.
- **Sekundär:**
  - je Organ;
  - robuste Statistik (Rangtest innerhalb Organ, Mittel der Ränge);
  - ohne die 25 Top-Gene der Diagonale (Liste aus 20);
  - Skelett mit der MATRIX-Definition als Referenz (nicht unabhängig).

## Lesart
- **Mindestens ein Weg primär signifikant und in Missense repliziert:** Die Hauptprodukt-Regel verallgemeinert.
- **Robustheit:** Hält das Ergebnis ohne die Top-25 Gene, trägt es nicht nur an Lehrbuchgenen.
- **Alle Wege negativ:** Die Regel bleibt skelettspezifisch (Abschluss dieser Linie).
