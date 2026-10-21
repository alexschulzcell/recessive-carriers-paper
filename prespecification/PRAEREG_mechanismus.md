# Vorab-Festlegung: Mechanismus-Regel über Organe (StilleTraeger, Stufe 4)

Festgelegt am 19.10.2026, **bevor** eine Klassen-Analyse außerhalb des Skeletts gerechnet wurde. Danach wird nichts an Hypothese, Klassen, Merkmalen oder Test geändert;
Abweichungen werden als solche markiert.

**Was vorher bekannt war (ehrlich):** Beim Skelett ist die Regel aus den Daten abgeleitet (Matrix/Sekretion dosisempfindlich, Zilie/Lysosom still; 05/10/11).
Die 25 stärksten Gene der Organ-Diagonale habe ich gesehen (u. a. LCAT, LIPC, ABCG5, TMPRSS6, TFR2, CLDN10, SLC22A12, ABCC2). Klassenzugehörigkeit oder
Klassenvergleiche außerhalb des Skeletts habe ich nicht berechnet.

## Hypothese
Rezessivität „hält“ (Träger still), wenn das Genprodukt in einem gepufferten Fluss mit Reservekapazität arbeitet: katabole/lysosomale Enzyme, Zilien-Maschinerie.
Rezessivität „bricht“ (Träger zeigen die Krankheit abgeschwächt), wenn das Genprodukt in großer Menge hergestellt, sekretiert oder strukturell verbaut wird:
sekretorischer und struktureller Hochdurchsatz.

## Klassen (GO, Nachfahren eingeschlossen, keine Handauswahl; Reihenfolge = Vorrang bei Mehrfachzugehörigkeit)
1. **SEKRET** = extrazelluläre Region (GO:0005576, CC) ODER strukturelle Molekülaktivität (GO:0005198, MF) ODER ER→Golgi-Transport (GO:0006888)
   ODER ECM-Organisation (GO:0030198) ODER Kollagen-/GAG-/Proteoglykan-Biosynthese (GO:0032964, GO:0006024, GO:0030166)
2. **KATABOL** = Lysosom (GO:0005764, CC) ODER katabolischer Prozess (GO:0009056, BP)
3. **ZILIE** = Zilie (GO:0005929, CC)
4. **andere**

## Gene, Organe, Merkmale
Rein rezessive Mendel-Gene mit Organ-Profil (HPO-Anteil ≥ 0,2; wie 14/20), Merkmale je Organ wie in 14 (alle zugeordneten Merkmale).
Genebass UKB, Träger ≥ 10, keep-Flags. Organ-passender Gen-Wert = Mittel des doppelt zentrierten, **locus-bereinigten** z² (V3 aus 20) über die passenden Merkmale.

## Primärtest
- **Datensatz:** pLoF, **ohne das Skelett**. Skelett-Organ und Skelett-Merkmale werden entfernt, damit der Test unabhängig von der Herleitung ist.
- **Statistik:** Mittelwert des organ-passenden Gen-Werts in SEKRET minus KATABOL.
- **Signifikanz:** Permutation der Klassen-Etiketten innerhalb des Hauptorgans (Organ mit größtem HPO-Anteil), B = 10.000, einseitig (SEKRET > KATABOL).
- **Vorhersage:** Differenz > 0, p < 0,05.

## Sekundär
- Missense|LC, gleicher Test.
- SEKRET gegen ZILIE und SEKRET gegen andere.
- Je Organ (Blut, Immun, Niere, Leber, Lipide, Lunge, Muskel, Nerven, Kreislauf), wo je Klasse ≥ 5 Gene vorhanden sind.
- Regression: Gen-Wert ~ Klasse + log10 Träger + LOEUF + Hauptorgan.
- Skelett allein als Bestätigung der Herleitung (nicht unabhängig).
- Kacser-Burns direkt: Enzym (EC laut HGNC) gegen Nicht-Enzym über alle Organe (Erwartung nach 09: kein Unterschied).

## Lesart
- **Primärtest positiv:** Die Regel verallgemeinert über Organe.
- **Primärtest negativ:** Die Regel ist skelettspezifisch (Wachstumsfuge als Hochdurchsatz-Matrixorgan). Das ist ebenfalls eine Aussage, nur eine kleinere.
