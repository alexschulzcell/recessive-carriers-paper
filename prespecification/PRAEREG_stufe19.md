> **Excerpt.** This is an excerpt of the pre-specification document of October 20, 2026, limited to the parts used in the
> paper (P2 and P3, for the cohorts analyzed there). Omitted text is marked [...]; the remaining text is verbatim.
> SHA-256 of the complete original document: d7be89f9ac8ed9721bc76cd187a0865ba4f3f6b3fd8c986b87a091b7862e219e
> (available to editors and reviewers on request).

# Vorab festgelegt: Stufe 19 – fünf Härtungen für Nature Genetics (`work/45_*`)

Festgelegt am 20.10.2026, **bevor** die jeweiligen Zahlen berechnet wurden. [...]

## P1 [...]

## P2 – Konsistenz über Kohorten verschiedener Abstammung [...]
Gen-Wert = Mittel von s_erwartet · z über die Zellen der Richtungskarte (31), pLoF, je Kohorte: UKB (Genebass, europäisch), Genes & Health
(REGENIE Maske A, südasiatisch), [...]. Spearman der Gen-Werte [...] UKB–G&H [...] über Gene mit Daten in beiden;
Null: Gen-Etiketten in einer Kohorte permutieren (B = 2.000), einseitig. Kriterium: alle drei ρ > 0, mindestens zwei mit p < 0,05.

## P3 – Hefe und Mensch auf Gen-Ebene
Orthologe Mensch–Hefe aus Alliance (ORTHOLOGY-ALLIANCE_COMBINED, IsBestScore = Yes und IsBestRevScore = Yes), SGD-ID → ORF über SGD_features.
Hefe-Wert: Gen-Diagonale aus `work/40_hefe_gene.tsv` (Hoepfner, primär); Hillenmeyer als zweite Zeile (Gen-Diagonale aus 41).
Mensch-Wert: Gen-Wert aus P2 (UKB pLoF) [...]. Nur rein rezessive menschliche Gene mit Richtungskarte.
Test: Spearman Hefe-Wert ~ Mensch-Wert, einseitig positiv, Permutation (B = 5.000). Kriterium p < 0,05 mit UKB [...].
Lesart: Positiv = dieselben Gene sind in beiden Arten „halb empfindlich“. Null wird als „Analogie, aber keine Gen-Konkordanz“ berichtet.

## P4 [...]

## P5 [...]
