# Vorab festgelegt: Sagt die Regel den Selektionsdruck voraus – ohne ein einziges Phänotypdatum? (Stufe 11, `work/34_*`)

Festgelegt am 20.10.2026, **bevor** gerechnet wurde. Interne Festlegung, kein öffentlicher Zeitstempel.

## Frage
Die Regel aus Stufe 5/9A lautet: Halbe Dosis zeigt sich beim Träger dort, wo das Gen **organspezifisch** in seinem Krankheitsorgan arbeitet.
Wenn das stimmt, sollten genau diese rezessiven Gene in der Population auch **stärker gegen heterozygote LoF selektiert** sein. Das lässt sich
aus der Allelhäufigkeit in gnomAD prüfen, ganz ohne UKB-Merkmale – eine dritte, unabhängige Messebene (Populationsgenetik).

## Einheit und Daten
- Gene: rein rezessive Gene aus `organ_anteile_AR.tsv.gz` (wie Stufe 5), **Hauptorgan** = Organ mit dem größten HPO-Anteil, Anteil ≥ 0,2.
- Ohne Skelett (HPA hat keinen Knorpel/Knochen; wie 22/27). Organ→Gewebe und Organ→Produzentenzelle exakt wie in `22_hauptprodukt.py` (TIS) und `27_zelltyp_W4.py` (PROD).
- Spezifität im Krankheitsorgan: `spez` = log2((nTPM im Organgewebe + 1)/(Median über alle Gewebe + 1)) wie 22; Zelltyp `zspez` wie 27.
- Allgemeine Gewebespezifität als Störgröße: `tau` (Yanai) über alle HPA-Gewebe, plus log10 Median-Expression.
- Nur Gene mit exp_lof ≥ 10 (LOEUF sonst unzuverlässig).

## Zielgrößen
- **Primär: LOEUF** (gnomAD v2.1.1, `oe_lof_upper`). Rein aus LoF-Häufigkeit, enthält keine Expressionsmerkmale. gnomAD v2 enthält keine UKB.
- Sekundär: s_het (GeneBayes, Zeng 2024). **Achtung Zirkularität:** Der GeneBayes-Prior nutzt Genmerkmale, darunter Expression. Nur als Sekundärgröße.
- Sekundär: Punktschätzer obs/exp LoF.

## Statistik
- Regression: rank(−LOEUF) ~ spez + tau + log10 Median-Expression + log10 CDS-Länge + Organ (fest). Teststatistik: Koeffizient von `spez`.
- Null: Organ-Etiketten zwischen Genen permutieren (B = 10.000), `spez` neu berechnen, Koeffizient neu schätzen. p einseitig (Erwartung positiv).
- Primärkriterium: Gewebe-`spez`, LOEUF, p < 0,05.
- Sekundär: Zelltyp-`zspez` (gleiche Logik); s_het; und als Brücke zum Phänotyp: Korrelation des Trägereffekts (locus-bereinigter Organwert aus 22) mit −LOEUF innerhalb der Organe.

## Lesart vorab
- Positiv: Die Regel sagt aus der Expression allein voraus, welche rezessiven Gene von der Selektion „als halb-dominant“ behandelt werden. Verbindet Stufe 5 mit Judd 2026/Fridman 2021.
- Null: Die Regel beschreibt den Phänotyp, aber die Selektion folgt anderer Pleiotropie (vereinbar mit Stufe 2, Test 18). Auch das ist eine Antwort.
- Negativ (spezifische Gene weniger constrained): spricht für die bekannte Regel „gewebespezifisch = weniger essenziell“, die wir dann nicht ausreichend herausadjustiert hätten.

---

# Nachtrag 20.10.2026: Richtung × Spezifität (Stufe 12, `work/36`, `work/37`)

**Post-hoc-Befund in der UKB (aus dem Atlas, `36`):** Bei Zellen mit |z| > 2 zeigen Trägereffekte in Gewebe-Spezifitäts-Tertil 3
zu 0,677 in Krankheitsrichtung (n 257), in Tertil 1 zu 0,513 (n 199); Permutation der Tertile innerhalb Organ p 0,0002; synonym flach.
Das war **nicht** vorab festgelegt.

**Vorab festgelegte Replikation in Genes & Health 44k (vor dem Rechnen):** dieselbe Richtungskarte (31), dieselben Spezifitäts-Tertile
(aus `ATLAS_traeger.tsv`, unverändert), signierte REGENIE-z wie 33. Statistik: Anteil Krankheitsrichtung in Tertil 3 minus Tertil 1,
Zellen mit |z| > 2; Null durch Permutation der Tertile innerhalb Organ (B = 5000), einseitig. Primär Maske C (pLoF+Missense, mehr Zellen),
Maske A zweite Zeile, Maske D (synonym) Negativkontrolle. Kriterium: Maske C p < 0,05 und synonym nicht signifikant.
Erwartung: Richtung gleich, wegen 33k deutlich weniger Zellen; „nicht repliziert" wird so berichtet.
