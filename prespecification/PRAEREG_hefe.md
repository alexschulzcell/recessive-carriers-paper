# Vorab festgelegt: Bedingungs-Diagonale in der Hefe – zeigt der heterozygote Stamm seinen Defekt dort, wo der homozygote leidet? (Stufe 16, `work/40_*`)

Festgelegt am 20.10.2026, **bevor** gerechnet wurde. Bekannt: Format der Dateien, Vorzeichen (negativ = empfindlich; TUB2 unter Benomyl Rang 5/6.681),
Matrixgrößen. Nicht bekannt: jede Zahl dieses Tests.

## Analogie zum Menschen
| Mensch | Hefe |
|---|---|
| biallelischer Patient | homozygote Deletion (HOP) |
| Träger | heterozygote Deletion (HIP) |
| Krankheitsorgan | Bedingungen (Substanz × Konzentration), in denen die homozygote Deletion stark empfindlich ist |
| Organ-Diagonale | heterozygote Empfindlichkeit gerade in diesen Bedingungen |

Vorteil: Experiment mit definierten Nullallelen, isogen, keine LD-, Stratifikations- oder Ascertainment-Probleme.

## Daten
Hoepfner et al. 2014 (Microbiol Res; Dryad doi:10.5061/dryad.v5m8v), `HIP_scores.txt`, `HOP_scores.txt` (Gene × Experimente, adjustierte Scores,
negativ = Wachstumsdefekt). Experimente gepaart über den Schlüssel „Exp. <Substanz>_<Konzentration>“ (HIP- und HOP-Spalte derselben Substanz und Dosis).
Gene: nicht-essentiell, mit HIP- und HOP-Daten.

## Definitionen
- Empfindlichkeit s = −Score (positiv = empfindlich).
- **„Krankheitsbedingungen“ eines Gens:** gepaarte Experimente mit HOP-Score ≤ −3. Gene mit mindestens 3 solchen Experimenten.
- Heterozygoter Wert: HIP-Empfindlichkeit, doppelt zentriert (Experiment-Mittel, dann Gen-Mittel), wie beim Menschen.

## Statistik
- **Primär (T16-1):** mittlere zentrierte HIP-Empfindlichkeit in den Krankheitsbedingungen des Gens. Null: Bedingungsmengen zwischen Genen permutieren
  (B = 1.000), einseitig. Kriterium p < 0,05.
- **Richtung (T16-2):** Anteil der Krankheitsbedingungen mit HIP-Score < 0 (heterozygot ebenfalls empfindlich), gegen dieselbe Permutation.
- **Dosis (T16-3, beschreibend):** „Hefe-h“ = Median über Gene von (mittlere HIP-Empfindlichkeit in Krankheitsbedingungen − sonst) /
  (mittlere HOP-Empfindlichkeit dort − sonst). Vergleich mit dem Menschen (h ≈ 0,01 für Größe, Untergrenze).
- **Regel (T16-4):** Die menschliche Regel „halbe Dosis zeigt sich, wo das Gen spezifisch arbeitet“, übertragen: Gene mit **wenigen** Krankheitsbedingungen
  (bedingungsspezifisch) gegen Gene mit vielen (pleiotrop). Tertile der Anzahl HOP-Treffer; Vorhersage: Diagonale (Gen-Wert aus T16-1) im spezifischen
  Tertil größer als im pleiotropen. Test: Spearman zwischen Gen-Diagonale und Trefferzahl (Erwartung negativ) plus Tertil-Vergleich (Mann-Whitney).
- Sensitivität: Schwelle −2 und −4; nur „Ad.“-Scores; Ausschluss von Experimenten mit > 5 % HOP-Treffern (unspezifisch toxisch).

## Lesart vorab
- T16-1 positiv: Die Diagonale gilt artübergreifend, im kausalen Experiment. Spezieskonvergenz Hefe–Mensch.
- T16-4 positiv: Auch die Regel überträgt sich.
- Null: Die Bedingungs-Diagonale der Hefe ist keine Analogie zum Organ (Bedingungen sind keine Gewebe); als solches zu berichten.

---

# Nachtrag 20.10.2026: Unabhängige Replikation mit Hillenmeyer et al. 2008 (Stanford, Giaever/Nislow; T16-R, `work/41`)

Festgelegt **vor** dem Rechnen. Bekannt: Format (z-Scores, **positiv = empfindlich**; ERG11 het unter Fluconazol Rang 2/5.915, TUB3 hom unter Benomyl Rang 8/4.742),
het 5.915 Gene × 726 Experimente (332 Bedingungen), hom 4.742 × 418 (161 Bedingungen). Quelle: Web-Supplement, erhalten im Internet Archive
(Snapshot 07.12.2015, `het.z_result_nm.pub`, `hom.z_result_nm.pub`), auf E unter `hillenmeyer_roh/`.

Anders als Hoepfner wurden HIP und HOP **getrennt** gemessen (eigene Experimente, eigene Tage). Paarung daher über die **Substanz** (Bedingungsname),
nicht über das Experiment. Je Gen × Substanz: Mittel der z über die Experimente der Substanz.

- „Krankheitssubstanzen“ eines Gens: hom-Mittel z ≥ 3. Gene mit ≥ 2 solchen Substanzen, mit het- und hom-Daten.
- T16-R1: mittlere doppelt zentrierte het-Empfindlichkeit in den Krankheitssubstanzen; Null: Substanzmengen zwischen Genen permutieren (B = 1.000), einseitig. Kriterium p < 0,05.
- T16-R2: Anteil het z > 0 in Krankheitssubstanzen gegen dieselbe Permutation.
- T16-R3 (beschreibend): h = mittleres Δhet / mittleres Δhom (z-Skala; wie 40c nicht als Fitness-Dominanz interpretierbar).
- Sensitivität: Schwelle 2 und 4; Kontrolle Messtag: het-Experimente vom selben Tag wie ein hom-Treffer-Experiment, aber andere Substanz (Erwartung null).
- Gen-Konkordanz mit Hoepfner: Spearman der Gen-Diagonalen (Gene in beiden Datensätzen), Erwartung positiv.
