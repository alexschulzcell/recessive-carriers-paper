# Vorab festgelegt: Organ-Diagonale auf klinischen Diagnosen (Stufe 15, `work/39_*`)

Festgelegt am 20.10.2026, **bevor** Diagnose-Daten geladen wurden. Bekannt: Stufen 1–13 (Laborwerte/Messgrößen), Stufe 1 (Arthrose/Osteoporose/Fraktur null).

## Frage
Haben Träger rezessiver Krankheitsgene ein erhöhtes Risiko für **klinische Diagnosen im Organ ihrer rezessiven Krankheit**? Bisher waren alle
Endpunkte Messgrößen. Diagnosen machen aus „milde Form der Krankheit“ eine klinisch sichtbare Aussage (biologisch gerahmt, keine klinische Validierung).

## Daten
Genebass (UKB 394k Exome) ICD-10 first occurrence (UKB-Feld 130000ff.), dreistellige Codes mit ≥ 1.000 Fällen, Sets pLoF, missense|LC, synonymous.
Organ-Zuordnung der ICD-Blöcke a priori (vor dem Laden): Blut D50–D77; Immun D80–D89; Niere N00–N29; Leber K70–K83; Lipide E78; Glukose E10–E16;
Lunge J40–J99; Muskel G70–G73, M60–M63; Nerven G00–G69, G74–G99, F00–F03, F70–F79; Skelett M15–M19, M40–M43, M80–M94, Q65–Q79; Kreislauf I00–I99.
Ergibt 122 Codes (`work/39_icd_auswahl.tsv`). Gene und Organ-Etiketten wie Stufe 2 (rein rezessiv, HPO-Anteil ≥ 0,2).

## Statistik
- **Gerichtet** (Diagnose = Krankheitsrichtung, Erwartung: Träger häufiger krank): signiertes z = sign(BETA_Burden)·|z|, Gene ≥ 10 Träger, |z| ≤ √50.
- Zentrierung nach Code und Gen (wie 30/31). Gen-Wert je Organ = Mittel der zentrierten z über die Codes des Organs.
- Statistik: Mittel der passenden Gen×Organ-Zellen; Null: Organ-Etiketten zwischen Genen permutieren (B = 10.000), einseitig.
- **Primär:** pLoF, p < 0,05. Missense zweite Zeile, synonym Negativkontrolle (Erwartung null).
- Sekundär: richtungsfrei (z²) wie Stufe 2; leave-one-organ-out; ohne Kreislauf (größter Block, schwächster Bezug).
- Beschreibend: Odds Ratio je Träger für die stärksten Gen-Diagnose-Paare.

## Lesart vorab
- Positiv: Trägerschaft zeigt sich bis auf Diagnose-Ebene im Krankheitsorgan.
- Null: Effekte bleiben subklinisch (passt zu h ≈ 1 %), ehrliche Aussage mit Power-Angabe.
