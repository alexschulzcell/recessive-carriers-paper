# Vorab festgelegt: Unabhängige Replikation in BioBank Japan (Stufe 21, Nachtrag zu `PRAEREG_finngen.md`, `work/49_*`)

Festgelegt am 28.10.2026, **bevor** eine BBJ-Merkmalsdatei geladen wurde.
Bekannt: Download-Liste von PheWeb.jp (Merkmale, n, Dateinamen), die NBDC-Seiten hum0014/hum0197 und die Methoden von Sakaue et al.
(Nat Genet 2021; Imputation Minimac3 mit BBJ-WGS n = 1.037 + 1KGP3, Rsq ≥ 0,7, 13,5 Mio. Varianten, GRCh37).
Zum Zeitpunkt der Festlegung läuft der FinnGen-Download (`49a`); **kein FinnGen-Ergebnis ist gerechnet oder angesehen.**
Nicht bekannt: jeder Assoziationswert aus BBJ.

Zweck: eine zweite große, frei verfügbare, UKB-unabhängige Kohorte anderer Abstammung (ostasiatisch, Japan). Sie prüft zugleich, ob ein
FinnGen-Ergebnis an finnisch angereicherten Gründervarianten hängt: BBJ trägt andere seltene Varianten in denselben Genen.

## Daten
BioBank Japan, GWAS-Summary-Statistiken aus hum0197 (PheWeb.jp „Downloads“, nur die **BBJ-Spalte**, nie die Meta-Analysen mit UKB oder
FinnGen), offen ohne Einschränkung (Sakaue 2021: „openly shared without any restrictions“; NBDC-Richtlinie). Zitat Sakaue et al. 2021.
Imputiert, GRCh37, lineare Regression; keine UKB-, FinnGen- oder G&H-Überlappung.

| Merkmal | Datei | n |
|---|---|---|
| height | hum0197.v3.BBJ.Hei.v1.zip | 165.056 |
| hb | hum0197.v3.BBJ.Hb.v1.zip | 152.447 |
| mcv | hum0197.v3.BBJ.MCV.v1.zip | 129.832 |
| rbc | hum0197.v3.BBJ.RBC.v1.zip | 153.512 |
| plt | hum0197.v3.BBJ.PLT.v1.zip | 148.623 |
| wbc | hum0197.v3.BBJ.WBC.v1.zip | 154.355 |
| lymph | hum0197.v3.BBJ.LYM.v1.zip | 95.717 |
| creat | hum0197.v3.BBJ.sCr.v1.zip | 150.266 |
| urate | hum0197.v3.BBJ.UA.v1.zip | 129.405 |
| urea | hum0197.v3.BBJ.BUN.v1.zip | 148.767 |
| phosphate | hum0197.v5.gwas.P.v1.txt.gz | 46.848 |
| alt | hum0197.v3.BBJ.ALT.v1.zip | 150.545 |
| ast | hum0197.v3.BBJ.AST.v1.zip | 150.068 |
| ggt | hum0197.v3.BBJ.GGT.v1.zip | 133.471 |
| bili | hum0197.v3.BBJ.TBil.v1.zip | 124.341 |
| albumin | hum0197.v3.BBJ.Alb.v1.zip | 120.539 |
| ldl | hum0197.v3.BBJ.LDLC.v1.zip | 72.866 |
| hdl | hum0197.v3.BBJ.HDLC.v1.zip | 74.970 |
| tg | hum0197.v3.BBJ.TG.v1.zip | 111.667 |
| chol | hum0197.v3.BBJ.TC.v1.zip | 135.808 |
| glucose | hum0197.v3.BBJ.Glucose.v1.zip | 133.336 |
| hba1c | hum0197.v3.BBJ.HbA1c.v1.zip | 71.221 |
| sbp | hum0197.v3.BBJ.SBP.v1.zip | 145.505 |

Nicht vorhanden: cysc, sitheight, ebmd, fev1, fvc, grip, fluidint, reaction. Phosphat stammt aus einer späteren Freigabe (v5); es wird
genommen wie geliefert.

Laden (`work/49b_bbj_download.py`): Archiv laden, SHA-256 protokollieren, Spalten nach Bedeutung zuordnen (Chromosom, Position, REF, ALT,
Frequenz des ALT-Allels, BETA, SE; Name im Log), gespeichert werden nur Varianten mit MAF < 1 %, die pLoF (GRCh37-Position) sind oder in einem
kodierenden Exon ± 2 bp (refGene hg19) liegen; danach wird das Archiv gelöscht. Chromosom 23 wird als X gelesen. Fehlt eine Spalte,
bricht das Skript ab, statt zu raten.

## Varianten, Gen-Statistiken, Tests
**Alles wie in `PRAEREG_finngen.md`**, mit GRCh37-Abgleich (`lof19`, `refGene_hg19`, wie FinMetSeq in 25/35). Tests B1–B8 entsprechen F1–F8
mit denselben Kriterien und Mindestgrößen. Organzuordnung zusätzlich sbp → Kreislauf (wie 14/19).
**Primärkriterium für BBJ: B4 (gen-eigene Richtung) mit pLoF P < 0,05.**
B8 (h) nur, wenn die Quelle (README im Archiv oder Sakaue 2021) Größe als rangnormalisiert bzw. in SD-Einheiten ausweist; sonst „nicht prüfbar“.

## Abstammungsübergreifende Lesart („Gründer-Bias“)
1. F4 (FinnGen) und B4 (BBJ) werden **getrennt** berichtet.
2. **X1:** Spearman der Gen-Werte (Mittel s · z über die Richtungszellen, pLoF) FinnGen gegen BBJ, Permutation der Gen-Etiketten (B = 2.000,
   einseitig), Kriterium ρ > 0, P < 0,05. Zusätzlich BBJ–UKB und BBJ–G&H wie F6.
3. **X2:** Fisher-Kombination der F4/B4-P-Werte mit Genes & Health (35, variantenbasiert) über die drei überlappungsfreien Kohorten
   FinnGen, BBJ, G&H, je Varianten-Set. Zweite Zeile.
4. Lesart: Repliziert die Richtung in FinnGen **und** BBJ (beide pLoF oder beide kodierend P < 0,05), hängt das Ergebnis nicht an einer
   Abstammung. Repliziert sie nur in FinnGen, wird das als mögliche Folge finnischer Gründervarianten berichtet.

## Vorab benannte Grenzen
- Imputation mit Rsq ≥ 0,7: Sehr seltene Varianten fehlen weitgehend; in BBJ sind weniger pLoF-Gene zu erwarten als in FinnGen.
- Quantitative Merkmale in Sakaue 2021 nicht durchgehend als rangnormalisiert belegt; z = BETA/SE ist davon unabhängig, β_het (B8) nicht.
- Lipide und HbA1c mit n ≈ 71.000–75.000 sind deutlich kleiner.
