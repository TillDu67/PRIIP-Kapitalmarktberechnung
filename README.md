# PRIIP Kapitalmarktsimulation — Masterarbeit 2026

Dieses Repository enthält den im Rahmen der Masterarbeit

> **„Analyse und Implementierung von Kapitalmarktsimulationen unter Verwendung des Standardverfahrens für PRIIP der Kategorie 4"**  
> Till Duchaczek, Otto-von-Guericke Universität Magdeburg, 2026

entwickelten R-Code. Die Arbeit entstand in Kooperation mit den **Öffentlichen Versicherungen Sachsen-Anhalt (ÖSA)**.

---

## Inhalt

| Skript | Beschreibung |
|---|---|
| `main.R` | Zentrales Steuerungsskript — koordiniert den vollständigen Berechnungsablauf |
| `parameter_einlesen.R` | Liest Modellparameter aus der Excel-Eingabedatei |
| `vorberechnungen.R` | Berechnet pfadunabhängige Größen: ψ(t), V(t,T), kumulierte ψ-Integrale |
| `simulation.R` | Vektorisierte Monte-Carlo-Simulation der Basisprozesse x(t), y(t), S(t) |
| `renditen.R` | Berechnung von Swaprate K_d(t), Mischrendite und Fondsentwicklung |
| `output.R` | Schreibt Ergebnisse in Excel-Ausgabedatei |
| `szenario_matrix.R` | Implementierung des Szenario-Matrix-Ansatzes nach Günther & Hieber (2024) |
| `vergleich_mc_sm.R` | Vergleich Monte-Carlo vs. Szenario-Matrix: Konvergenz, Laufzeit, Genauigkeit |

---

## Modell

Das implementierte Kapitalmarktmodell entspricht dem **PIA-Modell** (Produktinformationsstelle Altersvorsorge) für die PRIIP-Kategorie 4:

- **Zinsdynamik:** Zwei-Faktor-Vasicek-Modell (G2++ Modell) mit Nelson-Siegel-Svensson-Zinskurve
- **Aktiendynamik:** Black-Scholes-Modell
- **Simulation:** Euler-Maruyama-Schema, 10.000 Pfade, 480 Monatschritte

Zusätzlich wurde der **Szenario-Matrix-Ansatz** nach Günther und Hieber (2024) implementiert und mit der Monte-Carlo-Simulation verglichen.

---

## Voraussetzungen

- R ≥ 4.3.3
- Pakete: `readxl`, `openxlsx`, `mvtnorm`

Installation der Pakete:
```r
install.packages(c("readxl", "openxlsx", "mvtnorm"))
```

---

## Ausführung

Die Eingabedaten (PIA-Parameter, Modellparameter) werden aus einer Excel-Datei gelesen, die aus Vertraulichkeitsgründen nicht im Repository enthalten ist. Der Pfad zur Eingabedatei wird in `main.R` gesetzt:

```r
pfad_input <- "Pfad/zur/Eingabedatei.xlsx"
art_override <- "Deckungsstock"  # oder "Fonds"
```

Anschließend wird der vollständige Durchlauf gestartet:

```r
source("main.R")
```

### Nur Szenario-Matrix-Ansatz

Der Szenario-Matrix-Ansatz kann nach einem vollständigen `main.R`-Durchlauf separat aufgerufen werden:

```r
source("szenario_matrix.R")
sm <- sm_berechnung(params, vb, N = 10, delta_grid = 8.25,
                    T_jahre = 25, g_garantie = 0.005, alpha = 0.465)
```

### MC vs. SM Vergleich

```r
source("szenario_matrix.R")
source("vergleich_mc_sm.R")
vgl <- vergleiche_alle(params, vb, ap, res,
                       N_werte = c(5, 8, 10, 15, 20),
                       delta_g = 8.25)
```

---

## Laufzeiten (Intel Core i7-7700K, 32 GB RAM)

### VBA vs. R (Kapitalmarktsimulation)

| Schritt | VBA | R | Faktor |
|---|---|---|---|
| Kapitalmarktsimulation (Deckungsstock, d=11) | ~20 Min. | ~85 Sek. | ~14× |
| Kapitalmarktsimulation (Fonds, d=6) | ~21 Min. | ~73 Sek. | ~17× |
| Vertragsberechnung (44 Varianten) | ~350 Min. | ~12 Min. | ~29× |

### Monte-Carlo vs. Szenario-Matrix (Basisprozesse)

| Verfahren | Laufzeit | Rel. Fehler |
|---|---|---|
| MC Euler-Maruyama (480 Monatschritte, VBA-Standard) | ~20 Min. | Referenz |
| MC Euler-Maruyama vektorisiert (480 Monatschritte, R) | ~1,8 Sek. | Referenz |
| MC exakt, Jahresgitter (Lemma 5.1.1) | ~2 Sek. | Referenz |
| SM N=5 | 0,4 Sek. | 20,2 % |
| SM N=8 | 1,5 Sek. | 38,5 % |
| SM N=10 | 3,0 Sek. | 3,1 % |
| SM N=15 | 13 Sek. | 15,1 % |
| SM N=20 | 41 Sek. | 16,2 % |

> **Hinweis:** Der Szenario-Matrix-Ansatz ersetzt nur die Simulation der
> Basisprozesse (~2 Sek. in R). Der dominante Rechenaufwand liegt in der
> Swapratenberechnung K_d(t) (~73 Sek.), die vom SM-Ansatz nicht berührt wird.


---

## Literatur

- Günther, S. & Hieber, P. (2024): *Efficient simulation and valuation of equity-indexed annuities under a two-factor G2++ model.* European Actuarial Journal.
- Graf, S. & Korn, R. (2020): *A guide to Monte Carlo simulation concepts for assessment of risk-return profiles for regulatory purposes.* European Actuarial Journal.
- Produktinformationsstelle Altersvorsorge GmbH (2024): *PIA-Kapitalmarktmodell.*

---

## Lizenz

Die zugrundeliegenden Modelldaten und Parameter unterliegen den Vertraulichkeitsbestimmungen der ÖSA Versicherungen Sachsen-Anhalt.
