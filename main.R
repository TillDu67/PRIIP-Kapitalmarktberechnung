# main.R
# Hauptskript -- PRIIP Kapitalmarktberechnung 

library(readxl)
library(openxlsx)

source("parameter_einlesen.R")
source("vorberechnungen.R")
source("simulation.R")
source("renditen.R")
source("output.R")

# --- Pfad zur Eingabedatei ---
pfad_input <- "C:/Users/TillD/OneDrive/Desktop/PRIIP-20260610T153054Z-3-001/PRIIP/2024/2024_PRIIP_Kapitalmarktberechnung ohne Makros.xlsx"

# art_override:
#   NULL            -> Art wird aus Zelle "Art" (C31) in Excel gelesen
#   "Deckungsstock" -> erzwingt Deckungsstock-Parameter (H-Spalte in Input-Datei)
#   "Fonds"         -> erzwingt Fonds-Parameter (I-Spalte in Input-Datei)
art_override <- "Deckungsstock"

# --- Parameter einlesen ---
cat("Lese Parameter ein...\n")
params <- lese_parameter(pfad_input, art_override = art_override)
params$Pfad_Output <- "C:/Users/TillD/OneDrive/Desktop/PRIIP-20260610T153054Z-3-001/PRIIP/R Output"
art <- params$Art
cat("Art:", art, "| Pfade:", params$n_pfade, "| Duration:", params$Duration, "\n\n")

# --- Vorberechnungen (psi, V_t_T, int_psi_t_T) ---
cat("Starte Vorberechnungen...\n")
start <- Sys.time()

vb <- berechne_vorberechnungen(
  a_x = params$a_x, a_y = params$a_y,
  sigma_x = params$sigma_x, sigma_y = params$sigma_y,
  rho_xy = params$rho_xy,
  beta_0 = params$beta_0, beta_1 = params$beta_1,
  beta_2 = params$beta_2, beta_3 = params$beta_3,
  tau_1 = params$tau_1, tau_2 = params$tau_2,
  delta = params$delta,
  linear_beginn = params$LinearBeginn,
  linear_wert = params$Linear_Wert
)

cat("Vorberechnungen fertig:", round(difftime(Sys.time(), start, units = "secs"), 1), "Sekunden\n\n")

# --- Simulation (x, y, r, I_N_oK, I_S_oK, [I_N, I_S]) ---
cat("Starte Simulation...\n")
start <- Sys.time()

ap <- simuliere_pfade(params, vb, art = art)

cat("Simulation fertig:", round(difftime(Sys.time(), start, units = "secs"), 1), "Sekunden\n\n")

# --- Renditen (K_swap, R_B_d, Mischrendite / Fondsentwicklung) ---
cat("Starte Renditenberechnung...\n")
start <- Sys.time()

res <- berechne_renditen(params, vb, ap, art = art)

cat("Renditen fertig:", round(difftime(Sys.time(), start, units = "secs"), 1), "Sekunden\n\n")

# --- Output ---
schreibe_output(res$output, art = art, pfad_output = params$Pfad_Output)

cat("Alles fertig!\n")

# ============================================================
# --- Szenario-Matrix vs. Monte-Carlo Vergleich ---
# Drei Analysen:
#   1. SM-Konvergenz (N=5,8,10,15,20)
#   2. MC auf EIA-Basis als sauberer Referenzwert
#   3. Bezug zur bestehenden PRIIP-MC
# ============================================================
source("szenario_matrix.R")
source("vergleich_mc_sm.R")

vgl <- vergleiche_alle(
  params     = params,
  vb         = vb,
  ap         = ap,
  res        = res,
  N_werte    = c(5, 8, 10, 15, 20,25,30,35),
  delta_g    = 8.25,     # Standardwert aus Günther & Hieber (2024)
  T_jahre    = 25,
  g_garantie = 0.005,
  alpha      = 0.465,
  rho1       = -0.15,
  rho2       =  0.15,
  n_mc       = 10000
)