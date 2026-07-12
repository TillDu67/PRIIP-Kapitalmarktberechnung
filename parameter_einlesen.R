# parameter_einlesen.R
# Liest alle Modellparameter aus dem "Input"-Sheet der Excel-Datei.
#
# Wichtig: Deckungsstock und Fonds haben unterschiedliche Werte fuer
# sigma_I, Ueberrendite, Anteil_Aktie, Duration, F_0/F_1/F_2, ...
#
# readxl liefert für solche Formelzellen nur den gecachten Wert (so wie
# Excel die Datei zuletzt gespeichert hat) - das entspricht aber nur einem
# der beiden "Art"-Faelle. Deshalb lesen wir hier direkt die statischen
# H- bzw. I-Spalten (in der Excel-Input Datei) und wählen selbst anhand von Art (Zelle C31) aus.
#
# art_override kann "Deckungsstock" oder "Fonds" sein, um die Auswahl zu
# erzwingen, unabhaengig davon was in C31 steht (z.B. um beide Fälle aus
# derselben Datei zu rechnen).

lese_parameter <- function(pfad_input, art_override = NULL) {
  
  input <- readxl::read_excel(pfad_input, sheet = "Input", col_names = FALSE)
  
  # Hilfsfunktionen: Zelle(Zeile, Spalte) -- Spalten als Zahl (A=1, B=2, ..., H=8, I=9)
  zelle <- function(row, col) {
    val <- input[row, col][[1]]
    as.numeric(val)
  }
  zelle_chr <- function(row, col) {
    val <- input[row, col][[1]]
    as.character(val)
  }
  
  # Art bestimmen
  art <- if (!is.null(art_override)) art_override else zelle_chr(31, 3)  # C31
  
  if (!(art %in% c("Deckungsstock", "Fonds"))) {
    warning(paste0("Unerwarteter Art-Wert: '", art, "' - erwarte 'Deckungsstock' oder 'Fonds'"))
  }
  
  # IF(Art="Fonds", I<row>, H<row>) -- H=Spalte 8, I=Spalte 9
  art_abhaengig <- function(row) {
    if (art == "Fonds") zelle(row, 9) else zelle(row, 8)
  }
  
  sigma_I_N <- art_abhaengig(32)  # \sigma_I
  
  params <- list(
    # Allgemein
    Art         = art,
    Pfade       = zelle(1, 3),
    Pfad_Output = zelle_chr(2, 3),
    Jahr        = zelle(3, 3),
    
    # Zins-/Aktienmodell (Art-unabhaengig)
    l_x     = zelle(10, 3),
    l_y     = zelle(11, 3),
    d_x     = zelle(12, 3),
    d_y     = zelle(13, 3),
    tau     = zelle(14, 3),
    a_x     = zelle(15, 3),
    sigma_x = zelle(16, 3),
    a_y     = zelle(17, 3),
    sigma_y = zelle(18, 3),
    rho_xy  = zelle(19, 3),
    beta_0  = zelle(20, 3),
    beta_1  = zelle(21, 3),
    beta_2  = zelle(22, 3),
    beta_3  = zelle(23, 3),
    tau_1   = zelle(24, 3),
    tau_2   = zelle(25, 3),
    
    # Nelson-Siegel-Svensson / lineare Extrapolation
    LinearBeginn = zelle(27, 3),
    Linear_Wert  = zelle(28, 3),
    delta        = zelle(29, 3),
    
    # Art-abhängige Parameter (Deckungsstock vs. Fonds)
    sigma_s   = art_abhaengig(9),
    sigma_I_N = sigma_I_N,
    sigma_I_S = 1.5 * sigma_I_N,   # Excel: sigma_I_Stress = 1.5 * sigma_I
    lambda_N  = art_abhaengig(34), # Überrendite
    lambda_S  = 0,                 # Überrendite_Stress (immer 0)
    K_f       = art_abhaengig(39), # Portfoliotransaktionskosten/Fondskosten
    Ant_Akt   = art_abhaengig(40), # Anteil der Aktien am Deckungsstock
    Duration  = art_abhaengig(36), # d
    
    F_0_oK = art_abhaengig(41),
    F_0    = art_abhaengig(42),
    F_1_oK = art_abhaengig(43),
    F_1    = art_abhaengig(44),
    F_2_oK = art_abhaengig(45),
    F_2    = art_abhaengig(46),
    
    Initialer_Fondskurs = zelle(47, 3)
  )
  
  params$pfad_input <- pfad_input
  params$n_pfade     <- params$Pfade
  
  return(params)
}