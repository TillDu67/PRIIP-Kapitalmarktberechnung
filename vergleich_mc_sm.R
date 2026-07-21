# # vergleich_mc_sm.R
# # Vollständiger Vergleich: Szenario-Matrix vs. Monte-Carlo
# #
# # Drei Analysen:
# #   Teil 1 — Konvergenz des SM (E[Z_T] als Funktion von N)
# #   Teil 2 — MC auf gemeinsamer EIA-Basis (sauberer Referenzwert)
# #   Teil 3 — Bezug zur bestehenden PRIIP-MC (näherungsweise)
# #
# # Voraussetzung: params, vb, ap, res aus main.R vorhanden.
# # Aufruf:
# #   source("szenario_matrix.R")
# #   source("vergleich_mc_sm.R")
# #   vgl <- vergleiche_alle(params, vb, ap, res)
# 
# # ============================================================
# # HAUPTFUNKTION
# # ============================================================
# vergleiche_alle <- function(params, vb, ap, res,
#                             N_werte    = c(5, 8, 10, 15, 20),
#                             delta_g    = 8.25,
#                             T_jahre    = 25,
#                             g_garantie = 0.005,
#                             alpha      = 0.465,
#                             rho1       = -0.15,
#                             rho2       =  0.15,
#                             n_mc       = 10000) {
#   
#   cat("============================================================\n")
#   cat("Vergleich: Szenario-Matrix vs. Monte-Carlo\n")
#   cat(sprintf("T=%d Jahre, g=%.3f, alpha=%.3f, delta=%.2f\n",
#               T_jahre, g_garantie, alpha, delta_g))
#   cat("============================================================\n\n")
#   
#   # Teil 2: MC auf EIA-Basis (Referenzwert)
#   cat("--- Teil 2: Monte-Carlo Referenz (EIA-Struktur) ---\n")
#   mc <- mc_eia(params, vb, T_jahre, g_garantie, alpha,
#                rho1, rho2, n_mc = n_mc)
#   cat(sprintf("  E[Z_T]_MC  = %.6f\n",  mc$E_ZT))
#   cat(sprintf("  Std.err.   = %.6f\n",  mc$se))
#   cat(sprintf("  95%%-KI     = [%.6f, %.6f]\n", mc$CI_low, mc$CI_high))
#   cat(sprintf("  Laufzeit   = %.1f Sek.\n\n", mc$laufzeit))
#   
#   # Teil 1: SM Konvergenzsweep
#   cat("--- Teil 1: SM Konvergenz über N ---\n")
#   sm_ergebnisse <- vector("list", length(N_werte))
#   for (k in seq_along(N_werte)) {
#     N <- N_werte[k]
#     cat(sprintf("  N=%d ...", N))
#     sm <- sm_berechnung(params, vb,
#                         N          = N,
#                         delta_grid = delta_g,
#                         T_jahre    = T_jahre,
#                         g_garantie = g_garantie,
#                         alpha      = alpha,
#                         rho1       = rho1,
#                         rho2       = rho2,
#                         verbose    = FALSE)
#     fehler_pct <- abs(sm$E_ZT - mc$E_ZT) / abs(mc$E_ZT) * 100
#     im_KI      <- (sm$E_ZT >= mc$CI_low && sm$E_ZT <= mc$CI_high)
#     sm_ergebnisse[[k]] <- list(
#       N          = N,
#       E_ZT       = sm$E_ZT,
#       fehler_pct = fehler_pct,
#       laufzeit   = sm$laufzeit,
#       im_KI      = im_KI
#     )
#     cat(sprintf(" E[Z_T]=%.5f  Fehler=%.3f%%  %.1fs  KI:%s\n",
#                 sm$E_ZT, fehler_pct, sm$laufzeit,
#                 ifelse(im_KI, "JA", "NEIN")))
#   }
#   
#   # Teil 3: Bezug zur bestehenden PRIIP-MC
#   cat("\n--- Teil 3: Bezug zur bestehenden PRIIP-MC ---\n")
#   priip <- priip_naehering(res, params, T_jahre)
#   cat(sprintf("  E[Mischrendite kumuliert, %d J.] = %.6f\n",
#               T_jahre, priip$E_kumuliert))
#   cat(sprintf("  Hinweis: nicht direkt mit E[Z_T] vergleichbar\n"))
#   cat(sprintf("  (Mischrendite ≠ EIA-Struktur; kein Diskontierungsfaktor)\n\n"))
#   
#   # --- Ergebnistabelle ---
#   drucke_tabelle(sm_ergebnisse, mc)
#   
#   # --- Plots ---
#   plot_konvergenz_sm(sm_ergebnisse, mc, delta_g)
#   
#   return(list(mc = mc, sm = sm_ergebnisse, priip = priip))
# }
# 
# 
# # ============================================================
# # TEIL 2: MC auf gemeinsamer EIA-Basis
# # ============================================================
# # Simuliert n_mc Pfade mit exakter Jahres-Simulation (Lemma 2.1).
# # Berechnet E[Z_T] = E[prod_t max(e^g, S_t^alpha) * e^{-int_r_t}]
# #
# # 3-Schritt-Simulation:
# #   Schritt 1: (x_t, y_t) ~ N(mu_12, Sigma_2x2)
# #   Schritt 2: int_r | (x_t,y_t) ~ N(mu_bar, sigma_bar2)
# #              [sigma_bar2 ≈ 0 bei rho≈-1 -> näherungsweise deterministisch]
# #   Schritt 3: ln(S_t/S_{t-1}) | (x_t,y_t) ~ N(mu_tilde, sigma2_tilde)
# #              [unabhängig von int_r bedingt auf (x_t,y_t)]
# 
# mc_eia <- function(params, vb, T_jahre, g_garantie, alpha,
#                    rho1, rho2, n_mc = 10000) {
#   
#   start <- Sys.time()
#   
#   a  <- params$a_x;  b  <- params$a_y
#   nu <- params$sigma_x; eta <- params$sigma_y
#   rho_r    <- params$rho_xy
#   sigma    <- params$sigma_s
#   lambda_S <- params$lambda_S   # = 0 (Q-Maß)
#   d_x <- params$d_x; d_y <- params$d_y
#   l_x <- params$l_x; l_y <- params$l_y
#   tau <- params$tau
#   
#   # Sigma-Komponenten (Lemma 2.1)
#   sig <- berechne_sigma(a, b, nu, eta, rho_r, sigma, rho1, rho2)
#   S   <- sig$Sigma
#   S11 <- S[1,1]; S12 <- S[1,2]; S13 <- S[1,3]; S14 <- S[1,4]
#   S22 <- S[2,2]; S23 <- S[2,3]; S24 <- S[2,4]
#   S33 <- S[3,3]; S34 <- S[3,4]; S44 <- S[4,4]
#   ga1 <- sig$ga1; gb1 <- sig$gb1
#   
#   # 2x2-Cholesky für (x_t, y_t) — immer positiv definit
#   Sigma_xy <- matrix(c(S11, S12, S12, S22), 2, 2)
#   L_xy     <- t(chol(Sigma_xy))   # Sigma_xy = L %*% t(L)
#   
#   # Bedingte Verteilung von int_r und ln_S gegeben (x_t, y_t):
#   det_xy   <- S11 * S22 - S12^2
#   A_x      <- (S22*S13  - S12*S23)  / det_xy   # für int_r
#   A_y      <- (S11*S23  - S12*S13)  / det_xy
#   At_x     <- (S22*S14  - S12*S24)  / det_xy   # für ln_S (tilted)
#   At_y     <- (S11*S24  - S12*S14)  / det_xy
#   
#   # Bedingte Varianzen
#   sigma_bar2   <- max(S33 - (S22*S13^2 - 2*S12*S13*S23 + S11*S23^2)/det_xy, 0)
#   sigma2_tilde <- max(S44 - (S22*S14^2 - 2*S12*S14*S24 + S11*S24^2)/det_xy, 0)
#   sd_bar       <- sqrt(sigma_bar2)
#   sd_tilde     <- sqrt(sigma2_tilde)
#   
#   # Jährliche int_psi-Werte aus vb
#   int_psi_jahres <- numeric(T_jahre)
#   for (t in 1:T_jahre)
#     int_psi_jahres[t] <- vb$int_psi_t_T[1, 12*t+1] -
#     vb$int_psi_t_T[1, 12*(t-1)+1]
#   
#   exp_a <- exp(-a); exp_b <- exp(-b)
#   
#   # Simulation
#   set.seed(123)
#   ZT_vec <- numeric(n_mc)
#   
#   for (k in 1:n_mc) {
#     x_t <- 0.0; y_t <- 0.0
#     log_z <- 0.0
#     
#     for (t in 1:T_jahre) {
#       # Risikoprämien
#       t_mitte <- 12 * t - 6
#       if (t_mitte <= tau) { lx <- d_x; ly <- d_y
#       } else               { lx <- l_x; ly <- l_y }
#       
#       # Schritt 1: (x_t, y_t) simulieren
#       mu1 <- exp_a * x_t + lx * (1 - exp_a)
#       mu2 <- exp_b * y_t + ly * (1 - exp_b)
#       z2  <- L_xy %*% rnorm(2)
#       x_t1 <- mu1 + z2[1]
#       y_t1 <- mu2 + z2[2]
#       
#       # Schritt 2: int_r | (x_t1, y_t1)
#       mu3    <- int_psi_jahres[t] + ga1*x_t + gb1*y_t +
#         lx*(1-ga1) + ly*(1-gb1)
#       mu_bar <- mu3 + A_x*(x_t1 - mu1) + A_y*(y_t1 - mu2)
#       int_r  <- mu_bar + sd_bar * rnorm(1)   
#       
#       # Schritt 3: ln(S_t/S_{t-1}) | (x_t1, y_t1)
#       # Vereinfacht mit lambda_S=0:
#       mu1t     <- mu1 - S13  
#       mu2t     <- mu2 - S23
#       mu3t     <- mu3 + lambda_S - 0.5*sigma^2 -
#         (S33 + sigma*(sig$Delta))
#       mu_tilde <- mu3t + At_x*(x_t1 - mu1t) + At_y*(y_t1 - mu2t)
#       ln_S     <- mu_tilde + sd_tilde * rnorm(1)
#       
#       # Jährlicher Beitrag L_t = max(g, alpha*ln_S) - int_r
#       log_z <- log_z + max(g_garantie, alpha * ln_S) - int_r
#       
#       x_t <- x_t1; y_t <- y_t1
#     }
#     ZT_vec[k] <- exp(log_z)
#   }
#   
#   E_ZT <- mean(ZT_vec)
#   se   <- sd(ZT_vec) / sqrt(n_mc)
#   list(
#     E_ZT     = E_ZT,
#     se       = se,
#     CI_low   = E_ZT - 1.96*se,
#     CI_high  = E_ZT + 1.96*se,
#     ZT_vec   = ZT_vec,
#     laufzeit = as.numeric(difftime(Sys.time(), start, units="secs"))
#   )
# }
# 
# 
# # ============================================================
# # TEIL 3: Näherung aus bestehender PRIIP-MC
# # ============================================================
# priip_naehering <- function(res, params, T_jahre) {
#   
#   n_jahre   <- min(T_jahre, 40)
#   n_pfade   <- params$n_pfade
#   
#   if (params$Art == "Deckungsstock") {
#     # R_N_mat: n_pfade x 480 (monatlich, gleiche Werte pro Jahr)
#     # Jährliche Renditen: Spalten 12, 24, ..., 12*n_jahre
#     R_mat <- res$R_N_mat
#     idx   <- seq(12, 12*n_jahre, by=12)
#     R_j   <- R_mat[, idx, drop=FALSE]   # n_pfade x n_jahre
#     
#     # Kumuliertes Produkt der jährlichen Faktoren (1+R_t)
#     kum <- apply(1 + R_j, 1, prod)   # pro Pfad
#     E_kum <- mean(kum)
#     se_kum <- sd(kum) / sqrt(n_pfade)
#   } else {
#     F_mat <- res$F_N_mat
#     start <- params$Initialer_Fondskurs
#     kum   <- F_mat[, 12*n_jahre] / start
#     E_kum <- mean(kum)
#     se_kum <- sd(kum) / sqrt(n_pfade)
#   }
#   
#   list(
#     E_kumuliert = E_kum,
#     se          = se_kum,
#     CI_low      = E_kum - 1.96*se_kum,
#     CI_high     = E_kum + 1.96*se_kum,
#     art         = params$Art,
#     n_jahre     = n_jahre
#   )
# }
# 
# 
# # ============================================================
# # AUSGABE: Ergebnistabelle
# # ============================================================
# drucke_tabelle <- function(sm_ergebnisse, mc) {
#   cat("============================================================\n")
#   cat("Ergebnistabelle: SM vs. MC\n")
#   cat(sprintf("MC-Referenz: E[Z_T] = %.6f  95%%-KI [%.6f, %.6f]\n",
#               mc$E_ZT, mc$CI_low, mc$CI_high))
#   cat(sprintf("%-5s  %-12s  %-12s  %-10s  %-8s\n",
#               "N", "E[Z_T]_SM", "Fehler [%]", "Zeit [s]", "im KI"))
#   cat(strrep("-", 55), "\n")
#   for (e in sm_ergebnisse) {
#     cat(sprintf("%-5d  %-12.6f  %-12.4f  %-10.1f  %-8s\n",
#                 e$N, e$E_ZT, e$fehler_pct, e$laufzeit,
#                 ifelse(e$im_KI, "JA", "NEIN")))
#   }
#   cat("============================================================\n\n")
# }
# 
# 
# # ============================================================
# # PLOTS: Konvergenz
# # ============================================================
# plot_konvergenz_sm <- function(sm_ergebnisse, mc, delta_g) {
#   
#   N_vec   <- sapply(sm_ergebnisse, `[[`, "N")
#   E_vec   <- sapply(sm_ergebnisse, `[[`, "E_ZT")
#   Err_vec <- sapply(sm_ergebnisse, `[[`, "fehler_pct")
#   T_vec   <- sapply(sm_ergebnisse, `[[`, "laufzeit")
#   
#   par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
#   
#   # Plot 1: E[Z_T] vs N mit MC-Konfidenzband
#   plot(N_vec, E_vec, type = "b", pch = 16, col = "steelblue",
#        xlab = "Gitterpunkte N", ylab = "E[Z_T]",
#        main = sprintf("Konvergenz E[Z_T]\n(delta=%.2f)", delta_g),
#        ylim = range(c(E_vec, mc$CI_low, mc$CI_high)) * c(0.99, 1.01))
#   polygon(c(min(N_vec), max(N_vec), max(N_vec), min(N_vec)),
#           c(mc$CI_low, mc$CI_low, mc$CI_high, mc$CI_high),
#           col = rgb(0.2, 0.8, 0.2, 0.2), border = NA)
#   abline(h = mc$E_ZT, lwd = 2, lty = 1)
#   legend("bottomright", c("SM", "MC", "95%-KI"),
#          col = c("steelblue", "black", "lightgreen"),
#          lwd = c(2, 2, 8), pch = c(16, NA, NA), bty = "n")
#   
#   # Plot 2: Relativer Fehler vs N (log-log)
#   plot(N_vec, Err_vec, type = "b", pch = 16, col = "firebrick",
#        log = "xy",
#        xlab = "Gitterpunkte N (log)", ylab = "Relativer Fehler [%] (log)",
#        main = "Relativer Fehler vs. N\n(log-log)")
#   abline(h = 0.1, lty = 2, col = "grey50")
#   text(min(N_vec), 0.12, "0.1%-Schwelle", adj = 0, col = "grey50", cex = 0.8)
#   
#   # Quadratische Konvergenz einzeichnen (Referenzlinie)
#   if (length(N_vec) >= 2) {
#     N_ref <- N_vec[1]
#     E_ref <- Err_vec[1]
#     N_fit <- seq(min(N_vec), max(N_vec), length.out = 50)
#     lines(N_fit, E_ref * (N_fit/N_ref)^(-2), lty = 3, col = "grey30")
#     text(max(N_fit)*0.9, E_ref*(max(N_fit)/N_ref)^(-2)*1.5,
#          "O(N⁻²)", col = "grey30", cex = 0.8)
#   }
#   
#   # Plot 3: Laufzeit vs N
#   plot(N_vec, T_vec, type = "b", pch = 16, col = "darkorange",
#        xlab = "Gitterpunkte N", ylab = "Laufzeit [Sek.]",
#        main = "Laufzeit SM\n(O(N⁴))")
#   # O(N^4) Referenzlinie
#   if (length(N_vec) >= 2) {
#     N_ref <- N_vec[1]; T_ref <- T_vec[1]
#     N_fit <- seq(min(N_vec), max(N_vec), length.out = 50)
#     lines(N_fit, T_ref * (N_fit/N_ref)^4, lty = 3, col = "grey30")
#     text(max(N_fit)*0.85, T_ref*(max(N_fit)/N_ref)^4*1.3,
#          "O(N⁴)", col = "grey30", cex = 0.8)
#   }
#   # MC-Laufzeit als horizontale Referenz
#   abline(h = 85, lty = 2, col = "steelblue")
#   text(min(N_vec), 85*1.1, "MC (~85s)", col = "steelblue", cex = 0.8, adj = 0)
#   
#   par(mfrow = c(1, 1))
#   cat("Plots angezeigt.\n")
# }

# vergleich_mc_sm.R
# Vollständiger Vergleich: Szenario-Matrix vs. Monte-Carlo
#
# Drei Analysen:
#   Teil 1 — Konvergenz des SM (E[Z_T] als Funktion von N)
#   Teil 2 — MC auf gemeinsamer EIA-Basis (sauberer Referenzwert)
#   Teil 3 — Bezug zur bestehenden PRIIP-MC (näherungsweise)
#
# Voraussetzung: params, vb, ap, res aus main.R vorhanden.
# Aufruf:
#   source("szenario_matrix.R")
#   source("vergleich_mc_sm.R")
#   vgl <- vergleiche_alle(params, vb, ap, res)

# ============================================================
# HAUPTFUNKTION
# ============================================================
vergleiche_alle <- function(params, vb, ap, res,
                            N_werte    = c(5, 8, 10, 15, 20,25,30,35),
                            delta_g    = 8.25,
                            T_jahre    = 25,
                            g_garantie = 0.005,
                            alpha      = 0.465,
                            rho1       = -0.15,
                            rho2       =  0.15,
                            n_mc       = 10000) {
  
  cat("============================================================\n")
  cat("Vergleich: Szenario-Matrix vs. Monte-Carlo\n")
  cat(sprintf("T=%d Jahre, g=%.3f, alpha=%.3f, delta=%.2f\n",
              T_jahre, g_garantie, alpha, delta_g))
  cat("============================================================\n\n")
  
  # Teil 2: MC auf EIA-Basis (Referenzwert)
  cat("--- Teil 2: Monte-Carlo Referenz (EIA-Struktur) ---\n")
  mc <- mc_eia(params, vb, T_jahre, g_garantie, alpha,
               rho1, rho2, n_mc = n_mc)
  cat(sprintf("  E[Z_T]_MC  = %.6f\n",  mc$E_ZT))
  cat(sprintf("  Std.err.   = %.6f\n",  mc$se))
  cat(sprintf("  95%%-KI     = [%.6f, %.6f]\n", mc$CI_low, mc$CI_high))
  cat(sprintf("  Laufzeit   = %.1f Sek.\n\n", mc$laufzeit))
  
  # Teil 1: SM Konvergenzsweep
  cat("--- Teil 1: SM Konvergenz über N ---\n")
  sm_ergebnisse <- vector("list", length(N_werte))
  for (k in seq_along(N_werte)) {
    N <- N_werte[k]
    cat(sprintf("  N=%d ...", N))
    sm <- sm_berechnung(params, vb,
                        N          = N,
                        delta_grid = delta_g,
                        T_jahre    = T_jahre,
                        g_garantie = g_garantie,
                        alpha      = alpha,
                        rho1       = rho1,
                        rho2       = rho2,
                        verbose    = FALSE)
    fehler_pct <- abs(sm$E_ZT - mc$E_ZT) / abs(mc$E_ZT) * 100
    im_KI      <- (sm$E_ZT >= mc$CI_low && sm$E_ZT <= mc$CI_high)
    sm_ergebnisse[[k]] <- list(
      N          = N,
      E_ZT       = sm$E_ZT,
      fehler_pct = fehler_pct,
      laufzeit   = sm$laufzeit,
      im_KI      = im_KI
    )
    cat(sprintf(" E[Z_T]=%.5f  Fehler=%.3f%%  %.1fs  KI:%s\n",
                sm$E_ZT, fehler_pct, sm$laufzeit,
                ifelse(im_KI, "JA", "NEIN")))
  }
  
  # Teil 3: Bezug zur bestehenden PRIIP-MC
  cat("\n--- Teil 3: Bezug zur bestehenden PRIIP-MC ---\n")
  priip <- priip_naehering(res, params, T_jahre)
  cat(sprintf("  E[Mischrendite kumuliert, %d J.] = %.6f\n",
              T_jahre, priip$E_kumuliert))
  cat(sprintf("  Hinweis: nicht direkt mit E[Z_T] vergleichbar\n"))
  cat(sprintf("  (Mischrendite ≠ EIA-Struktur; kein Diskontierungsfaktor)\n\n"))
  
  # --- Ergebnistabelle ---
  drucke_tabelle(sm_ergebnisse, mc)
  
  # --- Plots ---
  plot_konvergenz_sm(sm_ergebnisse, mc, delta_g)
  
  return(list(mc = mc, sm = sm_ergebnisse, priip = priip))
}


# ============================================================
# TEIL 2: MC auf gemeinsamer EIA-Basis
# ============================================================
# Simuliert n_mc Pfade mit exakter Jahres-Simulation (Lemma 2.1).
# Berechnet E[Z_T] = E[prod_t max(e^g, S_t^alpha) * e^{-int_r_t}]
#
# Der Zustandsvektor (x_t, y_t, int_r, ln_S) wird direkt als
# gemeinsam 4-dimensional normalverteilt gezogen (volle Cholesky-
# Zerlegung von Sigma, vgl. Anhang C.0.9), nicht mehr ueber eine
# bedingte 3-Schritt-Zerlegung. Der Grund: int_r und ln_S sind
# beide vom selben Zinsintegral getrieben (siehe Gl. 4-6) und daher
# bedingt auf (x_t,y_t) nicht unabhaengig -- die alte Fassung hat
# genau das angenommen. Ausserdem wurden dort fuer ln_S die
# getilteten Parameter aus Satz 3.1 verwendet; die sind fuer die
# momenterzeugende Funktion des SM-Ansatzes richtig, nicht aber fuer
# eine direkte Simulation von Z_T unter dem Bewertungsmass selbst.

mc_eia <- function(params, vb, T_jahre, g_garantie, alpha,
                   rho1, rho2, n_mc = 10000) {
  
  start <- Sys.time()
  
  a  <- params$a_x;  b  <- params$a_y
  nu <- params$sigma_x; eta <- params$sigma_y
  rho_r <- params$rho_xy
  
  # sigma_s ist auch hier (wie in sm_berechnung, szenario_matrix.R) nur
  # der Normierungsfaktor aus der VBA-Referenz, nicht die tatsaechliche
  # Aktienvolatilitaet -- die ist sigma_I_N. lambda_N ist entsprechend
  # erst nach Normierung mit sigma_I_N/sigma_s die tatsaechliche
  # Ueberrendite. Beide Groessen muessen mit dem SM-Ansatz uebereinstimmen,
  # sonst vergleicht man zwei verschieden parametrisierte Modelle.
  sigma    <- params$sigma_I_N
  lambda_S <- params$lambda_N * params$sigma_I_N / params$sigma_s
  
  d_x <- params$d_x; d_y <- params$d_y
  l_x <- params$l_x; l_y <- params$l_y
  tau <- params$tau
  
  # Sigma-Komponenten (Lemma 2.1)
  sig <- berechne_sigma(a, b, nu, eta, rho_r, sigma, rho1, rho2)
  Sigma4 <- sig$Sigma
  ga1 <- sig$ga1; gb1 <- sig$gb1
  
  # Volle 4x4-Cholesky-Zerlegung. Sigma4 ist bei den ÖSA-Parametern
  # streng positiv definit (vgl. Anhang C.0.9) -- ein Sonderfall
  # ueber eine bedingte Zerlegung ist nicht mehr noetig.
  L <- t(chol(Sigma4))   # Sigma4 = L %*% t(L)
  
  # Jährliche int_psi-Werte aus vb
  int_psi_jahres <- numeric(T_jahre)
  for (t in 1:T_jahre)
    int_psi_jahres[t] <- vb$int_psi_t_T[1, 12*t+1] -
    vb$int_psi_t_T[1, 12*(t-1)+1]
  
  exp_a <- exp(-a); exp_b <- exp(-b)
  
  # Simulation -- ueber alle n_mc Pfade gleichzeitig, Schleife nur
  # noch ueber die T Simulationsjahre.
  set.seed(123)
  x     <- numeric(n_mc)
  y     <- numeric(n_mc)
  log_z <- numeric(n_mc)
  
  for (t in 1:T_jahre) {
    
    # Risikoprämien
    t_mitte <- 12 * t - 6
    if (t_mitte <= tau) { lx <- d_x; ly <- d_y
    } else               { lx <- l_x; ly <- l_y }
    
    # Bedingte Erwartungswerte nach Lemma 2.1, gegeben (x_{t-1}, y_{t-1})
    mu1 <- exp_a * x + lx * (1 - exp_a)
    mu2 <- exp_b * y + ly * (1 - exp_b)
    mu3 <- int_psi_jahres[t] + ga1*x + gb1*y +
      lx*(1-ga1) + ly*(1-gb1)
    mu4 <- mu3 + lambda_S - 0.5*sigma^2
    
    # Ein gemeinsamer 4-dim. Zug pro Pfad
    E <- L %*% matrix(rnorm(4 * n_mc), nrow = 4)
    
    x_neu <- mu1 + E[1, ]
    y_neu <- mu2 + E[2, ]
    int_r <- mu3 + E[3, ]
    ln_S  <- mu4 + E[4, ]
    
    # Jährlicher Beitrag L_t = max(g, alpha*ln_S) - int_r
    log_z <- log_z + pmax(g_garantie, alpha * ln_S) - int_r
    
    x <- x_neu; y <- y_neu
  }
  
  ZT_vec <- exp(log_z)
  E_ZT <- mean(ZT_vec)
  se   <- sd(ZT_vec) / sqrt(n_mc)
  list(
    E_ZT     = E_ZT,
    se       = se,
    CI_low   = E_ZT - 1.96*se,
    CI_high  = E_ZT + 1.96*se,
    ZT_vec   = ZT_vec,
    laufzeit = as.numeric(difftime(Sys.time(), start, units="secs"))
  )
}


# ============================================================
# TEIL 3: Näherung aus bestehender PRIIP-MC
# ============================================================
priip_naehering <- function(res, params, T_jahre) {
  
  n_jahre   <- min(T_jahre, 40)
  n_pfade   <- params$n_pfade
  
  if (params$Art == "Deckungsstock") {
    # R_N_mat: n_pfade x 480 (monatlich, gleiche Werte pro Jahr)
    # Jährliche Renditen: Spalten 12, 24, ..., 12*n_jahre
    R_mat <- res$R_N_mat
    idx   <- seq(12, 12*n_jahre, by=12)
    R_j   <- R_mat[, idx, drop=FALSE]   # n_pfade x n_jahre
    
    # Kumuliertes Produkt der jährlichen Faktoren (1+R_t)
    kum <- apply(1 + R_j, 1, prod)   # pro Pfad
    E_kum <- mean(kum)
    se_kum <- sd(kum) / sqrt(n_pfade)
  } else {
    F_mat <- res$F_N_mat
    start <- params$Initialer_Fondskurs
    kum   <- F_mat[, 12*n_jahre] / start
    E_kum <- mean(kum)
    se_kum <- sd(kum) / sqrt(n_pfade)
  }
  
  list(
    E_kumuliert = E_kum,
    se          = se_kum,
    CI_low      = E_kum - 1.96*se_kum,
    CI_high     = E_kum + 1.96*se_kum,
    art         = params$Art,
    n_jahre     = n_jahre
  )
}


# ============================================================
# AUSGABE: Ergebnistabelle
# ============================================================
drucke_tabelle <- function(sm_ergebnisse, mc) {
  cat("============================================================\n")
  cat("Ergebnistabelle: SM vs. MC\n")
  cat(sprintf("MC-Referenz: E[Z_T] = %.6f  95%%-KI [%.6f, %.6f]\n",
              mc$E_ZT, mc$CI_low, mc$CI_high))
  cat(sprintf("%-5s  %-12s  %-12s  %-10s  %-8s\n",
              "N", "E[Z_T]_SM", "Fehler [%]", "Zeit [s]", "im KI"))
  cat(strrep("-", 55), "\n")
  for (e in sm_ergebnisse) {
    cat(sprintf("%-5d  %-12.6f  %-12.4f  %-10.1f  %-8s\n",
                e$N, e$E_ZT, e$fehler_pct, e$laufzeit,
                ifelse(e$im_KI, "JA", "NEIN")))
  }
  cat("============================================================\n\n")
}


# ============================================================
# PLOTS: Konvergenz
# ============================================================
plot_konvergenz_sm <- function(sm_ergebnisse, mc, delta_g) {
  
  N_vec   <- sapply(sm_ergebnisse, `[[`, "N")
  E_vec   <- sapply(sm_ergebnisse, `[[`, "E_ZT")
  Err_vec <- sapply(sm_ergebnisse, `[[`, "fehler_pct")
  T_vec   <- sapply(sm_ergebnisse, `[[`, "laufzeit")
  
  par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
  
  # Plot 1: E[Z_T] vs N mit MC-Konfidenzband
  plot(N_vec, E_vec, type = "b", pch = 16, col = "steelblue",
       xlab = "Gitterpunkte N", ylab = "E[Z_T]",
       main = sprintf("Konvergenz E[Z_T]\n(delta=%.2f)", delta_g),
       ylim = range(c(E_vec, mc$CI_low, mc$CI_high)) * c(0.99, 1.01))
  polygon(c(min(N_vec), max(N_vec), max(N_vec), min(N_vec)),
          c(mc$CI_low, mc$CI_low, mc$CI_high, mc$CI_high),
          col = rgb(0.2, 0.8, 0.2, 0.2), border = NA)
  abline(h = mc$E_ZT, lwd = 2, lty = 1)
  legend("bottomright", c("SM", "MC", "95%-KI"),
         col = c("steelblue", "black", "lightgreen"),
         lwd = c(2, 2, 8), pch = c(16, NA, NA), bty = "n")
  
  # Plot 2: Relativer Fehler vs N (log-log)
  plot(N_vec, Err_vec, type = "b", pch = 16, col = "firebrick",
       log = "xy",
       xlab = "Gitterpunkte N (log)", ylab = "Relativer Fehler [%] (log)",
       main = "Relativer Fehler vs. N\n(log-log)")
  abline(h = 0.1, lty = 2, col = "grey50")
  text(min(N_vec), 0.12, "0.1%-Schwelle", adj = 0, col = "grey50", cex = 0.8)
  
  # Quadratische Konvergenz einzeichnen (Referenzlinie)
  if (length(N_vec) >= 2) {
    N_ref <- N_vec[1]
    E_ref <- Err_vec[1]
    N_fit <- seq(min(N_vec), max(N_vec), length.out = 50)
    lines(N_fit, E_ref * (N_fit/N_ref)^(-2), lty = 3, col = "grey30")
    text(max(N_fit)*0.9, E_ref*(max(N_fit)/N_ref)^(-2)*1.5,
         "O(N⁻²)", col = "grey30", cex = 0.8)
  }
  
  # Plot 3: Laufzeit vs N
  plot(N_vec, T_vec, type = "b", pch = 16, col = "darkorange",
       xlab = "Gitterpunkte N", ylab = "Laufzeit [Sek.]",
       main = "Laufzeit SM\n(O(N⁴))")
  # O(N^4) Referenzlinie
  if (length(N_vec) >= 2) {
    N_ref <- N_vec[1]; T_ref <- T_vec[1]
    N_fit <- seq(min(N_vec), max(N_vec), length.out = 50)
    lines(N_fit, T_ref * (N_fit/N_ref)^4, lty = 3, col = "grey30")
    text(max(N_fit)*0.85, T_ref*(max(N_fit)/N_ref)^4*1.3,
         "O(N⁴)", col = "grey30", cex = 0.8)
  }
  # Hinweis: die MC-Referenzlaufzeit wird bewusst nicht als horizontale
  # Linie eingezeichnet -- sie liegt bei der vektorisierten 4D-Simulation
  # (siehe mc_eia oben) deutlich unterhalb der Achsenskala aller SM-Laufzeiten
  # und wird stattdessen im Textoutput von vergleiche_alle() ausgegeben.
  
  par(mfrow = c(1, 1))
  cat("Plots angezeigt.\n")
}