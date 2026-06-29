# # simulation.R
# # Pfadsimulation: x(t), y(t), r(t), I_N_oK, I_S_oK, (I_N, I_S für Fonds)
# 
# simuliere_pfade <- function(params, vb, art = "Deckungsstock") {
#   
#   # Parameter auspacken
#   a_x      <- params$a_x
#   a_y      <- params$a_y
#   sigma_x  <- params$sigma_x
#   sigma_y  <- params$sigma_y
#   rho_xy   <- params$rho_xy
#   sigma_I_N <- params$sigma_I_N
#   sigma_I_S <- params$sigma_I_S
#   sigma_s  <- params$sigma_s
#   lambda_N <- params$lambda_N
#   lambda_S <- params$lambda_S
#   l_x      <- params$l_x
#   l_y      <- params$l_y
#   d_x      <- params$d_x
#   d_y      <- params$d_y
#   tau      <- params$tau
#   K_f      <- ifelse(art == "Fonds", params$K_f, 0)
#   delta    <- params$delta
#   n_pfade  <- params$n_pfade
#   
#   F_2_oK  <- params$F_2_oK
#   F_1_oK  <- params$F_1_oK
#   F_0_oK  <- params$F_0_oK
#   F_2     <- params$F_2
#   F_1     <- params$F_1
#   F_0     <- params$F_0
#   
#   psi         <- vb$psi
#   V_t_T       <- vb$V_t_T
#   int_psi_t_T <- vb$int_psi_t_T
#   
#   # Output vorbereiten
#   if (art == "Deckungsstock") {
#     output <- matrix("", nrow = n_pfade * 2, ncol = 1)
#   } else {
#     output <- matrix("", nrow = n_pfade * 4, ncol = 1)
#   }
#   
#   for (sznr in 1:n_pfade) {
#     
#     if (sznr %% 500 == 0) cat("Simuliere Pfad", sznr, "von", n_pfade, "\n")
#     
#     set.seed(sznr)
#     
#     # Initialisierung
#     x <- numeric(481)
#     y <- numeric(481)
#     r <- numeric(481)
#     
#     I_N_oK <- numeric(505)
#     I_S_oK <- numeric(505)
#     I_N    <- numeric(505)
#     I_S    <- numeric(505)
#     
#     I_N_oK[1] <- F_2_oK;  I_N_oK[13] <- F_1_oK;  I_N_oK[25] <- F_0_oK
#     I_S_oK[1] <- F_2_oK;  I_S_oK[13] <- F_1_oK;  I_S_oK[25] <- F_0_oK
#     I_N[1]    <- F_2;     I_N[13]    <- F_1;     I_N[25]    <- F_0
#     I_S[1]    <- F_2;     I_S[13]    <- F_1;     I_S[25]    <- F_0
#     
#     # Output-Zeilen initialisieren
#     if (art == "Deckungsstock") {
#       output[(sznr - 1) * 2 + 1, 1] <- paste0(sznr, ";R(t)_ohne_K_f_Normal;1")
#       output[(sznr - 1) * 2 + 2, 1] <- paste0(sznr, ";R(t)_ohne_K_f_Stress;1")
#     } else {
#       Fonds_initial <- params$Fonds_initial
#       output[(sznr - 1) * 4 + 1, 1] <- paste0(sznr, ";R(t)_ohne_K_f_Normal;", Fonds_initial)
#       output[(sznr - 1) * 4 + 2, 1] <- paste0(sznr, ";R(t)_Normal;",          Fonds_initial)
#       output[(sznr - 1) * 4 + 3, 1] <- paste0(sznr, ";R(t)_ohne_K_f_Stress;", Fonds_initial)
#       output[(sznr - 1) * 4 + 4, 1] <- paste0(sznr, ";R(t)_Stress;",          Fonds_initial)
#       
#       F_N_oK <- Fonds_initial; F_S_oK <- Fonds_initial
#       F_N    <- Fonds_initial; F_S    <- Fonds_initial
#       R_N    <- 1; R_N_oK <- 1
#       R_S    <- 1; R_S_oK <- 1
#       P_akt  <- 1
#     }
#     
#     # Pfadschleife
#     for (t in 1:480) {
#       
#       # Zeitabhängige Risikoprämie
#       if (t <= tau) {
#         lambda_x <- d_x
#         lambda_y <- d_y
#       } else {
#         lambda_x <- l_x
#         lambda_y <- l_y
#       }
#       
#       # Zufallszahlen
#       U_x <- runif(1)
#       U_y <- runif(1)
#       U_f <- runif(1)
#       
#       Z_x <- qnorm(U_x)
#       Z_y <- qnorm(U_y)
#       Z_f <- qnorm(U_f)
#       
#       # Vasicek x(t) und y(t)
#       x[t + 1] <- x[t] * exp(-a_x * delta) +
#         Z_x * sigma_x * sqrt((1 - exp(-2 * a_x * delta)) / (2 * a_x))
#       
#       y[t + 1] <- y[t] * exp(-a_y * delta) +
#         Z_x * rho_xy   * sigma_y * sqrt((1 - exp(-2 * a_y * delta)) / (2 * a_y)) +
#         Z_y * sqrt(1 - rho_xy^2) * sigma_y * sqrt((1 - exp(-2 * a_y * delta)) / (2 * a_y))
#       
#       # Short Rate
#       r[t + 1] <- x[t + 1] + lambda_x * (1 - exp(-a_x * t / 12)) +
#         y[t + 1] + lambda_y * (1 - exp(-a_y * t / 12)) + psi[t + 1]
#       
#       # Aktienindex Normal ohne Kosten
#       I_N_oK[t + 25] <- I_N_oK[t + 24] * exp(
#         delta * (r[t + 1] + r[t]) / 2 +
#           (lambda_N * sigma_I_N / sigma_s - 0.5 * sigma_I_N^2) * delta +
#           Z_f * sigma_I_N * sqrt(delta)
#       )
#       
#       # Aktienindex Stress ohne Kosten
#       I_S_oK[t + 25] <- I_S_oK[t + 24] * exp(
#         (lambda_S * sigma_I_S / sigma_s - 0.5 * sigma_I_S^2) * delta +
#           Z_f * sigma_I_S * sqrt(delta)
#       )
#       
#       if (art == "Fonds") {
#         # Normal mit Kosten
#         #I_N[t + 25] <- I_N[t + 24] * exp(
#        #   delta * (r[t + 1] + r[t]) / 2 +
#       #      (lambda_N * sigma_I_N / sigma_s - 0.5 * sigma_I_N^2 - K_f) * delta +
#      #       Z_f * sigma_I_N * sqrt(delta)
#     #    )
#         
#         # Stress mit Kosten
#    #     I_S[t + 25] <- I_S[t + 24] * exp(
#   #        (lambda_S * sigma_I_S / sigma_s - 0.5 * sigma_I_S^2 - K_f) * delta +
#  #           Z_f * sigma_I_S * sqrt(delta)
# #        )
#    #   }
#       
#   #  } # Ende t-Schleife
#     
#     # Ergebnisse speichern — wird in renditen.R weiterverarbeitet
#     #if (sznr == 1) {
#     #  alle_pfade <- list(
#    #     x      = matrix(0, nrow = n_pfade, ncol = 481),
#   #      y      = matrix(0, nrow = n_pfade, ncol = 481),
#  #       r      = matrix(0, nrow = n_pfade, ncol = 481),
# #        I_N_oK = matrix(0, nrow = n_pfade, ncol = 505),
#         #I_S_oK = matrix(0, nrow = n_pfade, ncol = 505),
#        # I_N    = matrix(0, nrow = n_pfade, ncol = 505),
#       #  I_S    = matrix(0, nrow = n_pfade, ncol = 505)
#      # )
#     #}
#     
#    # alle_pfade$x[sznr, ]      <- x
#   #  alle_pfade$y[sznr, ]      <- y
#  #   alle_pfade$r[sznr, ]      <- r
# #    alle_pfade$I_N_oK[sznr, ] <- I_N_oK
#     #alle_pfade$I_S_oK[sznr, ] <- I_S_oK
#     #alle_pfade$I_N[sznr, ]    <- I_N
#    # alle_pfade$I_S[sznr, ]    <- I_S
#     
#   #} # Ende sznr-Schleife
#   
#  # return(alle_pfade)
# #}
# 

# simulation.R
# Pfadsimulation vektorisiert: alle 10.000 Pfade gleichzeitig

simuliere_pfade <- function(params, vb, art = "Deckungsstock") {
  
  a_x      <- params$a_x
  a_y      <- params$a_y
  sigma_x  <- params$sigma_x
  sigma_y  <- params$sigma_y
  rho_xy   <- params$rho_xy
  sigma_I_N <- params$sigma_I_N
  sigma_I_S <- params$sigma_I_S
  sigma_s  <- params$sigma_s
  lambda_N <- params$lambda_N
  lambda_S <- params$lambda_S
  l_x      <- params$l_x
  l_y      <- params$l_y
  d_x      <- params$d_x
  d_y      <- params$d_y
  tau      <- params$tau
  K_f      <- ifelse(art == "Fonds", params$K_f, 0)
  delta    <- params$delta
  n_pfade  <- params$n_pfade
  
  psi <- vb$psi
  
  # Zufallszahlen für alle Pfade auf einmal ziehen
  # Jede Matrix: n_pfade x 480
  set.seed(42)
  Z_x <- matrix(qnorm(matrix(runif(n_pfade * 480), n_pfade, 480)), n_pfade, 480)
  Z_y <- matrix(qnorm(matrix(runif(n_pfade * 480), n_pfade, 480)), n_pfade, 480)
  Z_f <- matrix(qnorm(matrix(runif(n_pfade * 480), n_pfade, 480)), n_pfade, 480)
  
  # Vorberechnete Konstanten
  exp_ax <- exp(-a_x * delta)
  exp_ay <- exp(-a_y * delta)
  sx <- sigma_x * sqrt((1 - exp(-2 * a_x * delta)) / (2 * a_x))
  sy <- sigma_y * sqrt((1 - exp(-2 * a_y * delta)) / (2 * a_y))
  
  # Matrizen für x, y, r (n_pfade x 481, Spalte 1 = t=0)
  x <- matrix(0, n_pfade, 481)
  y <- matrix(0, n_pfade, 481)
  r <- matrix(0, n_pfade, 481)
  
  # Zeitabhängige Risikoprämien (Vektor der Länge 480)
  lambda_x_t <- ifelse(1:480 <= tau, d_x, l_x)
  lambda_y_t <- ifelse(1:480 <= tau, d_y, l_y)
  
  # Vasicek Schleife — nur noch über Zeit (480 Schritte), nicht über Pfade
  for (t in 1:480) {
    x[, t + 1] <- x[, t] * exp_ax + Z_x[, t] * sx
    
    y[, t + 1] <- y[, t] * exp_ay +
      Z_x[, t] * rho_xy   * sy +
      Z_y[, t] * sqrt(1 - rho_xy^2) * sy
    
    r[, t + 1] <- x[, t + 1] +
      lambda_x_t[t] * (1 - exp(-a_x * t / 12)) +
      y[, t + 1] +
      lambda_y_t[t] * (1 - exp(-a_y * t / 12)) +
      psi[t + 1]
  }
  
  # Aktienindex (n_pfade x 505)
  I_N_oK <- matrix(0, n_pfade, 505)
  I_S_oK <- matrix(0, n_pfade, 505)
  I_N    <- matrix(0, n_pfade, 505)
  I_S    <- matrix(0, n_pfade, 505)
  
  # Startwerte (Vorperioden)
  I_N_oK[, 1]  <- params$F_2_oK;  I_N_oK[, 13] <- params$F_1_oK;  I_N_oK[, 25] <- params$F_0_oK
  I_S_oK[, 1]  <- params$F_2_oK;  I_S_oK[, 13] <- params$F_1_oK;  I_S_oK[, 25] <- params$F_0_oK
  I_N[, 1]     <- params$F_2;     I_N[, 13]    <- params$F_1;     I_N[, 25]    <- params$F_0
  I_S[, 1]     <- params$F_2;     I_S[, 13]    <- params$F_1;     I_S[, 25]    <- params$F_0
  
  # GBM Schleife — nur über Zeit
  drift_N    <- (lambda_N * sigma_I_N / sigma_s - 0.5 * sigma_I_N^2) * delta
  drift_S    <- (lambda_S * sigma_I_S / sigma_s - 0.5 * sigma_I_S^2) * delta
  drift_N_Kf <- drift_N - K_f * delta
  drift_S_Kf <- drift_S - K_f * delta
  vol_N      <- sigma_I_N * sqrt(delta)
  vol_S      <- sigma_I_S * sqrt(delta)
  
  for (t in 1:480) {
    zins_N <- delta * (r[, t + 1] + r[, t]) / 2
    
    I_N_oK[, t + 25] <- I_N_oK[, t + 24] * exp(zins_N + drift_N    + Z_f[, t] * vol_N)
    I_S_oK[, t + 25] <- I_S_oK[, t + 24] * exp(         drift_S    + Z_f[, t] * vol_S)
    
    if (art == "Fonds") {
      I_N[, t + 25] <- I_N[, t + 24] * exp(zins_N + drift_N_Kf + Z_f[, t] * vol_N)
      I_S[, t + 25] <- I_S[, t + 24] * exp(         drift_S_Kf + Z_f[, t] * vol_S)
    }
  }
  
  return(list(
    x      = x,
    y      = y,
    r      = r,
    I_N_oK = I_N_oK,
    I_S_oK = I_S_oK,
    I_N    = I_N,
    I_S    = I_S
  ))
}