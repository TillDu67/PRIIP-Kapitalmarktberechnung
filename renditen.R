# # renditen.R
# # K_swap, R_B_d und Mischrendite vektorisiert
# 
# berechne_renditen <- function(params, vb, alle_pfade, art = "Deckungsstock") {
#   
#   a_x      <- params$a_x
#   a_y      <- params$a_y
#   delta    <- params$delta
#   l_x      <- params$l_x
#   l_y      <- params$l_y
#   d_x      <- params$d_x
#   d_y      <- params$d_y
#   tau      <- params$tau
#   K_f      <- ifelse(art == "Fonds", params$K_f, 0)
#   Ant_Akt  <- params$Ant_Akt
#   Duration <- params$Duration
#   n_pfade  <- params$n_pfade
#   Jahr     <- params$Jahr
#   
#   int_psi_t_T <- vb$int_psi_t_T
#   V_t_T       <- vb$V_t_T
#   
#   x <- alle_pfade$x
#   y <- alle_pfade$y
#   
#   I_N_oK <- alle_pfade$I_N_oK
#   I_S_oK <- alle_pfade$I_S_oK
#   I_N    <- alle_pfade$I_N
#   I_S    <- alle_pfade$I_S
#   
#   n_swap <- 12 * 2 * Duration  # = 144
#   
#   # --- Historische Swaprates einlesen ---
#   swap_hist <- readxl::read_excel(
#     params$pfad_input,
#     sheet = "HistorischeSwaprates",
#     col_names = FALSE
#   )
#   zeile_start <- (Jahr - 2012) * 12 + 3
#   K_swap_hist <- as.numeric(swap_hist[zeile_start:(zeile_start + 479), 11][[1]])
#   
#   # K_swap: n_pfade x 961
#   K_swap <- matrix(0, n_pfade, 961)
#   K_swap[, 1:480] <- matrix(rep(K_swap_hist, each = n_pfade), n_pfade, 480)
#   
#   # K_swap(480) deterministisch
#   TT_vec <- 1:(n_swap + 1)
#   P_vec  <- exp(-int_psi_t_T[1, TT_vec + 1] + 0.5 * V_t_T[1, TT_vec + 1])
#   K_swap[, 481] <- (1 - P_vec[n_swap + 1]) / sum(P_vec)
#   
#   # --- R_B_d Initialisierung (t=0..24) ---
#   R_B_d <- matrix(0, n_pfade, 505)
#   sum2 <- sum(K_swap_hist[(480 - n_swap):(479)])
#   R_B_d[, 1] <- sum2 / n_swap
#   for (t in 1:24) {
#     sum2 <- sum2 - K_swap_hist[480 + t - n_swap] + K_swap[1, 480 + t]
#     R_B_d[, t + 1] <- sum2 / n_swap
#   }
#   
#   # Zeitabhängige Risikoprämien
#   lambda_x_t <- ifelse(1:480 <= tau, d_x, l_x)
#   lambda_y_t <- ifelse(1:480 <= tau, d_y, l_y)
#   
#   # sum2 pro Pfad initialisieren
#   sum2_vec <- rep(n_swap * R_B_d[1, 25], n_pfade)
#   
#   # Output-Matrizen
#   if (art == "Deckungsstock") {
#     R_N_mat <- matrix(0, n_pfade, 480)
#     R_S_mat <- matrix(0, n_pfade, 480)
#   } else {
#     Fonds_initial <- params$Fonds_initial
#     F_N_oK_mat <- matrix(0, n_pfade, 480)
#     F_N_mat    <- matrix(0, n_pfade, 480)
#     F_S_oK_mat <- matrix(0, n_pfade, 480)
#     F_S_mat    <- matrix(0, n_pfade, 480)
#     F_N_oK <- rep(Fonds_initial, n_pfade)
#     F_N    <- rep(Fonds_initial, n_pfade)
#     F_S_oK <- rep(Fonds_initial, n_pfade)
#     F_S    <- rep(Fonds_initial, n_pfade)
#     R_N    <- rep(1, n_pfade)
#     R_N_oK <- rep(1, n_pfade)
#     R_S    <- rep(1, n_pfade)
#     R_S_oK <- rep(1, n_pfade)
#   }
#   
#   # --- Hauptschleife ---
#   cat("Starte Hauptschleife...\n")
#   start_loop <- Sys.time()
#   
#   for (t in 1:480) {
#     
#     lx <- lambda_x_t[t]
#     ly <- lambda_y_t[t]
#     xt <- x[, t + 1]
#     yt <- y[, t + 1]
#     
#     # K_swap vektorisiert
#     TT_vec <- (t + 1):(t + n_swap + 1)
#     h_vec  <- TT_vec - t
#     Bx <- (1 - exp(-a_x * h_vec * delta)) / a_x
#     By <- (1 - exp(-a_y * h_vec * delta)) / a_y
#     int_psi_vec <- int_psi_t_T[t + 1, TT_vec + 1]
#     V_vec       <- 0.5 * V_t_T[t + 1, TT_vec + 1]
#     x_term <- xt + lx * (1 - exp(-a_x * t * delta))
#     y_term <- yt + ly * (1 - exp(-a_y * t * delta))
#     log_P <- outer(-x_term, Bx) + outer(-y_term, By)
#     log_P <- log_P + matrix(rep(-int_psi_vec + V_vec, each = n_pfade), n_pfade, n_swap + 1)
#     P_mat  <- exp(log_P)
#     P_last <- P_mat[, n_swap + 1]
#     sum_P  <- rowSums(P_mat)
#     K_swap[, 480 + t + 1] <- (1 - P_last) / sum_P
#     
#     # R_B_d
#     sum2_vec <- sum2_vec - K_swap[, 480 + t - n_swap + 1] + K_swap[, 480 + t + 1]
#     R_B_d[, 24 + t + 1] <- sum2_vec / n_swap
#     
#     # Mischrendite Deckungsstock (jährlich)
#     if (art == "Deckungsstock" && t %% 12 == 0) {
#       R_N_oK_t <- (
#         (Ant_Akt * I_N_oK[, t+25] / I_N_oK[, t+13] + (1-Ant_Akt) * (1 + R_B_d[, t+25])) *
#           (Ant_Akt * I_N_oK[, t+13] / I_N_oK[, t+ 1] + (1-Ant_Akt) * (1 + R_B_d[, t+13])) *
#           (Ant_Akt * I_N_oK[, t+ 1] / I_N_oK[, t-11] + (1-Ant_Akt) * (1 + R_B_d[, t+ 1]))
#       )^(1/3) - 1
#       R_S_oK_t <- (
#         (Ant_Akt * I_S_oK[, t+25] / I_S_oK[, t+13] + (1-Ant_Akt) * (1 + R_B_d[, t+25])) *
#           (Ant_Akt * I_S_oK[, t+13] / I_S_oK[, t+ 1] + (1-Ant_Akt) * (1 + R_B_d[, t+13])) *
#           (Ant_Akt * I_S_oK[, t+ 1] / I_S_oK[, t-11] + (1-Ant_Akt) * (1 + R_B_d[, t+ 1]))
#       )^(1/3) - 1
#       for (h in 1:12) {
#         R_N_mat[, t - 12 + h] <- R_N_oK_t
#         R_S_mat[, t - 12 + h] <- R_S_oK_t
#       }
#     }
#     
#     # Fondsentwicklung (monatlich)
#     if (art == "Fonds") {
#       P_Temp <- exp(
#         -int_psi_t_T[t,     t + 12*Duration + 1]
#         - (1-exp(-a_x * 12*Duration * delta)) / a_x * (x[, t]   + lx * (1-exp(-a_x*(t-1)*delta)))
#         - (1-exp(-a_y * 12*Duration * delta)) / a_y * (y[, t]   + ly * (1-exp(-a_y*(t-1)*delta)))
#         + 0.5 * V_t_T[t, t + 12*Duration + 1]
#       )
#       P_akt <- exp(
#         -int_psi_t_T[t+1,   t + 12*Duration + 1]
#         - (1-exp(-a_x * 12*Duration * delta)) / a_x * (x[, t+1] + lx * (1-exp(-a_x*t*delta)))
#         - (1-exp(-a_y * 12*Duration * delta)) / a_y * (y[, t+1] + ly * (1-exp(-a_y*t*delta)))
#         + 0.5 * V_t_T[t+1, t + 12*Duration + 1]
#       )
#       R_N_Temp <- R_N;   R_N   <- R_N   * P_akt/P_Temp * exp(-K_f*delta)
#       R_N_oK_T <- R_N_oK; R_N_oK <- R_N_oK * P_akt/P_Temp
#       R_S_Temp <- R_S;   R_S   <- R_S   * P_akt/P_Temp * exp(-K_f*delta)
#       R_S_oK_T <- R_S_oK; R_S_oK <- R_S_oK * P_akt/P_Temp
#       F_N    <- Ant_Akt*F_N    * I_N[,   t+25]/I_N[,   t+24] + (1-Ant_Akt)*F_N    * R_N   /R_N_Temp
#       F_N_oK <- Ant_Akt*F_N_oK * I_N_oK[,t+25]/I_N_oK[,t+24] + (1-Ant_Akt)*F_N_oK * R_N_oK/R_N_oK_T
#       F_S    <- Ant_Akt*F_S    * I_S[,   t+25]/I_S[,   t+24] + (1-Ant_Akt)*F_S    * R_S   /R_S_Temp
#       F_S_oK <- Ant_Akt*F_S_oK * I_S_oK[,t+25]/I_S_oK[,t+24] + (1-Ant_Akt)*F_S_oK * R_S_oK/R_S_oK_T
#       F_N_mat[,    t] <- F_N
#       F_N_oK_mat[, t] <- F_N_oK
#       F_S_mat[,    t] <- F_S
#       F_S_oK_mat[, t] <- F_S_oK
#     }
#     
#   } # Ende t
#   
#   cat("Hauptschleife:", round(difftime(Sys.time(), start_loop, units="secs"), 1), "Sekunden\n")
#   
#   # --- Output zusammenbauen ---
#   cat("Baue Output zusammen...\n")
#   start_output <- Sys.time()
#   
#   if (art == "Deckungsstock") {
#     output <- vector("list", n_pfade * 2)
#     for (sznr in 1:n_pfade) {
#       output[[(sznr-1)*2+1]] <- paste0(sznr, ";R(t)_ohne_K_f_Normal;1;",
#                                        paste(sprintf("%.6f", R_N_mat[sznr,]), collapse=";"))
#       output[[(sznr-1)*2+2]] <- paste0(sznr, ";R(t)_ohne_K_f_Stress;1;",
#                                        paste(sprintf("%.6f", R_S_mat[sznr,]), collapse=";"))
#     }
#   } else {
#     output <- vector("list", n_pfade * 4)
#     for (sznr in 1:n_pfade) {
#       output[[(sznr-1)*4+1]] <- paste0(sznr, ";R(t)_ohne_K_f_Normal;", Fonds_initial, ";",
#                                        paste(sprintf("%#010.6f", F_N_oK_mat[sznr,]), collapse=";"))
#       output[[(sznr-1)*4+2]] <- paste0(sznr, ";R(t)_Normal;",          Fonds_initial, ";",
#                                        paste(sprintf("%#010.6f", F_N_mat[sznr,]),    collapse=";"))
#       output[[(sznr-1)*4+3]] <- paste0(sznr, ";R(t)_ohne_K_f_Stress;", Fonds_initial, ";",
#                                        paste(sprintf("%#010.6f", F_S_oK_mat[sznr,]), collapse=";"))
#       output[[(sznr-1)*4+4]] <- paste0(sznr, ";R(t)_Stress;",          Fonds_initial, ";",
#                                        paste(sprintf("%#010.6f", F_S_mat[sznr,]),    collapse=";"))
#     }
#   }
#   
#   cat("Output-Zusammenbau:", round(difftime(Sys.time(), start_output, units="secs"), 1), "Sekunden\n")
#   
#   return(output)
# }

# renditen.R
# K_swap, R_B_d und Mischrendite vektorisiert (ohne outer())

# berechne_renditen <- function(params, vb, alle_pfade, art = "Deckungsstock") {
#   
#   a_x      <- params$a_x
#   a_y      <- params$a_y
#   delta    <- params$delta
#   l_x      <- params$l_x
#   l_y      <- params$l_y
#   d_x      <- params$d_x
#   d_y      <- params$d_y
#   tau      <- params$tau
#   K_f      <- ifelse(art == "Fonds", params$K_f, 0)
#   Ant_Akt  <- params$Ant_Akt
#   Duration <- params$Duration
#   n_pfade  <- params$n_pfade
#   Jahr     <- params$Jahr
#   
#   int_psi_t_T <- vb$int_psi_t_T
#   V_t_T       <- vb$V_t_T
#   
#   x <- alle_pfade$x
#   y <- alle_pfade$y
#   
#   I_N_oK <- alle_pfade$I_N_oK
#   I_S_oK <- alle_pfade$I_S_oK
#   I_N    <- alle_pfade$I_N
#   I_S    <- alle_pfade$I_S
#   
#   n_swap <- 12 * 2 * Duration  # = 144
#   n_h    <- n_swap + 1         # = 145
#   
#   # --- Historische Swaprates einlesen ---
#   swap_hist <- readxl::read_excel(
#     params$pfad_input,
#     sheet = "HistorischeSwaprates",
#     col_names = FALSE
#   )
#   zeile_start <- (Jahr - 2012) * 12 + 3
#   K_swap_hist <- as.numeric(swap_hist[zeile_start:(zeile_start + 479), 11][[1]])
#   
#   # K_swap: n_pfade x 961
#   K_swap <- matrix(0, n_pfade, 961)
#   K_swap[, 1:480] <- matrix(rep(K_swap_hist, each = n_pfade), n_pfade, 480)
#   
#   # K_swap(480) deterministisch
#   TT_vec0 <- 1:n_h
#   P_vec0  <- exp(-int_psi_t_T[1, TT_vec0 + 1] + 0.5 * V_t_T[1, TT_vec0 + 1])
#   K_swap[, 481] <- (1 - P_vec0[n_h]) / sum(P_vec0)
#   
#   # --- R_B_d Initialisierung (t=0..24) ---
#   R_B_d <- matrix(0, n_pfade, 505)
#   sum2 <- sum(K_swap_hist[(480 - n_swap):(479)])
#   R_B_d[, 1] <- sum2 / n_swap
#   for (t in 1:24) {
#     sum2 <- sum2 - K_swap_hist[480 + t - n_swap] + K_swap[1, 480 + t]
#     R_B_d[, t + 1] <- sum2 / n_swap
#   }
#   
#   # Zeitabhängige Risikoprämien
#   lambda_x_t <- ifelse(1:480 <= tau, d_x, l_x)
#   lambda_y_t <- ifelse(1:480 <= tau, d_y, l_y)
#   
#   # sum2 pro Pfad initialisieren
#   sum2_vec <- rep(n_swap * R_B_d[1, 25], n_pfade)
#   
#   # Output-Matrizen
#   if (art == "Deckungsstock") {
#     R_N_mat <- matrix(0, n_pfade, 480)
#     R_S_mat <- matrix(0, n_pfade, 480)
#   } else {
#     Fonds_initial <- params$Fonds_initial
#     F_N_oK_mat <- matrix(0, n_pfade, 480)
#     F_N_mat    <- matrix(0, n_pfade, 480)
#     F_S_oK_mat <- matrix(0, n_pfade, 480)
#     F_S_mat    <- matrix(0, n_pfade, 480)
#     F_N_oK <- rep(Fonds_initial, n_pfade)
#     F_N    <- rep(Fonds_initial, n_pfade)
#     F_S_oK <- rep(Fonds_initial, n_pfade)
#     F_S    <- rep(Fonds_initial, n_pfade)
#     R_N    <- rep(1, n_pfade)
#     R_N_oK <- rep(1, n_pfade)
#     R_S    <- rep(1, n_pfade)
#     R_S_oK <- rep(1, n_pfade)
#   }
#   
#   # --- Hauptschleife ---
#   cat("Starte Hauptschleife...\n")
#   start_loop <- Sys.time()
#   
#   for (t in 1:480) {
#     
#     lx <- lambda_x_t[t]
#     ly <- lambda_y_t[t]
#     xt <- x[, t + 1]
#     yt <- y[, t + 1]
#     
#     # K_swap vektorisiert OHNE outer():
#     # P(pfad, h) = exp(-int_psi[h] + V[h]) * exp(-Bx[h]*x_term[pfad]) * exp(-By[h]*y_term[pfad])
#     # sum_P = sum_h C[h] * exp(-Bx[h]*x_term) * exp(-By[h]*y_term)
#     # -> für jedes h einzeln akkumulieren (Vektor der Länge n_pfade), keine n_pfade x n_h Matrix
#     TT_vec <- (t + 1):(t + n_h)
#     h_vec  <- TT_vec - t
#     Bx <- (1 - exp(-a_x * h_vec * delta)) / a_x
#     By <- (1 - exp(-a_y * h_vec * delta)) / a_y
#     int_psi_vec <- int_psi_t_T[t + 1, TT_vec + 1]
#     V_vec       <- 0.5 * V_t_T[t + 1, TT_vec + 1]
#     C <- exp(-int_psi_vec + V_vec)  # Länge n_h, konstant über Pfade
#     
#     x_term <- xt + lx * (1 - exp(-a_x * t * delta))  # Länge n_pfade
#     y_term <- yt + ly * (1 - exp(-a_y * t * delta))
#     
#     sum_P <- numeric(n_pfade)
#     for (h in 1:n_h) {
#       Ph <- C[h] * exp(-Bx[h] * x_term - By[h] * y_term)
#       sum_P <- sum_P + Ph
#       if (h == n_h) P_last <- Ph
#     }
#     
#     K_swap[, 480 + t + 1] <- (1 - P_last) / sum_P
#     
#     # R_B_d
#     sum2_vec <- sum2_vec - K_swap[, 480 + t - n_swap + 1] + K_swap[, 480 + t + 1]
#     R_B_d[, 24 + t + 1] <- sum2_vec / n_swap
#     
#     # Mischrendite Deckungsstock (jährlich)
#     if (art == "Deckungsstock" && t %% 12 == 0) {
#       R_N_oK_t <- (
#         (Ant_Akt * I_N_oK[, t+25] / I_N_oK[, t+13] + (1-Ant_Akt) * (1 + R_B_d[, t+25])) *
#           (Ant_Akt * I_N_oK[, t+13] / I_N_oK[, t+ 1] + (1-Ant_Akt) * (1 + R_B_d[, t+13])) *
#           (Ant_Akt * I_N_oK[, t+ 1] / I_N_oK[, t-11] + (1-Ant_Akt) * (1 + R_B_d[, t+ 1]))
#       )^(1/3) - 1
#       R_S_oK_t <- (
#         (Ant_Akt * I_S_oK[, t+25] / I_S_oK[, t+13] + (1-Ant_Akt) * (1 + R_B_d[, t+25])) *
#           (Ant_Akt * I_S_oK[, t+13] / I_S_oK[, t+ 1] + (1-Ant_Akt) * (1 + R_B_d[, t+13])) *
#           (Ant_Akt * I_S_oK[, t+ 1] / I_S_oK[, t-11] + (1-Ant_Akt) * (1 + R_B_d[, t+ 1]))
#       )^(1/3) - 1
#       for (h in 1:12) {
#         R_N_mat[, t - 12 + h] <- R_N_oK_t
#         R_S_mat[, t - 12 + h] <- R_S_oK_t
#       }
#     }
#     
#     # Fondsentwicklung (monatlich)
#     if (art == "Fonds") {
#       P_Temp <- exp(
#         -int_psi_t_T[t,     t + 12*Duration + 1]
#         - (1-exp(-a_x * 12*Duration * delta)) / a_x * (x[, t]   + lx * (1-exp(-a_x*(t-1)*delta)))
#         - (1-exp(-a_y * 12*Duration * delta)) / a_y * (y[, t]   + ly * (1-exp(-a_y*(t-1)*delta)))
#         + 0.5 * V_t_T[t, t + 12*Duration + 1]
#       )
#       P_akt <- exp(
#         -int_psi_t_T[t+1,   t + 12*Duration + 1]
#         - (1-exp(-a_x * 12*Duration * delta)) / a_x * (x[, t+1] + lx * (1-exp(-a_x*t*delta)))
#         - (1-exp(-a_y * 12*Duration * delta)) / a_y * (y[, t+1] + ly * (1-exp(-a_y*t*delta)))
#         + 0.5 * V_t_T[t+1, t + 12*Duration + 1]
#       )
#       R_N_Temp <- R_N;   R_N   <- R_N   * P_akt/P_Temp * exp(-K_f*delta)
#       R_N_oK_T <- R_N_oK; R_N_oK <- R_N_oK * P_akt/P_Temp
#       R_S_Temp <- R_S;   R_S   <- R_S   * P_akt/P_Temp * exp(-K_f*delta)
#       R_S_oK_T <- R_S_oK; R_S_oK <- R_S_oK * P_akt/P_Temp
#       F_N    <- Ant_Akt*F_N    * I_N[,   t+25]/I_N[,   t+24] + (1-Ant_Akt)*F_N    * R_N   /R_N_Temp
#       F_N_oK <- Ant_Akt*F_N_oK * I_N_oK[,t+25]/I_N_oK[,t+24] + (1-Ant_Akt)*F_N_oK * R_N_oK/R_N_oK_T
#       F_S    <- Ant_Akt*F_S    * I_S[,   t+25]/I_S[,   t+24] + (1-Ant_Akt)*F_S    * R_S   /R_S_Temp
#       F_S_oK <- Ant_Akt*F_S_oK * I_S_oK[,t+25]/I_S_oK[,t+24] + (1-Ant_Akt)*F_S_oK * R_S_oK/R_S_oK_T
#       F_N_mat[,    t] <- F_N
#       F_N_oK_mat[, t] <- F_N_oK
#       F_S_mat[,    t] <- F_S
#       F_S_oK_mat[, t] <- F_S_oK
#     }
#     
#   } # Ende t
#   
#   cat("Hauptschleife:", round(difftime(Sys.time(), start_loop, units="secs"), 1), "Sekunden\n")
#   
#   # --- Output zusammenbauen (vektorisiert mit apply) ---
#   cat("Baue Output zusammen...\n")
#   start_output <- Sys.time()
#   
#   if (art == "Deckungsstock") {
#     output <- vector("list", n_pfade * 2)
#     fmt_N <- matrix(sprintf("%.6f", R_N_mat), n_pfade, 480)
#     fmt_S <- matrix(sprintf("%.6f", R_S_mat), n_pfade, 480)
#     rows_N <- apply(fmt_N, 1, paste, collapse = ";")
#     rows_S <- apply(fmt_S, 1, paste, collapse = ";")
#     for (sznr in 1:n_pfade) {
#       output[[(sznr-1)*2+1]] <- paste0(sznr, ";R(t)_ohne_K_f_Normal;1;", rows_N[sznr])
#       output[[(sznr-1)*2+2]] <- paste0(sznr, ";R(t)_ohne_K_f_Stress;1;", rows_S[sznr])
#     }
#   } else {
#     output <- vector("list", n_pfade * 4)
#     fmt_N_oK <- matrix(sprintf("%#010.6f", F_N_oK_mat), n_pfade, 480)
#     fmt_N    <- matrix(sprintf("%#010.6f", F_N_mat),    n_pfade, 480)
#     fmt_S_oK <- matrix(sprintf("%#010.6f", F_S_oK_mat), n_pfade, 480)
#     fmt_S    <- matrix(sprintf("%#010.6f", F_S_mat),    n_pfade, 480)
#     rows_N_oK <- apply(fmt_N_oK, 1, paste, collapse = ";")
#     rows_N    <- apply(fmt_N,    1, paste, collapse = ";")
#     rows_S_oK <- apply(fmt_S_oK, 1, paste, collapse = ";")
#     rows_S    <- apply(fmt_S,    1, paste, collapse = ";")
#     for (sznr in 1:n_pfade) {
#       output[[(sznr-1)*4+1]] <- paste0(sznr, ";R(t)_ohne_K_f_Normal;", Fonds_initial, ";", rows_N_oK[sznr])
#       output[[(sznr-1)*4+2]] <- paste0(sznr, ";R(t)_Normal;",          Fonds_initial, ";", rows_N[sznr])
#       output[[(sznr-1)*4+3]] <- paste0(sznr, ";R(t)_ohne_K_f_Stress;", Fonds_initial, ";", rows_S_oK[sznr])
#       output[[(sznr-1)*4+4]] <- paste0(sznr, ";R(t)_Stress;",          Fonds_initial, ";", rows_S[sznr])
#     }
#   }
#   
#   cat("Output-Zusammenbau:", round(difftime(Sys.time(), start_output, units="secs"), 1), "Sekunden\n")
#   
#   return(output)
# }

# renditen.R
# K_swap, R_B_d und Mischrendite - voll vektorisiert (Schleife über h, nicht t x h)

# renditen.R
# K_swap, R_B_d und Mischrendite - voll vektorisiert (Schleife über h, nicht t x h)


####################derzeit bester Stand
# berechne_renditen <- function(params, vb, alle_pfade, art = "Deckungsstock") {
#   
#   a_x      <- params$a_x
#   a_y      <- params$a_y
#   delta    <- params$delta
#   l_x      <- params$l_x
#   l_y      <- params$l_y
#   d_x      <- params$d_x
#   d_y      <- params$d_y
#   tau      <- params$tau
#   K_f      <- ifelse(art == "Fonds", params$K_f, 0)
#   Ant_Akt  <- params$Ant_Akt
#   Duration <- params$Duration
#   n_pfade  <- params$n_pfade
#   Jahr     <- params$Jahr
#   
#   int_psi_t_T <- vb$int_psi_t_T
#   V_t_T       <- vb$V_t_T
#   
#   x <- alle_pfade$x  # n_pfade x 481
#   y <- alle_pfade$y
#   
#   I_N_oK <- alle_pfade$I_N_oK
#   I_S_oK <- alle_pfade$I_S_oK
#   I_N    <- alle_pfade$I_N
#   I_S    <- alle_pfade$I_S
#   
#   n_swap <- 12 * 2 * Duration  # = 144
#   n_h    <- n_swap + 1         # = 145
#   
#   # --- Historische Swaprates einlesen ---
#   swap_hist <- readxl::read_excel(
#     params$pfad_input,
#     sheet = "HistorischeSwaprates",
#     col_names = FALSE
#   )
#   zeile_start <- (Jahr - 2012) * 12 + 3
#   K_swap_hist <- as.numeric(swap_hist[zeile_start:(zeile_start + 479), 11][[1]])
#   
#   # K_swap: n_pfade x 961
#   K_swap <- matrix(0, n_pfade, 961)
#   K_swap[, 1:480] <- matrix(rep(K_swap_hist, each = n_pfade), n_pfade, 480)
#   
#   # K_swap(480) deterministisch (t=0)
#   TT_vec0 <- 1:n_h
#   P_vec0  <- exp(-int_psi_t_T[1, TT_vec0 + 1] + 0.5 * V_t_T[1, TT_vec0 + 1])
#   K_swap[, 481] <- (1 - P_vec0[n_h]) / sum(P_vec0)
#   
#   cat("Vorbereitung K_swap-Simulation (vektorisiert)...\n")
#   start_prep <- Sys.time()
#   
#   # --- Bx, By: konstant über t (h_vec ist immer 1:n_h) ---
#   h_vec <- 1:n_h
#   Bx <- (1 - exp(-a_x * h_vec * delta)) / a_x   # Länge n_h
#   By <- (1 - exp(-a_y * h_vec * delta)) / a_y   # Länge n_h
#   
#   # --- C[t,h] = exp(-int_psi_t_T[t+1,t+1+h] + 0.5*V_t_T[t+1,t+1+h]), t=1..480, h=1..n_h ---
#   C <- matrix(0, 480, n_h)
#   rows <- 2:481  # entspricht t+1 für t=1..480
#   for (h in 1:n_h) {
#     cols <- (1:480) + 1 + h
#     idx  <- cbind(rows, cols)
#     C[, h] <- exp(-int_psi_t_T[idx] + 0.5 * V_t_T[idx])
#   }
#   
#   # --- Zeitabhängige Risikoprämien ---
#   lambda_x_t <- ifelse(1:480 <= tau, d_x, l_x)
#   lambda_y_t <- ifelse(1:480 <= tau, d_y, l_y)
#   
#   # --- X_TERM_T, Y_TERM_T: 480 x n_pfade (transponiert: Zeile=t, Spalte=Pfad) ---
#   X_TERM_T <- t(x[, 2:481]) + (lambda_x_t * (1 - exp(-a_x * (1:480) * delta)))
#   Y_TERM_T <- t(y[, 2:481]) + (lambda_y_t * (1 - exp(-a_y * (1:480) * delta)))
#   
#   cat("Vorbereitung:", round(difftime(Sys.time(), start_prep, units="secs"), 1), "Sekunden\n")
#   
#   # --- Hauptschleife: nur 145 Iterationen über h, jeweils 480 x n_pfade Matrizen ---
#   cat("Starte K_swap-Simulation (Schleife über h)...\n")
#   start_h <- Sys.time()
#   
#   sum_P_T <- matrix(0, 480, n_pfade)
#   P_last_T <- NULL
#   
#   for (h in 1:n_h) {
#     # Ein einziger exp()-Aufruf statt zwei (EX*EY getrennt) - spart Hälfte der Allokationen
#     contrib <- C[, h] * exp(-(Bx[h] * X_TERM_T + By[h] * Y_TERM_T))  # 480 x n_pfade
#     sum_P_T <- sum_P_T + contrib
#     if (h == n_h) P_last_T <- contrib
#   }
#   
#   K_swap_sim_T <- (1 - P_last_T) / sum_P_T   # 480 x n_pfade
#   K_swap[, 482:961] <- t(K_swap_sim_T)        # zurück zu n_pfade x 480
#   
#   cat("K_swap-Simulation:", round(difftime(Sys.time(), start_h, units="secs"), 1), "Sekunden\n")
#   
#   # --- R_B_d Initialisierung (t=0..24) ---
#   R_B_d <- matrix(0, n_pfade, 505)
#   sum2 <- sum(K_swap_hist[(480 - n_swap):(479)])
#   R_B_d[, 1] <- sum2 / n_swap
#   for (t in 1:24) {
#     sum2 <- sum2 - K_swap_hist[480 + t - n_swap] + K_swap[1, 480 + t]
#     R_B_d[, t + 1] <- sum2 / n_swap
#   }
#   
#   sum2_vec <- rep(n_swap * R_B_d[1, 25], n_pfade)
#   
#   # Output-Matrizen
#   if (art == "Deckungsstock") {
#     R_N_mat <- matrix(0, n_pfade, 480)
#     R_S_mat <- matrix(0, n_pfade, 480)
#   } else {
#     Fonds_initial <- params$Fonds_initial
#     F_N_oK_mat <- matrix(0, n_pfade, 480)
#     F_N_mat    <- matrix(0, n_pfade, 480)
#     F_S_oK_mat <- matrix(0, n_pfade, 480)
#     F_S_mat    <- matrix(0, n_pfade, 480)
#     F_N_oK <- rep(Fonds_initial, n_pfade)
#     F_N    <- rep(Fonds_initial, n_pfade)
#     F_S_oK <- rep(Fonds_initial, n_pfade)
#     F_S    <- rep(Fonds_initial, n_pfade)
#     R_N    <- rep(1, n_pfade)
#     R_N_oK <- rep(1, n_pfade)
#     R_S    <- rep(1, n_pfade)
#     R_S_oK <- rep(1, n_pfade)
#   }
#   
#   # --- Restliche Schleife: nur noch O(n_pfade) pro Iteration, 480 Iterationen ---
#   cat("Starte R_B_d / Mischrendite-Schleife...\n")
#   start_loop <- Sys.time()
#   
#   for (t in 1:480) {
#     
#     lx <- lambda_x_t[t]
#     ly <- lambda_y_t[t]
#     
#     # R_B_d gleitender Durchschnitt (sliding window)
#     sum2_vec <- sum2_vec - K_swap[, 480 + t - n_swap + 1] + K_swap[, 480 + t + 1]
#     R_B_d[, 24 + t + 1] <- sum2_vec / n_swap
#     
#     # Mischrendite Deckungsstock (jährlich)
#     if (art == "Deckungsstock" && t %% 12 == 0) {
#       R_N_oK_t <- (
#         (Ant_Akt * I_N_oK[, t+25] / I_N_oK[, t+13] + (1-Ant_Akt) * (1 + R_B_d[, t+25])) *
#           (Ant_Akt * I_N_oK[, t+13] / I_N_oK[, t+ 1] + (1-Ant_Akt) * (1 + R_B_d[, t+13])) *
#           (Ant_Akt * I_N_oK[, t+ 1] / I_N_oK[, t-11] + (1-Ant_Akt) * (1 + R_B_d[, t+ 1]))
#       )^(1/3) - 1
#       R_S_oK_t <- (
#         (Ant_Akt * I_S_oK[, t+25] / I_S_oK[, t+13] + (1-Ant_Akt) * (1 + R_B_d[, t+25])) *
#           (Ant_Akt * I_S_oK[, t+13] / I_S_oK[, t+ 1] + (1-Ant_Akt) * (1 + R_B_d[, t+13])) *
#           (Ant_Akt * I_S_oK[, t+ 1] / I_S_oK[, t-11] + (1-Ant_Akt) * (1 + R_B_d[, t+ 1]))
#       )^(1/3) - 1
#       for (h in 1:12) {
#         R_N_mat[, t - 12 + h] <- R_N_oK_t
#         R_S_mat[, t - 12 + h] <- R_S_oK_t
#       }
#     }
#     
#     # Fondsentwicklung (monatlich)
#     if (art == "Fonds") {
#       P_Temp <- exp(
#         -int_psi_t_T[t,     t + 12*Duration + 1]
#         - (1-exp(-a_x * 12*Duration * delta)) / a_x * (x[, t]   + lx * (1-exp(-a_x*(t-1)*delta)))
#         - (1-exp(-a_y * 12*Duration * delta)) / a_y * (y[, t]   + ly * (1-exp(-a_y*(t-1)*delta)))
#         + 0.5 * V_t_T[t, t + 12*Duration + 1]
#       )
#       P_akt <- exp(
#         -int_psi_t_T[t+1,   t + 12*Duration + 1]
#         - (1-exp(-a_x * 12*Duration * delta)) / a_x * (x[, t+1] + lx * (1-exp(-a_x*t*delta)))
#         - (1-exp(-a_y * 12*Duration * delta)) / a_y * (y[, t+1] + ly * (1-exp(-a_y*t*delta)))
#         + 0.5 * V_t_T[t+1, t + 12*Duration + 1]
#       )
#       R_N_Temp <- R_N;   R_N   <- R_N   * P_akt/P_Temp * exp(-K_f*delta)
#       R_N_oK_T <- R_N_oK; R_N_oK <- R_N_oK * P_akt/P_Temp
#       R_S_Temp <- R_S;   R_S   <- R_S   * P_akt/P_Temp * exp(-K_f*delta)
#       R_S_oK_T <- R_S_oK; R_S_oK <- R_S_oK * P_akt/P_Temp
#       F_N    <- Ant_Akt*F_N    * I_N[,   t+25]/I_N[,   t+24] + (1-Ant_Akt)*F_N    * R_N   /R_N_Temp
#       F_N_oK <- Ant_Akt*F_N_oK * I_N_oK[,t+25]/I_N_oK[,t+24] + (1-Ant_Akt)*F_N_oK * R_N_oK/R_N_oK_T
#       F_S    <- Ant_Akt*F_S    * I_S[,   t+25]/I_S[,   t+24] + (1-Ant_Akt)*F_S    * R_S   /R_S_Temp
#       F_S_oK <- Ant_Akt*F_S_oK * I_S_oK[,t+25]/I_S_oK[,t+24] + (1-Ant_Akt)*F_S_oK * R_S_oK/R_S_oK_T
#       F_N_mat[,    t] <- F_N
#       F_N_oK_mat[, t] <- F_N_oK
#       F_S_mat[,    t] <- F_S
#       F_S_oK_mat[, t] <- F_S_oK
#     }
#     
#   } # Ende t
#   
#   cat("R_B_d/Mischrendite-Schleife:", round(difftime(Sys.time(), start_loop, units="secs"), 1), "Sekunden\n")
#   
#   # --- Output zusammenbauen ---
#   cat("Baue Output zusammen...\n")
#   start_output <- Sys.time()
#   
#   if (art == "Deckungsstock") {
#     output <- vector("list", n_pfade * 2)
#     fmt_N <- matrix(sprintf("%.6f", R_N_mat), n_pfade, 480)
#     fmt_S <- matrix(sprintf("%.6f", R_S_mat), n_pfade, 480)
#     rows_N <- do.call(paste, c(as.data.frame(fmt_N), sep = ";"))
#     rows_S <- do.call(paste, c(as.data.frame(fmt_S), sep = ";"))
#     for (sznr in 1:n_pfade) {
#       output[[(sznr-1)*2+1]] <- paste0(sznr, ";R(t)_ohne_K_f_Normal;1;", rows_N[sznr])
#       output[[(sznr-1)*2+2]] <- paste0(sznr, ";R(t)_ohne_K_f_Stress;1;", rows_S[sznr])
#     }
#   } else {
#     output <- vector("list", n_pfade * 4)
#     fmt_N_oK <- matrix(sprintf("%#010.6f", F_N_oK_mat), n_pfade, 480)
#     fmt_N    <- matrix(sprintf("%#010.6f", F_N_mat),    n_pfade, 480)
#     fmt_S_oK <- matrix(sprintf("%#010.6f", F_S_oK_mat), n_pfade, 480)
#     fmt_S    <- matrix(sprintf("%#010.6f", F_S_mat),    n_pfade, 480)
#     rows_N_oK <- do.call(paste, c(as.data.frame(fmt_N_oK), sep = ";"))
#     rows_N    <- do.call(paste, c(as.data.frame(fmt_N),    sep = ";"))
#     rows_S_oK <- do.call(paste, c(as.data.frame(fmt_S_oK), sep = ";"))
#     rows_S    <- do.call(paste, c(as.data.frame(fmt_S),    sep = ";"))
#     for (sznr in 1:n_pfade) {
#       output[[(sznr-1)*4+1]] <- paste0(sznr, ";R(t)_ohne_K_f_Normal;", Fonds_initial, ";", rows_N_oK[sznr])
#       output[[(sznr-1)*4+2]] <- paste0(sznr, ";R(t)_Normal;",          Fonds_initial, ";", rows_N[sznr])
#       output[[(sznr-1)*4+3]] <- paste0(sznr, ";R(t)_ohne_K_f_Stress;", Fonds_initial, ";", rows_S_oK[sznr])
#       output[[(sznr-1)*4+4]] <- paste0(sznr, ";R(t)_Stress;",          Fonds_initial, ";", rows_S[sznr])
#     }
#   }
#   
#   cat("Output-Zusammenbau:", round(difftime(Sys.time(), start_output, units="secs"), 1), "Sekunden\n")
#   
#   return(output)
# }

# renditen.R
# K_swap, R_B_d und Mischrendite - voll vektorisiert (Schleife über h, nicht t x h)

berechne_renditen <- function(params, vb, alle_pfade, art = "Deckungsstock") {
  
  a_x      <- params$a_x
  a_y      <- params$a_y
  delta    <- params$delta
  l_x      <- params$l_x
  l_y      <- params$l_y
  d_x      <- params$d_x
  d_y      <- params$d_y
  tau      <- params$tau
  K_f      <- ifelse(art == "Fonds", params$K_f, 0)
  Ant_Akt  <- params$Ant_Akt
  Duration <- params$Duration
  n_pfade  <- params$n_pfade
  Jahr     <- params$Jahr
  
  int_psi_t_T <- vb$int_psi_t_T
  V_t_T       <- vb$V_t_T
  
  x <- alle_pfade$x  # n_pfade x 481
  y <- alle_pfade$y
  
  I_N_oK <- alle_pfade$I_N_oK
  I_S_oK <- alle_pfade$I_S_oK
  I_N    <- alle_pfade$I_N
  I_S    <- alle_pfade$I_S
  
  n_swap <- 12 * 2 * Duration  # = 144
  n_h    <- n_swap + 1         # = 145
  
  # --- Historische Swaprates einlesen ---
  swap_hist <- readxl::read_excel(
    params$pfad_input,
    sheet = "HistorischeSwaprates",
    col_names = FALSE
  )
  zeile_start <- (Jahr - 2012) * 12 + 3
  K_swap_hist <- as.numeric(swap_hist[zeile_start:(zeile_start + 479), 11][[1]])
  
  # K_swap: n_pfade x 961
  K_swap <- matrix(0, n_pfade, 961)
  K_swap[, 1:480] <- matrix(rep(K_swap_hist, each = n_pfade), n_pfade, 480)
  
  # K_swap(480) deterministisch (t=0)
  TT_vec0 <- 1:n_h
  P_vec0  <- exp(-int_psi_t_T[1, TT_vec0 + 1] + 0.5 * V_t_T[1, TT_vec0 + 1])
  K_swap[, 481] <- (1 - P_vec0[n_h]) / sum(P_vec0)
  
  cat("Vorbereitung K_swap-Simulation (vektorisiert)...\n")
  start_prep <- Sys.time()
  
  # Bx, By: konstant über t (h_vec ist immer 1:n_h)
  h_vec <- 1:n_h
  Bx <- (1 - exp(-a_x * h_vec * delta)) / a_x   # Länge n_h
  By <- (1 - exp(-a_y * h_vec * delta)) / a_y   # Länge n_h
  
  # C[t,h] = exp(-int_psi_t_T[t+1,t+1+h] + 0.5*V_t_T[t+1,t+1+h]), t=1..480, h=1..n_h
  C <- matrix(0, 480, n_h)
  rows <- 2:481  # entspricht t+1 für t=1..480
  for (h in 1:n_h) {
    cols <- (1:480) + 1 + h
    idx  <- cbind(rows, cols)
    C[, h] <- exp(-int_psi_t_T[idx] + 0.5 * V_t_T[idx])
  }
  
  # Zeitabhängige Risikoprämien
  lambda_x_t <- ifelse(1:480 <= tau, d_x, l_x)
  lambda_y_t <- ifelse(1:480 <= tau, d_y, l_y)
  
  # X_TERM_T, Y_TERM_T: 480 x n_pfade (transponiert: Zeile=t, Spalte=Pfad)
  X_TERM_T <- t(x[, 2:481]) + (lambda_x_t * (1 - exp(-a_x * (1:480) * delta)))
  Y_TERM_T <- t(y[, 2:481]) + (lambda_y_t * (1 - exp(-a_y * (1:480) * delta)))
  
  cat("Vorbereitung:", round(difftime(Sys.time(), start_prep, units="secs"), 1), "Sekunden\n")
  
  # Hauptschleife: nur 145 Iterationen über h, jeweils 480 x n_pfade Matrizen
  cat("Starte K_swap-Simulation (Schleife über h)...\n")
  start_h <- Sys.time()
  
  sum_P_T <- matrix(0, 480, n_pfade)
  P_last_T <- NULL
  
  for (h in 1:n_h) {
    # Ein einziger exp()-Aufruf statt zwei (EX*EY getrennt) - spart Hälfte der Allokationen
    contrib <- C[, h] * exp(-(Bx[h] * X_TERM_T + By[h] * Y_TERM_T))  # 480 x n_pfade
    sum_P_T <- sum_P_T + contrib
    if (h == n_h) P_last_T <- contrib
  }
  
  K_swap_sim_T <- (1 - P_last_T) / sum_P_T   # 480 x n_pfade
  K_swap[, 482:961] <- t(K_swap_sim_T)        # zurück zu n_pfade x 480
  
  cat("K_swap-Simulation:", round(difftime(Sys.time(), start_h, units="secs"), 1), "Sekunden\n")
  
  # R_B_d Initialisierung (s=0..24)
  # VBA: R_B_d(s) = mean(K_swap(456-n_swap+s .. 455+s)), allgemein für beliebiges n_swap
  R_B_d <- matrix(0, n_pfade, 505)
  sum2 <- sum(K_swap_hist[(457 - n_swap):456])
  R_B_d[, 1] <- sum2 / n_swap
  for (t in 1:24) {
    sum2 <- sum2 - K_swap_hist[456 - n_swap + t] + K_swap_hist[456 + t]
    R_B_d[, t + 1] <- sum2 / n_swap
  }
  
  sum2_vec <- rep(n_swap * R_B_d[1, 25], n_pfade)
  
  # Output-Matrizen
  if (art == "Deckungsstock") {
    R_N_mat <- matrix(0, n_pfade, 480)
    R_S_mat <- matrix(0, n_pfade, 480)
  } else {
    Fonds_initial <- params$Initialer_Fondskurs
    F_N_oK_mat <- matrix(0, n_pfade, 480)
    F_N_mat    <- matrix(0, n_pfade, 480)
    F_S_oK_mat <- matrix(0, n_pfade, 480)
    F_S_mat    <- matrix(0, n_pfade, 480)
    F_N_oK <- rep(Fonds_initial, n_pfade)
    F_N    <- rep(Fonds_initial, n_pfade)
    F_S_oK <- rep(Fonds_initial, n_pfade)
    F_S    <- rep(Fonds_initial, n_pfade)
    R_N    <- rep(1, n_pfade)
    R_N_oK <- rep(1, n_pfade)
    R_S    <- rep(1, n_pfade)
    R_S_oK <- rep(1, n_pfade)
  }
  
  # Restliche Schleife: nur noch O(n_pfade) pro Iteration, 480 Iterationen
  cat("Starte R_B_d / Mischrendite-Schleife...\n")
  start_loop <- Sys.time()
  
  for (t in 1:480) {
    
    lx <- lambda_x_t[t]
    ly <- lambda_y_t[t]
    
    # R_B_d gleitender Durchschnitt (sliding window)
    # VBA: sum_2 -= K_swap(335+t) ; sum_2 += K_swap(479+t)
    sum2_vec <- sum2_vec - K_swap[, 480 - n_swap + t] + K_swap[, 480 + t]
    R_B_d[, 24 + t + 1] <- sum2_vec / n_swap
    
    # Mischrendite Deckungsstock (jährlich)
    if (art == "Deckungsstock" && t %% 12 == 0) {
      R_N_oK_t <- (
        (Ant_Akt * I_N_oK[, t+25] / I_N_oK[, t+13] + (1-Ant_Akt) * (1 + R_B_d[, t+25])) *
          (Ant_Akt * I_N_oK[, t+13] / I_N_oK[, t+ 1] + (1-Ant_Akt) * (1 + R_B_d[, t+13])) *
          (Ant_Akt * I_N_oK[, t+ 1] / I_N_oK[, t-11] + (1-Ant_Akt) * (1 + R_B_d[, t+ 1]))
      )^(1/3) - 1
      R_S_oK_t <- (
        (Ant_Akt * I_S_oK[, t+25] / I_S_oK[, t+13] + (1-Ant_Akt) * (1 + R_B_d[, t+25])) *
          (Ant_Akt * I_S_oK[, t+13] / I_S_oK[, t+ 1] + (1-Ant_Akt) * (1 + R_B_d[, t+13])) *
          (Ant_Akt * I_S_oK[, t+ 1] / I_S_oK[, t-11] + (1-Ant_Akt) * (1 + R_B_d[, t+ 1]))
      )^(1/3) - 1
      for (h in 1:12) {
        R_N_mat[, t - 12 + h] <- R_N_oK_t
        R_S_mat[, t - 12 + h] <- R_S_oK_t
      }
    }
    
    # Fondsentwicklung (monatlich)
    if (art == "Fonds") {
      P_Temp <- exp(
        -int_psi_t_T[t,     t + 12*Duration + 1]
        - (1-exp(-a_x * 12*Duration * delta)) / a_x * (x[, t]   + lx * (1-exp(-a_x*(t-1)*delta)))
        - (1-exp(-a_y * 12*Duration * delta)) / a_y * (y[, t]   + ly * (1-exp(-a_y*(t-1)*delta)))
        + 0.5 * V_t_T[t, t + 12*Duration + 1]
      )
      P_akt <- exp(
        -int_psi_t_T[t+1,   t + 12*Duration + 1]
        - (1-exp(-a_x * 12*Duration * delta)) / a_x * (x[, t+1] + lx * (1-exp(-a_x*t*delta)))
        - (1-exp(-a_y * 12*Duration * delta)) / a_y * (y[, t+1] + ly * (1-exp(-a_y*t*delta)))
        + 0.5 * V_t_T[t+1, t + 12*Duration + 1]
      )
      R_N_Temp <- R_N;   R_N   <- R_N   * P_akt/P_Temp * exp(-K_f*delta)
      R_N_oK_T <- R_N_oK; R_N_oK <- R_N_oK * P_akt/P_Temp
      R_S_Temp <- R_S;   R_S   <- R_S   * P_akt/P_Temp * exp(-K_f*delta)
      R_S_oK_T <- R_S_oK; R_S_oK <- R_S_oK * P_akt/P_Temp
      F_N    <- Ant_Akt*F_N    * I_N[,   t+25]/I_N[,   t+24] + (1-Ant_Akt)*F_N    * R_N   /R_N_Temp
      F_N_oK <- Ant_Akt*F_N_oK * I_N_oK[,t+25]/I_N_oK[,t+24] + (1-Ant_Akt)*F_N_oK * R_N_oK/R_N_oK_T
      F_S    <- Ant_Akt*F_S    * I_S[,   t+25]/I_S[,   t+24] + (1-Ant_Akt)*F_S    * R_S   /R_S_Temp
      F_S_oK <- Ant_Akt*F_S_oK * I_S_oK[,t+25]/I_S_oK[,t+24] + (1-Ant_Akt)*F_S_oK * R_S_oK/R_S_oK_T
      F_N_mat[,    t] <- F_N
      F_N_oK_mat[, t] <- F_N_oK
      F_S_mat[,    t] <- F_S
      F_S_oK_mat[, t] <- F_S_oK
    }
    
  } # Ende t
  
  cat("R_B_d/Mischrendite-Schleife:", round(difftime(Sys.time(), start_loop, units="secs"), 1), "Sekunden\n")
  
  # Output zusammenbauen
  cat("Baue Output zusammen...\n")
  start_output <- Sys.time()
  
  if (art == "Deckungsstock") {
    output <- vector("list", n_pfade * 2)
    fmt_N <- matrix(sprintf("%.6f", R_N_mat), n_pfade, 480)
    fmt_S <- matrix(sprintf("%.6f", R_S_mat), n_pfade, 480)
    rows_N <- do.call(paste, c(as.data.frame(fmt_N), sep = ";"))
    rows_S <- do.call(paste, c(as.data.frame(fmt_S), sep = ";"))
    for (sznr in 1:n_pfade) {
      output[[(sznr-1)*2+1]] <- paste0(sznr, ";R(t)_ohne_K_f_Normal;1;", rows_N[sznr])
      output[[(sznr-1)*2+2]] <- paste0(sznr, ";R(t)_ohne_K_f_Stress;1;", rows_S[sznr])
    }
  } else {
    output <- vector("list", n_pfade * 4)
    fmt_N_oK <- matrix(sprintf("%#010.6f", F_N_oK_mat), n_pfade, 480)
    fmt_N    <- matrix(sprintf("%#010.6f", F_N_mat),    n_pfade, 480)
    fmt_S_oK <- matrix(sprintf("%#010.6f", F_S_oK_mat), n_pfade, 480)
    fmt_S    <- matrix(sprintf("%#010.6f", F_S_mat),    n_pfade, 480)
    rows_N_oK <- do.call(paste, c(as.data.frame(fmt_N_oK), sep = ";"))
    rows_N    <- do.call(paste, c(as.data.frame(fmt_N),    sep = ";"))
    rows_S_oK <- do.call(paste, c(as.data.frame(fmt_S_oK), sep = ";"))
    rows_S    <- do.call(paste, c(as.data.frame(fmt_S),    sep = ";"))
    for (sznr in 1:n_pfade) {
      output[[(sznr-1)*4+1]] <- paste0(sznr, ";R(t)_ohne_K_f_Normal;", Fonds_initial, ";", rows_N_oK[sznr])
      output[[(sznr-1)*4+2]] <- paste0(sznr, ";R(t)_Normal;",          Fonds_initial, ";", rows_N[sznr])
      output[[(sznr-1)*4+3]] <- paste0(sznr, ";R(t)_ohne_K_f_Stress;", Fonds_initial, ";", rows_S_oK[sznr])
      output[[(sznr-1)*4+4]] <- paste0(sznr, ";R(t)_Stress;",          Fonds_initial, ";", rows_S[sznr])
    }
  }
  
  cat("Output-Zusammenbau:", round(difftime(Sys.time(), start_output, units="secs"), 1), "Sekunden\n")
  
  if (art == "Deckungsstock") {
    return(list(output = output, R_N_mat = R_N_mat, R_S_mat = R_S_mat, K_swap = K_swap, R_B_d = R_B_d))
  } else {
    return(list(output = output, F_N_mat = F_N_mat, F_S_mat = F_S_mat,
                F_N_oK_mat = F_N_oK_mat, F_S_oK_mat = F_S_oK_mat, K_swap = K_swap, R_B_d = R_B_d))
  }
}