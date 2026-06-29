# szenario_matrix.R
# Szenario-Matrix-Ansatz nach Günther & Hieber (2024)
# "Efficient simulation and valuation of equity-indexed annuities
#  under a two-factor G2++ model", European Actuarial Journal
#
# Wiederverwendet: vb$psi, vb$int_psi_t_T, vb$V_t_T aus vorberechnungen.R
# Benötigt: install.packages("mvtnorm")
#
# Öffentliche Einstiegsfunktion: sm_berechnung()
# Aufruf aus main.R:
#   source("szenario_matrix.R")
#   sm <- sm_berechnung(params, vb, N = 10, delta_grid = 8.25,
#                       T_jahre = 25, g = 0.005, alpha = 0.465)

# ============================================================
# MODUL 1: Kovarianzmatrix Sigma (Lemma 2.1)
# ============================================================
# Berechnet die 4x4-Kovarianzmatrix der gemeinsamen Verteilung
# (x_t, y_t, integral_r, ln(S_t/S_{t-1})) | F_{t-1}
# nach Lemma 2.1 im Paper.
#
# Parameter entsprechen deiner Namenskonvention:
#   a     = a_x   (mean reversion speed x)
#   b     = a_y   (mean reversion speed y)
#   nu    = sigma_x (Volatilität x)
#   eta   = sigma_y (Volatilität y)
#   rho_r = rho_xy  (Korrelation x-y)
#   sigma = sigma_s (Aktienvolatilität)
#   rho1  = Korrelation S mit x-Brownian
#   rho2  = Korrelation S mit y-Brownian (NICHT rho_xy*rho1 + ...)
#
# Hinweis: Im Paper gilt rho2_tilde = rho1*rho_r + rho2*sqrt(1-rho_r^2)
# das ist die Korrelation zwischen dem y-Brownschen und dem Aktien-
# Brownschen. Dein params$rho_xy ist rho_r aus dem Paper.

berechne_sigma <- function(a, b, nu, eta, rho_r, sigma, rho1, rho2) {
  
  # Hilfsfunktion g_k(s) = (1 - exp(-k*s)) / k
  g <- function(k, s) (1 - exp(-k * s)) / k
  
  ga1  <- g(a,     1)   # g_a(1)   = (1-exp(-a))/a
  gb1  <- g(b,     1)   # g_b(1)
  gab1 <- g(a + b, 1)   # g_{a+b}(1)
  
  # Exakte OU-Varianzen: Var(x_t | x_{t-1}) = nu^2 * g_{2a}(1)
  # Das Paper schreibt "nu^2 * g_a^2(1)" — das ist die Notation
  # g_a^2 = g_{2a}  (quadratisches Argument, nicht quadrierter Wert).
  # Korrekte Herleitung: Ito-Isometrie aus Appendix C, Gl. (C3):
  #   Var(x_t) = nu^2 * integral_0^1 exp(-2a(1-s)) ds = nu^2 * g_{2a}(1)
  g2a1 <- g(2*a, 1)   # g_{2a}(1) = (1-exp(-2a))/(2a)
  g2b1 <- g(2*b, 1)   # g_{2b}(1)
  
  # rho2_tilde: effektive Korrelation Aktie-y (Formel unter Gl. 6)
  rho2_tilde <- rho1 * rho_r + rho2 * sqrt(1 - rho_r^2)
  
  # --- Kovarianzeinträge nach Lemma 2.1 (mit korrigierter Diagonale) ---
  
  # Sigma_11 = Var(x_t) — exakte OU-Formel
  S11 <- nu^2 * g2a1
  
  # Sigma_22 = Var(y_t) — exakte OU-Formel
  S22 <- eta^2 * g2b1
  
  # Sigma_33 = Var(integral r ds)
  # Formel aus Appendix C (unveraendert, nutzt ga1/gb1/gab1):
  S33 <- (nu/a)^2  * (1 - 2*ga1 + g2a1) +
    (eta/b)^2 * (1 - 2*gb1 + g2b1) +
    2 * rho_r * nu * eta / (a * b) * (1 - ga1 - gb1 + gab1)
  
  # Sigma_12 = Cov(x_t, y_t)
  S12 <- rho_r * nu * eta * gab1
  
  # Sigma_13 = Cov(x_t, integral r)
  # = nu^2/2 * g_{2a}(1) + rho_r*nu*eta/b * (ga(1) - g_{a+b}(1))
  # Hinweis: Paper schreibt nu^2/2 * ga^2(1); korrekt: nu^2/2 * g_{2a}(1)
  S13 <- nu^2 / 2 * g2a1 + rho_r * nu * eta / b * (ga1 - gab1)
  
  # Sigma_23 = Cov(y_t, integral r)
  S23 <- eta^2 / 2 * g2b1 + rho_r * nu * eta / a * (gb1 - gab1)
  
  # Delta = rho1*sigma*nu/a*(1-ga1) + rho2_tilde*sigma*eta/b*(1-gb1)
  Delta <- rho1 * sigma * nu / a * (1 - ga1) +
    rho2_tilde * sigma * eta / b * (1 - gb1)
  
  # Sigma_14 = Cov(x_t, ln S_t/S_{t-1})
  S14 <- S13 + rho1 * sigma * nu * ga1
  
  # Sigma_24 = Cov(y_t, ln S_t/S_{t-1})
  S24 <- S23 + rho2_tilde * sigma * eta * gb1
  
  # Sigma_34 = Cov(integral r, ln S_t/S_{t-1})
  S34 <- S33 + Delta
  
  # Sigma_44 = Var(ln S_t/S_{t-1})
  S44 <- S34 + sigma^2
  
  # 4x4 Kovarianzmatrix (symmetrisch)
  Sigma <- matrix(c(
    S11, S12, S13, S14,
    S12, S22, S23, S24,
    S13, S23, S33, S34,
    S14, S24, S34, S44
  ), nrow = 4, ncol = 4, byrow = TRUE)
  
  return(list(
    Sigma      = Sigma,
    ga1        = ga1,
    gb1        = gb1,
    gab1       = gab1,
    g2a1       = g2a1,
    g2b1       = g2b1,
    Delta      = Delta,
    rho2_tilde = rho2_tilde
  ))
}


# ============================================================
# MODUL 2: Gitter für (x, y) nach Gl. (10) und (11)
# ============================================================
# Erstellt N^2 Gitterpunkte für (x, y) als gleichmäßiges
# 2D-Gitter auf der GEMEINSAMEN (marginalen) Verteilung von (x_t, y_t).
#
# Hintergrund: Das Paper (Gl. 10/11) konstruiert das y-Gitter bedingt
# auf x (sd_y|x). Das funktioniert gut wenn |rho_xy| klein ist.
# Bei starker Korrelation (|rho| > ~0.9, wie hier rho=-0.992) kollabiert
# sd_y|x auf nahezu null und das Gitter deckt nur einen Bruchteil der
# marginalen y-Verteilung ab — die Zeilensummen der P-Matrix werden
# weit unter 1. Daher: Gitter auf MARGINALE Verteilung aufbauen.
# Die pmvnorm-Berechnung in Modul 3 ist davon unabhängig korrekt.

erstelle_gitter <- function(N, delta_g, Sigma) {
  
  S11 <- Sigma[1, 1]
  S12 <- Sigma[1, 2]
  S22 <- Sigma[2, 2]
  
  sd_x <- sqrt(S11)
  sd_y <- sqrt(S22)
  
  # Bedingte Streuung y|x (für Information, nicht für Gitterbau)
  var_y_given_x <- S22 - S12^2 / S11
  sd_y_given_x  <- sqrt(max(var_y_given_x, 0))  # Schutz gegen num. Fehler
  rho_xy        <- S12 / (sd_x * sd_y)
  
  # Gitterbau-Strategie: bedingte Konstruktion (Paper) wenn |rho| < 0.9,
  # sonst marginale Konstruktion (robuster bei extremer Korrelation).
  use_conditional <- abs(rho_xy) < 0.9
  
  # x-Gitter: immer auf marginaler x-Verteilung
  x_grid <- seq(-delta_g * sd_x, delta_g * sd_x, length.out = N)
  dx <- (x_grid[N] - x_grid[1]) / (2 * max(N - 1, 1))
  
  n2 <- N^2
  x_pairs <- numeric(n2)
  y_pairs <- numeric(n2)
  ki_idx  <- integer(n2)
  li_idx  <- integer(n2)
  
  if (use_conditional) {
    # --- Originale Paper-Methode (Gl. 11): y-Gitter bedingt auf x_i ---
    dy <- (2 * delta_g * sd_y_given_x) / (2 * max(N - 1, 1))
    idx <- 1
    for (ki in 1:N) {
      xi            <- x_grid[ki]
      mu_y_given_xi <- S12 / S11 * xi
      y_grid_i <- seq(
        mu_y_given_xi - delta_g * sd_y_given_x,
        mu_y_given_xi + delta_g * sd_y_given_x,
        length.out = N
      )
      for (li in 1:N) {
        x_pairs[idx] <- xi
        y_pairs[idx] <- y_grid_i[li]
        ki_idx[idx]  <- ki
        li_idx[idx]  <- li
        idx <- idx + 1
      }
    }
    cat(sprintf(
      "  Gitter: Paper-Methode (|rho|=%.3f < 0.9), sd_y|x=%.5f\n",
      abs(rho_xy), sd_y_given_x
    ))
  } else {
    # --- Marginale Methode (robust bei |rho| nahe 1): ---
    # y-Gitter gleichmäßig auf [-delta_g*sd_y, +delta_g*sd_y]
    # für jedes x_i gleich. Die bivariate Struktur wird allein
    # über pmvnorm in Modul 3 erfasst.
    y_grid <- seq(-delta_g * sd_y, delta_g * sd_y, length.out = N)
    dy <- (y_grid[N] - y_grid[1]) / (2 * max(N - 1, 1))
    idx <- 1
    for (ki in 1:N) {
      for (li in 1:N) {
        x_pairs[idx] <- x_grid[ki]
        y_pairs[idx] <- y_grid[li]
        ki_idx[idx]  <- ki
        li_idx[idx]  <- li
        idx <- idx + 1
      }
    }
    cat(sprintf(
      "  Gitter: Marginale Methode (|rho|=%.3f >= 0.9), sd_y_marginal=%.5f\n",
      abs(rho_xy), sd_y
    ))
  }
  
  return(list(
    x             = x_pairs,
    y             = y_pairs,
    ki            = ki_idx,
    li            = li_idx,
    x_grid        = x_grid,
    dx            = dx,
    dy            = dy,
    sd_y_gx       = sd_y_given_x,
    var_y_gx      = var_y_given_x,
    sd_y_marginal = sd_y,
    rho_xy        = rho_xy,
    N             = N,
    delta_g       = delta_g
  ))
}


# ============================================================
# MODUL 3: Übergangswahrscheinlichkeiten P (Gl. 13)
# ============================================================
# Berechnet die N^2 x N^2 Matrix P der Übergangswahrscheinlichkeiten.
# P[i,j] = P(x_t in Intervall j, y_t in Intervall j | x_{t-1}=x_i, y_{t-1}=y_i)
#
# Zeitinvariant! Wird einmalig berechnet.
# Verwendet bivariate Normalverteilung aus dem mvtnorm-Paket.
#
# Laufzeit: O(N^4) — für N=10 ca. 1-2 Sek., N=20 ca. 20-30 Sek.

berechne_P_matrix <- function(gitter, sig, a, b) {
  
  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    stop("Paket 'mvtnorm' wird benoetigt: install.packages('mvtnorm')")
  }
  
  N   <- gitter$N
  n2  <- N^2
  dx  <- gitter$dx
  dy  <- gitter$dy
  
  # Aus Sigma: Kovarianz von (x_t, y_t) | (x_{t-1}, y_{t-1})
  # Die BEDINGTE Kovarianzmatrix von (x_t, y_t) gegeben F_{t-1}
  # ist der 2x2-Block Sigma[1:2, 1:2] aus Lemma 2.1
  # (unabhängig von x_{t-1}, y_{t-1} — Markov-Eigenschaft)
  S11 <- sig$Sigma[1, 1]
  S12 <- sig$Sigma[1, 2]
  S22 <- sig$Sigma[2, 2]
  Sigma_xy <- matrix(c(S11, S12, S12, S22), 2, 2)
  
  exp_a <- exp(-a)
  exp_b <- exp(-b)
  
  P <- matrix(0, nrow = n2, ncol = n2)
  
  cat("Berechne P-Matrix (", n2, "x", n2, "Szenarien)...\n")
  pb <- txtProgressBar(min = 0, max = n2, style = 3)
  
  for (i in 1:n2) {
    # Ausgangszustand: x_{t-1} = gitter$x[i], y_{t-1} = gitter$y[i]
    xi_prev <- gitter$x[i]
    yi_prev <- gitter$y[i]
    
    # Bedingte Erwartung von (x_t, y_t) | (x_{t-1}, y_{t-1}) nach Lemma 2.1
    # mu1 = e^{-a} * x_{t-1}  (lambda_x = 0 hier, da wir die INNOVATIONS-
    # Verteilung meinen; lambda_x wird im phi-Term berücksichtigt)
    # Korrekt: mu_xt = e^{-a}*x_{t-1}, mu_yt = e^{-b}*y_{t-1}
    mu_xt <- exp_a * xi_prev
    mu_yt <- exp_b * yi_prev
    mu_cond <- c(mu_xt, mu_yt)
    
    for (j in 1:n2) {
      # Zielzustand: x_t in [x_j - dx, x_j + dx], y_t in [y_j - dy, y_j + dy]
      xj <- gitter$x[j]
      yj <- gitter$y[j]
      
      lower <- c(xj - dx, yj - dy)
      upper <- c(xj + dx, yj + dy)
      
      # Bivariate Normalwahrscheinlichkeit P(lower < (x_t,y_t) < upper | F_{t-1})
      P[i, j] <- mvtnorm::pmvnorm(
        lower = lower,
        upper = upper,
        mean  = mu_cond,
        sigma = Sigma_xy
      )[1]
    }
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Zeilensummen-Check (jede Zeile sollte ≈ 1 sein, Abweichung = Randmasse)
  row_sums <- rowSums(P)
  cat(sprintf(
    "P-Matrix: min Zeilensumme = %.4f, max = %.4f (ideal: 1.0)\n",
    min(row_sums), max(row_sums)
  ))
  
  return(P)
}


# ============================================================
# MODUL 4: Bedingte MGF phi nach Theorem 3.1
# ============================================================
# Berechnet für jedes Szenario-Paar (i,j) den Wert
#   phi_{x_{t-1}, y_{t-1}, x_t, y_t}(u=1, t)
# = E[e^{L_{t-1,t}} | x_{t-1}, y_{t-1}, x_t, y_t]
#
# Das ergibt direkt den bedingten Erwartungswert des diskontierten
# jährlichen Beitrags zur EIA-Ablaufleistung.
#
# Zeitabhängig wegen int_psi(t-1,t)!

berechne_phi_matrix <- function(t_jahr, gitter, sig,
                                a, b, nu, eta, rho_r,
                                sigma, rho1, rho2_tilde,
                                lambda_x, lambda_y, lambda_S,
                                int_psi_t,
                                g_garantie, alpha) {
  
  N   <- gitter$N
  n2  <- N^2
  
  # Abkürzungen aus der Sigma-Struktur
  S   <- sig$Sigma      # 4x4
  S11 <- S[1,1]; S12 <- S[1,2]; S13 <- S[1,3]
  S22 <- S[2,2]; S23 <- S[2,3]; S33 <- S[3,3]
  ga1 <- sig$ga1
  gb1 <- sig$gb1
  Delta <- sig$Delta
  
  exp_a <- exp(-a)
  exp_b <- exp(-b)
  
  # Determinante des 2x2-Blocks [[S11, S12],[S12, S22]]
  det_xy <- S11 * S22 - S12^2
  
  # Vektoren für alle n2 Szenarien (Ausgangszustand i)
  # Erwartungswerte mu1(i), mu2(i) für gegebenes (x_{t-1}, y_{t-1})
  mu1_vec <- exp_a * gitter$x + lambda_x * (1 - exp_a)
  mu2_vec <- exp_b * gitter$y + lambda_y * (1 - exp_b)
  
  # mu3: bedingte Erwartung des Integrals ∫r ds | F_{t-1}
  # mu3 = int_psi(t-1,t) + ga1*x_{t-1} + gb1*y_{t-1} + lambda_x*(1-ga1) + lambda_y*(1-gb1)
  mu3_vec <- int_psi_t +
    sig$ga1 * gitter$x +
    sig$gb1 * gitter$y +
    lambda_x * (1 - sig$ga1) +
    lambda_y * (1 - sig$gb1)
  
  # mu4: bedingte Erwartung von ln(S_t/S_{t-1}) | F_{t-1}
  # = mu3 + lambda_S - 0.5*sigma^2
  mu4_vec <- mu3_vec + lambda_S - 0.5 * sigma^2
  
  # --- Für u=1: Tilted parameters (Theorem 3.1) ---
  # mu1_tilde(u=1) = mu1 - sigma13
  # mu2_tilde(u=1) = mu2 - sigma23
  # mu3_tilde(u=1) = mu3 + lambda_S - 0.5*sigma^2
  #                        - (S33 + sigma*(rho1*nu/a*(1-ga1) + rho2_tilde*eta/b*(1-gb1)))
  sigma_sum <- sigma * (rho1 * nu / a * (1 - ga1) + rho2_tilde * eta / b * (1 - gb1))
  
  mu1t_vec <- mu1_vec - S13      # = mu1 - u*S13, u=1
  mu2t_vec <- mu2_vec - S23      # = mu2 - u*S23, u=1
  mu3t_base <- lambda_S - 0.5 * sigma^2 - (S33 + sigma_sum)
  # mu3_tilde = int_psi + ga1*x + gb1*y + lambda_x*(1-ga1) + lambda_y*(1-gb1) + mu3t_base
  
  # Sigma_tilde (gleiche Kovarianzstruktur der ersten 3 Komponenten,
  # aber mit "~"-Einträgen für die Aktienkomponente):
  S13t <- S13 + sigma * rho1 * nu * ga1
  S23t <- S23 + sigma * rho2_tilde * eta * gb1
  S33t <- sigma^2 + S33 + sigma * (rho1 * nu / a * (1 - ga1) + rho2_tilde * eta / b * (1 - gb1))
  # S11t=S11, S22t=S22, S12t=S12 (unveraendert)
  
  det_xy_t <- S11 * S22 - S12^2  # = det_xy (S11,S12,S22 unveraendert)
  
  # Bedingte Varianz des Integrals ∫r ds | x_t, y_t (Formel für sigma_bar^2):
  # sigma_bar^2 = S33 - [S13, S23] * inv([[S11,S12],[S12,S22]]) * [S13,S23]'
  sigma_bar2 <- S33 -
    (S22 * S13^2 - 2 * S12 * S13 * S23 + S11 * S23^2) / det_xy
  
  # Bedingte Varianz von ln(S) | x_t, y_t (unter Qt_u, Theorem 3.1):
  # sigma2_tilde = S33t - [S13t,S23t]*inv([S11,S12;S12,S22])*[S13t,S23t]'
  sigma2_tilde <- S33t -
    (S22 * S13t^2 - 2 * S12 * S13t * S23t + S11 * S23t^2) / det_xy_t
  
  # Matrix C(t) aufbauen: Zeile i = Ausgangszustand, Spalte j = Zielzustand
  C <- matrix(0, nrow = n2, ncol = n2)
  
  # Hilfskonstante für Inverse 2x2:
  # inv([[S11,S12],[S12,S22]]) = 1/det * [[S22,-S12],[-S12,S11]]
  # [S13,S23] * inv * [dx,dy]' = (S22*S13 - S12*S23)*dx/det + (S11*S23 - S12*S13)*dy/det
  
  A_x <- (S22 * S13 - S12 * S23) / det_xy    # Koeffizient für (x_j - mu1_i)
  A_y <- (S11 * S23 - S12 * S13) / det_xy    # Koeffizient für (y_j - mu2_i)
  
  # Tilted Koeffizienten:
  At_x <- (S22 * S13t - S12 * S23t) / det_xy_t
  At_y <- (S11 * S23t - S12 * S13t) / det_xy_t
  
  for (i in 1:n2) {
    
    # Bedingte Erwartung von ∫r | x_t=x_j, y_t=y_j
    # mu_bar(i,j) = mu3_vec[i] + A_x*(x_j - mu1_vec[i]) + A_y*(y_j - mu2_vec[i])
    mu_bar_vec <- mu3_vec[i] +
      A_x * (gitter$x - mu1_vec[i]) +
      A_y * (gitter$y - mu2_vec[i])
    
    # Bedingte Erwartung von ln(S) | x_t=x_j, y_t=y_j (unter Qt_{u=1})
    mu3t_vec_i <- int_psi_t + mu3t_base +
      sig$ga1 * gitter$x[i] +
      sig$gb1 * gitter$y[i] +
      lambda_x * (1 - sig$ga1) +
      lambda_y * (1 - sig$gb1)
    
    mu_tilde_vec <- mu3t_vec_i +
      At_x * (gitter$x - mu1t_vec[i]) +
      At_y * (gitter$y - mu2t_vec[i])
    
    # Theorem 3.1 (Gl. 16), u=1:
    # phi = exp(-mu_bar + sigma_bar^2/2) *
    #   [ exp(g)*Phi((g/alpha - mu_tilde)/sigma_tilde)
    #     + exp(alpha*mu_tilde + 0.5*alpha^2*sigma2_tilde)
    #       * Phi((-g/alpha + mu_tilde + alpha*sigma2_tilde)/sigma_tilde) ]
    
    sigma_tilde <- sqrt(sigma2_tilde)
    
    term1 <- exp(g_garantie) *
      pnorm((g_garantie / alpha - mu_tilde_vec) / sigma_tilde)
    
    term2 <- exp(alpha * mu_tilde_vec + 0.5 * alpha^2 * sigma2_tilde) *
      pnorm((-g_garantie / alpha + mu_tilde_vec + alpha * sigma2_tilde) / sigma_tilde)
    
    C[i, ] <- exp(-mu_bar_vec + sigma_bar2 / 2) * (term1 + term2)
  }
  
  return(C)
}


# ============================================================
# MODUL 5: Hauptfunktion — Iteration nach Algorithm 3.1
# ============================================================
# Berechnet E[Z_T] via Szenario-Matrix-Ansatz.
#
# Argumente:
#   params      — Liste aus parameter_einlesen.R
#   vb          — Liste aus vorberechnungen.R (psi, int_psi_t_T, V_t_T)
#   N           — Anzahl Gitterpunkte pro Dimension (default: 10)
#   delta_grid  — Breite des Diskretisierungsintervalls (default: 8.25)
#   T_jahre     — Zeithorizont in Jahren (default: 25)
#   g_garantie  — Jährliche Mindestrendite (stetig, default: 0.005 = 0.5%)
#   alpha       — Partizipationsrate (default: 0.465)
#   rho1        — Korrelation Aktie mit x-Brownschen (default: -0.15)
#   rho2        — Korrelation Aktie mit y-Brownschen (default: +0.15)
#   verbose     — Fortschrittsausgabe (default: TRUE)
#
# Rückgabe: Liste mit
#   E_ZT        — E[diskontierte Ablaufleistung Z_T]
#   laufzeit    — Gesamtlaufzeit in Sekunden
#   N, delta_grid — verwendete Parameter
#   A_verlauf   — Vektor der akkumulierten Erwartungswerte nach jedem Jahr

sm_berechnung <- function(params, vb,
                          N           = 10,
                          delta_grid  = 8.25,
                          T_jahre     = 25,
                          g_garantie  = 0.005,
                          alpha       = 0.465,
                          rho1        = -0.15,
                          rho2        =  0.15,
                          verbose     = TRUE) {
  
  start_gesamt <- Sys.time()
  
  # Parameter auspacken
  a       <- params$a_x
  b       <- params$a_y
  nu      <- params$sigma_x
  eta     <- params$sigma_y
  rho_r   <- params$rho_xy
  sigma   <- params$sigma_s
  lambda_S <- params$lambda_N   # Risikoprämie Aktie (Normal-Szenario)
  delta   <- params$delta       # Monatlicher Zeitschritt (= 1/12)
  
  # Zeitabhängige Risikoprämien: vor/nach Tau
  # Für den SM-Ansatz auf Jahresbasis vereinfachen wir:
  # lambda_x/y = d_x/d_y im ersten Jahr (t=1..tau/12 Jahre),
  # danach l_x/l_y
  # (tau ist in Monaten angegeben)
  tau_jahre <- params$tau / 12
  
  # int_psi(t-1, t) für t=1..T_jahre (jährliche Integrale)
  # In vb$int_psi_t_T ist das kumulierte Integral gespeichert.
  # int_psi(t-1,t) = int_psi_t_T[1, 12*(t-1)+1 .. 12*t+1] aufsummiert
  # Da int_psi_t_T monatliche Trapezintegrale enthält:
  # int_psi_jahres[t] = sum_{m=1}^{12} int_psi_neu[12*(t-1) + m]
  #
  # Einfacher: int_psi_t_T[1, 12*t+1] - int_psi_t_T[1, 12*(t-1)+1]
  # = int_psi_t_T[Zeile 1 = t=0, kumuliert von 0 bis 12*t Monate]
  int_psi_jahres <- numeric(T_jahre)
  for (t in 1:T_jahre) {
    # Monatsspalten: t=0 -> Spalte 1, t=k -> Spalte k+1
    col_start <- 12 * (t - 1) + 1 + 1   # Spaltenindex für Monat 12*(t-1)
    col_end   <- 12 * t + 1              # Spaltenindex für Monat 12*t
    # kumuliertes Integral von 0 bis 12*t minus kumuliertes bis 12*(t-1)
    int_psi_jahres[t] <- vb$int_psi_t_T[1, col_end] - vb$int_psi_t_T[1, col_start - 1 + 1]
    # Alternativformel (direkter):
    # int_psi_jahres[t] <- vb$int_psi_t_T[1, 12*t + 1] - vb$int_psi_t_T[1, 12*(t-1) + 1]
  }
  
  if (verbose) cat("=== Szenario-Matrix-Ansatz (N =", N, ", delta =", delta_grid, ") ===\n")
  
  # --- Modul 1: Sigma ---
  if (verbose) cat("Berechne Kovarianzmatrix Sigma (jaehrlich)...\n")
  # Achtung: Sigma ist auf JAHRESBASIS (Zeitschritt = 1 Jahr)
  # aber unsere Parameter a, b, nu, eta sind monatlich kalibriert!
  # Wir müssen auf Jahresbasis skalieren: ga(1) mit a_jaehrlich = a_monatlich * (1/12)? 
  # NEIN: Im G2++/Vasicek-Modell sind a,b dimensionslos-pro-Zeiteinheit.
  # Wenn delta=1/12 und a in 1/Monat: dann ist ga(1 Jahr) = (1-exp(-a*12))/a
  # Das Paper arbeitet mit Jahresschritten, also ga(1) = (1-exp(-a))/a mit a=Jahresrate.
  # Dein params$a_x ist die MONATLICHE mean-reversion-speed (in simulation.R: exp(-a_x*delta))
  # -> a_jaehrlich = a_x (da delta bereits herausmultipliziert ist in den Formeln)
  # Aber: g_a(1) im Paper = (1-exp(-a*1))/a mit Jahresschritt
  # -> Im Code: g_a(1 Jahr) = (1 - exp(-a_x * 12 * delta)) / a_x = (1-exp(-a_x))/a_x
  #    da 12*delta = 1 Jahr
  # Fazit: Wenn a_x in deinen params die monatliche Geschwindigkeit ist,
  # dann verwenden wir a_x direkt (sie bezieht sich auf 1 Monat als Zeiteinheit).
  # Für Jahresschritte: ersetze a -> a*12 in allen ga(1)-Ausdrücken? 
  # -> NEIN! Im Paper gilt: ga(s) = (1-exp(-a*s))/a, und s=1 bedeutet 1 JAHR.
  # -> Wir müssen die jährliche mean-reversion angeben: a_jahr = a_monat * 12.
  # Hinweis: Das wird im Kommentar explizit ausgewiesen.
  
  # Konversion: Modellparameter von monatlich auf jährlich
  # (da SM-Ansatz auf jährlichem Gitter operiert)
  a_jahr <- a      # params$a_x ist bereits auf Monatsbasis; da delta=1/12:
  b_jahr <- b      # exp(-a_x * delta) = exp(-a_x/12) pro Monat
  # => pro Jahr: exp(-a_x) = exp(-a_jahr * 1)
  # => a_jahr = a_x (gleicher numerischer Wert, andere Zeiteinheit)
  
  # Wichtiger Unterschied: sigma muss ebenfalls auf Jahresbasis skaliert werden!
  # In simulation.R: Z_x * sigma_x * sqrt((1 - exp(-2*a_x*delta)) / (2*a_x))
  # Das ist die MONATLICHE Innovations-Std.abw.
  # Für Lemma 2.1 (Jahresschritt): nu_jahr = sigma_x (entspricht nu im Paper)
  # da die Paper-Formeln für g_a(1) mit a = a_monat * 12 kalibriert sein müssen.
  #
  # PRAGMATISCHE LÖSUNG: Wir übergeben a_x direkt und setzen s=1 in g_k(s).
  # ga(1) = (1-exp(-a_x))/a_x, das entspricht genau der Jahres-Übergangsdynamik
  # in deiner Simulation (exp(-a_x * delta) mit delta=1 Jahr = 12*1/12).
  
  sig <- berechne_sigma(
    a     = a_jahr,
    b     = b_jahr,
    nu    = nu,
    eta   = eta,
    rho_r = rho_r,
    sigma = sigma,
    rho1  = rho1,
    rho2  = rho2
  )
  
  if (verbose) {
    cat(sprintf("  Sigma[1,1] = %.6f (Var x_t)\n",  sig$Sigma[1,1]))
    cat(sprintf("  Sigma[2,2] = %.6f (Var y_t)\n",  sig$Sigma[2,2]))
    cat(sprintf("  Sigma[3,3] = %.6f (Var int r)\n", sig$Sigma[3,3]))
    cat(sprintf("  Sigma[4,4] = %.6f (Var ln S)\n",  sig$Sigma[4,4]))
    cat(sprintf("  ga(1)=%.6f, gb(1)=%.6f\n", sig$ga1, sig$gb1))
  }
  
  # --- Modul 2: Gitter ---
  if (verbose) cat(sprintf("Erstelle Gitter (%d x %d = %d Szenarien)...\n", N, N, N^2))
  gitter <- erstelle_gitter(N, delta_grid, sig$Sigma)
  
  # Plausibilitätsprüfung: implizierter effektiver Zinsbereich an den Gitterrändern.
  # Da x und y stark korreliert sind (rho=-0.992), kompensieren sie sich fast
  # vollständig. Der relevante Zinsbereich basiert auf sd(x+y), nicht x_max+y_max.
  if (verbose) {
    psi_1j   <- int_psi_jahres[1]
    S11      <- sig$Sigma[1,1]; S12 <- sig$Sigma[1,2]; S22 <- sig$Sigma[2,2]
    sd_r_eff <- sqrt(max(S11 + S22 + 2*S12, 0))   # effektive sd(x_t + y_t)
    r_min    <- psi_1j - delta_grid * sd_r_eff
    r_max    <- psi_1j + delta_grid * sd_r_eff
    cat(sprintf(
      "  Effektiver Zinsbereich (delta*sd(x+y)): r_min=%.1f%%, r_max=%.1f%%\n",
      r_min * 100, r_max * 100
    ))
    if (r_min < -0.20 || r_max > 0.35) {
      warning(sprintf(paste0(
        "delta_grid=%.2f erzeugt einen sehr weiten Zinsbereich (%.1f%% .. %.1f%%).\n",
        "  Empfehlung: delta_grid zwischen 3.0 und 8.25 fuer diese Parameter pruefen."
      ), delta_grid, r_min * 100, r_max * 100))
    }
  }
  
  if (verbose) cat("Berechne P-Matrix...\n")
  start_P <- Sys.time()
  P <- berechne_P_matrix(gitter, sig, a_jahr, b_jahr)
  if (verbose) cat(sprintf("  P-Matrix: %.1f Sekunden\n",
                           as.numeric(difftime(Sys.time(), start_P, units="secs"))))
  
  # --- Anfangszustand: (x_0, y_0) = (0, 0) ---
  # Index i_0: Welcher Gitterpunkt ist am nächsten an (0,0)?
  abstand <- (gitter$x)^2 + (gitter$y)^2
  i0 <- which.min(abstand)
  if (verbose) cat(sprintf("  Anfangszustand i0=%d: x=%.4f, y=%.4f\n",
                           i0, gitter$x[i0], gitter$y[i0]))
  
  # --- Initialer Zustandsvektor A_0 ---
  # A_0 ist ein N^2-Vektor mit 1 in Position i0 (= Startbedingung x_0=y_0=0)
  n2 <- N^2
  A <- numeric(n2)
  A[i0] <- 1.0
  
  A_verlauf <- numeric(T_jahre)
  
  # --- Modul 4 & 5: Iteration über t = 1..T_jahre ---
  if (verbose) cat(sprintf("Starte Iteration (T=%d Jahre)...\n", T_jahre))
  start_iter <- Sys.time()
  
  for (t in 1:T_jahre) {
    
    # Zeitabhängige Risikoprämien (Monat 12*(t-1)+1 .. 12*t)
    # Vereinfachung: Wert in der Mitte des Jahres entscheidet
    t_mitte <- 12 * t - 6  # Monat in der Mitte des t-ten Jahres
    if (t_mitte <= params$tau) {
      lambda_x <- params$d_x
      lambda_y <- params$d_y
    } else {
      lambda_x <- params$l_x
      lambda_y <- params$l_y
    }
    
    # Jährliches Integral int_psi(t-1, t)
    int_psi_t <- int_psi_jahres[t]
    
    # Phi-Matrix C(t)
    C_t <- berechne_phi_matrix(
      t_jahr      = t,
      gitter      = gitter,
      sig         = sig,
      a           = a_jahr,
      b           = b_jahr,
      nu          = nu,
      eta         = eta,
      rho_r       = rho_r,
      sigma       = sigma,
      rho1        = rho1,
      rho2_tilde  = sig$rho2_tilde,
      lambda_x    = lambda_x,
      lambda_y    = lambda_y,
      lambda_S    = lambda_S,
      int_psi_t   = int_psi_t,
      g_garantie  = g_garantie,
      alpha       = alpha
    )
    
    # Szenariomatrix Q(t) = P elementweise * C(t)  [Gl. 14]
    # Q[i,j] = P(Zustand i -> Zustand j) * phi(i->j)
    Q_t <- P * C_t
    
    # Update: A_t[j] = sum_i Q[i,j] * A_{t-1}[i]  [Gl. 15]
    # = (t(Q) %*% A)[j] — Transposition nötig, da Q Zeile=von, Spalte=nach
    # crossprod(Q, A) = t(Q) %*% A, effizienter als explizites t()
    A <- as.vector(crossprod(Q_t, A))
    
    # E[Z_t] = Summe aller Einträge von A (da 1' * A_T = E[Z_T])
    A_verlauf[t] <- sum(A)
    
    if (verbose) {
      cat(sprintf("  Jahr %2d/%d: E[Z_t] = %.6f\n", t, T_jahre, A_verlauf[t]))
    }
  }
  
  laufzeit <- as.numeric(difftime(Sys.time(), start_gesamt, units = "secs"))
  
  # E[Z_T] = sum(A_T) gemäß Algorithm 3.1
  E_ZT <- A_verlauf[T_jahre]
  
  if (verbose) {
    cat(sprintf("\n=== Ergebnis ===\n"))
    cat(sprintf("  E[Z_T] = %.6f\n", E_ZT))
    cat(sprintf("  Laufzeit gesamt: %.1f Sekunden\n", laufzeit))
  }
  
  return(list(
    E_ZT       = E_ZT,
    A_verlauf  = A_verlauf,
    laufzeit   = laufzeit,
    N          = N,
    delta_grid = delta_grid,
    T_jahre    = T_jahre,
    g_garantie = g_garantie,
    alpha      = alpha,
    gitter     = gitter,
    P          = P,
    sig        = sig
  ))
}


# ============================================================
# HILFSFUNKTIONEN: Validierung & Debugging
# ============================================================

# Prüft die Sigma-Matrix gegen bekannte Paper-Werte (Table 1)
validiere_sigma_paper <- function() {
  cat("=== Validierung Sigma gegen Paper-Parameter (Table 1) ===\n")
  
  # Parameter aus Table 1 des Papers
  a_p    <- 0.3912;  b_p <- 0.0785
  nu_p   <- 0.0201;  eta_p <- 0.0135
  rho_r_p <- -0.6450
  sigma_p <- 0.1
  rho1_p  <- -0.15;  rho2_p <- 0.15
  
  sig <- berechne_sigma(a_p, b_p, nu_p, eta_p, rho_r_p, sigma_p, rho1_p, rho2_p)
  
  cat("Kovarianzmatrix Sigma:\n")
  print(round(sig$Sigma, 8))
  cat(sprintf("ga(1) = %.6f\n", sig$ga1))
  cat(sprintf("gb(1) = %.6f\n", sig$gb1))
  cat(sprintf("rho2_tilde = %.6f\n", sig$rho2_tilde))
  
  # Positivitätsprüfung
  ev <- eigen(sig$Sigma)$values
  cat(sprintf("Eigenwerte Sigma: %s\n", paste(round(ev, 8), collapse=", ")))
  cat(sprintf("Sigma positiv definit: %s\n", all(ev > 0)))
  
  return(invisible(sig))
}

# Prüft das Gitter: Visualisierung analog Fig. 2 im Paper
validiere_gitter <- function(N = 5, delta_g = 2.0) {
  cat("=== Validierung Gitter (N=", N, ", delta=", delta_g, ") ===\n")
  
  # Paper-Parameter
  a_p <- 0.3912; b_p <- 0.0785
  nu_p <- 0.0201; eta_p <- 0.0135
  rho_r_p <- -0.6450
  sigma_p <- 0.1; rho1_p <- -0.15; rho2_p <- 0.15
  
  sig <- berechne_sigma(a_p, b_p, nu_p, eta_p, rho_r_p, sigma_p, rho1_p, rho2_p)
  g   <- erstelle_gitter(N, delta_g, sig$Sigma)
  
  cat(sprintf("N^2 = %d Gitterpunkte\n", N^2))
  cat(sprintf("x in [%.4f, %.4f]\n", min(g$x), max(g$x)))
  cat(sprintf("y in [%.4f, %.4f]\n", min(g$y), max(g$y)))
  cat("Erste 5 (x,y)-Paare:\n")
  print(data.frame(x=round(g$x[1:5],5), y=round(g$y[1:5],5)))
  
  return(invisible(g))
}