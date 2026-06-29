# vorberechnungen.R
# Deterministische Vorberechnungen: psi, V_t_T, int_psi_t_T

berechne_vorberechnungen <- function(a_x, a_y, sigma_x, sigma_y, rho_xy,
                                     beta_0, beta_1, beta_2, beta_3,
                                     tau_1, tau_2, delta,
                                     linear_beginn, linear_wert) {
  
  T_max <- 960  # 80 Jahre * 12 Monate
  
  # --- f_m und psi ---
  f_m <- numeric(T_max + 1)  # Index 0..960
  psi <- numeric(T_max + 2)  # Index 0..961
  f_m[1] <- 0
  psi[1] <- 0
  
  for (t in 1:T_max) {
    # NSS-Forwardkurve
    if (t <= linear_beginn) {
      f_m[t + 1] <- 1/100 * (
        beta_0
        + beta_1 * ((1 - exp(-t * delta / tau_1)) / (t * delta / tau_1))
        + beta_2 * ((1 - exp(-t * delta / tau_1)) / (t * delta / tau_1) - exp(-t * delta / tau_1))
        + beta_3 * ((1 - exp(-t * delta / tau_2)) / (t * delta / tau_2) - exp(-t * delta / tau_2))
      )
    } else {
      f_m[t + 1] <- linear_wert
    }
    
    # psi(t), vgl. Formel (1.3)
    psi[t + 1] <- f_m[t + 1] +
      ((sigma_x / a_x)^2) / 2 * (1 - exp(-a_x * t * delta))^2 +
      ((sigma_y / a_y)^2) / 2 * (1 - exp(-a_y * t * delta))^2 +
      rho_xy * sigma_x * sigma_y / (a_x * a_y) *
      (1 - exp(-a_x * t * delta)) * (1 - exp(-a_y * t * delta))
  }
  
  # --- V_t_T (482 x 961) ---
  # VBA: V_t_T(t+1, TT+1), t=0..481, TT=t+1..960
  # R:   V_t_T[t+1, TT+1], gleiche Indizierung
  V_t_T <- matrix(0, nrow = 482, ncol = 961)
  
  for (t in 0:481) {
    for (TT in (t + 1):960) {
      h <- TT - t  # Laufzeit in Monaten
      V_t_T[t + 1, TT + 1] <-
        (sigma_x / a_x)^2 * (
          h * delta + 2/a_x * exp(-a_x * h * delta)
          - 1/(2*a_x) * exp(-2*a_x * h * delta) - 3/(2*a_x)
        ) +
        (sigma_y / a_y)^2 * (
          h * delta + 2/a_y * exp(-a_y * h * delta)
          - 1/(2*a_y) * exp(-2*a_y * h * delta) - 3/(2*a_y)
        ) +
        2 * rho_xy * sigma_x * sigma_y / (a_x * a_y) * (
          h * delta
          + (exp(-a_x * h * delta) - 1) / a_x
          + (exp(-a_y * h * delta) - 1) / a_y
          - (exp(-(a_x + a_y) * h * delta) - 1) / (a_x + a_y)
        )
    }
  }
  
  # --- int_psi_neu (Trapezregel, delta-Intervalle) ---
  int_psi_neu <- numeric(T_max)  # Index 1..960
  for (t in 0:(T_max - 1)) {
    int_psi_neu[t + 1] <- delta * (psi[t + 2] + psi[t + 3]) / 2
  }
  
  # --- int_psi_t_T (482 x 961, kumulierte Integrale) ---
  int_psi_t_T <- matrix(0, nrow = 482, ncol = 961)
  for (t in 0:481) {
    for (TT in (t + 1):960) {
      int_psi_t_T[t + 1, TT + 1] <- int_psi_t_T[t + 1, TT] + int_psi_neu[TT]
    }
  }
  
  return(list(
    psi        = psi,
    V_t_T      = V_t_T,
    int_psi_t_T = int_psi_t_T
  ))
}