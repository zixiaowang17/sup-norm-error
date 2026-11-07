# =============================================================================
# Configuration for strong-signal experiment (Theorem 3.2)
#
# Two-panel figure: gamma = 0.5 (left) and gamma = 2.0 (right).
# x-axis = alpha, fixed lambda = 1.
# =============================================================================

N_MC         <- 500
N_WORKERS    <- 48
N_VAL        <- 1000
GAMMA_VALS   <- c(0.5, 2.0)
SIGMA_EPS    <- 1.0
LAMBDA_FIXED <- 1.0
ALPHA_GRID   <- seq(0.12, 0.20, by = 0.02)
S0           <- 1
C_STAR       <- 1

# =============================================================================
# Beta construction
#
# beta = ( c* * r_n,  1/sqrt(n), ..., 1/sqrt(n) )
#         |-- s0 --|  |--------- p - s0 ---------|
#
# noise_scale = sqrt(2 log p / n)
# r_n         = noise_scale * n^alpha
# =============================================================================

make_beta <- function(n, p, alpha, s0 = S0) {
  noise_scale <- sqrt(2 * log(p) / n)
  r_n         <- noise_scale * n^alpha

  beta       <- rep(1 / sqrt(n), p)
  beta[1:s0] <- C_STAR * r_n

  list(
    beta        = beta,
    s0          = s0,
    alpha       = alpha,
    noise_scale = noise_scale,
    r_n         = r_n,
    beta_inf    = max(abs(beta)),
    B_n         = sum(beta^2)
  )
}
