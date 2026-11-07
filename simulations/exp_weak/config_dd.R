# =============================================================================
# Configuration for double-descent experiment (weak-signal regime)
#
# x-axis = gamma (aspect ratio p/n)
# y-axis = ||beta_hat - beta_0||_inf / r_n  (normalized error)
# =============================================================================

N_MC         <- 500
N_WORKERS    <- max(1, parallelly::availableCores() - 1)
N_VAL        <- 1000
SIGMA_EPS    <- 1.0
RHO          <- 0
B            <- 0.5
A_VAL        <- 2.0
LAMBDA_DD    <- 1e-5

# Gamma grid: ~20 points, denser near gamma=1
GAMMA_UNDER  <- c(0.2, 0.4, 0.6, 0.75, 0.8, 0.85, 0.9, 0.93, 0.95, 0.97)
GAMMA_OVER   <- c(1.03, 1.05, 1.07, 1.1, 1.15, 1.25, 1.5, 2.0, 2.5, 3.0)
GAMMA_GRID   <- c(GAMMA_UNDER, GAMMA_OVER)

make_beta <- function(n, p, a) {
  r_n <- sqrt(2 * log(p) / n)

  lead <- a * r_n
  residual_sq <- B - lead^2
  if (residual_sq < 0) {
    stop(sprintf("B - (a*r_n)^2 = %.4f < 0: a=%.2f is too large for B=%.2f, n=%d, p=%d",
                 residual_sq, a, B, n, p))
  }

  tail_val <- sqrt(residual_sq / (p - 1))

  beta       <- rep(tail_val, p)
  beta[1]    <- lead

  list(
    beta = beta,
    a    = a,
    r_n  = r_n,
    B    = sum(beta^2),
    beta_inf = max(abs(beta))
  )
}
