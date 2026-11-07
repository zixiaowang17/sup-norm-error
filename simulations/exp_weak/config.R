# =============================================================================
# Configuration for weak-signal experiment (Theorem: weak_signal)
#
# Two-panel figure: gamma = 0.5 (left) and gamma = 2.0 (right).
# x-axis = lambda (ridge parameter)
# y-axis = ||beta_hat - beta_0||_inf / r_n  (normalized error)
# =============================================================================

N_MC         <- 500
N_WORKERS    <- max(1, parallel::detectCores() - 1)
N_VAL        <- 1000
SIGMA_EPS    <- 1.0
RHO          <- 0
B            <- 1.2
GAMMA_VALS   <- c(0.5, 2.0)
A_VAL        <- 3.5
LAMBDA_GRID  <- exp(seq(log(1e-3), log(100), length.out = 12))

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
