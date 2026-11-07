# =============================================================================
# Theory functions for strong-signal regime (Theorem 3.2)
# All return RAW ||beta_hat - beta_0||_inf
# =============================================================================

source("../R/00_functions.R")  # m_neg_lambda

# Ridge (lambda > 0 fixed): lambda * m_gamma(-lambda) * c* * r_n
strong_ridge_theory <- function(lambda, gamma, c_star, r_n) {
  m <- m_neg_lambda(lambda, gamma)
  lambda * m * c_star * r_n
}

# Null: c* * r_n
strong_null_theory <- function(c_star, r_n) {
  c_star * r_n
}

# MoM: sqrt(B + sigma_eps^2) * noise_scale
strong_mom_theory <- function(noise_scale, B, sigma_eps2) {
  sqrt(B + sigma_eps2) * noise_scale
}

# OLS (gamma < 1 only): sigma_eps / sqrt(1 - gamma) * noise_scale
strong_ols_theory <- function(noise_scale, gamma, sigma_eps) {
  sigma_eps / sqrt(1 - gamma) * noise_scale
}

# Ridgeless / min-norm (gamma > 1): (1 - 1/gamma) * c* * r_n
strong_ridgeless_theory <- function(gamma, c_star, r_n) {
  (1 - 1 / gamma) * c_star * r_n
}
