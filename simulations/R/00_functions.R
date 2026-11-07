
# =============================================================================
# Basic Formula Functions
# =============================================================================

#' m(-lambda): Stieltjes transform evaluated at -lambda
#' @param lambda Regularization parameter (positive)
#' @param gamma Aspect ratio p/n
m_neg_lambda <- function(lambda, gamma) {
 ((lambda - gamma + 1) - sqrt((lambda + gamma + 1)^2 - 4 * gamma)) / (2 * gamma * (-lambda))
}


#' sqrt_psi(lambda) = sqrt(m(-lambda) - lambda * m'(-lambda) + sigma^2 * lambda^2 * m'(-lambda))
#' @param lambda Regularization parameter (positive)
#' @param gamma Aspect ratio p/n
#' @param sigma Signal strength sigma_beta (default 0 gives the noise-only term)
sqrt_psi <- function(lambda, gamma, sigma) {
  m <- m_neg_lambda(lambda, gamma)
  m_prime <- m_neg_lambda_deriv(lambda, gamma)
  return(sqrt(m - lambda * m_prime + sigma^2 * lambda^2 * m_prime))
}

#' sqrt_psi with separate sigma_beta and sigma_epsilon
#' @param lambda Regularization parameter (positive)
#' @param gamma Aspect ratio p/n
#' @param sigma_beta Signal strength
#' @param sigma_epsilon Noise standard deviation
sqrt_psi_beta_epsilon <- function(lambda, gamma, sigma_beta, sigma_epsilon) {
  m <- m_neg_lambda(lambda, gamma)
  m_prime <- m_neg_lambda_deriv(lambda, gamma)
  return(sqrt(sigma_epsilon^2 * (m - lambda * m_prime) + sigma_beta^2 * lambda^2 * m_prime))
}


# =============================================================================
# Frechet Distribution Functions
# =============================================================================

#' Frechet CDF: F(t) = exp(-(t/scale)^{-alpha}) for t > 0
#' @param t Quantile value(s)
#' @param alpha Shape parameter (tail index)
#' @param scale Scale parameter = lambda * m(-lambda) * sigma_beta
frechet_cdf <- function(t, alpha, scale) {
  ifelse(t <= 0, 0, exp(-(t / scale)^(-alpha)))
}

#' Frechet quantile function (inverse CDF)
#' F^{-1}(u) = scale * (-log(u))^{-1/alpha}
#' @param u Probability value(s) in (0, 1)
#' @param alpha Shape parameter (tail index)
#' @param scale Scale parameter
frechet_quantile <- function(u, alpha, scale) {
  scale * (-log(u))^(-1/alpha)
}


m_neg_lambda_deriv <- function(lambda, gamma) {
  sqrt_term <- sqrt((gamma + lambda + 1)^2 - 4 * gamma)
  term1 <- gamma + lambda - sqrt_term + 1
  term2 <- gamma - lambda + sqrt_term - 1
  numerator <- -lambda * term1 + sqrt_term * term2
  denominator <- 2 * gamma * lambda^2 * sqrt_term
  return(numerator / denominator)
}

# =============================================================================
# Hill Estimator
# =============================================================================

#' Hill estimator of tail index alpha
#' @param x Vector of positive values (e.g. |X'y/n|)
#' @param k Number of top order statistics to use
#' @return Estimated tail index alpha_hat
hill_estimator <- function(x, k) {
  stopifnot(k >= 2, k < length(x))
  x_sorted <- sort(x, decreasing = TRUE)
  1 / mean(log(x_sorted[1:k] / x_sorted[k + 1]))
}

# =============================================================================
# CI Half-Width (Frechet quantile based)
# =============================================================================

#' Compute simultaneous CI half-width Delta
#' @param alpha_val Tail index
#' @param lambda Ridge parameter
#' @param gamma_val Aspect ratio p/n
#' @param sigma_beta Signal strength
#' @param n Sample size
#' @param q Significance level (e.g. 0.05 for 95% CI)
compute_delta <- function(alpha_val, lambda, gamma_val, sigma_beta, n, q) {
  c_alpha <- 2 * gamma((alpha_val + 1) / 2) /
    (sqrt(alpha_val * pi) * gamma(alpha_val / 2)) *
    alpha_val^((alpha_val - 1) / 2)
  ell_n <- c_alpha^(1 / alpha_val)
  m_val <- m_neg_lambda(lambda, gamma_val)
  S <- ell_n * lambda * m_val * sigma_beta
  t_1q <- S * (-log(1 - q) / gamma_val)^(-1 / alpha_val)
  n^(-(0.5 - 1 / alpha_val)) * t_1q
}
