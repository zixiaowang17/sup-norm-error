
source("../R/00_functions.R")  # m_neg_lambda, m_neg_lambda_deriv

mp_density <- function(x, gamma) {
  a <- (1 - sqrt(gamma))^2
  b <- (1 + sqrt(gamma))^2
  ifelse(x > a & x < b,
         sqrt((b - x) * (x - a)) / (2 * pi * x),
         0)
}

mp_int_lambda_sq <- function(lambda, gamma) {
  a_mp <- (1 - sqrt(gamma))^2
  b_mp <- (1 + sqrt(gamma))^2

  # Continuous part: integrate lambda^2/(x+lambda)^2 * f(x) dx
  integrand <- function(x) {
    lambda^2 / (x + lambda)^2 * mp_density(x, gamma)
  }
  cont <- integrate(integrand, lower = a_mp, upper = b_mp,
                    rel.tol = 1e-10, abs.tol = 1e-14)$value

  # Atom at 0 with mass max(0, 1 - gamma): contributes lambda^2/(0+lambda)^2 * mass = mass

  atom_mass <- max(0, 1 - gamma)
  atom_contrib <- atom_mass  # lambda^2/lambda^2 = 1

  cont + atom_contrib
}

#' h(lambda, gamma) = m(-lambda) - lambda * m'(-lambda)
#' This equals int 1/(x+lambda)^2 d underline{F}(x) ... actually let's verify.
#' We have m'(-lambda) = d/d(lambda) m(-lambda) (derivative w.r.t. lambda).
#' h = m(-lambda) - lambda * m'(-lambda)
mp_h <- function(lambda, gamma) {
  m <- m_neg_lambda(lambda, gamma)
  m_prime <- m_neg_lambda_deriv(lambda, gamma)
  m - lambda * m_prime
}

weak_tau <- function(lambda, gamma, B, sigma_eps2) {
  m <- m_neg_lambda(lambda, gamma)
  h <- mp_h(lambda, gamma)

  tau_sq <- h * (sigma_eps2 + B * (lambda * m)^2)

  if (tau_sq < 0) {
    warning(sprintf("tau_sq = %.6e < 0 at lambda=%.4f, gamma=%.4f",
                    tau_sq, lambda, gamma))
    return(0)
  }

  sqrt(tau_sq)
}

#' Null: R_null = a
weak_null_theory <- function(a) {
  a
}

#' MoM: R_MoM = sqrt(B + sigma_eps^2)
weak_mom_theory <- function(B, sigma_eps2) {
  sqrt(B + sigma_eps2)
}

#' OLS (gamma < 1 only): R_OLS = sigma_eps / sqrt(1 - gamma)
weak_ols_theory <- function(gamma, sigma_eps) {
  if (gamma >= 1) return(NA_real_)
  sigma_eps / sqrt(1 - gamma)
}

#' Ridge: R_ridge(lambda) = max{ lambda*m(-lambda)*a + tau*sqrt(rho), tau }
#' With rho = 0: R_ridge(lambda) = max{ lambda*m(-lambda)*a, tau }
weak_ridge_theory <- function(lambda, gamma, a, B, sigma_eps2, rho = 0) {
  m <- m_neg_lambda(lambda, gamma)
  tau <- weak_tau(lambda, gamma, B, sigma_eps2)

  bias_term <- lambda * m * a + tau * sqrt(rho)
  R <- max(bias_term, tau)
  R
}

weak_ridge_opt <- function(lambda_grid, gamma, a, B, sigma_eps2, rho = 0) {
  R_vals <- sapply(lambda_grid, function(lam) {
    weak_ridge_theory(lam, gamma, a, B, sigma_eps2, rho)
  })
  idx <- which.min(R_vals)
  list(
    lambda_opt = lambda_grid[idx],
    R_opt      = R_vals[idx],
    R_all      = R_vals
  )
}
