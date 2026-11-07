# =============================================================================
# Strong-Signal Simulation (Theorem 3.2)
#
# Usage: Rscript 01_run_sim.R <gamma>
#   gamma: 0.5 or 2.0
#
# Sweeps alpha grid with fixed lambda. One job per gamma.
# All errors are RAW ||beta_hat - beta_0||_inf.
# =============================================================================

library(future.apply)
source("config.R")
source("00_functions.R")

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 1)
gamma_val <- as.numeric(args[1])

n <- N_VAL
p <- round(gamma_val * n)

message(sprintf("=== gamma = %.1f (p = %d) ===", gamma_val, p))

simulate_one_rep <- function(rep_id, beta, n, p) {
  X   <- matrix(rnorm(n * p), n, p)
  eps <- rnorm(n, sd = SIGMA_EPS)
  y   <- X %*% beta + eps

  # MoM
  beta_MoM <- as.vector(crossprod(X, y)) / n
  err_MoM  <- max(abs(beta_MoM - beta))

  # Null
  err_null <- max(abs(beta))

  # OLS (gamma < 1 only)
  if (p < n) {
    beta_OLS <- as.vector(solve(crossprod(X), crossprod(X, y)))
    err_OLS  <- max(abs(beta_OLS - beta))
  } else {
    err_OLS  <- NA
  }

  # Ridgeless / min-norm (gamma > 1 only)
  if (p > n) {
    beta_rl <- as.vector(crossprod(X, solve(tcrossprod(X), y)))
    err_rl  <- max(abs(beta_rl - beta))
  } else {
    err_rl  <- NA
  }

  # Ridge with fixed lambda
  if (p <= n) {
    XtX_n    <- crossprod(X) / n
    beta_hat <- as.vector(solve(XtX_n + LAMBDA_FIXED * diag(p),
                                crossprod(X, y) / n))
  } else {
    XXt_n    <- tcrossprod(X) / n
    beta_hat <- as.vector(crossprod(X, solve(XXt_n + LAMBDA_FIXED * diag(n),
                                            y)) / n)
  }
  err_ridge <- max(abs(beta_hat - beta))

  data.frame(
    err_ridge     = err_ridge,
    err_MoM       = err_MoM,
    err_OLS       = err_OLS,
    err_null      = err_null,
    err_ridgeless = err_rl,
    rep           = rep_id
  )
}

# --- Run ---
n_workers <- min(N_WORKERS, max(1, parallel::detectCores() - 1))
plan(multisession, workers = n_workers)
dir.create("results", showWarnings = FALSE)

for (alpha in ALPHA_GRID) {
  bp   <- make_beta(n, p, alpha, S0)
  beta <- bp$beta

  message(sprintf("  alpha=%.2f: r_n=%.4f, noise=%.4f, ratio=%.2f, B_n=%.4f",
                  alpha, bp$r_n, bp$noise_scale, bp$r_n / bp$noise_scale, bp$B_n))

  all_results <- future_lapply(1:N_MC, function(i) {
    simulate_one_rep(i, beta, n, p)
  }, future.seed = TRUE)

  final_df       <- do.call(rbind, all_results)
  final_df$alpha <- alpha
  final_df$gamma <- gamma_val

  outfile <- sprintf("results/sim_gamma%.1f_alpha%.2f.rds", gamma_val, alpha)
  saveRDS(final_df, outfile)
  message(sprintf("    Saved: %s (%d rows)", outfile, nrow(final_df)))
}

message("Done.")
