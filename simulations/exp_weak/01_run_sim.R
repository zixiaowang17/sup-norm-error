library(future.apply)
source("config.R")
source("00_functions.R")

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 1)
gamma_val <- as.numeric(args[1])

n <- N_VAL
p <- round(gamma_val * n)
r_n <- sqrt(2 * log(p) / n)
a <- A_VAL

bp   <- make_beta(n, p, a)
beta <- bp$beta

# Verification
stopifnot(abs(sum(beta^2) - B) < 1e-10)
stopifnot(abs(beta[1] / r_n - a) < 1e-10)

message(sprintf("=== gamma=%.1f (p=%d, r_n=%.4f, a=%g, B=%.2f) ===",
                gamma_val, p, r_n, a, B))
message(sprintf("  beta_inf/r_n=%.4f, ||beta||_2^2=%.6f", bp$beta_inf / r_n, bp$B))

# --- Simulation function ---
simulate_one_rep <- function(rep_id, beta, n, p, r_n) {
  X   <- matrix(rnorm(n * p), n, p)
  eps <- rnorm(n, sd = SIGMA_EPS)
  y   <- X %*% beta + eps

  Xty_n <- as.vector(crossprod(X, y)) / n

  # MoM
  err_MoM  <- max(abs(Xty_n - beta)) / r_n

  # Null
  err_null <- max(abs(beta)) / r_n

  # OLS (gamma < 1 only)
  if (p < n) {
    beta_OLS <- as.vector(solve(crossprod(X), crossprod(X, y)))
    err_OLS  <- max(abs(beta_OLS - beta)) / r_n
  } else {
    err_OLS  <- NA
  }

  # Ridge over full lambda grid
  if (p <= n) {
    XtX_n <- crossprod(X) / n
    ridge_errs <- sapply(LAMBDA_GRID, function(lam) {
      beta_hat <- as.vector(solve(XtX_n + lam * diag(p), Xty_n))
      max(abs(beta_hat - beta)) / r_n
    })
  } else {
    XXt_n <- tcrossprod(X) / n
    ridge_errs <- sapply(LAMBDA_GRID, function(lam) {
      beta_hat <- as.vector(crossprod(X, solve(XXt_n + lam * diag(n), y)) / n)
      max(abs(beta_hat - beta)) / r_n
    })
  }

  c(ridge_errs, err_MoM, err_null, err_OLS)
}

# --- Run ---
n_workers <- min(N_WORKERS, max(1, parallel::detectCores() - 1))
plan(multisession, workers = n_workers)
dir.create("results", showWarnings = FALSE)

all_results <- future_lapply(1:N_MC, function(i) {
  simulate_one_rep(i, beta, n, p, r_n)
}, future.seed = TRUE)

# Stack into matrix: N_MC x (n_lambda + 3)
mat <- do.call(rbind, all_results)
n_lam <- length(LAMBDA_GRID)

result <- list(
  ridge_errs = mat[, 1:n_lam],          # N_MC x n_lambda
  mom_errs   = mat[, n_lam + 1],        # N_MC
  null_errs  = mat[, n_lam + 2],        # N_MC (all same)
  ols_errs   = mat[, n_lam + 3],        # N_MC (NA for gamma > 1)
  lambda_grid = LAMBDA_GRID,
  gamma      = gamma_val,
  a          = a,
  n          = n,
  p          = p
)

outfile <- sprintf("results/sim_gamma%.1f.rds", gamma_val)
saveRDS(result, outfile)
message(sprintf("Saved: %s (%d reps, %d lambdas)", outfile, N_MC, n_lam))
message("Done.")
