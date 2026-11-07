library(future.apply)
source("config_dd.R")
source("00_functions.R")

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 1)
gamma_val <- as.numeric(args[1])

n <- N_VAL
p <- floor(gamma_val * n)
r_n <- sqrt(2 * log(p) / n)
a <- A_VAL
lam <- LAMBDA_DD

bp   <- make_beta(n, p, a)
beta <- bp$beta

stopifnot(abs(sum(beta^2) - B) < 1e-10)
stopifnot(abs(beta[1] / r_n - a) < 1e-10)

message(sprintf("=== gamma=%.2f (p=%d, r_n=%.4f, a=%g, B=%.2f, lambda=%.1e) ===",
                gamma_val, p, r_n, a, B, lam))

simulate_one_rep <- function(rep_id, beta, n, p, r_n, lam) {
  X   <- matrix(rnorm(n * p), n, p)
  eps <- rnorm(n, sd = SIGMA_EPS)
  y   <- X %*% beta + eps

  Xty_n <- as.vector(crossprod(X, y)) / n

  if (p <= n) {
    XtX_n <- crossprod(X) / n
    beta_hat <- as.vector(solve(XtX_n + lam * diag(p), Xty_n))
  } else {
    XXt_n <- tcrossprod(X) / n
    beta_hat <- as.vector(crossprod(X, solve(XXt_n + lam * diag(n), y)) / n)
  }

  max(abs(beta_hat - beta)) / r_n
}

n_workers <- min(N_WORKERS, max(1, parallel::detectCores() - 1))
plan(multisession, workers = n_workers)
dir.create("results", showWarnings = FALSE)

errs <- future_sapply(1:N_MC, function(i) {
  simulate_one_rep(i, beta, n, p, r_n, lam)
}, future.seed = TRUE)

result <- list(
  errs       = errs,
  gamma      = gamma_val,
  a          = a,
  n          = n,
  p          = p,
  lambda     = lam
)

outfile <- sprintf("results/dd_gamma%.2f.rds", gamma_val)
saveRDS(result, outfile)
message(sprintf("Saved: %s (%d reps, mean=%.4f, sd=%.4f)",
                outfile, N_MC, mean(errs), sd(errs)))
