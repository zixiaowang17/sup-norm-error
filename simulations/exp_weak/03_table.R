# =============================================================================
# Generate LaTeX table: Empirical mean / theory ratio for weak-signal experiment
# =============================================================================

source("config.R")
source("00_functions.R")

n <- N_VAL
a <- A_VAL

# Select ~5 lambda values from the grid (evenly spaced in log scale)
idx_sel <- round(seq(1, length(LAMBDA_GRID), length.out = 5))
lam_sel <- LAMBDA_GRID[idx_sel]

rows <- list()

for (gamma_val in GAMMA_VALS) {
  p   <- round(gamma_val * n)
  r_n <- sqrt(2 * log(p) / n)

  f <- sprintf("results/sim_gamma%.1f.rds", gamma_val)
  if (!file.exists(f)) stop("Missing: ", f)
  res <- readRDS(f)

  for (i in seq_along(lam_sel)) {
    lam <- lam_sel[i]
    j   <- idx_sel[i]  # column index in ridge_errs

    # Empirical means
    emp_ridge <- mean(res$ridge_errs[, j])
    emp_null  <- mean(res$null_errs)
    emp_mom   <- mean(res$mom_errs)
    emp_ols   <- if (gamma_val < 1) mean(res$ols_errs, na.rm = TRUE) else NA

    # Theory
    thy_ridge <- weak_ridge_theory(lam, gamma_val, a, B, SIGMA_EPS^2, RHO)
    thy_null  <- weak_null_theory(a)
    thy_mom   <- weak_mom_theory(B, SIGMA_EPS^2)
    thy_ols   <- weak_ols_theory(gamma_val, SIGMA_EPS)

    # Ratios
    r_ridge <- emp_ridge / thy_ridge
    r_null  <- emp_null  / thy_null
    r_mom   <- emp_mom   / thy_mom
    r_ols   <- if (!is.na(emp_ols) && !is.na(thy_ols)) emp_ols / thy_ols else NA

    rows[[length(rows) + 1]] <- list(
      gamma = gamma_val, lam = lam,
      r_ridge = r_ridge, r_null = r_null, r_mom = r_mom, r_ols = r_ols,
      first = (i == 1)
    )
  }
}

# --- Print LaTeX ---
cat("\\begin{table}[ht]\n")
cat("  \\centering\n")
cat(sprintf("  \\caption{Empirical mean / theory ratio ($a = %.1f$, $B = %.1f$, $n = %d$, $n_{\\mathrm{MC}} = %d$).}\n",
            A_VAL, B, N_VAL, N_MC))
cat("  \\label{tab:weak-ratio}\n")
cat("  \\small\n")
cat("  \\begin{tabular}{cc cccc}\n")
cat("  \\toprule\n")
cat("  $\\gamma$ & $\\lambda$ & Ridge & Null & MoM & OLS \\\\\n")
cat("  \\midrule\n")

prev_gamma <- NULL
for (r in rows) {
  if (!is.null(prev_gamma) && r$gamma != prev_gamma) {
    cat("  \\addlinespace\n")
  }
  gam_str <- if (r$first) sprintf("%.1f", r$gamma) else ""
  ols_str <- if (is.na(r$r_ols)) "---" else sprintf("%.3f", r$r_ols)
  cat(sprintf("  %s & %.3f & %.3f & %.3f & %.3f & %s \\\\\n",
              gam_str, r$lam, r$r_ridge, r$r_null, r$r_mom, ols_str))
  prev_gamma <- r$gamma
}

cat("  \\bottomrule\n")
cat("  \\end{tabular}\n")
cat("\\end{table}\n")

message("Done.")
