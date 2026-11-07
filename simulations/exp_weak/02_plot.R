# =============================================================================
# Plotting: Weak-signal regime (Theorem: weak_signal)
#
# Two-panel figure:
#   Left:  gamma = 0.5 — Ridge path + Null, MoM, OLS horizontal lines
#   Right: gamma = 2.0 — Ridge path + Null, MoM horizontal lines
#
# x-axis: lambda (log scale)
# y-axis: ||beta_hat - beta_0||_inf / r_n  (normalized error)
# =============================================================================

library(ggplot2)
library(dplyr)
library(gridExtra)
source("config.R")
source("00_functions.R")

n <- N_VAL
a <- A_VAL

# --- Colors (same as strong-signal) ---
col_ridge <- "#0072B2"   # blue
col_null  <- "#009E73"   # green
col_mom   <- "#D55E00"   # vermillion
col_ols   <- "#CC79A7"   # pink

# --- Theme (same as strong-signal) ---
theme_paper <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, size = base_size + 1),
      axis.line  = element_line(linewidth = 0.3),
      axis.ticks = element_line(linewidth = 0.3),
      axis.title = element_text(size = base_size + 1),
      axis.text  = element_text(size = base_size, color = "black"),
      plot.margin = margin(5, 10, 5, 5)
    )
}

# --- Build one panel ---
build_panel <- function(gamma_val) {
  p   <- round(gamma_val * n)
  r_n <- sqrt(2 * log(p) / n)

  f <- sprintf("results/sim_gamma%.1f.rds", gamma_val)
  if (!file.exists(f)) stop("Missing: ", f)
  res <- readRDS(f)

  # --- Empirical ridge path ---
  ridge_mean <- colMeans(res$ridge_errs)
  ridge_sd   <- apply(res$ridge_errs, 2, sd)

  # --- Empirical horizontal estimators ---
  emp_null <- mean(res$null_errs)
  emp_mom  <- mean(res$mom_errs)
  emp_ols  <- if (gamma_val < 1) mean(res$ols_errs, na.rm = TRUE) else NA

  # --- Theory ---
  thy_ridge <- sapply(LAMBDA_GRID, function(lam) {
    weak_ridge_theory(lam, gamma_val, a, B, SIGMA_EPS^2, RHO)
  })
  thy_null <- weak_null_theory(a)
  thy_mom  <- weak_mom_theory(B, SIGMA_EPS^2)
  thy_ols  <- weak_ols_theory(gamma_val, SIGMA_EPS)

  # --- Build long-format data for unified legend ---
  # Empirical (solid lines)
  emp_long <- bind_rows(
    data.frame(lambda = LAMBDA_GRID, value = ridge_mean,
               lo = ridge_mean - 2 * ridge_sd, hi = ridge_mean + 2 * ridge_sd,
               Estimator = "Ridge"),
    data.frame(lambda = LAMBDA_GRID, value = emp_null,
               lo = emp_null - 2 * sd(res$null_errs),
               hi = emp_null + 2 * sd(res$null_errs),
               Estimator = "Null"),
    data.frame(lambda = LAMBDA_GRID, value = emp_mom,
               lo = emp_mom - 2 * sd(res$mom_errs),
               hi = emp_mom + 2 * sd(res$mom_errs),
               Estimator = "MoM")
  )
  if (gamma_val < 1 && !is.na(emp_ols)) {
    emp_long <- bind_rows(emp_long,
      data.frame(lambda = LAMBDA_GRID, value = emp_ols,
                 lo = emp_ols - 2 * sd(res$ols_errs, na.rm = TRUE),
                 hi = emp_ols + 2 * sd(res$ols_errs, na.rm = TRUE),
                 Estimator = "OLS"))
  }

  # Theory (dashed lines)
  thy_long <- bind_rows(
    data.frame(lambda = LAMBDA_GRID, value = thy_ridge, Estimator = "Ridge"),
    data.frame(lambda = LAMBDA_GRID, value = thy_null,  Estimator = "Null"),
    data.frame(lambda = LAMBDA_GRID, value = thy_mom,   Estimator = "MoM")
  )
  if (gamma_val < 1 && !is.na(thy_ols)) {
    thy_long <- bind_rows(thy_long,
      data.frame(lambda = LAMBDA_GRID, value = thy_ols, Estimator = "OLS"))
  }

  # Factor levels and maps (same order/style as strong-signal)
  est_levels <- if (gamma_val < 1) c("Null", "Ridge", "MoM", "OLS") else c("Null", "Ridge", "MoM")
  emp_long$Estimator <- factor(emp_long$Estimator, levels = est_levels)
  thy_long$Estimator <- factor(thy_long$Estimator, levels = est_levels)

  color_map <- c("Ridge" = col_ridge, "Null" = col_null, "MoM" = col_mom, "OLS" = col_ols)
  shape_map <- c("Ridge" = 16, "Null" = 17, "MoM" = 15, "OLS" = 18)

  pp <- ggplot() +
    geom_ribbon(data = emp_long, aes(x = lambda, ymin = lo, ymax = hi, fill = Estimator),
                alpha = 0.25) +
    geom_line(data = emp_long, aes(x = lambda, y = value, color = Estimator), linewidth = 0.6) +
    geom_point(data = emp_long, aes(x = lambda, y = value, color = Estimator, shape = Estimator),
               size = 2) +
    geom_line(data = thy_long, aes(x = lambda, y = value, color = Estimator),
              linetype = "dashed", linewidth = 0.7) +
    scale_color_manual(values = color_map) +
    scale_fill_manual(values = color_map, guide = "none") +
    scale_shape_manual(values = shape_map) +
    scale_x_log10() +
    guides(color = guide_legend(title = NULL), shape = guide_legend(title = NULL)) +
    labs(
      title = bquote(gamma == .(gamma_val)),
      x = expression(lambda),
      y = expression(group("||", hat(beta) - beta[0], "||")[infinity] / r[n])
    ) +
    coord_cartesian(ylim = c(0, NA)) +
    theme_paper() +
    theme(legend.position = c(0.05, 0.88), legend.justification = c(0, 1),
          legend.text = element_text(size = 9), legend.key.size = unit(0.5, "cm"))

  list(plot = pp, emp_long = emp_long, thy_long = thy_long)
}

# --- Build both panels ---
message("Building left panel (gamma = 0.5)...")
left  <- build_panel(0.5)
message("Building right panel (gamma = 2.0)...")
right <- build_panel(2.0)

# Shared y-axis
y_max <- max(c(left$emp_long$hi, left$thy_long$value,
               right$emp_long$hi, right$thy_long$value), na.rm = TRUE) * 1.05
left$plot  <- left$plot  + coord_cartesian(ylim = c(0, y_max))
right$plot <- right$plot + coord_cartesian(ylim = c(0, y_max))

dir.create("figures", showWarnings = FALSE)
p_combined <- arrangeGrob(left$plot, right$plot, ncol = 2, widths = c(1, 1))
ggsave("figures/weak_lambda_sweep_twopanel.pdf", p_combined,
       width = 10, height = 4.5, device = "pdf")
message("Saved: figures/weak_lambda_sweep_twopanel.pdf")

# --- LaTeX table: empirical / theory ratio ---
print_ratio_table <- function(gamma_val) {
  p   <- round(gamma_val * n)
  r_n <- sqrt(2 * log(p) / n)

  f <- sprintf("results/sim_gamma%.1f.rds", gamma_val)
  res <- readRDS(f)

  ridge_mean <- colMeans(res$ridge_errs)
  emp_null   <- mean(res$null_errs)
  emp_mom    <- mean(res$mom_errs)
  emp_ols    <- if (gamma_val < 1) mean(res$ols_errs, na.rm = TRUE) else NA

  thy_ridge <- sapply(LAMBDA_GRID, function(lam) {
    weak_ridge_theory(lam, gamma_val, a, B, SIGMA_EPS^2, RHO)
  })
  thy_null <- weak_null_theory(a)
  thy_mom  <- weak_mom_theory(B, SIGMA_EPS^2)
  thy_ols  <- weak_ols_theory(gamma_val, SIGMA_EPS)

  rows <- character(length(LAMBDA_GRID))
  for (i in seq_along(LAMBDA_GRID)) {
    lam_str <- if (i == 1) sprintf("%.1f", gamma_val) else ""
    r_ridge <- sprintf("%.3f", ridge_mean[i] / thy_ridge[i])
    r_null  <- sprintf("%.3f", emp_null / thy_null)
    r_mom   <- sprintf("%.3f", emp_mom / thy_mom)
    r_ols   <- if (!is.na(emp_ols) && !is.na(thy_ols)) sprintf("%.3f", emp_ols / thy_ols) else "---"
    rows[i] <- sprintf(" %s & %.3f & %s & %s & %s & %s \\\\",
                        lam_str, LAMBDA_GRID[i], r_ridge, r_null, r_mom, r_ols)
  }
  rows
}

cat("\n\\begin{table}[ht]\n\\centering\n")
cat(sprintf("\\caption{Empirical mean / theory ratio ($a = %.1f$, $B = %.1f$, $n_{\\mathrm{MC}} = %d$).}\n", a, B, N_MC))
cat("\\label{tab:weak-ratio}\n\\small\n")
cat("\\begin{tabular}{cc cccc}\n\\toprule\n")
cat("$\\gamma$ & $\\lambda$ & Ridge & Null & MoM & OLS \\\\\n\\midrule\n")
cat(paste(print_ratio_table(0.5), collapse = "\n"), "\n")
cat("\\addlinespace\n")
cat(paste(print_ratio_table(2.0), collapse = "\n"), "\n")
cat("\\bottomrule\n\\end{tabular}\n\\end{table}\n")

message("Done.")
