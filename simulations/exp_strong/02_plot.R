# =============================================================================
# Plotting: Strong-signal regime (Theorem 3.2)
#
# Two-panel figure:
#   Left:  gamma = 0.5 — Ridge, Null, MoM, OLS
#   Right: gamma = 2.0 — Ridge, Null, MoM, Ridgeless
#
# x-axis: alpha
# y-axis: raw ||beta_hat - beta_0||_inf
# Solid lines + points: MC means
# Dashed lines: theory
# Shaded bands: mean ± 2 sd
# =============================================================================

library(ggplot2)
library(dplyr)
library(gridExtra)
source("config.R")
source("00_functions.R")

n <- N_VAL

# --- Colors ---
col_ridge <- "#0072B2"   # blue
col_null  <- "#009E73"   # green
col_mom   <- "#D55E00"   # vermillion
col_ols   <- "#CC79A7"   # pink
col_rl    <- "#E69F00"   # gold (ridgeless)

# --- Theme ---
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
  p <- round(gamma_val * n)

  # --- Load and aggregate ---
  all_alpha <- list()
  for (alpha in ALPHA_GRID) {
    f <- sprintf("results/sim_gamma%.1f_alpha%.2f.rds", gamma_val, alpha)
    if (!file.exists(f)) { warning("Missing: ", f); next }
    df <- readRDS(f)
    one_rep <- df %>% distinct(rep, .keep_all = TRUE)
    row <- data.frame(
      alpha       = alpha,
      ridge_mean  = mean(df$err_ridge),
      ridge_lo    = mean(df$err_ridge) - 2 * sd(df$err_ridge),
      ridge_hi    = mean(df$err_ridge) + 2 * sd(df$err_ridge),
      null_mean   = mean(one_rep$err_null),
      mom_mean    = mean(one_rep$err_MoM),
      mom_lo      = mean(one_rep$err_MoM) - 2 * sd(one_rep$err_MoM),
      mom_hi      = mean(one_rep$err_MoM) + 2 * sd(one_rep$err_MoM)
    )
    if (gamma_val < 1) {
      row$ols_mean <- mean(one_rep$err_OLS, na.rm = TRUE)
      row$ols_lo   <- mean(one_rep$err_OLS, na.rm = TRUE) - 2 * sd(one_rep$err_OLS, na.rm = TRUE)
      row$ols_hi   <- mean(one_rep$err_OLS, na.rm = TRUE) + 2 * sd(one_rep$err_OLS, na.rm = TRUE)
    }
    all_alpha[[length(all_alpha) + 1]] <- row
  }
  agg <- do.call(rbind, all_alpha) %>% arrange(alpha)

  # --- Theory ---
  alpha_fine <- seq(min(ALPHA_GRID), max(ALPHA_GRID), length.out = 100)
  theory <- data.frame(alpha = alpha_fine) %>%
    rowwise() %>%
    mutate(
      bp       = list(make_beta(n, p, alpha, S0)),
      r_n_val  = bp$r_n,
      noise_sc = bp$noise_scale,
      B_n      = bp$B_n,
      ridge    = strong_ridge_theory(LAMBDA_FIXED, gamma_val, C_STAR, r_n_val),
      null_val = strong_null_theory(C_STAR, r_n_val),
      mom      = strong_mom_theory(noise_sc, B_n, SIGMA_EPS^2)
    ) %>%
    ungroup()

  if (gamma_val < 1) {
    theory$ols <- strong_ols_theory(theory$noise_sc, gamma_val, SIGMA_EPS)
  }
  # --- Reshape for legend-friendly plotting ---
  # Empirical long format
  emp_long <- bind_rows(
    agg %>% transmute(alpha, value = ridge_mean, lo = ridge_lo, hi = ridge_hi, Estimator = "Ridge"),
    agg %>% transmute(alpha, value = null_mean,  lo = null_mean, hi = null_mean, Estimator = "Null"),
    agg %>% transmute(alpha, value = mom_mean,   lo = mom_lo,   hi = mom_hi,   Estimator = "MoM")
  )
  if (gamma_val < 1) {
    emp_long <- bind_rows(emp_long,
      agg %>% transmute(alpha, value = ols_mean, lo = ols_lo, hi = ols_hi, Estimator = "OLS"))
  }

  # Theory long format
  thy_long <- bind_rows(
    theory %>% transmute(alpha, value = ridge,    Estimator = "Ridge"),
    theory %>% transmute(alpha, value = null_val, Estimator = "Null"),
    theory %>% transmute(alpha, value = mom,      Estimator = "MoM")
  )
  if (gamma_val < 1) {
    thy_long <- bind_rows(thy_long,
      theory %>% transmute(alpha, value = ols, Estimator = "OLS"))
  }

  # Color and shape maps
  est_levels <- if (gamma_val < 1) c("Null", "Ridge", "MoM", "OLS") else c("Null", "Ridge", "MoM")
  emp_long$Estimator <- factor(emp_long$Estimator, levels = est_levels)
  thy_long$Estimator <- factor(thy_long$Estimator, levels = est_levels)

  color_map <- c("Ridge" = col_ridge, "Null" = col_null, "MoM" = col_mom, "OLS" = col_ols)
  shape_map <- c("Ridge" = 16, "Null" = 17, "MoM" = 15, "OLS" = 18)

  pp <- ggplot() +
    geom_ribbon(data = emp_long, aes(x = alpha, ymin = lo, ymax = hi, fill = Estimator), alpha = 0.15) +
    geom_line(data = emp_long, aes(x = alpha, y = value, color = Estimator), linewidth = 0.6) +
    geom_point(data = emp_long, aes(x = alpha, y = value, color = Estimator, shape = Estimator), size = 2) +
    geom_line(data = thy_long, aes(x = alpha, y = value, color = Estimator),
              linetype = "dashed", linewidth = 0.7) +
    scale_color_manual(values = color_map) +
    scale_fill_manual(values = color_map, guide = "none") +
    scale_shape_manual(values = shape_map) +
    guides(color = guide_legend(title = NULL), shape = guide_legend(title = NULL)) +
    labs(
      title = bquote(gamma == .(gamma_val)),
      x = expression(alpha),
      y = expression(group("||", hat(beta) - beta[0], "||")[infinity])
    ) +
    coord_cartesian(ylim = c(0, NA)) +
    theme_paper() +
    theme(legend.position = "top", legend.justification = "left")

  list(plot = pp, agg = agg, theory = theory)
}

# --- Build both panels with shared y-axis ---
left  <- build_panel(0.5)
right <- build_panel(2.0)

# Compute shared y-axis range
y_max <- max(c(left$agg$ridge_hi, left$agg$null_mean,
               right$agg$ridge_hi, right$agg$null_mean), na.rm = TRUE) * 1.1
left$plot  <- left$plot  + coord_cartesian(ylim = c(0, y_max))
right$plot <- right$plot + coord_cartesian(ylim = c(0, y_max))

p_combined <- arrangeGrob(left$plot, right$plot, ncol = 2, widths = c(1, 1))
ggsave("strong_alpha_sweep.pdf", p_combined, width = 10, height = 4.5, device = "pdf")
message("Saved: strong_alpha_sweep.pdf")

# =====================================================================
# Ratio table (both gammas combined)
# =====================================================================
ratio_rows <- list()
for (gamma_val in GAMMA_VALS) {
  p <- round(gamma_val * n)
  for (alpha in ALPHA_GRID) {
    f <- sprintf("results/sim_gamma%.1f_alpha%.2f.rds", gamma_val, alpha)
    if (!file.exists(f)) next
    df <- readRDS(f)
    one_rep <- df %>% distinct(rep, .keep_all = TRUE)
    bp <- make_beta(n, p, alpha, S0)

    row <- data.frame(
      gamma       = gamma_val,
      alpha       = alpha,
      ratio_rn    = round(bp$r_n / bp$noise_scale, 2),
      Ridge_ratio = round(mean(df$err_ridge) / strong_ridge_theory(LAMBDA_FIXED, gamma_val, C_STAR, bp$r_n), 3),
      Null_ratio  = round(mean(one_rep$err_null) / strong_null_theory(C_STAR, bp$r_n), 3),
      MoM_ratio   = round(mean(one_rep$err_MoM) / strong_mom_theory(bp$noise_scale, bp$B_n, SIGMA_EPS^2), 3)
    )
    if (gamma_val < 1) {
      row$OLS_ratio <- round(mean(one_rep$err_OLS, na.rm = TRUE) /
                               strong_ols_theory(bp$noise_scale, gamma_val, SIGMA_EPS), 3)
    }
    ratio_rows[[length(ratio_rows) + 1]] <- row
  }
}
ratio_tab <- dplyr::bind_rows(ratio_rows)

cat("\n=== Empirical Mean / Theory Ratio ===\n")
print(as.data.frame(ratio_tab), row.names = FALSE)

# --- Save LaTeX ---
tex <- character()
tex <- c(tex, "\\begin{table}[ht]")
tex <- c(tex, "\\centering")
tex <- c(tex, sprintf("\\caption{Empirical mean / theory ratio ($\\lambda = %g$, $s_0 = %d$, $n_{\\mathrm{MC}} = %d$).}", LAMBDA_FIXED, S0, N_MC))
tex <- c(tex, "\\label{tab:strong-ratio}")
tex <- c(tex, "\\small")
tex <- c(tex, "\\begin{tabular}{cc c ccc c}")
tex <- c(tex, "\\toprule")
tex <- c(tex, "$\\gamma$ & $\\alpha$ & $n^\\alpha$ & Ridge & Null & MoM & OLS \\\\")
tex <- c(tex, "\\midrule")

prev_gamma <- -1
for (i in seq_len(nrow(ratio_tab))) {
  r <- ratio_tab[i, ]
  if (r$gamma != prev_gamma && prev_gamma > 0) tex <- c(tex, "\\addlinespace")
  g_str <- if (r$gamma != prev_gamma) sprintf("%.1f", r$gamma) else ""
  last_col <- if ("OLS_ratio" %in% names(r) && !is.na(r$OLS_ratio)) sprintf("%.3f", r$OLS_ratio) else "---"
  tex <- c(tex, sprintf("%s & %.2f & %.2f & %.3f & %.3f & %.3f & %s \\\\",
                        g_str, r$alpha, r$ratio_rn,
                        r$Ridge_ratio, r$Null_ratio, r$MoM_ratio, last_col))
  prev_gamma <- r$gamma
}

tex <- c(tex, "\\bottomrule")
tex <- c(tex, "\\end{tabular}")
tex <- c(tex, "\\end{table}")
writeLines(tex, "results/strong_ratio_table.tex")
message("Saved: results/strong_ratio_table.tex")

message("Done.")
