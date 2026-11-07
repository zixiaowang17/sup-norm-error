# =============================================================================
# Plotting: Double descent in the weak-signal regime
#
# x-axis: gamma = p/n
# y-axis: ||beta_hat - beta_0||_inf / r_n  (normalized error)
# =============================================================================

library(ggplot2)
library(dplyr)
source("config_dd.R")
source("00_functions.R")

n <- N_VAL
a <- A_VAL
lam <- LAMBDA_DD

# --- Load results ---
emp_list <- list()
for (gam in GAMMA_GRID) {
  f <- sprintf("results/dd_gamma%.2f.rds", gam)
  if (!file.exists(f)) {
    message("Missing: ", f)
    next
  }
  res <- readRDS(f)
  emp_list[[length(emp_list) + 1]] <- data.frame(
    gamma = gam,
    mean  = mean(res$errs),
    sd    = sd(res$errs)
  )
}
emp_df <- bind_rows(emp_list)
emp_df$lo <- emp_df$mean - 2 * emp_df$sd
emp_df$hi <- emp_df$mean + 2 * emp_df$sd

# --- Theory: ridgeless limiting constant ---
# Under-parameterized (gamma < 1): sigma_eps / sqrt(1 - gamma)
# Over-parameterized (gamma > 1): weak_ridge_theory with lambda -> 0
gamma_thy_under <- seq(0.05, 0.99, by = 0.01)
gamma_thy_over  <- seq(1.01, 3.0,  by = 0.01)

thy_under <- data.frame(
  gamma = gamma_thy_under,
  value = SIGMA_EPS / sqrt(1 - gamma_thy_under)
)

thy_over <- data.frame(
  gamma = gamma_thy_over,
  value = sapply(gamma_thy_over, function(g) {
    weak_ridge_theory(LAMBDA_DD, g, a, B, SIGMA_EPS^2, RHO)
  })
)

thy_df <- bind_rows(thy_under, thy_over)

# --- Colors ---
col_emp <- "#0072B2"
col_thy <- "#D55E00"

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

# --- Plot ---
pp <- ggplot() +
  geom_ribbon(data = emp_df,
              aes(x = gamma, ymin = lo, ymax = hi),
              fill = col_emp, alpha = 0.2) +
  geom_line(data = emp_df,
            aes(x = gamma, y = mean, color = "Empirical"),
            linewidth = 0.7) +
  geom_point(data = emp_df,
             aes(x = gamma, y = mean, color = "Empirical"),
             size = 1.2) +
  geom_line(data = thy_df,
            aes(x = gamma, y = value, color = "Theory"),
            linetype = "dashed", linewidth = 0.7) +
  geom_vline(xintercept = 1, linetype = "dotted", color = "grey40") +
  scale_color_manual(values = c("Empirical" = col_emp, "Theory" = col_thy)) +
  labs(
    x = expression(gamma == p / n),
    y = expression(group("||", hat(beta)[lambda] - beta[0], "||")[infinity] / r[n])
  ) +
  coord_cartesian(ylim = c(0, NA)) +
  theme_paper() +
  theme(legend.position = c(0.80, 0.85),
        legend.title = element_blank(),
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.5, "cm"))

dir.create("figures", showWarnings = FALSE)
ggsave("figures/weak_double_descent.pdf", pp,
       width = 6, height = 4, device = "pdf")
message("Saved: figures/weak_double_descent.pdf")

# --- LaTeX table: empirical mean / theory ratio ---
cat("\n\\begin{table}[ht]\n\\centering\n")
cat(sprintf("\\caption{Empirical mean / theory ratio for nearly ridgeless regression ($\\lambda = 10^{-5}$, $a = %g$, $B = %g$, $n_{\\mathrm{MC}} = %d$).}\n",
            a, B, N_MC))
cat("\\label{tab:dd-ratio}\n\\small\n")
cat("\\begin{tabular}{ccccc}\n\\toprule\n")
cat("$\\gamma$ & Empirical & Theory & Ratio \\\\\n\\midrule\n")

for (i in seq_len(nrow(emp_df))) {
  gam <- emp_df$gamma[i]
  emp_val <- emp_df$mean[i]
  if (gam < 1) {
    thy_val <- SIGMA_EPS / sqrt(1 - gam)
  } else {
    thy_val <- weak_ridge_theory(lam, gam, a, B, SIGMA_EPS^2, RHO)
  }
  ratio <- emp_val / thy_val
  cat(sprintf("%.2f & %.3f & %.3f & %.3f \\\\\n", gam, emp_val, thy_val, ratio))
}

cat("\\bottomrule\n\\end{tabular}\n\\end{table}\n")
