library(dplyr)
library(ggplot2)
library(htmlwidgets)
library(stargazer)
library(pagedown)
library(haven)

# load quantile on quantiles functions (plot_qq_surface, plot_qq_slices, etc.)
source("R/00_2_functions.R")

#################################################################
## Shared parameters (must match the three bootstrap scripts)
#################################################################

z_min_max <- c(-7, 9)
y_min_max <- c(-7, 9)

taus <- round(seq(0.1, 0.9, by = 0.1), 1)

# rescaling to one st.dev.
shocks <- read_dta("data/czone_exposure_by_period_v5_gh.dta")
sd_d <- sd(shocks$d_tradeusch_p1_2000_2012, na.rm = TRUE)
#################################################################
## 1. Linear IV approach 
#################################################################

bootstrap_results <- readRDS(
  "plots/iv_qoq/linear_approach/bootstrap_data_linear.rds"
)
bootstrap_results_pct <- (exp(bootstrap_results * sd_d) - 1) * 100

coef_hat <- apply(bootstrap_results_pct, c(2, 3), mean, na.rm = TRUE)
se_hat   <- apply(bootstrap_results_pct, c(2, 3), sd, na.rm = TRUE)
ci_low   <- apply(bootstrap_results_pct, c(2, 3), quantile, probs = 0.025, na.rm = TRUE)
ci_high  <- apply(bootstrap_results_pct, c(2, 3), quantile, probs = 0.975, na.rm = TRUE)
mean_se <- mean(se_hat, na.rm = TRUE)
median_se <- median(se_hat, na.rm = TRUE)

result_linear <- list(
  coef = coef_hat,
  se   = se_hat,
  ci_low = ci_low,
  ci_high = ci_high
)

# ---- surface plot ----
boot_surface <- plot_qq_surface(
  coef_hat,
  taus = taus,
  tplot = 1,
  mycolors = "",
  z_label = "IV Coefficient",
  group = "CZ",
  zlim = z_min_max
)

htmlwidgets::saveWidget(
  boot_surface,
  file = "plots/iv_qoq/linear_approach/surface_boot_linear.html",
  selfcontained = TRUE
)

webshot2::webshot(
  "plots/iv_qoq/linear_approach/surface_boot_linear.html",
  file = "plots/iv_qoq/linear_approach/surface_boot_linear.png",
  vwidth = 1050, vheight = 900, zoom = 3
)

# ---- slices plot ----
boot_slices <- plot_qq_slices(
  coef_hat,
  taus = taus,
  z_label = "IV Coefficient",
  group = "CZ",
  ylim = y_min_max
) +
  ggtitle("China Shock Impact - Linear Approach")

ggsave(
  filename = "plots/iv_qoq/linear_approach/slices_boot_linear.png",
  plot = boot_slices,
  width = 8,
  height = 5,
  dpi = 300,
  bg = "white"
)

# ---- significance table ----
sig <- (ci_low > 0) | (ci_high < 0)
stars <- ifelse(sig, "**", "")
coef_star <- matrix(
  paste0(round(coef_hat, 3), stars),
  nrow = nrow(coef_hat),
  ncol = ncol(coef_hat),
  dimnames = dimnames(coef_hat)
)

se_formatted <- matrix(
  paste0("(", sprintf("%.3f", se_hat), ")"),
  nrow = nrow(se_hat)
)

final_table <- matrix("",
                      nrow = nrow(coef_hat) * 2,
                      ncol = ncol(coef_hat),
                      dimnames = list(
                        rep(rownames(coef_hat), each = 2),
                        colnames(coef_hat)
                      ))

final_table[seq(1, nrow(final_table), by = 2), ] <- coef_star
final_table[seq(2, nrow(final_table), by = 2), ] <- se_formatted

make_bootstrap_table(
  coef_hat,
  se_hat,
  ci_low,
  ci_high,
  file = "plots/iv_qoq/linear_approach/bootstrap_results_linear.tex"
)


#################################################################
## 2. IVQR / Chernozhukov-Hansen approach 
#################################################################

bootstrap_results <- readRDS("plots/iv_qoq/IVQR_approach/bootstrap_data_IVQR.rds")
bootstrap_results_pct <- (exp(bootstrap_results * sd_d) - 1) * 100


coef_hat <- apply(bootstrap_results_pct, c(2, 3), mean, na.rm = TRUE)
se_hat   <- apply(bootstrap_results_pct, c(2, 3), sd, na.rm = TRUE)
ci_low   <- apply(bootstrap_results_pct, c(2, 3), quantile, probs = 0.025, na.rm = TRUE)
ci_high  <- apply(bootstrap_results_pct, c(2, 3), quantile, probs = 0.975, na.rm = TRUE)

result_ivqr <- list(
  coef = coef_hat,
  se   = se_hat,
  ci_low = ci_low,
  ci_high = ci_high
)

# ---- surface plot ----
boot_surface <- plot_qq_surface(
  coef_hat,
  taus = taus,
  tplot = 1,
  mycolors = "",
  z_label = "IV Coefficient",
  group = "CZ",
  zlim = z_min_max
)

htmlwidgets::saveWidget(
  boot_surface,
  file = "plots/iv_qoq/IVQR_approach/surface_boot_IVQR.html",
  selfcontained = TRUE
)

webshot2::webshot(
  "plots/iv_qoq/IVQR_approach/surface_boot_IVQR.html",
  file = "plots/iv_qoq/IVQR_approach/surface_boot_IVQR.png",
  vwidth = 1200, vheight = 1000, zoom = 3
)

# ---- slices plot ----
boot_slices <- plot_qq_slices(
  coef_hat,
  taus = taus,
  z_label = "IV Coefficient",
  group = "CZ",
  ylim = y_min_max
) +
  ggtitle("China Shock Impact - IVQR Approach")

ggsave(
  filename = "plots/iv_qoq/IVQR_approach/slices_boot_IVQR.png",
  plot = boot_slices,
  width = 8,
  height = 5,
  dpi = 300,
  bg = "white"
)

# ---- significance table ----
sig <- (ci_low > 0) | (ci_high < 0)
stars <- ifelse(sig, "**", "")
coef_star <- matrix(
  paste0(round(coef_hat, 3), stars),
  nrow = nrow(coef_hat),
  ncol = ncol(coef_hat),
  dimnames = dimnames(coef_hat)
)

se_formatted <- matrix(
  paste0("(", sprintf("%.3f", se_hat), ")"),
  nrow = nrow(se_hat)
)

final_table <- matrix("",
                      nrow = nrow(coef_hat) * 2,
                      ncol = ncol(coef_hat),
                      dimnames = list(
                        rep(rownames(coef_hat), each = 2),
                        colnames(coef_hat)
                      ))

final_table[seq(1, nrow(final_table), by = 2), ] <- coef_star
final_table[seq(2, nrow(final_table), by = 2), ] <- se_formatted

make_bootstrap_table(
  coef_hat,
  se_hat,
  ci_low,
  ci_high,
  file = "plots/iv_qoq/IVQR_approach/bootstrap_results_IVQR.tex"
)


#################################################################
## 3. Control function approach (from 04_03_IV_CF_boot.R)
#################################################################


bootstrap_results <- readRDS(
  "plots/iv_qoq/control_function_approach/bootstrap_data_CF.rds"
)
bootstrap_results_pct <- (exp(bootstrap_results * sd_d) - 1) * 100


coef_hat <- apply(bootstrap_results_pct, c(2, 3), mean, na.rm = TRUE)
se_hat   <- apply(bootstrap_results_pct, c(2, 3), sd, na.rm = TRUE)
ci_low   <- apply(bootstrap_results_pct, c(2, 3), quantile, probs = 0.025, na.rm = TRUE)
ci_high  <- apply(bootstrap_results_pct, c(2, 3), quantile, probs = 0.975, na.rm = TRUE)

result_cf <- list(
  coef = coef_hat,
  se   = se_hat,
  ci_low = ci_low,
  ci_high = ci_high
)

# ---- surface plot ----
boot_surface <- plot_qq_surface(
  coef_hat,
  taus = taus,
  tplot = 1,
  mycolors = "",
  z_label = "IV Coefficient",
  group = "CZ",
  zlim = z_min_max
)

htmlwidgets::saveWidget(
  boot_surface,
  file = "plots/iv_qoq/control_function_approach/surface_boot_CF.html",
  selfcontained = TRUE
)

webshot2::webshot(
  "plots/iv_qoq/control_function_approach/surface_boot_CF.html",
  file = "plots/iv_qoq/control_function_approach/surface_boot_CF.png",
  vwidth = 1050, vheight = 900, zoom = 3
)

# ---- slices plot ----
boot_slices <- plot_qq_slices(
  coef_hat,
  taus = taus,
  z_label = "IV Coefficient",
  group = "CZ",
  ylim = y_min_max
) +
  ggtitle("China Shock Impact - Control Function Approach")

ggsave(
  filename = "plots/iv_qoq/control_function_approach/slices_boot_CF.png",
  plot = boot_slices,
  width = 8,
  height = 5,
  dpi = 300,
  bg = "white"
)

# ---- significance table ----
sig <- (ci_low > 0) | (ci_high < 0)
stars <- ifelse(sig, "**", "")
coef_star <- matrix(
  paste0(round(coef_hat, 3), stars),
  nrow = nrow(coef_hat),
  ncol = ncol(coef_hat),
  dimnames = dimnames(coef_hat)
)

se_formatted <- matrix(
  paste0("(", sprintf("%.3f", se_hat), ")"),
  nrow = nrow(se_hat)
)

final_table <- matrix("",
                      nrow = nrow(coef_hat) * 2,
                      ncol = ncol(coef_hat),
                      dimnames = list(
                        rep(rownames(coef_hat), each = 2),
                        colnames(coef_hat)
                      ))

final_table[seq(1, nrow(final_table), by = 2), ] <- coef_star
final_table[seq(2, nrow(final_table), by = 2), ] <- se_formatted

make_bootstrap_table(
  coef_hat,
  se_hat,
  ci_low,
  ci_high,
  file = "plots/iv_qoq/control_function_approach/bootstrap_results_CF.tex"
)

#################################################################
## 4. Control function approach robust (from 04_03b_IV_CF_boot_robust.R)
#################################################################


bootstrap_results <- readRDS(
  "plots/iv_qoq/control_function_approach/bootstrap_data_CF_robust.rds"
)
bootstrap_results_pct <- (exp(bootstrap_results * sd_d) - 1) * 100


coef_hat <- apply(bootstrap_results_pct, c(2, 3), mean, na.rm = TRUE)
se_hat   <- apply(bootstrap_results_pct, c(2, 3), sd, na.rm = TRUE)
ci_low   <- apply(bootstrap_results_pct, c(2, 3), quantile, probs = 0.025, na.rm = TRUE)
ci_high  <- apply(bootstrap_results_pct, c(2, 3), quantile, probs = 0.975, na.rm = TRUE)

result_cf <- list(
  coef = coef_hat,
  se   = se_hat,
  ci_low = ci_low,
  ci_high = ci_high
)

# ---- surface plot ----
boot_surface <- plot_qq_surface(
  coef_hat,
  taus = taus,
  tplot = 1,
  mycolors = "",
  z_label = "Shock Impact",
  group = "CZ",
  zlim = z_min_max
)

htmlwidgets::saveWidget(
  boot_surface,
  file = "plots/iv_qoq/control_function_approach/surface_boot_CF_robust.html",
  selfcontained = TRUE
)

webshot2::webshot(
  "plots/iv_qoq/control_function_approach/surface_boot_CF_robust.html",
  file = "plots/iv_qoq/control_function_approach/surface_boot_CF_robust.png",
  vwidth = 1050, vheight = 900, zoom = 3
)

# ---- slices plot ----
boot_slices <- plot_qq_slices(
  coef_hat,
  taus = taus,
  z_label = "IV Coefficient",
  group = "CZ",
  ylim = y_min_max
) +
  ggtitle("China Shock Impact - Control Function Approach")

ggsave(
  filename = "plots/iv_qoq/control_function_approach/slices_boot_CF_robust.png",
  plot = boot_slices,
  width = 8,
  height = 5,
  dpi = 300,
  bg = "white"
)

# ---- significance table ----
sig <- (ci_low > 0) | (ci_high < 0)
stars <- ifelse(sig, "**", "")
coef_star <- matrix(
  paste0(round(coef_hat, 3), stars),
  nrow = nrow(coef_hat),
  ncol = ncol(coef_hat),
  dimnames = dimnames(coef_hat)
)

se_formatted <- matrix(
  paste0("(", sprintf("%.3f", se_hat), ")"),
  nrow = nrow(se_hat)
)

final_table <- matrix("",
                      nrow = nrow(coef_hat) * 2,
                      ncol = ncol(coef_hat),
                      dimnames = list(
                        rep(rownames(coef_hat), each = 2),
                        colnames(coef_hat)
                      ))

final_table[seq(1, nrow(final_table), by = 2), ] <- coef_star
final_table[seq(2, nrow(final_table), by = 2), ] <- se_formatted

make_bootstrap_table(
  coef_hat,
  se_hat,
  ci_low,
  ci_high,
  file = "plots/iv_qoq/control_function_approach/bootstrap_results_CF_robust.tex"
)


#################################################################
## 4. Control function approach unweighted (from 04_03c_IV_CF_boot_unweighted.R)
#################################################################


bootstrap_results <- readRDS(
  "plots/iv_qoq/control_function_approach/bootstrap_data_CF_unweighted.rds"
)
bootstrap_results_pct <- (exp(bootstrap_results * sd_d) - 1) * 100


coef_hat <- apply(bootstrap_results_pct, c(2, 3), mean, na.rm = TRUE)
se_hat   <- apply(bootstrap_results_pct, c(2, 3), sd, na.rm = TRUE)
ci_low   <- apply(bootstrap_results_pct, c(2, 3), quantile, probs = 0.025, na.rm = TRUE)
ci_high  <- apply(bootstrap_results_pct, c(2, 3), quantile, probs = 0.975, na.rm = TRUE)

result_cf <- list(
  coef = coef_hat,
  se   = se_hat,
  ci_low = ci_low,
  ci_high = ci_high
)

# ---- surface plot ----
boot_surface <- plot_qq_surface(
  coef_hat,
  taus = taus,
  tplot = 1,
  mycolors = "",
  z_label = "IV Coefficient",
  group = "CZ",
  zlim = z_min_max
)

htmlwidgets::saveWidget(
  boot_surface,
  file = "plots/iv_qoq/control_function_approach/surface_boot_CF_unweighted.html",
  selfcontained = TRUE
)

webshot2::webshot(
  "plots/iv_qoq/control_function_approach/surface_boot_CF_unweighted.html",
  file = "plots/iv_qoq/control_function_approach/surface_boot_CF_unweighted.png",
  vwidth = 1050, vheight = 900, zoom = 3
)

# ---- slices plot ----
boot_slices <- plot_qq_slices(
  coef_hat,
  taus = taus,
  z_label = "IV Coefficient",
  group = "CZ",
  ylim = y_min_max
) +
  ggtitle("China Shock Impact - Control Function Approach")

ggsave(
  filename = "plots/iv_qoq/control_function_approach/slices_boot_CF_unweighted.png",
  plot = boot_slices,
  width = 8,
  height = 5,
  dpi = 300,
  bg = "white"
)

# ---- significance table ----
sig <- (ci_low > 0) | (ci_high < 0)
stars <- ifelse(sig, "**", "")
coef_star <- matrix(
  paste0(round(coef_hat, 3), stars),
  nrow = nrow(coef_hat),
  ncol = ncol(coef_hat),
  dimnames = dimnames(coef_hat)
)

se_formatted <- matrix(
  paste0("(", sprintf("%.3f", se_hat), ")"),
  nrow = nrow(se_hat)
)

final_table <- matrix("",
                      nrow = nrow(coef_hat) * 2,
                      ncol = ncol(coef_hat),
                      dimnames = list(
                        rep(rownames(coef_hat), each = 2),
                        colnames(coef_hat)
                      ))

final_table[seq(1, nrow(final_table), by = 2), ] <- coef_star
final_table[seq(2, nrow(final_table), by = 2), ] <- se_formatted

make_bootstrap_table(
  coef_hat,
  se_hat,
  ci_low,
  ci_high,
  file = "plots/iv_qoq/control_function_approach/bootstrap_results_CF_unweighted.tex"
)








