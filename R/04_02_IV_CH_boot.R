library(haven)
library(ipumsr)
library(data.table)
library(sf)
library(Hmisc)
library(dplyr)
library(stringr)
library(tmap)
library(qqreg)
library(quantreg)
library(car)
library(ivreg)
library(Formula)
#install ivqr package
#install.packages("remotes")
#remotes::install_github("yuchang0321/IVQR")
library(IVQR)
#load quantile on quantiles functions
source("R/00_2_functions.R")

shocks <- read_dta("data/czone_exposure_by_period_v5_gh.dta")
controls <- read_dta("data/ADH_control_vars.dta")
data <- readRDS("data/data.rds")
ACS_pop <- read_dta("data/ACS_pop_emp_inc.dta")

# select important variables
data <- data %>%
  select(czone, wage, lnwage, pweight, STATEFIP)

# set seed
set.seed(1234)

valid_cz <- unique(shocks$czone)

# controls for 1990
controls_1990 <- controls %>%
  filter(year == 1990)  %>%
  select(czone, l_shind_manuf_cbp, l_sh_popedu_c, l_sh_popfborn,
         l_sh_empl_f, l_sh_routine33, l_task_outsource,
         region, sh_65up_all, sh_4064_all, sh_0017_all, sh_00up_nw, statefip)

# weighting by 2000 working age population by CZ as in Autor 2021
pop_2000 <- ACS_pop %>%
  filter(year == 2000) %>%
  group_by(czone) %>%
  summarise(
    population_2000 = mean(population_bl_all, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    total_population_2000 = sum(population_2000, na.rm = TRUE),
    pop_share_2000 = population_2000 / total_population_2000
  )

pop_2000 <- pop_2000 %>%
  filter(czone %in% valid_cz)

# merge controls and shocks by CZ for the first stage IV 
shocks <- shocks %>%
  left_join(controls_1990, by = "czone")

# just do the most complex one: IV_2 in Autor 2021 
formula_first_stage_IV_2 <- d_tradeusch_p1_2000_2012 ~ d_tradeotch_p1_lag_2000_2012 + l_shind_manuf_cbp + l_sh_popedu_c + l_sh_popfborn + l_sh_empl_f + l_sh_routine33 + l_task_outsource + region +  sh_65up_all + sh_4064_all + sh_0017_all + sh_00up_nw

# parameters 
# common limits
z_min_max <- c(-7, 9)
y_min_max <- c(-7, 9)

years <- 2005:2019
grid <- seq(-10, 10, by = 0.1)

# quantiles
taus <- round(seq(0.1, 0.9, by = 0.1), 1)


# create folders
dir.create(
  "plots/iv_qoq",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "plots/iv_qoq/IVQR_approach",
  recursive = TRUE,
  showWarnings = FALSE
)
#################################################################
## QoQ-IVQR-IV appraoch for all years with bootstrap std. errors 
#################################################################

# bootstrap function
bootstrap_qoq_iv_CH <- function(data, shocks) {
  
  # -------------------------------------------------------
  # 1. Resample data by state clusters
  # -------------------------------------------------------
  state_czone <- data %>%
    select(STATEFIP, czone) %>%
    distinct()
  
  boot_czone_index <- boot_draws %>%
    inner_join(state_czone, by = "STATEFIP") %>%
    mutate(czone_boot = paste0(czone, "_", boot_id))
  
  boot_data <- boot_czone_index %>%
    left_join(data, by = c("STATEFIP", "czone"))

  boot_shocks <- boot_czone_index %>%
    inner_join(shocks, by = c("czone")) %>%
    left_join(pop_2000, by = "czone")
  
  # -------------------------------------------------------
  # 2. QQ first stage
  # -------------------------------------------------------
  qq_first <- qq_first_stage_collapsed(
    y = boot_data$lnwage,
    group_vec = boot_data$czone_boot,
    taus = taus,
    pweights = boot_data$pweight
  )
  
  cz <- qq_first$group_order
  
  X_df <- boot_shocks %>%
    filter(czone_boot %in% cz) %>%
    mutate(
      czone_boot = factor(
        czone_boot,
        levels = cz
      )
    ) %>%
    arrange(czone_boot)
  
  # align shocks
  X <- model.matrix(
    ~ l_shind_manuf_cbp +
      l_sh_popedu_c +
      l_sh_popfborn +
      l_sh_empl_f +
      l_sh_routine33 +
      l_task_outsource +
      region +
      sh_65up_all +
      sh_4064_all +
      sh_0017_all +
      sh_00up_nw - 1,
    data = X_df
  )
  
  D <- X_df$d_tradeusch_p1_2000_2012
  Z <- X_df$d_tradeotch_p1_lag_2000_2012
  
  pop_w <- X_df$pop_share_2000  
  # -------------------------------------------------------
  # 3. Second stage
  # -------------------------------------------------------
  coeffs <- qq_second_stage_ivqr_group(
    fitted_first = qq_first$fitted,
    D = D,
    Z = Z,
    X = X,
    taus = taus,
    grid = grid
  )
  
  iv_coef <- coeffs$beta_D
  
  return(iv_coef)
}

# do the bootstrapping 

states <- unique(shocks$statefip)

B <- 100 # increase later (500–1000 in production)


bootstrap_results <- array(
  NA,
  dim = c(B, length(taus), length(taus)),
  dimnames = list(
    paste0("b_", 1:B),
    paste0("u_", taus),
    paste0("v_", taus)
  )
)


for (b in 1:B) {
  
  # sample states WITH replacement
  boot_states <- sample(states, length(states), replace = TRUE)
  
  state_map <- data.frame(
    STATEFIP = states,
    boot_id = seq_along(states)
  )
  
  boot_draws <- data.frame(
    boot_id = seq_along(boot_states),
    STATEFIP = boot_states
  )  
  boot_draws$statefip <- boot_draws$STATEFIP
  # run bootstrap
  bootstrap_results[b, , ] <- bootstrap_qoq_iv_CH(
    data = data, 
    shocks = shocks
  )
  
  cat("Bootstrap", b, "done\n")
}

saveRDS(
  bootstrap_results,
  file = "plots/iv_qoq/IVQR_approach/bootstrap_data_IVQR.rds"
)


B <- dim(bootstrap_results)[1]
U <- dim(bootstrap_results)[2]
V <- dim(bootstrap_results)[3]

coef_hat <- apply(bootstrap_results, c(2, 3), mean, na.rm = TRUE)
se_hat <- apply(bootstrap_results, c(2, 3), sd, na.rm = TRUE)
ci_low <- apply(bootstrap_results, c(2, 3), quantile, probs = 0.025, na.rm = TRUE)
ci_high <- apply(bootstrap_results, c(2, 3), quantile, probs = 0.975, na.rm = TRUE)

result <- list(
  coef = coef_hat,
  se   = se_hat,
  ci_low = ci_low,
  ci_high = ci_high
)

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
  file = "plots/iv_qoq/IVQR_approach/surface_boot_IVQR.html",
  selfcontained = TRUE
)


boot_slices <- plot_qq_slices(
  coef_hat,
  taus = taus,
  z_label = "Shock Impact",
  group = "CZ",
  ylim = y_min_max
) +
  ggtitle("China Shock Impact - Bootstrapped")

# save png
ggsave(
  filename = "plots/iv_qoq/IVQR_approach/slices_boot_IVQR.png",
  plot = boot_slices,
  width = 8,
  height = 5,
  dpi = 300,
  bg = "white"
)

# significance matrix
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

stargazer::stargazer(
  as.data.frame(final_table),
  summary = FALSE,
  rownames = TRUE,
  type = "latex",
  out = "plots/iv_qoq/IVQR_approach/bootstrap_results_IVQR.tex"
)


#################################################################
## QoQ-IVQR-IV approach: full-sample point estimates (no bootstrap)
#################################################################

# -------------------------------------------------------
# 1. QQ first stage (full sample)
# -------------------------------------------------------
qq_first_full <- qq_first_stage_collapsed(
  y = data$lnwage,
  group_vec = data$czone,
  taus = taus,
  pweights = data$pweight
)

cz_full <- qq_first_full$group_order

X_df_full <- shocks %>%
  left_join(pop_2000, by = "czone") %>%
  filter(czone %in% cz_full) %>%
  mutate(
    czone = factor(
      czone,
      levels = cz_full
    )
  ) %>%
  arrange(czone)

# align shocks
X_full <- model.matrix(
  ~ l_shind_manuf_cbp +
    l_sh_popedu_c +
    l_sh_popfborn +
    l_sh_empl_f +
    l_sh_routine33 +
    l_task_outsource +
    region +
    sh_65up_all +
    sh_4064_all +
    sh_0017_all +
    sh_00up_nw - 1,
  data = X_df_full
)

D_full <- X_df_full$d_tradeusch_p1_2000_2012
Z_full <- X_df_full$d_tradeotch_p1_lag_2000_2012

# -------------------------------------------------------
# 2. Second stage (full sample)
# -------------------------------------------------------
coeffs_full <- qq_second_stage_ivqr_group(
  fitted_first = qq_first_full$fitted,
  D = D_full,
  Z = Z_full,
  X = X_full,
  taus = taus,
  grid = grid
)

point_estimates_IVQR <- coeffs_full$beta_D

dimnames(point_estimates_IVQR) <- list(
  paste0("u_", taus),
  paste0("v_", taus)
)

saveRDS(
  point_estimates_IVQR,
  file = "plots/iv_qoq/IVQR_approach/point_estimates_IVQR.rds"
)

