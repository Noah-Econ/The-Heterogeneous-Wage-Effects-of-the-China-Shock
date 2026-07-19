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

# quantiles
taus <- round(seq(0.1, 0.9, by = 0.1), 1)

# controls for CF
control_vars <- c(
  "l_shind_manuf_cbp",
  "l_sh_popedu_c",
  "l_sh_popfborn",
  "l_sh_empl_f",
  "l_sh_routine33",
  "l_task_outsource",
  "region",
  "sh_65up_all",
  "sh_4064_all",
  "sh_0017_all",
  "sh_00up_nw"
)

# create folders
dir.create(
  "plots/iv_qoq",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "plots/iv_qoq/control_function_approach",
  recursive = TRUE,
  showWarnings = FALSE
)
#################################################################
## QoQ-CF-IV appraoch for all years with bootstrap std. errors 
#################################################################

# bootstrap function
bootstrap_qoq_iv_CF <- function(data, shocks) {
  
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
  # 2. Estimate Control Function
  # -------------------------------------------------------
  
  boot_shocks <- estimate_cf_distribution(
    data = boot_shocks,
    x_var = "d_tradeusch_p1_2000_2012",
    z_var = "d_tradeotch_p1_lag_2000_2012",
    control_vars = control_vars
  )
  
  # -------------------------------------------------------
  # 3. QQ first stage
  # -------------------------------------------------------
  qq_first <- qq_first_stage_collapsed(
    y = boot_data$lnwage,
    group_vec = boot_data$czone_boot,
    taus = taus,
    pweights = boot_data$pweight
  )
  
  cz <- qq_first$group_order
  
  # align shocks
  X_df <- boot_shocks %>%
    filter(czone_boot %in% cz) %>%
    mutate(
      czone_boot = factor(
        czone_boot,
        levels = cz
      )
    ) %>%
    arrange(czone_boot)
  
  X <- model.matrix(
    ~ d_tradeusch_p1_2000_2012 +
      l_shind_manuf_cbp +
      l_sh_popedu_c +
      l_sh_popfborn +
      l_sh_empl_f +
      l_sh_routine33 +
      l_task_outsource +
      region +
      sh_65up_all +
      sh_4064_all +
      sh_0017_all +
      sh_00up_nw +
      V_hat,
    data = X_df
  )
  
  pop_w <- X_df$pop_share_2000  
  # -------------------------------------------------------
  # 4. Second stage
  # -------------------------------------------------------
  coeffs <- qq_second_stage_group(
    fitted_first = qq_first$fitted,
    X = X,
    taus = taus
  )
  coefficients <- coeffs$coef
  
  
  iv_coef <- coefficients[
    "d_tradeusch_p1_2000_2012",
    ,
  ]
  
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

t0 <- Sys.time()

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
  bootstrap_results[b, , ] <- bootstrap_qoq_iv_CF(
    data = data,
    shocks = shocks)
  
  cat("Bootstrap", b, "done\n")
}

t1 <- Sys.time()
t <- t1 - t0

saveRDS(
  bootstrap_results,
  file = "plots/iv_qoq/control_function_approach/bootstrap_data_CF_unweighted.rds"
)


#################################################################
## QoQ-CF-IV approach: full-sample point estimates (no bootstrap)
#################################################################

# -------------------------------------------------------
# 1. Estimate Control Function (full sample, no resampling)
# -------------------------------------------------------
shocks_full <- shocks %>%
  left_join(pop_2000, by = "czone")

shocks_full <- estimate_cf_distribution(
  data = shocks_full,
  x_var = "d_tradeusch_p1_2000_2012",
  z_var = "d_tradeotch_p1_lag_2000_2012",
  control_vars = control_vars
)

# -------------------------------------------------------
# 2. QQ first stage (full sample)
# -------------------------------------------------------
qq_first_full <- qq_first_stage_collapsed(
  y = data$lnwage,
  group_vec = data$czone,
  taus = taus,
  pweights = data$pweight
)

cz_full <- qq_first_full$group_order

# align shocks
X_df_full <- shocks_full %>%
  filter(czone %in% cz_full) %>%
  mutate(
    czone = factor(
      czone,
      levels = cz_full
    )
  ) %>%
  arrange(czone)

X_full <- model.matrix(
  ~ d_tradeusch_p1_2000_2012 +
    l_shind_manuf_cbp +
    l_sh_popedu_c +
    l_sh_popfborn +
    l_sh_empl_f +
    l_sh_routine33 +
    l_task_outsource +
    region +
    sh_65up_all +
    sh_4064_all +
    sh_0017_all +
    sh_00up_nw +
    V_hat,
  data = X_df_full
)

pop_w_full <- X_df_full$pop_share_2000

# -------------------------------------------------------
# 3. Second stage (full sample)
# -------------------------------------------------------
coeffs_full <- qq_second_stage_group(
  fitted_first = qq_first_full$fitted,
  X = X_full,
  taus = taus
)

coefficients_full <- coeffs_full$coef

point_estimates_CF <- coefficients_full[
  "d_tradeusch_p1_2000_2012",
  ,
]

dimnames(point_estimates_CF) <- list(
  paste0("u_", taus),
  paste0("v_", taus)
)

saveRDS(
  point_estimates_CF,
  file = "plots/iv_qoq/control_function_approach/point_estimates_CF_unweighted.rds"
)

