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
library(AER)
library(purrr)
library(ggplot2)
library(sandwich)
library(lmtest)
library(stargazer)


shocks <- read_dta("data/czone_exposure_by_period_v5_gh.dta")
controls <- read_dta("data/ADH_control_vars.dta")
data <- readRDS("data/data.rds")
ACS_pop <- read_dta("data/ACS_pop_emp_inc.dta")

####### get rid of outliers 
data <- data %>%
  filter(wage > 5) %>%
  filter(wage < 500)

# for interpretation: 
summary(shocks$d_tradeusch_p1_2000_2012)

shocks

####################################################################################################
## 1. IV estimation (Approach analog to Autor et al. 2021)
####################################################################################################
dir.create("plots/autor_replication", showWarnings = FALSE)

controls_1990 <- controls %>%
  filter(year == 1990)

shocks <- shocks %>%
  left_join(controls_1990, by = "czone")

controls_1990$ba_pop

median_wage_cz <- data %>%
  filter(YEAR >= 2005, YEAR <= 2019) %>%
  filter(!is.na(czone)) %>%
  group_by(YEAR, czone) %>%
  summarise(
    median_wage = wtd.quantile(wage, weights = pweight, probs = 0.5, na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  )

IV_med_wage <- median_wage_cz %>%
  left_join(shocks, by = "czone")

# weighting by 2000 working age population by CZ
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

IV_med_wage <- IV_med_wage %>%
  left_join(
    pop_2000 %>% select(czone, pop_share_2000),
    by = "czone"
  )

# 1.1 IV estimation for every year

years <- 2005:2019
results_year <- list()
vcov_year <- list()

for (yr in years) {
  
  df_year <- IV_med_wage %>%
    filter(YEAR == yr)
  
  model <- ivreg(
    log(median_wage) ~ 
      d_tradeusch_p1_2000_2012 +
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
      sh_00up_nw
    |
      d_tradeotch_p1_lag_2000_2012 +
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
      sh_00up_nw,
    data = df_year,
    weights = pop_share_2000
  )
  
  results_year[[as.character(yr)]] <- model
  vcov_year[[as.character(yr)]] <- vcovCL(model, cluster = df_year$statefip)
  
}

coef_table_year <- do.call(rbind, lapply(years, function(yr) {
  
  model <- results_year[[as.character(yr)]]
  vcov_m <- vcov_year[[as.character(yr)]]
  
  est <- coef(model)["d_tradeusch_p1_2000_2012"]
  se  <- sqrt(vcov_m["d_tradeusch_p1_2000_2012",
                     "d_tradeusch_p1_2000_2012"])
  
  data.frame(
    year = yr,
    estimate = est,
    se = se,
    ci_low = est - 1.96 * se,
    ci_high = est + 1.96 * se
  )
}))

autor_plot_year <- ggplot(coef_table_year, aes(x = year, y = estimate)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "2000-2012 shock impact on pooled log Median Wage (2005-2019)",
    x = "Year",
    y = "Effect of trade exposure"
  )

ggsave(
  filename = "plots/autor_replication/autor_year_effect_median.pdf",
  plot = autor_plot_year,
  width = 6,
  height = 6,
  dpi = 300
)



# 1.2 pooled IV all years

model_pool_all_median <- ivreg(
  log(median_wage) ~ 
    d_tradeusch_p1_2000_2012 +
    l_shind_manuf_cbp +
    l_sh_popedu_c +
    l_sh_popfborn +
    l_sh_empl_f +
    l_task_outsource +
    region +
    sh_65up_all +
    sh_4064_all +
    sh_0017_all +
    sh_00up_nw
  |
    d_tradeotch_p1_lag_2000_2012 +
    l_shind_manuf_cbp +
    l_sh_popedu_c +
    l_sh_popfborn +
    l_sh_empl_f +
    l_task_outsource +
    region +
    sh_65up_all +
    sh_4064_all +
    sh_0017_all +
    sh_00up_nw,
  data = IV_med_wage,
  weights = pop_share_2000
)

# cluster by state
vcov_pool_all_median <- vcovCL(
  model_pool_all_median,
  cluster = IV_med_wage$statefip
)

se_pool_all_median <- sqrt(diag(vcov_pool_all_median))



# 1.3 pooled IV 5-year groups

IV_med_wage <- IV_med_wage %>%
  mutate(
    period = case_when(
      YEAR %in% 2005:2009 ~ "2005-2009",
      YEAR %in% 2010:2014 ~ "2010-2014",
      YEAR %in% 2015:2019 ~ "2015-2019",
      TRUE ~ NA_character_
    )
  )

periods_five <- unique(IV_med_wage$period)
periods_five <- periods_five[!is.na(periods_five)]

results_pool_five <- list()
vcov_five <- list()


for (p in periods_five) {
  
  df_period <- IV_med_wage %>%
    filter(period == p)
  
  model <- ivreg(
    log(median_wage) ~ 
      d_tradeusch_p1_2000_2012 +
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
      sh_00up_nw
    |
      d_tradeotch_p1_lag_2000_2012 +
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
      sh_00up_nw,
    data = df_period,
    weights = pop_share_2000
  )
  
  results_pool_five[[p]] <- model
  vcov_five[[as.character(p)]] <- vcovCL(model, cluster = df_period$statefip)
  
}

coef_table_five <- do.call(rbind, lapply(periods_five, function(p) {
  
  model <- results_pool_five[[as.character(p)]]
  vcov_m <- vcov_five[[as.character(p)]]
  
  est <- coef(model)["d_tradeusch_p1_2000_2012"]
  se  <- sqrt(vcov_m["d_tradeusch_p1_2000_2012",
                     "d_tradeusch_p1_2000_2012"])
  
  data.frame(
    period = p,
    estimate = est,
    se = se,
    ci_low = est - 1.96 * se,
    ci_high = est + 1.96 * se
  )
}))

autor_plot_five <- ggplot(coef_table_five, aes(x = period, y = estimate, group = 1)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high),
              alpha = 0.2)+
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "2000-2012 shock impact on pooled log Median Wage (2005-2019)",
    x = "Period",
    y = "Effect of trade exposure"
  )

ggsave(
  filename = "plots/autor_replication/autor_five_year.pdf",
  plot = autor_plot_five,
  width = 6,
  height = 6,
  dpi = 300
)




####################################################################################################
## 2. IV estimation mean wage (Approach analog to Autor et al. 2021)
####################################################################################################


mean_wage_cz <- data %>%
  filter(YEAR >= 2005, YEAR <= 2019) %>%
  filter(!is.na(czone)) %>%
  group_by(YEAR, czone) %>%
  summarise(
    mean_wage = weighted.mean(wage, w = pweight, na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  )

IV_mean_wage <- mean_wage_cz %>%
  left_join(shocks, by = "czone")



IV_mean_wage <- IV_mean_wage %>%
  left_join(
    pop_2000 %>% select(czone, pop_share_2000),
    by = "czone"
  )

# 2.1 IV estimation for every year

years <- 2005:2019
results_year <- list()
vcov_year <- list()

for (yr in years) {
  
  df_year <- IV_mean_wage %>%
    filter(YEAR == yr)
  
  model <- ivreg(
    log(mean_wage) ~ 
      d_tradeusch_p1_2000_2012 +
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
      sh_00up_nw
    |
      d_tradeotch_p1_lag_2000_2012 +
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
      sh_00up_nw,
    data = df_year,
    weights = pop_share_2000
  )
  
  results_year[[as.character(yr)]] <- model
  vcov_year[[as.character(yr)]] <- vcovCL(model, cluster = df_year$statefip)
  
}

coef_table_year <- do.call(rbind, lapply(years, function(yr) {
  
  model <- results_year[[as.character(yr)]]
  vcov_m <- vcov_year[[as.character(yr)]]
  
  est <- coef(model)["d_tradeusch_p1_2000_2012"]
  se  <- sqrt(vcov_m["d_tradeusch_p1_2000_2012",
                     "d_tradeusch_p1_2000_2012"])
  
  data.frame(
    year = yr,
    estimate = est,
    se = se,
    ci_low = est - 1.96 * se,
    ci_high = est + 1.96 * se
  )
}))

autor_plot_year_mean <- ggplot(coef_table_year, aes(x = year, y = estimate)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "2000-2012 shock impact on pooled log Mean Wage (2005-2019)",
    x = "Year",
    y = "Effect of trade exposure"
  )

ggsave(
  filename = "plots/autor_replication/autor_year_effect_mean.pdf",
  plot = autor_plot_year_mean,
  width = 6,
  height = 6,
  dpi = 300
)



# 2.2 pooled IV all years

model_pool_all_mean <- ivreg(
  log(mean_wage) ~ 
    d_tradeusch_p1_2000_2012 +
    l_shind_manuf_cbp +
    l_sh_popedu_c +
    l_sh_popfborn +
    l_sh_empl_f +
    l_task_outsource +
    region +
    sh_65up_all +
    sh_4064_all +
    sh_0017_all +
    sh_00up_nw
  |
    d_tradeotch_p1_lag_2000_2012 +
    l_shind_manuf_cbp +
    l_sh_popedu_c +
    l_sh_popfborn +
    l_sh_empl_f +
    l_task_outsource +
    region +
    sh_65up_all +
    sh_4064_all +
    sh_0017_all +
    sh_00up_nw,
  data = IV_mean_wage,
  weights = pop_share_2000
)

# cluster by state
vcov_pool_all_mean <- vcovCL(
  model_pool_all_mean,
  cluster = IV_mean_wage$statefip
)

se_pool_all_mean <- sqrt(diag(vcov_pool_all_mean))



# 2.3 pooled IV 5-year groups

IV_mean_wage <- IV_mean_wage %>%
  mutate(
    period = case_when(
      YEAR %in% 2005:2009 ~ "2005-2009",
      YEAR %in% 2010:2014 ~ "2010-2014",
      YEAR %in% 2015:2019 ~ "2015-2019",
      TRUE ~ NA_character_
    )
  )

periods_five <- unique(IV_mean_wage$period)
periods_five <- periods_five[!is.na(periods_five)]

results_pool_five <- list()
vcov_five <- list()


for (p in periods_five) {
  
  df_period <- IV_mean_wage %>%
    filter(period == p)
  
  model <- ivreg(
    log(mean_wage) ~ 
      d_tradeusch_p1_2000_2012 +
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
      sh_00up_nw
    |
      d_tradeotch_p1_lag_2000_2012 +
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
      sh_00up_nw,
    data = df_period,
    weights = pop_share_2000
  )
  
  results_pool_five[[p]] <- model
  vcov_five[[as.character(p)]] <- vcovCL(model, cluster = df_period$statefip)
  
}

coef_table_five <- do.call(rbind, lapply(periods_five, function(p) {
  
  model <- results_pool_five[[as.character(p)]]
  vcov_m <- vcov_five[[as.character(p)]]
  
  est <- coef(model)["d_tradeusch_p1_2000_2012"]
  se  <- sqrt(vcov_m["d_tradeusch_p1_2000_2012",
                     "d_tradeusch_p1_2000_2012"])
  
  data.frame(
    period = p,
    estimate = est,
    se = se,
    ci_low = est - 1.96 * se,
    ci_high = est + 1.96 * se
  )
}))

autor_plot_five_mean <- ggplot(coef_table_five, aes(x = period, y = estimate, group = 1)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high),
              alpha = 0.2)+
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "2000-2012 shock impact on pooled log mean Wage (2005-2019)",
    x = "Period",
    y = "Effect of trade exposure"
  )

ggsave(
  filename = "plots/autor_replication/autor_five_year_mean.pdf",
  plot = autor_plot_five_mean,
  width = 6,
  height = 6,
  dpi = 300
)

####################################################################################################
## 3. IV estimation sum of wage (Approach analog to Autor et al. 2021)
####################################################################################################
#data$INCWAGE_CPIU_2010 = total wage income per person 

sum_wage_cz <- data %>%
  filter(YEAR >= 2005, YEAR <= 2019) %>%
  filter(!is.na(czone)) %>%
  group_by(YEAR, czone) %>%
  summarise(
    sum_wage = sum(INCWAGE_CPIU_2010 * pweight, na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  )

IV_sum_wage <- sum_wage_cz %>%
  left_join(shocks, by = "czone")



IV_sum_wage <- IV_sum_wage %>%
  left_join(
    pop_2000 %>% select(czone, pop_share_2000),
    by = "czone"
  )

# 3.1 IV estimation for every year

years <- 2005:2019
results_year <- list()
vcov_year <- list()

for (yr in years) {
  
  df_year <- IV_sum_wage %>%
    filter(YEAR == yr)
  
  model <- ivreg(
    log(sum_wage) ~ 
      d_tradeusch_p1_2000_2012 +
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
      sh_00up_nw
    |
      d_tradeotch_p1_lag_2000_2012 +
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
      sh_00up_nw,
    data = df_year,
    weights = pop_share_2000
  )
  
  results_year[[as.character(yr)]] <- model
  vcov_year[[as.character(yr)]] <- vcovCL(model, cluster = df_year$statefip)
  
}

coef_table_year <- do.call(rbind, lapply(years, function(yr) {
  
  model <- results_year[[as.character(yr)]]
  vcov_m <- vcov_year[[as.character(yr)]]
  
  est <- coef(model)["d_tradeusch_p1_2000_2012"]
  se  <- sqrt(vcov_m["d_tradeusch_p1_2000_2012",
                     "d_tradeusch_p1_2000_2012"])
  
  data.frame(
    year = yr,
    estimate = est,
    se = se,
    ci_low = est - 1.96 * se,
    ci_high = est + 1.96 * se
  )
}))

autor_plot_year_sum <- ggplot(coef_table_year, aes(x = year, y = estimate)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "2000-2012 shock impact on pooled log sum Wage (2005-2019)",
    x = "Year",
    y = "Effect of trade exposure"
  )

ggsave(
  filename = "plots/autor_replication/autor_year_effect_sum.pdf",
  plot = autor_plot_year_sum,
  width = 6,
  height = 6,
  dpi = 300
)



# 3.2 pooled IV all years

model_pool_all_sum <- ivreg(
  log(sum_wage) ~ 
    d_tradeusch_p1_2000_2012 +
    l_shind_manuf_cbp +
    l_sh_popedu_c +
    l_sh_popfborn +
    l_sh_empl_f +
    l_task_outsource +
    region +
    sh_65up_all +
    sh_4064_all +
    sh_0017_all +
    sh_00up_nw
  |
    d_tradeotch_p1_lag_2000_2012 +
    l_shind_manuf_cbp +
    l_sh_popedu_c +
    l_sh_popfborn +
    l_sh_empl_f +
    l_task_outsource +
    region +
    sh_65up_all +
    sh_4064_all +
    sh_0017_all +
    sh_00up_nw,
  data = IV_sum_wage,
  weights = pop_share_2000
)

# cluster by state
vcov_pool_all_sum <- vcovCL(
  model_pool_all_sum,
  cluster = IV_sum_wage$statefip
)

se_pool_all_sum <- sqrt(diag(vcov_pool_all_sum))



# 3.3 pooled IV 5-year groups

IV_sum_wage <- IV_sum_wage %>%
  mutate(
    period = case_when(
      YEAR %in% 2005:2009 ~ "2005-2009",
      YEAR %in% 2010:2014 ~ "2010-2014",
      YEAR %in% 2015:2019 ~ "2015-2019",
      TRUE ~ NA_character_
    )
  )

periods_five <- unique(IV_sum_wage$period)
periods_five <- periods_five[!is.na(periods_five)]

results_pool_five <- list()
vcov_five <- list()


for (p in periods_five) {
  
  df_period <- IV_sum_wage %>%
    filter(period == p)
  
  model <- ivreg(
    log(sum_wage) ~ 
      d_tradeusch_p1_2000_2012 +
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
      sh_00up_nw
    |
      d_tradeotch_p1_lag_2000_2012 +
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
      sh_00up_nw,
    data = df_period,
    weights = pop_share_2000
  )
  
  results_pool_five[[p]] <- model
  vcov_five[[as.character(p)]] <- vcovCL(model, cluster = df_period$statefip)
  
}

coef_table_five <- do.call(rbind, lapply(periods_five, function(p) {
  
  model <- results_pool_five[[as.character(p)]]
  vcov_m <- vcov_five[[as.character(p)]]
  
  est <- coef(model)["d_tradeusch_p1_2000_2012"]
  se  <- sqrt(vcov_m["d_tradeusch_p1_2000_2012",
                     "d_tradeusch_p1_2000_2012"])
  
  data.frame(
    period = p,
    estimate = est,
    se = se,
    ci_low = est - 1.96 * se,
    ci_high = est + 1.96 * se
  )
}))

autor_plot_five_sum <- ggplot(coef_table_five, aes(x = period, y = estimate, group = 1)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high),
              alpha = 0.2)+
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "2000-2012 shock impact on pooled log sum Wage (2005-2019)",
    x = "Period",
    y = "Effect of trade exposure"
  )

ggsave(
  filename = "plots/autor_replication/autor_five_year_sum.pdf",
  plot = autor_plot_five_sum,
  width = 6,
  height = 6,
  dpi = 300
)




# 4.0 all models together for all years 
stargazer(
  model_pool_all_sum,
  model_pool_all_mean,
  model_pool_all_median,
  type = "latex",
  title = "Pooled IV Estimates of Trade Exposure on Wages",
  column.labels = c(
    "Log Sum Wage",
    "Log Mean Wage",
    "Log Median Wage"
  ),
  
  dep.var.labels.include = FALSE,
  
  se = list(
    se_pool_all_sum,
    se_pool_all_mean,
    se_pool_all_median
  ),
  
  keep = "d_tradeusch_p1_2000_2012",
  
  covariate.labels = c(
    "Trade Exposure"
  ),
  
  keep.stat = c("n"),
  
  digits = 3,
  
  notes = "Clustered standard errors at the state level in parentheses.",
  
  out = "plots/autor_replication/pooled_iv_results.tex"
)
