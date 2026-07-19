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
#load quantile on quantiles functions
source("R/00_2_functions.R")

shocks <- read_dta("data/czone_exposure_by_period_v5_gh.dta")
controls <- read_dta("data/ADH_control_vars.dta")
data <- readRDS("data/data.rds")
ACS_pop <- read_dta("data/ACS_pop_emp_inc.dta")

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


####### get rid of outliers 
data <- data %>%
  filter(wage > 5) %>%
  filter(wage < 500)

######################################################################################################################################################################
## 1. IV estimation (Approach analog to Autor et al. 2021)
######################################################################################################################################################################
# use controls 1990 Autor
controls_1990 <- controls %>%
  filter(year == 1990)
# merge controls and shocks by CZ for the first stage IV 
shocks <- shocks %>%
  left_join(controls_1990, by = "czone")

# just do the most complex one: IV_2 in Autor 2021 
formula_first_stage_IV_2 <- d_tradeusch_p1_2000_2012 ~ d_tradeotch_p1_lag_2000_2012 + l_shind_manuf_cbp + l_sh_popedu_c + l_sh_popfborn + l_sh_empl_f + l_sh_routine33 + l_task_outsource + region +  sh_65up_all + sh_4064_all + sh_0017_all + sh_00up_nw

# run first stage IV
first_stage <- lm(formula_first_stage_IV_2, data = shocks)
summary(first_stage) # high significance of instrument 
shocks$IV_d_tradeusch_p1_2000_2012 <- predict(first_stage) 


# merge shocks + controls
cz_master <- unique(shocks$czone)

data <- data %>%
  filter(czone %in% cz_master)
df <- data %>%
  left_join(
    shocks %>% select(czone, IV_d_tradeusch_p1_2000_2012),
    by = "czone"
  )


df <- df[!is.na(df$pweight),]
#################################################################
## Chetverikov cross sectionally over the whole time period
#################################################################

# create folders
dir.create(
  "plots/chetverikov_replication",
  recursive = TRUE,
  showWarnings = FALSE
)

# quantiles
taus <- round(seq(0.1, 0.9, by = 0.1), 1)

qq_first <- qq_first_stage_collapsed(
  y = df$lnwage,
  group_vec = df$czone,
  taus = taus, 
  pweights = df$pweight #weighted by sample weights
)

# -------------------------------------------------------------
# build X matrix
# -------------------------------------------------------------

# master order
cz <- qq_first$group_order

# X aligned
X_df <- shocks %>%
  filter(czone %in% cz) %>%
  mutate(czone = factor(czone, levels = cz)) %>%
  arrange(czone)

X <- model.matrix(
  ~ IV_d_tradeusch_p1_2000_2012 +
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
  data = X_df
)

# weights aligned
pop_w <- data.frame(czone = qq_first$group_order) %>%
  left_join(pop_2000, by = "czone") %>%
  pull(pop_share_2000)
  

# -------------------------------------------------------------
# second stage
# -------------------------------------------------------------
  
coeffs <- qq_second_stage_chetverikov(
  fitted_first = qq_first$fitted,
  X = X,
  taus = taus,
  weights = pop_w
)
  
beta <- coeffs$coef["IV_d_tradeusch_p1_2000_2012", ]
se   <- coeffs$se["IV_d_tradeusch_p1_2000_2012", ]

ci_low  <- beta - 1.96 * se
ci_high <- beta + 1.96 * se
# -------------------------------------------------------------
# coefficients plot
# -------------------------------------------------------------
coeff_plot <- ggplot(
  data.frame(
    u = taus,
    beta = beta,
    ci_low = ci_low,
    ci_high = ci_high
  ),
  aes(x = u, y = beta)
) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Quantile (u)",
    y = "IV coefficient"
  )
  
# save html widget
ggsave("plots/chetverikov_replication/quant_coeffs.png",
  plot = coeff_plot,
  width = 8,
  height = 5,
  dpi = 300, 
  bg= "white"
)


#################################################################
## Chetverikov from 2005 to 2019
#################################################################

# create folders
dir.create(
  "plots/chetverikov_replication",
  recursive = TRUE,
  showWarnings = FALSE
)

# quantiles
taus <- round(seq(0.1, 0.9, by = 0.1), 1)

df_2005 <- data %>% filter(YEAR == 2005)
df_2019 <- data %>% filter(YEAR == 2019)

qq_2005 <- qq_first_stage_collapsed(
  y = df_2005$lnwage,
  group_vec = df_2005$czone,
  taus = taus,
  pweights = df_2005$pweight
)

qq_2019 <- qq_first_stage_collapsed(
  y = df_2019$lnwage,
  group_vec = df_2019$czone,
  taus = taus,
  pweights = df_2019$pweight
)

cz <- intersect(qq_2005$group_order, qq_2019$group_order)
f2005 <- qq_2005$fitted[qq_2005$group_order %in% cz, , drop = FALSE]
f2019 <- qq_2019$fitted[qq_2019$group_order %in% cz, , drop = FALSE]

f2005 <- f2005[order(rownames(f2005)), , drop = FALSE]
f2019 <- f2019[order(rownames(f2019)), , drop = FALSE]

delta_fitted <- f2019 - f2005

# -------------------------------------------------------------
# build X matrix
# -------------------------------------------------------------


# X aligned
X_df <- shocks %>%
  filter(czone %in% cz) %>%
  mutate(czone = factor(czone, levels = cz)) %>%
  arrange(czone)

X <- model.matrix(
  ~ IV_d_tradeusch_p1_2000_2012 +
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
  data = X_df
)

# weights aligned
pop_w <- pop_2000 %>%
  filter(czone %in% cz) %>%
  mutate(czone = factor(czone, levels = cz)) %>%
  arrange(czone) %>%
  pull(pop_share_2000)

# -------------------------------------------------------------
# second stage
# -------------------------------------------------------------
coeffs <- qq_second_stage_chetverikov(
  fitted_first = delta_fitted,
  X = X,
  taus = taus,
  weights = pop_w
)

beta <- coeffs$coef["IV_d_tradeusch_p1_2000_2012", ]
se   <- coeffs$se["IV_d_tradeusch_p1_2000_2012", ]

ci_low  <- beta - 1.96 * se
ci_high <- beta + 1.96 * se
# -------------------------------------------------------------
# coefficients plot
# -------------------------------------------------------------
coeff_plot_2005_2019 <- ggplot(
  data.frame(
    u = taus,
    beta = beta,
    ci_low = ci_low,
    ci_high = ci_high
  ),
  aes(x = u, y = beta)
) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Quantile (u)",
    y = "IV coefficient"
  )

# save html widget
ggsave(
  filename = "plots/chetverikov_replication/quant_coeffs_2005_2019.png",
  plot = coeff_plot_2005_2019,
  width = 8,
  height = 5,
  dpi = 300, 
  bg= "white"
)



