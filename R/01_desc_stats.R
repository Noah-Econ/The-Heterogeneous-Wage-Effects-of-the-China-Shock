library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(sf)
library(tmap)
library(xtable)
library(stargazer)
library(haven)
library(stringr)
library(readr)




cz_geo <- st_read("data/cz_geo.shp")
data <- readRDS("data/data.rds")
controls <- read_dta("data/ADH_control_vars.dta")
shocks <- read_dta("data/czone_exposure_by_period_v5_gh.dta")
ACS_pop <- read_dta("data/ACS_pop_emp_inc.dta")
cz_names_raw <- read_tsv(
  "data/CZ_Names.txt",
  col_names = c("LMA_CZ", "FIPS", "CountyName", "TotalPopulation", "LaborForce"),
  skip = 1,
  col_types = cols(.default = "c")
)

# controls for 1990
controls_1990 <- controls %>%
  filter(year == 1990)  %>%
  select(czone, l_shind_manuf_cbp, l_sh_popedu_c, l_sh_popfborn,
         l_sh_empl_f, l_sh_routine33, l_task_outsource,
         region, sh_65up_all, sh_4064_all, sh_0017_all, sh_00up_nw, statefip)
shocks <- shocks %>%
  left_join(controls_1990, by = "czone")

plot_dir <- "plots/descriptive_statistics"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

summary_wage <- data %>%
  summarise(
    mean_wage = mean(wage, na.rm = TRUE),
    p1_wage   = quantile(wage, 0.01, na.rm = TRUE),
    p99_wage  = quantile(wage, 0.99, na.rm = TRUE),
    min_wage  = min(wage, na.rm = TRUE),
    max_wage  = max(wage, na.rm = TRUE)
  )

summary_wage





# 1. obs per year 
obs_year_full <- data[, .(observations = .N), by = YEAR][order(YEAR)]

setnames(obs_year_full, "observations", "Observations")

obs_year_tex <- xtable(
  obs_year_full,
  caption = "Number of Observations per Year",
  label = "tab:obs_year"
)

print(
  obs_year_tex,
  type = "latex",
  file = file.path(plot_dir, "obs_year.tex"),
  include.rownames = FALSE,
  sanitize.text.function = identity
)

# 2. wage, education and age per year
desc_table_year <- data %>%
  group_by(YEAR) %>%
  summarise(
    wage = paste0(
      round(mean(wage, na.rm = TRUE),2),
      " (", round(sd(wage, na.rm = TRUE),2), ")"
    ),
    educ = paste0(
      round(mean(educ_years, na.rm = TRUE),2),
      " (", round(sd(educ_years, na.rm = TRUE),2), ")"
    ),
    age = paste0(
      round(mean(AGE, na.rm = TRUE),2),
      " (", round(sd(AGE, na.rm = TRUE),2), ")"
    )
  )

desc_table_year_tex <- xtable(
  desc_table_year,
  caption = "Descriptive Statistics by Year",
  label = "tab:desc_year"
)

print(
  desc_table_year_tex,
  type = "latex",
  file = file.path(plot_dir, "desc_table_year.tex"),
  include.rownames = FALSE,
  sanitize.text.function = identity
)


# map of shocks
cz_geo_shock <- cz_geo %>%
  left_join(shocks, by = "czone")

cz_geo_shock <- cz_geo_shock %>%
  filter(!is.na(d_tradeusch_p1_2000_2012))

cz_geo_shock$d_tradeusch_p1_2000_2012_prct <- cz_geo_shock$d_tradeusch_p1_2000_2012 * 100

breaks = quantile(
  cz_geo_shock$d_tradeusch_p1_2000_2012_prct,
  probs = seq(0, 1, 0.1),
  na.rm = TRUE
)


cz_geo_shock <- st_set_crs(cz_geo_shock, 4326)
cz_geo_shock <- st_transform(cz_geo_shock, 5070)
shocks_map <- tm_shape(cz_geo_shock) +
  tm_polygons(
    col = "d_tradeusch_p1_2000_2012_prct",
    palette= "YlOrRd",
    title = "Exposure", 
    breaks = breaks, 
    value.na = "no"
  ) +
  tm_layout(
    frame = FALSE,
    outer.margins = 0
  )

tmap_save(shocks_map, "plots/descriptive_statistics/shocks_map.png", width = 8, height = 6, units = "in", dpi = 300)

cz_vals <- cz_geo_shock %>%
  st_drop_geometry() %>%
  select(czone, d_tradeusch_p1_2000_2012_prct)

# highest exposure
top_cz <- cz_vals %>%
  slice_max(d_tradeusch_p1_2000_2012_prct, n = 1)

# lowest exposure
bottom_cz <- cz_vals %>%
  slice_min(d_tradeusch_p1_2000_2012_prct, n = 1)

top_cz #Adams County, MS
bottom_cz # Phillips County, KS,  Smith County, KS

#descriptive statistics plots of the variables used
cz_vars <- shocks %>%
  select(
    d_tradeusch_p1_2000_2012,
    d_tradeotch_p1_lag_2000_2012,
    l_shind_manuf_cbp,
    l_sh_popedu_c,
    l_sh_popfborn,
    l_sh_empl_f,
    l_sh_routine33,
    l_task_outsource,
    sh_65up_all,
    sh_4064_all,
    sh_0017_all,
    sh_00up_nw
  ) %>%
  as.data.frame()


ind_vars <- data %>%
  select(wage, lnwage) %>%
  filter(!is.na(wage), !is.na(lnwage)) %>%
  as.data.frame()

# CZ-level table
stargazer(
  cz_vars,
  type = "latex",
  out = "plots/descriptive_statistics/desc_stats_cz.tex",
  title = "Descriptive Statistics: Commuting Zone Variables",
  label = "tab:desc_cz",
  summary.stat = c("mean", "sd", "p25", "median", "p75"),
  digits = 3,
  covariate.labels = c(
    "Trade exposure (US)",
    "Trade exposure (instrument)",
    "Manufacturing employment share",
    "College education share",
    "Foreign-born share",
    "Female employment share",
    "Routine task intensity",
    "Outsourcing task intensity",
    "Population share 65+",
    "Population share 40-64",
    "Population share 0-17",
    "Non-white share"
  )
)

# Individual-level wage table
stargazer(
  ind_vars,
  type = "latex",
  out = "plots/descriptive_statistics/desc_stats_wages.tex",
  title = "Descriptive Statistics: Wage Variables",
  label = "tab:desc_wages",
  summary.stat = c("mean", "sd","p25", "median", "p75"),
  digits = 3,
  covariate.labels = c(
    "Hourly wage (USD, 2010)",
    "Log hourly wage"
  )
)


# population in 2000
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



# county-level rows are the ones with a non-missing FIPS
# (the "Market Area Total" rows have FIPS blank and use the 3-digit CZ code directly)
cz_counties <- cz_names_raw %>%
  filter(!is.na(FIPS) & FIPS != "") %>%
  mutate(
    czone           = as.integer(LMA_CZ),        # keep full code, e.g. "00301" -> 301
    CountyName      = str_remove_all(CountyName, '"'),
    TotalPopulation = as.numeric(TotalPopulation)
  )

cz_name_lookup <- cz_counties %>%
  group_by(czone) %>%
  slice_max(TotalPopulation, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(czone, CZ_Name = CountyName)

top10_pop <- pop_2000 %>%
  slice_max(population_2000, n = 10) %>%
  transmute(
    CZ = as.integer(czone),
    Population = round(population_2000),
    `Pop. Share (%)` = round(pop_share_2000 * 100, 2)
  )

top10_pop_named <- top10_pop %>%
  left_join(cz_name_lookup, by = c("CZ" = "czone")) %>%
  relocate(CZ_Name, .after = CZ)

top10_pop_named

top10_pop_stargazer <- top10_pop_named %>%
  mutate(
    Population = format(Population, big.mark = ",", scientific = FALSE)
  ) %>%
  rename(
    `CZ` = CZ,
    `Name` = CZ_Name,
    `Population` = Population,
    `Pop. Share (%)` = `Pop. Share (%)`
  ) %>%
  as.data.frame()

stargazer(
  top10_pop_stargazer,
  type = "latex",
  summary = FALSE,
  rownames = FALSE,
  out = "plots/descriptive_statistics/top10_pop.tex",
  title = "Top 10 Commuting Zones by Working-Age Population Weight (2000)",
  label = "tab:top10_pop"
)
