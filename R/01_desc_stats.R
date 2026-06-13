library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(sf)
library(tmap)
library(xtable)

cz_geo <- st_read("data/cz_geo.shp")
data <- readRDS("data/data.rds")

plot_dir <- "plots/descriptive_statistics"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# 0 get rid of outliers 
data <- data %>%
  filter(wage > 7.25) %>%
  filter(wage < 500)
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
data_w <- data[data$wage > 0]
cutoffs <- quantile(data_w$wage, probs = c(0.01, 0.99), na.rm = TRUE)
data_w <- data_w[wage >= cutoffs[1] & wage <= cutoffs[2]]
obs_year_wage <- data_w[, .(observations = .N), by = YEAR][order(YEAR)]
setnames(obs_year_full, "observations", "Observations")
setnames(obs_year_wage, "observations", "Observations with Wage")
obs_year <- merge(obs_year_full, obs_year_wage, by = "YEAR")

obs_year_tex <- xtable(
  obs_year,
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

