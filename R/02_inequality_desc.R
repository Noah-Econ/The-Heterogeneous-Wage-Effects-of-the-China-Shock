library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(plotly)
library(tmap)
library(DescTools)
library(Hmisc)
library(modi)
library(magick)
library(webshot2)
library(htmltools)
library(jsonlite)
#load quantile on quantiles functions
source("R/00_2_functions.R")

cz_geo <- st_read("data/cz_geo.shp")
data <- readRDS("data/data.rds")

# 0.0 Exclude Outliers
data <- data %>%
  filter(wage > 7.25) %>%
  filter(wage < 500)

######################################################################################################################################################################
## 1. Between Group Inequality (Median Wage)
######################################################################################################################################################################
dir.create("plots/evolution_wage_inequality/between", showWarnings = FALSE)

# 1.1 Median Wage Plot by Gender
med_wage <- data %>%
  group_by(YEAR, SEX) %>%
  summarise(
    median_wage = median(wage, na.rm = TRUE, weights = pweight),
    .groups = "drop"
  ) %>%
  mutate(
    SEX = factor(SEX, levels = c(1, 2), labels = c("Men", "Women"))
  )

med_wage_gender <- ggplot(
  med_wage,
  aes(x = YEAR, y = median_wage, color = SEX)
) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("cornflowerblue", "deeppink3")) +
  labs(
    title = "Median Hourly Wage by Gender",
    x = "Year",
    y = "Median hourly wage (2010 USD)",
    color = NULL
  ) +
  coord_cartesian(ylim = c(10, 25)) +
  theme_minimal() +
  geom_text_repel(
    data = med_wage %>%
      group_by(SEX) %>%
      filter(YEAR == min(YEAR) | YEAR == max(YEAR)),
    aes(label = round(median_wage, 1), color = SEX),
    size = 3,
    nudge_y = 0.5, 
    show.legend = FALSE
  )

ggsave(
  filename = "plots/evolution_wage_inequality/between/median_hourly_wage_gender.pdf",
  plot = med_wage_gender,
  width = 8,
  height = 6,
  dpi = 300
)

# 1.2 Median Wage Plot by Industry
med_ind <- data %>%
  group_by(YEAR, naics_sector) %>%
  summarise(
    median_wage = median(wage, na.rm = TRUE, weights = pweight),
    .groups = "drop"
  )

med_wage_ind <- ggplot(
  med_ind,
  aes(x = YEAR, y = median_wage)
) +
  geom_line(linewidth = 1.1, color = "cornflowerblue") + 
  facet_wrap(~ naics_sector, scales = "fixed") +    
  labs(
    title = "Median Hourly Wage by Industry",
    x = "Year",
    y = "Median hourly wage (2010 USD)"
  ) +
  coord_cartesian(ylim = c(10, 30)) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 9, face = "bold")
  )

ggsave(
  filename = "plots/evolution_wage_inequality/between/median_hourly_wage_ind.pdf",
  plot = med_wage_ind,
  width = 8,
  height = 6,
  dpi = 300
)

# 1.3 richest and poorest industry by median wage
ind_rank <- med_ind %>%
  group_by(naics_sector) %>%
  summarise(
    avg_median_wage = mean(median_wage, na.rm = TRUE),
    .groups = "drop"
  )

richest_ind <- ind_rank %>%
  slice_max(avg_median_wage, n = 1) %>%
  pull(naics_sector)

poorest_ind <- ind_rank %>%
  slice_min(avg_median_wage, n = 1) %>%
  pull(naics_sector)

extreme_ind <- med_ind %>%
  filter(naics_sector %in% c(richest_ind, poorest_ind)) %>%
  mutate(
    group = ifelse(
      naics_sector == richest_ind,
      paste0("Professional, Scientific, and Technical Services"),
      paste0("Accommodation and Food Services")
    )
  )

rich_poor_ind_med_wage <- ggplot(extreme_ind, aes(x = YEAR, y = median_wage, color = group)) +
  geom_line(linewidth = 1.3) +
  geom_point(size = 2) +
  scale_color_manual(values = c("deeppink3", "cornflowerblue")) +
  labs(
    title = "Median Wage: Richest vs. Poorest Industry",
    x = "Year",
    y = "Median hourly wage",
    color = ""
  ) +
  coord_cartesian(ylim = c(10, 30)) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  geom_text_repel(
    data = extreme_ind %>%
      group_by(group) %>%
      filter(YEAR == min(YEAR) | YEAR == max(YEAR)),
    aes(label = round(median_wage, 1), color = group),
    size = 3,
    nudge_y = 0.5,
    show.legend = FALSE
  )

ggsave(
  filename = "plots/evolution_wage_inequality/between/rich_poor_industry_med_wage.pdf",
  plot = rich_poor_ind_med_wage,
  width = 6,
  height = 6,
  dpi = 300
)

# 1.4 animated map median wage 
median_wage_cz <- data %>%
  filter(YEAR >= 2005, YEAR <= 2019) %>%
  filter(!is.na(czone)) %>%
  filter(!STATEFIP %in% c(2, 15, 72)) %>%
  group_by(YEAR, czone) %>%
  summarise(
    median_wage = wtd.quantile(wage, weights = pweight, probs = 0.5, na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  )

pmap_med_wage <- cz_geo %>%
  left_join(median_wage_cz, by = "czone") %>%
  filter(!is.na(median_wage))

st_crs(pmap_med_wage) <- 4326
pmap_med_wage <- st_transform(pmap_med_wage, 5070)

tmap_mode("plot")

wage_limits <- range(pmap_med_wage$median_wage, na.rm = TRUE)

pmap_med_wage_anim <- tm_shape(pmap_med_wage) +
  tm_polygons(
    col = "median_wage",
    palette = "plasma",
    style = "cont",
    limits = wage_limits,
    title = "Weighted Median Wage"
  ) +
  tm_facets_pagewise(
    by = "YEAR",
    labeller = function(x) paste("Year", x)
  ) +
  tm_layout(
    frame = FALSE,
    outer.margins = 0,
    legend.outside = TRUE
  )

tmap_animation(
  pmap_med_wage_anim,
  filename = "plots/evolution_wage_inequality/between/median_wage_cz_2005_2019.gif",
  width = 1200,
  height = 800,
  delay = 120
)

# 1.5 richest and poorest CZ by median wage
cz_rank <- pmap_med_wage %>%
  group_by(czone) %>%
  summarise(
    avg_median_wage = mean(median_wage, na.rm = TRUE),
    .groups = "drop"
  )

richest_cz <- cz_rank %>%
  slice_max(avg_median_wage, n = 1) %>%
  pull(czone)

poorest_cz <- cz_rank %>%
  slice_min(avg_median_wage, n = 1) %>%
  pull(czone)

extreme_cz <- pmap_med_wage %>%
  filter(czone %in% c(richest_cz, poorest_cz)) %>%
  mutate(
    group = ifelse(czone == richest_cz, "Arlington County CZ", "Howell County CZ")
  )

rich_poor_cz_med_wage <- ggplot(extreme_cz, aes(x = YEAR, y = median_wage, color = group)) +
  geom_line(linewidth = 1.3) +
  geom_point(size = 2) +
  scale_color_manual(values = c(
    "Arlington County CZ" = "cornflowerblue",
    "Howell County CZ" = "deeppink3"
  )) +
  labs(
    title = "Median Wage: Richest vs. Poorest Commuting Zone",
    x = "Year",
    y = "Median hourly wage",
    color = ""
  ) +
  coord_cartesian(ylim = c(10, 30)) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  geom_text_repel(
    data = extreme_cz %>%
      group_by(group) %>%
      filter(YEAR == min(YEAR) | YEAR == max(YEAR)),
    aes(label = round(median_wage, 1), color = group),
    size = 3,
    nudge_y = 0.5,
    show.legend = FALSE
  )

ggsave(
  filename = "plots/evolution_wage_inequality/between/rich_poor_cz_med_wage.pdf",
  plot = rich_poor_cz_med_wage,
  width = 6,
  height = 6,
  dpi = 300
)


######################################################################################################################################################################
## 2. Within Group Inequality (Gini Coefficient of Wage)
######################################################################################################################################################################
dir.create("plots/evolution_wage_inequality/within", showWarnings = FALSE)

# 2.1.1 Gini in the US over time 
gini_wage_us <- data %>%
  filter(YEAR >= 2005, YEAR <= 2019) %>%
  filter(!is.na(wage), !is.na(pweight)) %>%
  group_by(YEAR) %>%
  summarise(
    gini = Gini(wage, weights = pweight, unbiased = TRUE),
    n_obs = n(),
    .groups = "drop"
  )

gini_plot_US <- ggplot(gini_wage_us, aes(x = YEAR, y = gini)) +
  geom_line(linewidth = 1.2, color = "deeppink3") +
  geom_point(size = 2, color = "deeppink3") +
  labs(
    title = "Gini Coefficient of Wage",
    x = "Year",
    y = "Gini Coefficient"
  ) +
  coord_cartesian(ylim = c(0.3, 0.45)) +
  theme_minimal()

ggsave(
  filename = "plots/evolution_wage_inequality/within/gini_us_2005_2019.pdf",
  plot = gini_plot_US,
  width = 10,
  height = 6,
  dpi = 300
)

# 2.1.2. 90/10 Ratio in the US over time
p9010_wage_us <- data %>%
  filter(YEAR >= 2005, YEAR <= 2019) %>%
  filter(!STATEFIP %in% c(2, 15, 72)) %>%
  filter(!is.na(wage), !is.na(pweight)) %>%
  group_by(YEAR) %>%
  summarise(
    p10 = wtd.quantile(wage, weights = pweight, probs = 0.10, na.rm = TRUE),
    p90 = wtd.quantile(wage, weights = pweight, probs = 0.90, na.rm = TRUE),
    ratio_90_10 = p90 / p10,
    .groups = "drop"
  )

p9010_plot_US <- ggplot(p9010_wage_us, aes(x = YEAR, y = ratio_90_10)) +
  geom_line(linewidth = 1.2, color = "cornflowerblue") +
  geom_point(size = 2, color = "cornflowerblue") +
  labs(
    title = "P90 / P10 Ratio of Wage",
    x = "Year",
    y = "P90 / P10 Ratio"
  ) +
  coord_cartesian(ylim = c(3, 6.5)) +
  theme_minimal()

ggsave(
  filename = "plots/evolution_wage_inequality/within/p9010_us_2005_2019.pdf",
  plot = p9010_plot_US,
  width = 10,
  height = 6,
  dpi = 300
)
 
# 2.2 Regional Analysis
# 2.2.1 Gini Regional
gini_wage_cz <- data %>%
  filter(YEAR >= 2005, YEAR <= 2019) %>%
  filter(!is.na(czone)) %>%
  filter(!is.na(wage)) %>%
  filter(!STATEFIP %in% c(2, 15, 72)) %>% 
  group_by(YEAR, czone) %>%
  summarise(
    gini = Gini(wage, weights = pweight, unbiased = T, conf.level = NA),
    n_obs = n(),
    .groups = "drop"
  )


pmap_gini_wage <- cz_geo %>%
  left_join(gini_wage_cz, by = "czone") %>%
  filter(!is.na(gini))

st_crs(pmap_gini_wage) <- 4326
pmap_gini_wage <- st_transform(pmap_gini_wage, 5070)

tmap_mode("plot")

wage_limits <- range(pmap_gini_wage$gini, na.rm = TRUE)


pmap_gini_wage_anim <- tm_shape(pmap_gini_wage) +
  tm_polygons(
    col = "gini",
    palette = "plasma",
    style = "cont",
    limits = wage_limits,
    title = "Weighted Gini"
  ) +
  tm_facets_pagewise(
    by = "YEAR",
    labeller = function(x) paste("Year", x)
  ) +
  tm_layout(
    frame = FALSE,
    outer.margins = 0,
    legend.outside = TRUE
  )


tmap_animation(
  pmap_gini_wage_anim,
  filename = "plots/evolution_wage_inequality/within/gini_wage_cz_2005_2019.gif",
  width = 1200,
  height = 800,
  delay = 120
)

# 2008 and 2019 to show in thesis
pmap_gini_2008 <- pmap_gini_wage %>% filter(YEAR == 2008)
pmap_gini_2019 <- pmap_gini_wage %>% filter(YEAR == 2019)

gini_map_2008 <- tm_shape(pmap_gini_2008) +
  tm_polygons(
    col = "gini",
    palette = "plasma",
    style = "cont",
    limits = wage_limits,
    title = "Weighted Gini"
  ) +
  tm_layout(
    title = "Year 2008",
    title.size = 1,
    title.position = c("center", "top"),
    frame = FALSE,
    legend.outside = TRUE
  )

gini_map_2019 <- tm_shape(pmap_gini_2019) +
  tm_polygons(
    col = "gini",
    palette = "plasma",
    style = "cont",
    limits = wage_limits,
    title = "Weighted Gini"
  ) +
  tm_layout(
    title = "Year 2019",
    title.size = 1,
    title.position = c("center", "top"),
    frame = FALSE,
    legend.outside = TRUE
  )

tmap_save(gini_map_2008, "plots/evolution_wage_inequality/within/gini_wage_cz_2008.png", width = 8, height = 6, units = "in", dpi = 300)
tmap_save(gini_map_2019, "plots/evolution_wage_inequality/within/gini_wage_cz_2019.png", width = 8, height = 6, units = "in", dpi = 300)
#2.2.2 p9010 over time over regions

p9010_wage_cz <- data %>%
  filter(YEAR >= 2005, YEAR <= 2019) %>%
  filter(!is.na(czone)) %>%
  filter(!is.na(wage), !is.na(pweight)) %>%
  filter(!STATEFIP %in% c(2, 15, 72)) %>% 
  group_by(YEAR, czone) %>%
  summarise(
    p10 = wtd.quantile(wage, weights = pweight, probs = 0.10, na.rm = TRUE),
    p90 = wtd.quantile(wage, weights = pweight, probs = 0.90, na.rm = TRUE),
    ratio_90_10 = p90 / p10,
    n_obs = n(),
    .groups = "drop"
  )

pmap_p9010_wage <- cz_geo %>%
  left_join(p9010_wage_cz, by = "czone") %>%
  filter(!is.na(ratio_90_10))

st_crs(pmap_p9010_wage) <- 4326
pmap_p9010_wage <- st_transform(pmap_p9010_wage, 5070)

tmap_mode("plot")

ratio_limits <- range(pmap_p9010_wage$ratio_90_10, na.rm = TRUE)

pmap_p9010_anim <- tm_shape(pmap_p9010_wage) +
  tm_polygons(
    col = "ratio_90_10",
    palette = "plasma",
    style = "cont",
    limits = ratio_limits,
    title = "P90/P10 Ratio"
  ) +
  tm_facets_pagewise(
    by = "YEAR",
    labeller = function(x) paste("Year", x)
  ) +
  tm_layout(
    frame = FALSE,
    outer.margins = 0,
    legend.outside = TRUE
  )

tmap_animation(
  pmap_p9010_anim,
  filename = "plots/evolution_wage_inequality/within/p9010_wage_cz_2005_2019.gif",
  width = 1200,
  height = 800,
  delay = 120
)


# 2.3 Evolution of wage inequality over time across industries
# 2.3.1 Gini Coefficient across industries 
gini_wage_industry <- data %>%
  filter(YEAR >= 2005, YEAR <= 2019) %>%
  filter(!is.na(wage), !is.na(pweight)) %>%
  filter(!is.na(naics_sector), naics_sector != "Unemployed") %>%
  group_by(YEAR, naics_sector) %>%
  summarise(
    gini = Gini(wage, weights = pweight, unbiased = TRUE),
    n_obs = n(),
    .groups = "drop"
  )


gini_plot_ind_facet <- ggplot(gini_wage_industry, 
                              aes(x = YEAR, y = gini, color = naics_sector)) +
  geom_line() +
  facet_wrap(~ naics_sector) +
  scale_color_manual(values = rep("deeppink3", length(unique(gini_wage_industry$naics_sector)))) +
  coord_cartesian(ylim = c(0.3, 0.45)) +
  theme_minimal() +
  labs(
    y = "Gini Coefficient",
    x = "Year"
  )+
  theme(legend.position = "none")


ggsave(
  filename = "plots/evolution_wage_inequality/within/gini_industry_us_facet_2005_2019.pdf",
  plot = gini_plot_ind_facet,
  width = 10,
  height = 6,
  dpi = 300
)



#2.3.2 Robustness test with 90/10 ratio
p9010_wage_industry <- data %>%
  filter(YEAR >= 2005, YEAR <= 2019) %>%
  filter(!STATEFIP %in% c(2, 15, 72)) %>%
  filter(!is.na(wage), !is.na(pweight)) %>%
  filter(!is.na(naics_sector), naics_sector != "Unemployed") %>%
  group_by(YEAR, naics_sector) %>%
  summarise(
    p10 = wtd.quantile(wage, weights = pweight, probs = 0.10, na.rm = TRUE),
    p90 = wtd.quantile(wage, weights = pweight, probs = 0.90, na.rm = TRUE),
    ratio_90_10 = p90 / p10,
    n_obs = n(),
    .groups = "drop"
  )

p9010_plot_ind_facet <- ggplot(p9010_wage_industry, 
                               aes(x = YEAR, y = ratio_90_10, color = naics_sector)) +
  geom_line() +
  facet_wrap(~ naics_sector) +
  scale_color_manual(values = rep("deeppink3", length(unique(p9010_wage_industry$naics_sector)))) +
  coord_cartesian(ylim = c(3, 6.5)) +  
  labs(
    y = "P90 / P10 Ratio",
    x = "Year"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(
  filename = "plots/evolution_wage_inequality/within/p9010_industry_us_facet_2005_2019.pdf",
  plot = p9010_plot_ind_facet,
  width = 10,
  height = 6,
  dpi = 300
)


######################################################################################################################################################################
## 3. Between and Within Group Inequality (Quantile on Quantiles)
######################################################################################################################################################################
dir.create("plots/evolution_wage_inequality/between_within", showWarnings = FALSE)

# 3.0 define variables 
taus  <- round(seq(0.01, 0.99, by = 0.01),2)

z_min_max <- c(
  q01 = quantile(data$wage, 0.01, na.rm = TRUE),
  q99 = quantile(data$wage, 0.99, na.rm = TRUE)
)
y_min_max <- c(
  q01 = quantile(data$wage, 0.01, na.rm = TRUE),
  80
)
years <- 2005:2019

#### 3.1 Regional Analysis (CZ)

dir.create("plots/evolution_wage_inequality/between_within/qq_surface_frames_cz", showWarnings = FALSE)
dir.create("plots/evolution_wage_inequality/between_within/qq_slices_frames_cz", showWarnings = FALSE)

qq_results_cz <- list()

for (yr in years) {
  
  cat("Processing year:", yr, "\n")
  
  # --- subset data
  data_year <- data %>%
    filter(YEAR == yr)
  
  # --- compute QQ matrix
  qq_mat <- make_qq_matrix(
    data = data_year,
    y = "wage",
    group = "czone",
    weight = "pweight",
    taus = taus
  )
  qq_results_cz[[as.character(yr)]] <- qq_mat
  
  # --- plot surface 
  p_surface <- plot_qq_surface(
    qq_mat,
    taus = taus,
    tplot = 3,
    mycolors = NULL,
    z_label = "Wage",
    group = "CZ", 
    zlim = z_min_max, 
    year = yr
  )
  
  # save as HTML
  htmlwidgets::saveWidget(
    p_surface,
    file = paste0("plots/evolution_wage_inequality/between_within/qq_surface_frames_cz/surface_cz_", yr, ".html"),
    selfcontained = TRUE
  )
  
  # --- slices 
  p_slices <- plot_qq_slices(
    qq_mat,
    taus = taus,
    z_label = "Wage",
    group = "CZ",
    ylim = y_min_max
  ) +
    ggtitle(paste("Year:", yr))
  
  # --- save slices
  
  ggsave(
    filename = paste0("plots/evolution_wage_inequality/between_within/qq_slices_frames_cz/slices_cz_", yr, ".png"),
    plot = p_slices,
    width = 8,
    height = 5,
    dpi = 300,
    bg = "white"
  )
}

saveRDS(
  qq_results_cz,
  "plots/evolution_wage_inequality/between_within/all_qq_surfaces_cz.rds"
)

## produce gifs of slices
slice_files <- list.files("plots/evolution_wage_inequality/between_within/qq_slices_frames_cz", full.names = TRUE)

slice_gif <- image_read(slice_files) %>%
  image_animate(fps = 0.25)

image_write(slice_gif, "plots/evolution_wage_inequality/between_within/qq_slices_cz_2005_2019.gif")

slice_files <- sort(list.files("plots/evolution_wage_inequality/between_within/qq_slices_frames_cz", full.names = TRUE))

years <- gsub(".*_(\\d+)\\.png", "\\1", slice_files)

html <- tags$html(
  tags$body(
    tags$div(
      tags$input(
        type = "range",
        min = 1,
        max = length(slice_files),
        value = 1,
        id = "slider"
      ),
      
      tags$h3(id = "year"),
      
      tags$img(id = "img", src = slice_files[1], style = "width:800px;")
    ),
    
    tags$script(HTML(sprintf("
      const images = %s;
      const years = %s;

      const slider = document.getElementById('slider');
      const img = document.getElementById('img');
      const year = document.getElementById('year');

      function update() {
        let i = slider.value - 1;
        img.src = images[i];
        year.innerHTML = 'Year: ' + years[i];
      }

      slider.oninput = update;
      update();
    ",
                             toJSON(slice_files),
                             toJSON(years)
    )))
  )
)

save_html(html, "plots/evolution_wage_inequality/between_within/qq_slices_cz_2005_2019.html")

# Animated Plotly Surface from SAVED QQ MATRICES

years <- 2005:2019

# load saved results
qq_results_cz <- readRDS(
  "plots/evolution_wage_inequality/between_within/all_qq_surfaces_cz.rds"
)

z_list <- list()
x_list <- list()
y_list <- list()

# thinner grid for speed
idx <- seq(1, length(taus), by = 3)

for (yr in years) {
  
  cat("Preparing year:", yr, "\n")
  
  qq_mat <- qq_results_cz[[as.character(yr)]]
  
  z_list[[as.character(yr)]] <-
    qq_mat[idx, idx, drop = FALSE]
  
  x_list[[as.character(yr)]] <- taus[idx]
  y_list[[as.character(yr)]] <- taus[idx]
}

# initial surface
first_year <- as.character(years[1])

p <- plot_ly(
  type = "surface",
  
  x = x_list[[first_year]],
  y = y_list[[first_year]],
  z = z_list[[first_year]],
  
  cmin = z_min_max[1],
  cmax = z_min_max[2],
  
  colorscale = NULL
)

# frames
frames <- lapply(years, function(yr) {
  
  list(
    name = as.character(yr),
    
    data = list(
      list(
        type = "surface",
        
        x = x_list[[as.character(yr)]],
        y = y_list[[as.character(yr)]],
        z = z_list[[as.character(yr)]],
        
        cmin = z_min_max[1],
        cmax = z_min_max[2],
        
        colorscale = NULL
      )
    )
  )
})

p$x$frames <- frames

# layout
p <- p %>%
  layout(
    
    title = "",
    
    scene = list(
      xaxis = list(
        title = "Between CZ",
        autorange = "reversed"
      ),
      
      yaxis = list(
        title = "Within CZ",
        autorange = "reversed"
      ),
      
      zaxis = list(
        title = "Wage",
        range = z_min_max
      )
    ),
    
    sliders = list(
      list(
        active = 0,
        
        currentvalue = list(
          prefix = "Year: "
        ),
        
        steps = lapply(years, function(yr) {
          
          list(
            label = as.character(yr),
            
            method = "animate",
            
            args = list(
              list(as.character(yr)),
              
              list(
                mode = "immediate",
                
                frame = list(
                  duration = 0,
                  redraw = TRUE
                ),
                
                transition = list(
                  duration = 0
                )
              )
            )
          )
        })
      )
    )
  )

# save
htmlwidgets::saveWidget(
  p,
  "plots/evolution_wage_inequality/between_within/qq_surface_cz_animation_2005_2019.html",
  selfcontained = TRUE
)


# 3.2 Pons analysis on industry level 

dir.create("plots/evolution_wage_inequality/between_within/qq_surface_frames_ind", showWarnings = FALSE)
dir.create("plots/evolution_wage_inequality/between_within/qq_slices_frames_ind", showWarnings = FALSE)

years <- 2005:2019

qq_results_ind <- list()

for (yr in years) {
  
  cat("Processing year:", yr, "\n")
  
  # --- subset data
  data_year <- data %>%
    filter(YEAR == yr)
  
  # --- compute QQ matrix
  qq_mat <- make_qq_matrix(
    data = data_year,
    y = "wage",
    group = "naics_sector",
    weight = "pweight",
    taus = taus
  )
  
  qq_results_ind[[as.character(yr)]] <- qq_mat
  
  # --- plot surface 
  p_surface <- plot_qq_surface(
    qq_mat,
    taus = taus,
    tplot = 3,
    mycolors = NULL,
    z_label = "Wage",
    group = "Industries", 
    zlim = z_min_max, 
    year = yr
  )
  
  # save as HTML
  htmlwidgets::saveWidget(
    p_surface,
    file = paste0("plots/evolution_wage_inequality/between_within/qq_surface_frames_ind/surface_ind_", yr, ".html"),
    selfcontained = TRUE
  )
  
  # --- slices 
  p_slices <- plot_qq_slices(
    qq_mat,
    taus = taus,
    z_label = "Wage",
    group = "Industries",
    ylim = y_min_max
  ) +
    ggtitle(paste("Year:", yr))
  
  # --- save slices
  
  ggsave(
    filename = paste0("plots/evolution_wage_inequality/between_within/qq_slices_frames_ind/slices_ind_", yr, ".png"),
    plot = p_slices,
    width = 8,
    height = 5,
    dpi = 300,
    bg = "white"
  )
}

saveRDS(
  qq_results_ind,
  "plots/evolution_wage_inequality/between_within/all_qq_surfaces_ind.rds"
)

## produce gifs of slices
slice_files <- list.files("plots/evolution_wage_inequality/between_within/qq_slices_frames_ind", full.names = TRUE)

slice_gif <- image_read(slice_files) %>%
  image_animate(fps = 0.25)

image_write(slice_gif, "plots/evolution_wage_inequality/between_within/qq_slices_ind_2005_2019.gif")

slice_files <- sort(list.files("plots/evolution_wage_inequality/between_within/qq_slices_frames_ind", full.names = TRUE))

years <- gsub(".*_(\\d+)\\.png", "\\1", slice_files)

html <- tags$html(
  tags$body(
    tags$div(
      tags$input(
        type = "range",
        min = 1,
        max = length(slice_files),
        value = 1,
        id = "slider"
      ),
      
      tags$h3(id = "year"),
      
      tags$img(id = "img", src = slice_files[1], style = "width:800px;")
    ),
    
    tags$script(HTML(sprintf("
      const images = %s;
      const years = %s;

      const slider = document.getElementById('slider');
      const img = document.getElementById('img');
      const year = document.getElementById('year');

      function update() {
        let i = slider.value - 1;
        img.src = images[i];
        year.innerHTML = 'Year: ' + years[i];
      }

      slider.oninput = update;
      update();
    ",
                             toJSON(slice_files),
                             toJSON(years)
    )))
  )
)

save_html(html, "plots/evolution_wage_inequality/between_within/qq_slices_ind_2005_2019.html")


# Animated Plotly Surface from SAVED QQ MATRICES

years <- 2005:2019

# load saved results
qq_results_ind <- readRDS(
  "plots/evolution_wage_inequality/between_within/all_qq_surfaces_cz.rds"
)

z_list <- list()
x_list <- list()
y_list <- list()

# thinner grid for speed
idx <- seq(1, length(taus), by = 3)

for (yr in years) {
  
  cat("Preparing year:", yr, "\n")
  
  qq_mat <- qq_results_ind[[as.character(yr)]]
  
  z_list[[as.character(yr)]] <-
    qq_mat[idx, idx, drop = FALSE]
  
  x_list[[as.character(yr)]] <- taus[idx]
  y_list[[as.character(yr)]] <- taus[idx]
}

# initial surface
first_year <- as.character(years[1])

p <- plot_ly(
  type = "surface",
  
  x = x_list[[first_year]],
  y = y_list[[first_year]],
  z = z_list[[first_year]],
  
  cmin = z_min_max[1],
  cmax = z_min_max[2],
  
  colorscale = NULL
)

# frames
frames <- lapply(years, function(yr) {
  
  list(
    name = as.character(yr),
    
    data = list(
      list(
        type = "surface",
        
        x = x_list[[as.character(yr)]],
        y = y_list[[as.character(yr)]],
        z = z_list[[as.character(yr)]],
        
        cmin = z_min_max[1],
        cmax = z_min_max[2],
        
        colorscale = NULL
      )
    )
  )
})

p$x$frames <- frames

# layout
p <- p %>%
  layout(
    
    title = "",
    
    scene = list(
      xaxis = list(
        title = "Between Industries",
        autorange = "reversed"
      ),
      
      yaxis = list(
        title = "Within Industries",
        autorange = "reversed"
      ),
      
      zaxis = list(
        title = "Wage",
        range = z_min_max
      )
    ),
    
    sliders = list(
      list(
        active = 0,
        
        currentvalue = list(
          prefix = "Year: "
        ),
        
        steps = lapply(years, function(yr) {
          
          list(
            label = as.character(yr),
            
            method = "animate",
            
            args = list(
              list(as.character(yr)),
              
              list(
                mode = "immediate",
                
                frame = list(
                  duration = 0,
                  redraw = TRUE
                ),
                
                transition = list(
                  duration = 0
                )
              )
            )
          )
        })
      )
    )
  )

# save
htmlwidgets::saveWidget(
  p,
  "plots/evolution_wage_inequality/between_within/qq_surface_ind_animation_2005_2019.html",
  selfcontained = TRUE
)

