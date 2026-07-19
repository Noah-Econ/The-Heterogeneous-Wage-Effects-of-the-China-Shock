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



######################################################################################################################################################################
## 1. Control function estimation following Imbens & Newey 2009 
######################################################################################################################################################################
# use controls 1990 Autor
controls_1990 <- controls %>%
  filter(year == 1990)
# merge controls and shocks by CZ for the first stage IV 
shocks <- shocks %>%
  left_join(controls_1990, by = "czone")

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

shocks <- estimate_cf_distribution(
  data = shocks,
  x_var = "d_tradeusch_p1_2000_2012",
  z_var = "d_tradeotch_p1_lag_2000_2012",
  control_vars = control_vars
)

# check estimated V 
V_hat_hist <- ggplot(data.frame(V_hat = shocks$V_hat), aes(V_hat)) +
  geom_histogram(binwidth = 0.05, color = "black", fill = "grey80") +
  theme_minimal()
ggsave(
  filename = "plots/descriptive_statistics/V_hat_hist.png",
  plot = V_hat_hist,
  width = 8,
  height = 5,
  dpi = 300,
  bg = "white"
)

formula_first_stage_IV_2 <- d_tradeusch_p1_2000_2012 ~ d_tradeotch_p1_lag_2000_2012 + l_shind_manuf_cbp + l_sh_popedu_c + l_sh_popfborn + l_sh_empl_f + l_sh_routine33 + l_task_outsource + region +  sh_65up_all + sh_4064_all + sh_0017_all + sh_00up_nw

# run first stage IV
first_stage <- lm(formula_first_stage_IV_2, data = shocks)
cor(
  shocks$V_hat,
  first_stage$residuals
)
plot(
  shocks$V_hat,
  first_stage$residuals,
  pch = 19,
  cex = 0.6,
  xlab = "First-stage residual",
  ylab = "Estimated control function"
)
formula_cf <- update(
  formula_first_stage_IV_2,
  . ~ . + V_hat
)
m2 <- lm(formula_cf,
  data = shocks
)
summary(first_stage)
summary(m2)
m2_clean <- lm(formula(m2), data = model.frame(m2))

stargazer(
  m2_clean,
  type = "latex",
  out = "plots/descriptive_statistics/first_stage_IV_CF.tex",
  title = "First-Stage Regression with Control Function",
  keep = c("d_tradeotch_p1_lag_2000_2012", "V_hat"),
  digits = 3,
  star.cutoffs = c(0.10, 0.05, 0.01),
  omit.stat = c("f", "ser"),
  dep.var.labels = c("d\\_tradeusch\\_p1\\_2000\\_2012")
)


#################################################################
## QoQ-control-function-IV appraoch for every year (2005–2019)
#################################################################

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

dir.create(
  "plots/iv_qoq/control_function_approach/slices",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "plots/iv_qoq/control_function_approach/surfaces",
  recursive = TRUE,
  showWarnings = FALSE
)

# parameters 
# common limits
z_min_max <- c(-7, 9)
y_min_max <- c(-7, 9)

years <- 2005:2019

# quantiles
taus <- round(seq(0.1, 0.9, by = 0.1), 1)

# results for every year 

results_iv <- list()

for (yr in years) {
  
  cat("Processing year:", yr, "\n")
  
  # -------------------------------------------------------------
  # subset year
  # -------------------------------------------------------------
  
  df_year <- data %>%
    filter(YEAR == yr)
  
  # merge shocks + controls
  df_year <- shocks %>%
    left_join(df_year, by = "czone")
  
  # -------------------------------------------------------------
  # first stage collapsed QQ
  # -------------------------------------------------------------
  
  qq_first <- qq_first_stage_collapsed(
    y = df_year$lnwage,
    group_vec = df_year$czone,
    taus = taus, 
    pweights = df_year$pweight #weighted by sample weights
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
  
  # weights aligned
  pop_w <- data.frame(czone = qq_first$group_order) %>%
    left_join(pop_2000, by = "czone") %>%
    pull(pop_share_2000)
  
  # -------------------------------------------------------------
  # second stage
  # -------------------------------------------------------------
  
  coeffs <- qq_second_stage_group(
    fitted_first = qq_first$fitted,
    X = X,
    taus = taus,
    weights = pop_w
  )
  
  coefficients <- coeffs$coef
  
  iv_coef <- coefficients[
    "d_tradeusch_p1_2000_2012",
    ,
  ]
  
  # save matrix
  results_iv[[as.character(yr)]] <- iv_coef
  
  # -------------------------------------------------------------
  # SURFACE PLOT
  # -------------------------------------------------------------
  
  p_surface <- plot_qq_surface(
    iv_coef,
    taus = taus,
    tplot = 1,
    mycolors = "",
    z_label = "Shock Impact",
    group = "CZ",
    zlim = z_min_max,
    year = yr
  )
  
  # save html widget
  htmlwidgets::saveWidget(
    p_surface,
    file = paste0(
      "plots/iv_qoq/control_function_approach/surfaces/ivl_surface_cz_",
      yr,
      ".html"
    ),
    selfcontained = TRUE
  )
  
  # -------------------------------------------------------------
  # SLICES PLOT
  # -------------------------------------------------------------
  
  p_slices <- plot_qq_slices(
    iv_coef,
    taus = taus,
    z_label = "Shock Impact",
    group = "CZ",
    ylim = y_min_max
  ) +
    ggtitle(
      paste(
        "China Shock Impact - Year",
        yr
      )
    )
  
  # save png
  ggsave(
    filename = paste0(
      "plots/iv_qoq/control_function_approach/slices/ivl_slices_cz_",
      yr,
      ".png"
    ),
    plot = p_slices,
    width = 8,
    height = 5,
    dpi = 300,
    bg = "white"
  )
}


saveRDS(
  results_iv,
  "plots/iv_qoq/control_function_approach/surfaces/all_ivl_surfaces.rds"
)




# Animated Plotly Surface for IV QoQ Coefficients


years <- 2005:2019


z_list <- list()
x_list <- list()
y_list <- list()

idx <- seq(1, length(taus), by = 3)

for (yr in years) {
  
  cat("Preparing year:", yr, "\n")
  
  iv_mat <- results_iv[[as.character(yr)]]
  
  z_list[[as.character(yr)]] <- iv_mat[idx, idx, drop = FALSE]
  
  x_list[[as.character(yr)]] <- taus[idx]
  y_list[[as.character(yr)]] <- taus[idx]
}


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

# -------------------------------------------------------------
# LAYOUT + SLIDER
# -------------------------------------------------------------

p <- p %>%
  layout(
    
    title = paste(
      "China Shock Impact - Year:",
      first_year
    ),
    
    scene = list(
      
      xaxis = list(
        title = "Between CZ Quantile",
        autorange = "reversed"
      ),
      
      yaxis = list(
        title = "Within CZ Quantile",
        autorange = "reversed"
      ),
      
      zaxis = list(
        title = "Shock Impact",
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

# -------------------------------------------------------------
# SAVE
# -------------------------------------------------------------

htmlwidgets::saveWidget(
  p,
  
  "plots/iv_qoq/control_function_approach/iv_qoq_animation_2005_2019.html",
  
  selfcontained = TRUE
)




#################################################################
## 2. QoQ-control-function-IV approach pooled in 5-year periods
#################################################################

# -------------------------------------------------------------
# define 5-year pooled periods
# -------------------------------------------------------------

data <- data %>%
  mutate(
    period5 = case_when(
      YEAR %in% 2005:2009 ~ "2005-2009",
      YEAR %in% 2010:2014 ~ "2010-2014",
      YEAR %in% 2015:2019 ~ "2015-2019",
      TRUE ~ NA_character_
    )
  )

periods5 <- c(
  "2005-2009",
  "2010-2014",
  "2015-2019"
)

# -------------------------------------------------------------
# folders
# -------------------------------------------------------------

dir.create(
  "plots/iv_qoq/control_function_approach/surfaces_5year",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "plots/iv_qoq/control_function_approach/slices_5year",
  recursive = TRUE,
  showWarnings = FALSE
)


# -------------------------------------------------------------
# store results
# -------------------------------------------------------------

results_iv_5year <- list()

# -------------------------------------------------------------
# LOOP OVER 5-YEAR PERIODS
# -------------------------------------------------------------

for (p in periods5) {
  
  cat("Processing period:", p, "\n")
  
  # -----------------------------------------------------------
  # subset pooled years
  # -----------------------------------------------------------
  
  df_period <- data %>%
    filter(period5 == p)
  
  # merge shocks + controls
  df_period <- shocks %>%
    left_join(df_period, by = "czone")
  
  # -----------------------------------------------------------
  # FIRST STAGE QQ
  # -----------------------------------------------------------
  
  qq_first <- qq_first_stage_collapsed(
    y = df_period$lnwage,
    group_vec = df_period$czone,
    taus = taus, 
    pweights = df_period$pweight #weighted by sample weights
  )
  
  # -----------------------------------------------------------
  # BUILD X MATRIX
  # -----------------------------------------------------------
  
  # master order
  cz <- qq_first$group_order
  
  # X aligned
  X_df <- shocks %>%
    filter(czone %in% cz) %>%
    mutate(czone = factor(czone, levels = cz)) %>%
    arrange(czone)
  
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
  
  # weights aligned
  pop_w <- data.frame(czone = qq_first$group_order) %>%
    left_join(pop_2000, by = "czone") %>%
    pull(pop_share_2000)
  
  # -------------------------------------------------------------
  # second stage
  # -------------------------------------------------------------
  
  coeffs <- qq_second_stage_group(
    fitted_first = qq_first$fitted,
    X = X,
    taus = taus,
    weights = pop_w
  )
  
  coefficients <- coeffs$coef
  
  iv_coef <- coefficients[
    "d_tradeusch_p1_2000_2012",
    ,
  ]
  
  # save matrix
  results_iv_5year[[p]] <- iv_coef
  
  # -----------------------------------------------------------
  # SURFACE PLOT
  # -----------------------------------------------------------
  
  p_surface <- plot_qq_surface(
    iv_coef,
    taus = taus,
    tplot = 1,
    mycolors = NULL,
    z_label = "Shock Impact",
    group = "CZ",
    zlim = z_min_max,
    year = p
  )
  
  htmlwidgets::saveWidget(
    p_surface,
    file = paste0(
      "plots/iv_qoq/control_function_approach/surfaces_5year/iv_surface_",
      p,
      ".html"
    ),
    selfcontained = TRUE
  )
  
  # -----------------------------------------------------------
  # SLICES
  # -----------------------------------------------------------
  
  p_slices <- plot_qq_slices(
    iv_coef,
    taus = taus,
    z_label = "Shock Impact",
    group = "CZ",
    ylim = y_min_max
  ) +
    ggtitle(
      paste(
        "China Shock Impact - 5 Years",
        p
      )
    )
  
  ggsave(
    filename = paste0(
      "plots/iv_qoq/control_function_approach/slices_5year/iv_slices_",
      p,
      ".png"
    ),
    plot = p_slices,
    width = 8,
    height = 5,
    dpi = 300,
    bg = "white"
  )
}

# -------------------------------------------------------------
# SAVE ALL RESULTS
# -------------------------------------------------------------

saveRDS(
  results_iv_5year,
  "plots/iv_qoq/control_function_approach/surfaces_5year/all_iv_surfaces_5year.rds"
)

#################################################################
## Animated plotly surface for pooled 5-year periods
#################################################################

periods5 <- c(
  "2005-2009",
  "2010-2014",
  "2015-2019"
)

z_list <- list()
x_list <- list()
y_list <- list()

idx <- seq(1, length(taus), by = 1)

for (p in periods5) {
  
  cat("Preparing:", p, "\n")
  
  iv_mat <- results_iv_5year[[p]]
  
  z_list[[p]] <- iv_mat[idx, idx, drop = FALSE]
  
  x_list[[p]] <- taus[idx]
  y_list[[p]] <- taus[idx]
}

# -------------------------------------------------------------
# initial plot
# -------------------------------------------------------------

first_period <- periods5[1]

p <- plot_ly(
  type = "surface",
  
  x = x_list[[first_period]],
  y = y_list[[first_period]],
  z = z_list[[first_period]],
  
  cmin = z_min_max[1],
  cmax = z_min_max[2],
  
  colorscale = "RdBu"
)

# -------------------------------------------------------------
# frames
# -------------------------------------------------------------

frames <- lapply(periods5, function(pname) {
  
  list(
    name = pname,
    
    data = list(
      list(
        type = "surface",
        
        x = x_list[[pname]],
        y = y_list[[pname]],
        z = z_list[[pname]],
        
        cmin = z_min_max[1],
        cmax = z_min_max[2],
        
        colorscale = NULL
      )
    )
  )
})

p$x$frames <- frames

# -------------------------------------------------------------
# layout + slider
# -------------------------------------------------------------

p <- p %>%
  layout(
    
    title = paste(
      "China Shock Impact",
      first_period
    ),
    
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
        title = "Shock Impact",
        range = z_min_max
      )
    ),
    
    sliders = list(
      list(
        active = 0,
        
        currentvalue = list(
          prefix = "Period: "
        ),
        
        steps = lapply(periods5, function(pname) {
          
          list(
            label = pname,
            
            method = "animate",
            
            args = list(
              list(pname),
              
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

# -------------------------------------------------------------
# save animation
# -------------------------------------------------------------

htmlwidgets::saveWidget(
  p,
  
  "plots/iv_qoq/control_function_approach/iv_qoq_animation_5year.html",
  
  selfcontained = TRUE
)

