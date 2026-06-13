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

####### get rid of outliers 
data <- data %>%
  filter(wage > 5) %>%
  filter(wage < 500)


######################################################################################################################################################################
##  Variable merging
######################################################################################################################################################################
# use controls 1990 Autor
controls_1990 <- controls %>%
  filter(year == 1990)
# merge controls and shocks by CZ for the first stage IV 
shocks <- shocks %>%
  left_join(controls_1990, by = "czone")


# define parameters
grid <- seq(-10, 10, by = 0.1)
taus <- round(seq(0.1, 0.9, by = 0.1), 1)
z_min_max <- c(-7, 9)
y_min_max <- c(-7, 9)

################################################################
## QoQ-CH-IVQR appraoch for every year (2005–2019)
#################################################################

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

dir.create(
  "plots/iv_qoq/IVQR_approach/slices",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "plots/iv_qoq/IVQR_approach/surfaces",
  recursive = TRUE,
  showWarnings = FALSE
)


years <- 2005:2019

results_iv <- list()

for (yr in years) {
  
  cat("Processing year:", yr, "\n")
  
  # -------------------------------------------------------------
  # subset year
  # -------------------------------------------------------------
  df_year <- data %>%
    filter(YEAR == yr)
  
  df_year <- shocks %>%
    left_join(df_year, by = "czone")
  
  # -------------------------------------------------------------
  # STEP 1: QQ FIRST STAGE
  # -------------------------------------------------------------
  qq_first <- qq_first_stage_collapsed(
    y = df_year$lnwage,
    group_vec = df_year$czone,
    taus = taus,
    pweights = df_year$pweight
  )
  
  cz <- qq_first$group_order
  
  # -------------------------------------------------------------
  # STEP 2: ALIGN SHOCKS MATRIX
  # -------------------------------------------------------------
  X_df <- shocks %>%
    filter(czone %in% cz) %>%
    mutate(czone = factor(czone, levels = cz)) %>%
    arrange(czone)
  
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
  
  # -------------------------------------------------------------
  # STEP 3: IV-QoQ SECOND STAGE
  # -------------------------------------------------------------
  coeffs <- qq_second_stage_ivqr_group(
    fitted_first = qq_first$fitted,
    D = D,
    Z = Z,
    X = X,
    taus = taus,
    grid = grid
  )
  
  iv_coef <- coeffs$beta_D
  
  results_iv[[as.character(yr)]] <- iv_coef
  
  # -------------------------------------------------------------
  # PLOT SURFACE
  # -------------------------------------------------------------
  p_surface <- plot_qq_surface(
    qq_mat = iv_coef,
    taus = taus,
    tplot = 1,
    mycolors = "RdBu",
    z_label = "Shock Impact",
    group = "CZ",
    zlim = z_min_max,
    year = yr
  )
  
  htmlwidgets::saveWidget(
    p_surface,
    file = paste0(
      "plots/iv_qoq/IVQR_approach/surfaces/ivqr_surface_cz_",
      yr,
      ".html"
    ),
    selfcontained = TRUE
  )
  
  # -------------------------------------------------------------
  # PLOT SLICES
  # -------------------------------------------------------------
  p_slices <- plot_qq_slices(
    iv_coef,
    taus = taus,
    z_label = "Shock Impact",
    group = "CZ"
  ) +
    ggtitle(paste("IV QoQ Shock Impact - Year", yr))
  
  ggsave(
    filename = paste0(
      "plots/iv_qoq/IVQR_approach/slices/ivqr_slices_cz_",
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
  "plots/iv_qoq/IVQR_approach/surfaces/all_ivqr_surfaces.rds"
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
  
  "plots/iv_qoq/IVQR_approach/ivqr_qoq_animation_2005_2019.html",
  
  selfcontained = TRUE
)


#################################################################
## 2. QoQ-IVQR approach pooled in 5-year periods
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
  "plots/iv_qoq/IVQR_approach/surfaces_5year",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "plots/iv_qoq/IVQR_approach/slices_5year",
  recursive = TRUE,
  showWarnings = FALSE
)

# -------------------------------------------------------------
# store results
# -------------------------------------------------------------

results_iv_5year <- list()


for (p in periods5) {
  
  cat("Processing period:", p, "\n")
  
  # ---------------------------------------------------------
  # subset
  # ---------------------------------------------------------
  df_period <- data %>%
    filter(period5 == p)
  
  df_period <- shocks %>%
    left_join(df_period, by = "czone")
  
  # ---------------------------------------------------------
  # FIRST STAGE (QoQ)
  # ---------------------------------------------------------
  qq_first <- qq_first_stage_collapsed(
    y = df_period$lnwage,
    group_vec = df_period$czone,
    taus = taus,
    pweights = df_period$pweight
  )
  
  cz <- qq_first$group_order
  
  # ---------------------------------------------------------
  # ALIGN SHOCKS
  # ---------------------------------------------------------
  X_df <- shocks %>%
    filter(czone %in% cz) %>%
    mutate(czone = factor(czone, levels = cz)) %>%
    arrange(czone)
  
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
  
  # ---------------------------------------------------------
  # SECOND STAGE (IVQR QoQ)
  # ---------------------------------------------------------
  coeffs <- qq_second_stage_ivqr_group(
    fitted_first = qq_first$fitted,
    D = D,
    Z = Z,
    X = X,
    taus = taus,
    grid = grid
  )
  
  iv_coef <- coeffs$beta_D
  
  results_iv_5year[[p]] <- iv_coef
  
  # ---------------------------------------------------------
  # SURFACE PLOT
  # ---------------------------------------------------------
  p_surface <- plot_qq_surface(
    qq_mat = iv_coef,
    taus = taus,
    tplot = 1,
    mycolors = "RdBu",
    z_label = "Shock Impact",
    group = "CZ",
    zlim = z_min_max,
    year = p
  )
  
  htmlwidgets::saveWidget(
    p_surface,
    file = paste0(
      "plots/iv_qoq/IVQR_approach/surfaces_5year/ivqr_surface_",
      p,
      ".html"
    ),
    selfcontained = TRUE
  )
  
  # ---------------------------------------------------------
  # SLICES
  # ---------------------------------------------------------
  p_slices <- plot_qq_slices(
    iv_coef,
    taus = taus,
    z_label = "Shock Impact",
    group = "CZ"
  ) +
    ggtitle(paste("IVQR QoQ Shock Impact -", p))
  
  ggsave(
    filename = paste0(
      "plots/iv_qoq/IVQR_approach/slices_5year/ivqr_slices_",
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

#############################################################
## SAVE RESULTS
#############################################################

saveRDS(
  results_iv_5year,
  "plots/iv_qoq/IVQR_approach/surfaces/all_ivqr_surfaces_5year.rds"
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
  
  "plots/iv_qoq/IVQR_approach/iv_qoq_animation_5year.html",
  
  selfcontained = TRUE
)





