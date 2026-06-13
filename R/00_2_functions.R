library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
#install ivqr package
#install.packages("remotes")
#remotes::install_github("yuchang0321/IVQR")
library(IVQR)






######################################################################################################################################################################
## 1. Functions for Quantile on Quantiles Descriptive Analysis (Without Covariates)
######################################################################################################################################################################

# 1.1 making qoq surface without covariates 
make_qq_matrix <- function(data, y, group, weight = "pweight", taus) {
  df <- data %>%
    dplyr::select(all_of(c(y, group, weight)))
  names(df) <- c("y", "group", "w")
  df <- df[!is.na(df$w) & df$w > 0, ]
  
  # -----------------------------
  # 1. WITHIN-GROUP QUANTILES
  # -----------------------------
  within_list <- split(df, df$group)
  
  within_mat <- t(vapply(within_list, function(g) {
    
    if (nrow(g) == 0 || all(is.na(g$y))) {
      return(rep(NA_real_, length(taus)))
    }
    
    sapply(taus, function(u)
      modi::weighted.quantile(g$y, g$w, prob = u, plot = FALSE)
    )
    
  }, FUN.VALUE = numeric(length(taus))))
  
  # -----------------------------
  # 2. BETWEEN-GROUP QUANTILES
  # -----------------------------
  group_w <- tapply(df$w, df$group, sum, na.rm = TRUE)
  
  valid_groups <- names(group_w)[group_w > 0 & !is.na(group_w)]
  
  within_mat <- within_mat[valid_groups, , drop = FALSE]
  group_w <- group_w[valid_groups]
  
  qq_mat <- t(vapply(seq_len(ncol(within_mat)), function(j) {
    
    y_vec  <- within_mat[, j]
    ww_vec <- group_w
    
    ok <- !is.na(y_vec) & !is.na(ww_vec)
    
    if (sum(ok) == 0) {
      return(rep(NA_real_, length(taus)))
    }
    
    sapply(taus, function(v)
      modi::weighted.quantile(y_vec[ok], ww_vec[ok], prob = v, plot = FALSE)
    )
    
  }, FUN.VALUE = numeric(length(taus))))

  rownames(qq_mat) <- paste0("u_", sprintf("%.2f", taus))
  colnames(qq_mat) <- paste0("v_", sprintf("%.2f", taus))
  
  return(qq_mat)
}

# 1.2 plotting qq surface in three dimensional space 
plot_qq_surface <- function(qq_mat, taus, tplot = 5, mycolors = "viridis",
                            reverse_axes = TRUE, z_label = "Wage",
                            group = "group", zlim = NULL, year = NULL) {
  
  idx <- seq(1, length(taus), by = tplot)
  z <- qq_mat[idx, idx, drop = FALSE]
  
  x <- taus[idx]
  y <- taus[idx]
  
  title_txt <- if (!is.null(year)) paste("Year:", year) else NULL
  
  p <- plot_ly(
    x = x,
    y = y,
    z = z
  ) %>%
    add_surface(
      z = z,
      colorscale = mycolors,
      cmin = if (!is.null(zlim)) zlim[1] else NULL,
      cmax = if (!is.null(zlim)) zlim[2] else NULL,
      reversescale = FALSE,
      opacity = 0.99,
      contours = list(
        x = list(show = TRUE, color = "gray", width = 1),
        y = list(show = TRUE, color = "gray", width = 1),
        z = list(show = FALSE)
      ),
      showscale = TRUE,
      colorbar = list(title = z_label)
    )
  
  scene_axes <- list(
    xaxis = list(
      title = paste("Between", group),
      showline = FALSE,
      zeroline = FALSE
    ),
    yaxis = list(
      title = paste("Within", group),
      showline = FALSE,
      zeroline = FALSE
    ),
    zaxis = list(
      title = z_label,
      showline = FALSE,
      zeroline = FALSE,
      tickformat = ".0f",
      range = zlim   
    )
  )
  
  if (reverse_axes) {
    scene_axes$xaxis$autorange <- "reversed"
    scene_axes$yaxis$autorange <- "reversed"
  }
  
  p <- p %>% layout(scene = scene_axes)
  
  # add year as title (for GIF frames)
  if (!is.null(year)) {
    p <- p %>% layout(title = title_txt)
  }
  
  return(p)
}

# 1.3 plotting qoq surface in slices 
plot_qq_slices <- function(qq_mat,
                           taus = NULL,
                           u_values = seq(0.1, 0.9, by = 0.1), 
                           z_label = "Wage", 
                           group = "group", 
                           ylim = NULL) {
  
  # --- handle taus
  if (is.null(taus)) {
    taus <- seq(0, 1, length.out = ncol(qq_mat))
  }
  
  taus <- round(taus, 2)
  u_values <- round(u_values, 2)
  
  # --- find closest indices
  u_idx <- sapply(u_values, function(u) {
    which.min(abs(taus - u))
  })
  
  # --- subset matrix
  qq_slice <- qq_mat[u_idx, , drop = FALSE]
  
  # --- reshape
  plot_dat <- as.data.frame(qq_slice)
  plot_dat$u <- taus[u_idx]
  
  plot_dat_long <- plot_dat %>%
    tidyr::pivot_longer(
      cols = -u,
      names_to = "v_index",
      values_to = "value"
    ) %>%
    dplyr::group_by(u) %>%
    dplyr::mutate(v = taus) %>%
    dplyr::ungroup()
  
  # --- ensure ordering
  plot_dat_long <- plot_dat_long %>%
    dplyr::mutate(u = factor(u, levels = sort(unique(u))))
  
  # --- label positions
  label_dat <- plot_dat_long %>%
    dplyr::group_by(u) %>%
    dplyr::slice_tail(n = 1) %>%
    dplyr::ungroup()
  
  # --- build plot (IMPORTANT FIX)
  p <- ggplot2::ggplot(plot_dat_long,
                       ggplot2::aes(x = v, y = value, group = u)) +
    ggplot2::geom_line(color = "black", linewidth = 0.4) +
    ggplot2::geom_text(
      data = label_dat,
      ggplot2::aes(label = paste0("u = ", u)),
      hjust = -0.1,
      size = 3.3
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.01, 0.15))
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      x = paste("Between", group),
      y = z_label
    )
  
  # --- fixed y-axis (correct placement)
  if (!is.null(ylim)) {
    p <- p + ggplot2::coord_cartesian(ylim = ylim)
  }
  
  return(p)
}






######################################################################################################################################################################
## 2. Functions for Quantile on Quantiles Analysis in two stages (Functions copied from https://github.com/martinapons/qqreg/blob/main/R/ )
######################################################################################################################################################################
# 2.1 fist stage of qoq approach
qq_first_stage <- function(y, X1, group_vec, taus, var_types,
                           parallel = TRUE, ncores = NULL, pweights = NULL) {
  groups <- unique(group_vec)
  n <- length(y)
  n_tau <- length(taus)
  
  # Function to run first stage for a single group
  run_group <- function(g) {
    idx <- which(group_vec == g)
    y_g <- y[idx]
    n_g <- length(y_g)
    fitted_g <- matrix(NA, n_g, n_tau)
    
    # Get design matrix for this group
    if (is.null(X1) || ncol(X1) == 0) {
      X1_g <- NULL
    } else {
      X1_g <- X1[idx, , drop = FALSE]
      
      # Remove zero-variance columns for this group
      var_g <- apply(X1_g, 2, stats::var, na.rm = TRUE)
      var_g[is.na(var_g)] <- 0
      keep_cols <- var_g > 1e-10
      
      if (any(keep_cols)) {
        X1_g <- X1_g[, keep_cols, drop = FALSE]
        
        # Check rank and remove linearly dependent columns
        qr_X <- qr(X1_g, tol = 1e-9, LAPACK = FALSE)
        if (qr_X$rank < ncol(X1_g)) {
          X1_g <- X1_g[, qr_X$pivot[seq_len(qr_X$rank)], drop = FALSE]
        }
      } else {
        X1_g <- NULL
      }
    }
    
    # Run quantile regression for each tau
    w_g <- if (is.null(pweights)) rep(1, length(y_g)) else pweights[idx] 
    for (t_idx in seq_along(taus)) { 
      if (is.null(X1_g) || ncol(X1_g) == 0) { 
        fitted_g[, t_idx] <- modi::weighted.quantile(y_g, w = w_g, prob = taus[t_idx], plot = FALSE) 
      } 
      else { 
        qr_fit <- quantreg::rq(y_g ~ X1_g, tau = taus[t_idx], weights = w_g) 
        fitted_g[, t_idx] <- qr_fit$fitted.values 
        }
      }
    
    # Sort fitted values within group to ensure monotonicity in u
    # Each row: sort across tau columns
    fitted_g <- t(apply(fitted_g, 1, sort))
    
    list(idx = idx, fitted = fitted_g)
  }
  
  # Run across all groups
  if (parallel && length(groups) > 1) {
    # Set up parallel backend
    if (is.null(ncores)) {
      ncores <- max(1, floor(0.8 * future::availableCores()))
    }
    future::plan(future::multisession, workers = ncores)
    
    results <- furrr::future_map(groups, run_group,
                                 .options = furrr::furrr_options(seed = TRUE))
  } else {
    results <- lapply(groups, run_group)
  }
  
  # Combine results into single matrix
  fitted <- matrix(NA, n, n_tau)
  for (res in results) {
    fitted[res$idx, ] <- res$fitted
  }
  
  colnames(fitted) <- paste0("u_", taus)
  return(fitted)
}


#' First Stage for Collapsed Data
#'
#' When no within-group varying variables exist, computes group-level quantiles
#' and returns collapsed data with one observation per group per tau.
#'
#' @param y Outcome vector.
#' @param group_vec Group indicator vector.
#' @param taus Numeric vector of quantile levels.
#'
#' @return A list with:
#'   \item{fitted}{Matrix of group-level quantiles (n_groups x n_tau)}
#'   \item{weights}{Vector of group sizes for weighting}
#'   \item{group_order}{Order of groups in the output}
#'
#' @keywords internal

#2.2 First stage of qoq appraoch when no within group covariates 
qq_first_stage_collapsed <- function(y, group_vec, taus, pweights = NULL) {
  groups <- unique(group_vec)
  n_groups <- length(groups)
  n_tau <- length(taus)
  
  fitted <- matrix(NA, n_groups, n_tau)
  weights <- numeric(n_groups)
  
  for (i in seq_along(groups)) {
    g <- groups[i]
    idx <- which(group_vec == g)
    y_g <- y[idx]
    weights[i] <- length(y_g)
    
    w_g <- if (is.null(pweights)) rep(1, length(y_g)) else pweights[idx] 
    for (t_idx in seq_along(taus)) {
      fitted[i, t_idx] <- modi::weighted.quantile(y_g, w = w_g, prob = taus[t_idx], plot = FALSE) 
    }
  }
  
  colnames(fitted) <- paste0("u_", taus)
  rownames(fitted) <- as.character(groups)
  
  list(
    fitted = fitted,
    weights = weights,
    group_order = groups
  )
}

# second stage qoq appraoch using first stage fitted values 
qq_second_stage_group <- function(fitted_first, X, taus, weights = NULL) {
  n_tau <- length(taus)
  n <- nrow(X)
  k <- ncol(X)
  
  # Initialize arrays
  coef_array <- array(NA, dim = c(k, n_tau, n_tau))
  fitted_array <- array(NA, dim = c(n, n_tau, n_tau))
  resid_array <- array(NA, dim = c(n, n_tau, n_tau))
  
  # Loop over u (rows of first-stage fitted values)
  for (u_idx in seq_along(taus)) {
    y_u <- fitted_first[, u_idx]
    
    # Loop over v
    for (v_idx in seq_along(taus)) {
      # Run quantile regression
      if (is.null(weights)) {
        qr_fit <- quantreg::rq(y_u ~ X - 1, tau = taus[v_idx])
      } else {
        qr_fit <- quantreg::rq(y_u ~ X - 1, tau = taus[v_idx], weights = weights)
      }
      
      coef_array[, u_idx, v_idx] <- qr_fit$coefficients
      fitted_array[, u_idx, v_idx] <- qr_fit$fitted.values
      resid_array[, u_idx, v_idx] <- qr_fit$residuals
    }
  }
  
  # Set dimension names
  dimnames(coef_array) <- list(
    colnames(X),
    paste0("u_", taus),
    paste0("v_", taus)
  )
  dimnames(fitted_array) <- list(
    NULL,
    paste0("u_", taus),
    paste0("v_", taus)
  )
  dimnames(resid_array) <- dimnames(fitted_array)
  
  list(
    coef = coef_array,
    fitted = fitted_array,
    residuals = resid_array
  )
}

######################################################################################################################################################################
## 3. Functions for Quantile on Quantiles Analysis in two stages with IVQR 
######################################################################################################################################################################
# fix IVQR to not include standard errors as there is a bug there
ivqr_a <- function (formula, taus = 0.5, data, grid, gridMethod = "Default", 
          ivqrMethod = "iqr", qrMethod = "br") 
{
  formula <- Formula(formula)
  if (length(formula)[2] < 3) {
    stop("If there's no control variable, specify the model like y ~ d | z | 1")
  }
  eps <- .Machine$double.eps^(2/3)
  if (length(taus) > 1) {
    if (any(taus < 0) || any(taus > 1)) 
      stop("invalid taus:  taus should be >= 0 and <= 1")
  }
  if (any(taus == 0)) 
    taus[taus == 0] <- eps
  if (any(taus == 1)) 
    taus[taus == 1] <- 1 - eps
  if (!any(grepl(qrMethod, c("br", "fn", "fnb")))) {
    stop("Please specify one of the quantreg solution method: br, fn, or fnb")
  }
  XZ <- formula(formula, lhs = 1, rhs = c(2, 3), collapse = TRUE)
  XZ <- model.matrix(XZ, data)
  D <- update(formula(formula, lhs = 1, rhs = 1), . ~ . - 1)
  D <- model.matrix(D, data)
  D_hat <- D
  for (i in 1:dim(D)[2]) {
    beta <- solve(crossprod(XZ), crossprod(XZ, D[, i]))
    D_hat[, i] <- XZ %*% beta
  }
  if (any(grepl(".ivqr_dhat", colnames(data)))) {
    stop("No names of variables in the data set should include .ivqr_dhat")
  }
  dhat_formula <- c("~")
  copy_data <- data
  for (i in 1:dim(D)[2]) {
    data[, paste(".ivqr_dhat", i, sep = "")] <- D_hat[, i]
    dhat_formula <- cbind(dhat_formula, paste(".ivqr_dhat", 
                                              i, sep = ""))
    if (i < dim(D)[2]) 
      dhat_formula <- cbind(dhat_formula, "+")
  }
  dhat_formula <- formula(paste(dhat_formula, collapse = ""))
  iqr_formula <- as.Formula(formula(formula, lhs = 1, rhs = 1), 
                            dhat_formula, formula(formula, lhs = 0, rhs = 3))
  coef <- list()
  endg_varnames <- formula(iqr_formula, lhs = 1, rhs = 1)
  endg_varnames <- update(endg_varnames, . ~ . - 1)
  inst_varnames <- formula(iqr_formula, lhs = 1, rhs = 2)
  inst_varnames <- update(inst_varnames, . ~ . - 1)
  exog_varnames <- formula(iqr_formula, lhs = 1, rhs = 3)
  D <- model.matrix(endg_varnames, data)
  PHI <- model.matrix(inst_varnames, data)
  X <- model.matrix(exog_varnames, data)
  coef$endg_var <- matrix(NA, ncol(D), length(taus))
  coef$inst_var <- matrix(NA, ncol(PHI), length(taus))
  coef$exog_var <- matrix(NA, ncol(X), length(taus))
  if (ncol(D) != ncol(PHI)) {
    warning("Shouldn't ncol(D) always equal ncol(PHI)?")
  }
  if (ncol(D) == 1) {
    if (!is.vector(grid)) {
      stop("Dimension of the grid does not match the numbers of endogenous variables. Grid should be a vector for dim(D) = 1")
    }
  }
  else if (ncol(D) == 2) {
    if (!is.list(grid)) 
      stop("Dimension of the grid does not match the number of endogenous variables. Grid should be a matrix with 2 rows for dim(D) == 2")
  }
  else if (ncol(D) >= 2) {
    stop("Complexity grows exponentially in the number of endogenous variables. This version only deals with cases when dim(D) <= 2")
  }
  fitted <- matrix(NA, nrow(X), length(taus))
  residuals <- matrix(NA, nrow(X), length(taus))
  if (!is.list(grid)) {
    grid_value <- matrix(NA, length(grid), length(taus))
  }
  else {
    grid_value <- array(NA, dim = c(length(grid[[1]]), length(grid[[2]]), 
                                    length(taus)))
  }
  error_tau_flag <- !logical(length(taus))
  for (i in 1:length(taus)) {
    ivqr_est <- IVQR:::ivqr.fit(iqr_formula, tau = taus[i], data, 
                         grid, gridMethod, ivqrMethod, qrMethod)
    if (is.list(ivqr_est)) {
      coef$endg_var[, i] <- ivqr_est$coef_endg_var
      coef$inst_var[, i] <- ivqr_est$coef_inst_var
      coef$exog_var[, i] <- ivqr_est$coef_exog_var
      rownames(coef$endg_var) <- names(ivqr_est$coef_endg_var)
      rownames(coef$inst_var) <- names(ivqr_est$coef_inst_var)
      rownames(coef$exog_var) <- names(ivqr_est$coef_exog_var)
      residuals[, i] <- ivqr_est$residuals
      fitted[, i] <- ivqr_est$fitted
      if (!is.list(grid)) {
        grid_value[, i] <- ivqr_est$grid_value
      }
      else {
        grid_value[, , i] <- ivqr_est$grid_value
      }
      error_tau_flag[i] <- FALSE
    }
  }
  taulabs <- paste("tau=", format(round(taus, 3)))
  colnames(coef$endg_var) <- colnames(coef$inst_var) <- colnames(coef$exog_var) <- taulabs
  fit <- list()
  class(fit) <- "ivqr"
  fit$coef <- coef
  fit$fitted <- fitted
  fit$residuals <- residuals
  fit$formula <- formula
  fit$taus <- taus
  fit$error_tau_flag <- error_tau_flag
  fit$data <- data
  fit$dim_d_d_k <- c(ncol(D), ncol(D), ncol(X))
  fit$n <- nrow(X)
  fit$obj_fcn <- grid_value
  fit$grid <- grid
  fit$copy_data <- copy_data
  fit$gridMethod <- gridMethod
  fit$ivqrMethod <- ivqrMethod
  fit$qrMethod <- qrMethod
  PSI <- cbind(PHI, X)
  DX <- cbind(D, X)
  fit$DX <- DX
  fit$PSI <- PSI
  fit$se <- NULL #adjusted from IVQR
  fit$vc <- NULL #adjusted from IVQR
  return(fit)
}


# second stage qoq appraoch with IVQR using first stage fitted values 
qq_second_stage_ivqr_group <- function(fitted_first, D, Z, X, taus, grid) 
  {
  n_tau <- length(taus)
  n     <- nrow(X)
  
  beta_D <- matrix(
    NA,
    nrow = n_tau,
    ncol = n_tau,
    dimnames = list(
      paste0("u_", taus),
      paste0("v_", taus)
    )
  )
  
  # -----------------------------------------------------------
  # fitted + residuals
  # -----------------------------------------------------------
  
  fitted_array <- array(
    NA,
    dim = c(n, n_tau, n_tau)
  )
  
  resid_array <- array(
    NA,
    dim = c(n, n_tau, n_tau)
  )
  
  # -----------------------------------------------------------
  # store full IVQR objects
  # -----------------------------------------------------------
  
  ivqr_fits <- vector(
    "list",
    n_tau * n_tau
  )
  
  dim(ivqr_fits) <- c(n_tau, n_tau)
  
  dimnames(ivqr_fits) <- list(
    paste0("u_", taus),
    paste0("v_", taus)
  )
  
  # -----------------------------------------------------------
  # estimation loop
  # -----------------------------------------------------------
  
  for (u_idx in seq_along(taus)) {
    
    y_u <- fitted_first[, u_idx]
    
    for (v_idx in seq_along(taus)) {
      
      data_df <- data.frame(
        y = y_u,
        D = D,
        Z = Z
      )
      
      data_df <- cbind(
        data_df,
        as.data.frame(X)
      )
      
      fit <- ivqr_a(
        formula = y ~ D | Z | X,
        data    = data_df,
        taus    = taus[v_idx],
        grid    = grid
      )
      
      ivqr_fits[[u_idx, v_idx]] <- fit
      
      beta_D[u_idx, v_idx] <-
        as.numeric(fit$coef$endg_var)
      
      fitted_array[, u_idx, v_idx] <-
        fit$fitted[, 1]
      
      resid_array[, u_idx, v_idx] <-
        fit$residuals[, 1]
    }
  }
  
  list(
    beta_D     = beta_D,
    fitted     = fitted_array,
    residuals  = resid_array,
    ivqr_fits  = ivqr_fits
  )
}



# create cf distribution 
estimate_cf_distribution <- function(data, x_var, z_var, control_vars, taus = seq(0.01, 0.99, by = 0.01)) 
  {
  rhs <- paste(
    c(z_var, control_vars),
    collapse = " + "
  )
  
  fml <- as.formula(
    paste(x_var, "~", rhs)
  )
  
  # Estimate conditional quantile process
  rq_fit <- rq(
    fml,
    tau = taus,
    data = data
  )
  
  X_obs <- data[[x_var]]
  
  n <- nrow(data)
  
  V_hat <- numeric(n)
  
  # Predicted conditional quantiles
  pred_q <- predict(rq_fit)
  
  for(i in seq_len(n)) {
    
    # observed treatment
    x_i <- X_obs[i]
    
    # predicted quantiles for this observation
    q_i <- pred_q[i, ]
    
    # estimated conditional rank
    V_hat[i] <- approx(
      x = q_i,
      y = taus,
      xout = x_i,
      rule = 2
    )$y
  }
  
  data$V_hat <- V_hat
  
  return(data)
}

############################
# Chetverikov 2016 functions but in crosssectional way
# second stage qoq appraoch using first stage fitted values 
qq_second_stage_chetverikov <- function(fitted_first, X, taus, weights = NULL) {
  
  n_tau <- length(taus)
  n <- nrow(X)
  k <- ncol(X)
  
  coef_mat <- matrix(NA, nrow = k, ncol = n_tau,
                     dimnames = list(colnames(X), paste0("u_", taus)))
  
  se_mat <- matrix(
    NA,
    nrow = k,
    ncol = n_tau,
    dimnames = list(
      colnames(X),
      paste0("u_", taus)
    )
  )  
  fitted_mat <- matrix(NA, nrow = n, ncol = n_tau)
  resid_mat  <- matrix(NA, nrow = n, ncol = n_tau)
  
  for (u_idx in seq_along(taus)) {
    
    y_u <- fitted_first[, u_idx]
    
    fit <- if (is.null(weights)) {
      lm(y_u ~ X-1)
    } else {
      lm(y_u ~ X-1, weights = weights)
    }
    
    cf <- coef(fit)
    coef_mat[, u_idx] <- cf
    
    vc <- vcov(fit)
    se_mat[, u_idx] <- sqrt(diag(vc))
    
    fitted_mat[, u_idx] <- fitted(fit)
    resid_mat[, u_idx] <- resid(fit)
  }
  
  list(
    coef = coef_mat,
    se = se_mat,
    fitted = fitted_mat,
    residuals = resid_mat
  )
}