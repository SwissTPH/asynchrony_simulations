# ==============================================================================
#  Useful functions for simulations and plots
#
#                 Younes Iggidr & Clara Champagne
#                             2025
#
#    _____         _           _______ _____  _    _ 
#   / ____|       (_)         |__   __|  __ \| |  | |
#  | (_____      ___ ___ ___     | |  | |__) | |__| |
#   \___ \ \ /\ / / / __/ __|    | |  |  ___/|  __  |
#   ____) \ V  V /| \__ \__ \    | |  | |    | |  | |
#  |_____/ \_/\_/ |_|___/___/    |_|  |_|    |_|  |_|
#
# ==============================================================================

#=================================== MODEL ===================================
rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source('model_implementation.r')

#=================================== IMPORTS ===================================
library(odin)
library(tidyverse)
library(ggplot2)
library(plotly)
library(patchwork)
library(cowplot)
library(nleqslv)
library(viridis)
library(tidyr)
library(rlang)
library(sensitivity)
library(lhs)
library(dplyr)
library(epiR)
library(randomForest)
library(FactoMineR)
library(factoextra)
library(cluster)
library(rpart)
library(rpart.plot)
library(fmsb)
library(pheatmap)
library(NbClust)


simulate_SIR <- function(n_days, start_interv, length_intervention, nb_studied_cycles, N, vars) {
  
  p_1 = vars["p_1"]
  p_2 = vars["p_2"]
  p_matrix = matrix(c(p_1, 1-p_1, 1-p_2, p_2), nrow = 2, byrow = TRUE)
  
  ## Interventions
  interv_times=c(0, start_interv+seq(0, 16)*365*length_intervention,n_days+2 )
  
  start_area1 = floor(start_interv/365)
  start_area2 = floor(n_days/365) + floor(start_interv/365)
  
  studied_duration_intervention = ceiling(length_intervention*nb_studied_cycles)
  end_area1 = start_area1 + studied_duration_intervention + 1
  end_area2 = start_area2 + studied_duration_intervention + 1
  z0 = c(vars["z01"], vars["z02"])
  x0 = c(vars["x01"], vars["x02"])
  alpha = c(vars["alpha1"], vars["alpha2"])
  rho = c(vars["rho1"], vars["rho2"])
  R0 = c(vars["R01"], vars["R02"])
  
  ## vector control
  W = rep(1, 2)
  W_1 = vars["W_1"]
  W_2 = vars["W_2"]
  
  interv_WA=cbind(c(1, rep(c(W_1,1), 8),1),
                  c(1, rep(c(W_2,1), 8),1)
  )
  
  interv_WB=cbind(c(1, rep(c(W_1,1), 8),1),
                  c(1, rep(c(1,W_2), 8),1)
  )
  
  W_A=create_timevarying_interventions(interv_WA, interv_times)
  W_B=create_timevarying_interventions(interv_WB, interv_times)
  
  myparameters <- list(N = N,
                       W = W,
                       alpha = alpha,
                       rho = rho,
                       R0 = R0, # assuming R0 low in Honduras foci but higher in Nicaragua
                       p = p_matrix,
                       z0 = z0,
                       x0 = x0, 
                       r=vars["r"]) # 10% of initial prevalence rate seems a good assumptions but even lower probably should be more relevant
  myparameters_A=myparameters
  myparameters_A$W=W_A
  
  myparameters_B=myparameters
  myparameters_B$W=W_B
  
  mysimul_A=run_SISconnected_deterministic(parameters=myparameters_A, n_days=n_days)
  mysimul_B=run_SISconnected_deterministic(parameters=myparameters_B, n_days=n_days)
  
  return(list(mysimul_A, mysimul_B))
}


compute_equilibrium_prevalence <- function(R0, N, r, p_1, p_2) {
  
  N_1 = N[1]
  N_2 = N[2]
  a = p_1
  b = 1-a
  d = p_2
  c = 1-d
  
  R0_1 = R0[1]
  R0_2 = R0[2]
  
  alpha = r*R0_1*a*N_1/(a*N_1 + c*N_1)
  gamma = r*R0_2*b*N_1/(b*N_1 + d*N_2)
  beta = r*R0_1*c*N_2/(a*N_1 + c*N_1)
  eta = r*R0_2*d*N_2/(b*N_1 + d*N_2)
  A = (a*alpha + b*gamma - r)
  B = (a*alpha + b*gamma)
  C = (a*beta + b*eta)
  E = (c*alpha + d*gamma - r)
  F = (c*alpha + d*gamma)
  G = (c*beta + d*eta)
  
  
  
  system <- function(vars) {
    x <- vars[1]
    y <- vars[2]
    
    eq1 <- A * x - B * x^2 - C * x * y + C * y
    eq2 <- E * y - F * y^2 - G * x * y + G * x
    
    return(c(eq1, eq2))
  }
  
  jacobian <- function(vars) {
    x <- vars[1]
    y <- vars[2]
    j11 <- A - 2*B*x - C*y         # ∂eq1/∂x
    j12 <- -C*x + C                # ∂eq1/∂y
    j21 <- -G*y + G                # ∂eq2/∂x
    j22 <- E - 2*F*y - G*x         # ∂eq2/∂y
    return(matrix(c(j11, j12, j21, j22), nrow=2, byrow=TRUE))
  }
  
  start <- c(0.5, 0.5)
  
  sol <- nleqslv(start, 
                 system, 
                 jac = jacobian,
                 method = "Newton",
                 control = list(ftol=1e-20, xtol=1e-20, maxit=50))
  return(c(sol$x[1], sol$x[2]))
}

compute_R0 <- function(x0, N, r, p_1, p_2) {
  
  N_1 = N[1]
  N_2 = N[2]
  a = p_1
  b = 1-a
  d = p_2
  c = 1-d
  X_1 = x0[1]
  X_2 = x0[2]
  
  
  alpha = (a*N_1*X_1 + c*N_2*X_2)/(a*N_1 + c*N_1)
  beta = (b*N_1*X_1 + d*N_2*X_2)/(b*N_1 + d*N_2)
  
  
  
  system <- function(vars) {
    x <- vars[1]
    y <- vars[2]
    
    eq1 <- a*alpha*x + b*beta*y - X_1/(1-X_1)
    eq2 <- c*alpha*x + d*beta*y - X_2/(1-X_2)
    
    return(c(eq1, eq2))
  }
  
  
  start <- c(1, 1)
  
  sol <- nleqslv(start, 
                 system)
  return(c(sol$x[1], sol$x[2]))
}



compute_metrics <- function(incidence_async, incidence_sync, n_days, 
                            start_interv, length_intervention,
                            nb_studied_cycles) {
  
  start_area1 = floor(start_interv/365)
  start_area2 = floor(n_days/365) + floor(start_interv/365)
  studied_duration_intervention = ceiling(length_intervention*nb_studied_cycles)
  end_area1 = start_area1 + studied_duration_intervention + 1
  end_area2 = start_area2 + studied_duration_intervention + 1
  
  sum_async_area1 = sum(incidence_async[c(start_area1:end_area1)])
  sum_async_area2 = sum(incidence_async[c(start_area2:end_area2)])
  sum_async = sum_async_area1 + sum_async_area2
  
  
  sum_sync_area1 = sum(incidence_sync[c(start_area1:end_area1)])
  sum_sync_area2 = sum(incidence_sync[c(start_area2:end_area2)])
  sum_sync = sum_sync_area1 + sum_sync_area2
  
  
  #Asynchrony Induced Growth Rate
  AIGR <- (sum_async - sum_sync)/sum_sync
  AIGR_area1 <- (sum_async_area1 - sum_sync_area1)/sum_sync_area1
  AIGR_area2 <- (sum_async_area2 - sum_sync_area2)/sum_sync_area2
  
  #Asynchrony Induced Growth Rate by year
  AIGR_year <- AIGR / studied_duration_intervention
  AIGR_year_area1 <- AIGR_area1 / studied_duration_intervention
  AIGR_year_area2 <- AIGR_area2 / studied_duration_intervention
  
  
  # Increase in infected in absolute values
  AIG <- round(sum_async - sum_sync)
  AIG_area1 <- round(sum_async_area1 - sum_sync_area1)
  AIG_area2 <- round(sum_async_area2 - sum_sync_area2)
  
  # By Year
  AIG_year <- round(AIG / studied_duration_intervention)
  AIG_year_area1 <- round(AIG_area1 / studied_duration_intervention) 
  AIG_year_area2 <- round(AIG_area2 / studied_duration_intervention)
  
  return(c("AIGR" = AIGR, "AIGR_area1" = AIGR_area1, "AIGR_area2" = AIGR_area2, 
           "AIGR_year" = AIGR_year, "AIGR_year_area1" = AIGR_year_area1, 
           "AIGR_year_area2" = AIGR_year_area2, 
           "AIG" = AIG, "AIG_area1" = AIG_area1, "AIG_area2" = AIG_area2,
           "AIG_year" = AIG_year, "AIG_year_area1" = AIG_year_area1,
           "AIG_year_area2" = AIG_year_area2))
}



metrics_computation <- function(myvars) {
  df <- data.frame("x0_1" = c(0), "x0_2" = c(0), "z0_1" = c(0), "z0_2" = c(0), 
                   "p_1" = c(0), "p_2" = c(0), "r"=c(0), "rinv" = c(0), 
                   "time_intervention" = c(0), "W_1" = c(0), "W_2" = c(0),
                   "R0_1" = c(0), "R0_2" = c(0), "RC_1" = c(0), "RC_2" = c(0), 
                   "AIGR" = c(0), "AIGR_area1" = c(0), "AIGR_area2" = c(0), 
                   "AIGR_year" = c(0), "AIGR_year_area1" = c(0), 
                   "AIGR_year_area2" = c(0), 
                   "AIG" = c(0), "AIG_area1" = c(0), "AIG_area2" = c(0),
                   "AIG_year" = c(0), "AIG_year_area1" = c(0),
                   "AIG_year_area2" = c(0))
  t <- Sys.time()
  for (i in 1:nrow(myvars)) {
    if (i%%100 == 0) {
      print(i)
      print(Sys.time()-t)}
    n_foci=2
    n_days=30000
    
    ## population
    N=c(1000, 1000)
    
    r = myvars[i,]$r
    rinv = myvars[i,]$rinv
    start_interv = 365
    R0 = c(myvars[i,]$R0_1, myvars[i,]$R0_2)
    p_1 = myvars[i,]$p_1
    p_2 = myvars[i,]$p_2
    
    W_1 = myvars[i,]$W_1
    W_2 = myvars[i,]$W_2
    
    x0 = compute_equilibrium_prevalence(R0, N, r, p_1, p_2)
    z0 = x0 *r*365*1000
    
    
    vars <- c("R0" = R0, "p_1" = p_1, 
              "p_2" = p_2, "W_1" = W_1, "W_2" = W_2, 
              "alpha" = rep(0, 2), "rho" = rep(1, 2),
              "z0" = z0, "x0" = x0, "r" = r)
    
    RC_1 = R0[1] * W_1
    RC_2 = R0[2] * W_2
    
    length_intervention = myvars[i, ]$time_intervention
    nb_studied_cycles = 6
    simulations = simulate_SIR(n_days, start_interv, length_intervention, nb_studied_cycles, N, 
                               vars)
    
    mysimul_0=simulations[[1]]
    mysimul_A=simulations[[2]]
    mysimul_B=simulations[[3]]
    mysimul_0A=simulations[[4]]
    mysimul_0B=simulations[[5]]
    
    
    metrics = compute_metrics(mysimul_B$annual_incidence$value, mysimul_A$annual_incidence$value, 
                              n_days, start_interv, length_intervention,
                              nb_studied_cycles)
    
    df[i+1,] <- c(x0[1], x0[2], z0[1], z0[2], p_1, p_2, r, rinv, 
                  length_intervention, W_1, W_2, vars["x01"], vars["x02"],
                  R0[1], R0[2], RC_1, RC_2, metrics)
    if (i%%2000 == 0)
    {
      write.csv(df, "/scicore/home/pothin/iggidr0000/simul_incidence.csv")
    }
  }
  
  return(df)
}

AIG1_computation <- function(myvars) {
  n = nrow(myvars)
  metric <- c()
  n_foci=2
  n_days=35000
  nb_p <- 0
  ## population
  N=c(1000, 1000)
  for (i in 1:n) {
    if (i%%100 == 0) {
      print(i)
      print(Sys.time()-t)}
    r = 1/myvars[i,]$rinv
    start_interv = 730
    R0 = c(myvars[i,]$R0_1, myvars[i,]$R0_2)
    p_1 = 1-myvars[i,]$p_1
    p_2 = 1-myvars[i,]$p_2
    
    W_1 = 1-myvars[i,]$W_1
    W_2 = 1-myvars[i,]$W_2
    x0 = compute_equilibrium_prevalence(R0, N, r, p_1, p_2)
    z0 = x0 *r*365*1000
    if (x0[1] <= 0 | x0[2] <= 0) {
      print("Fail")
      nb_p <- nb_p + 1
      metric <- c(metric, 0)
    }
    else {
      
      
      vars <- c("R0" = R0, "p_1" = p_1, 
                "p_2" = p_2, "W_1" = W_1, "W_2" = W_2, 
                "alpha" = rep(0, 2), "rho" = rep(1, 2),
                "z0" = z0, "x0" = x0, "r" = r)
      
      RC_1 = R0[1] * W_1
      RC_2 = R0[2] * W_2
      
      length_intervention = myvars[i, ]$time_intervention
      nb_studied_cycles = 6
      simulations = simulate_SIR(n_days, start_interv, length_intervention, nb_studied_cycles, N, 
                                 vars)
      mysimul_A=simulations[[1]]
      mysimul_B=simulations[[2]]
      
      metrics = compute_metrics(mysimul_B$annual_incidence$value, mysimul_A$annual_incidence$value, 
                                n_days, start_interv, length_intervention,
                                nb_studied_cycles)
      
      metric <- c(metric, metrics["AIG_year_area1"])
    }
  }
  print("Nb of under 0 prevalences :")
  print(nb_p)
  return(metric)
}


###########################
# Visualization functions #
###########################


visualize_prevalence <- function(prevalence, start_interv, length_intervention) {
  
  master_plot <- ggplot(prevalence , aes(x = t, y = value)) + 
    geom_line(aes(color = scenario),linetype = 1, linewidth=0.5, alpha=1) + 
    theme_minimal() + 
    ylab("prevalence") + xlab('days') +
    facet_wrap(variable ~. , scales = "free")+
    # xlim(19900, 23500)
    xlim(0, start_interv + 365*6*length_intervention + 365*5) +
    ylim(0, 0.6)
  master_plot
}


visualize_incidence <- function(incidence, start_year, length_intervention, year_start_plot, year_end_plot, max_y, gap, metrics_labels) {
  

  make_segments <- function(xmin, xmax, y, length, space, 
                            colors) {
    seqs <- seq(xmin, xmax, by = length + space)
    segs <- do.call(rbind, lapply(seqs, function(start) {
      ends <- start + length
      if (ends > xmax) return(NULL)
      data.frame(x = start, xend = ends, y = y, yend = y)
    }))
    segs$col <- rep(colors, length.out = nrow(segs))
    segs
  }
  
  seg1_graph1 <- make_segments(1, year_end_plot+1, y = max_y, length = length_intervention, space = length_intervention, colors = "Asynchronous")
  seg2_graph1 <- make_segments(1, year_end_plot+1, y = max_y+gap, length = length_intervention, space = length_intervention, colors = "Synchronous")
  
  seg1_graph2 <- make_segments(1+length_intervention, length_intervention + year_end_plot+1, y = max_y, length = length_intervention, space = length_intervention, colors = "Asynchronous")
  seg2_graph2 <- make_segments(1, year_end_plot+1, y = max_y+gap, length = length_intervention, space = length_intervention, colors = "Synchronous")
  
  
  seg1_graph1$area <- "Area1"
  seg2_graph1$area <- "Area1"
  seg1_graph2$area <- "Area2"
  seg2_graph2$area <- "Area2"
  
  segments_all <- bind_rows(seg1_graph1, seg2_graph1, seg1_graph2, seg2_graph2)
  
  
  incidence %>%
    ggplot(aes(x = t/365, y = value)) +
    geom_ribbon(
      data = incidence %>%
        tidyr::pivot_wider(names_from = scenario, values_from = value),
      aes(x = t/365, ymin = Synchronous, ymax = Asynchronous),
      fill = "burlywood",
      alpha = 0.4,
      inherit.aes = FALSE
    ) +
    geom_segment(data = segments_all,
                 aes(x = x, xend = xend, y = y, yend = yend, color = col),
                 inherit.aes = FALSE, linewidth = 2) +
    geom_line(aes(color = scenario), size = 1.5, alpha = 1) +
    facet_wrap(area ~ ., scales = "free", labeller = labeller(area = c("Area1" = "Area 1", "Area2" = "Area 2"))) +
    labs(y = "Annual incidence", x = "Years", color = "Scenario") +
    xlim(year_start_plot, year_end_plot) +
    ylim(0, NA) +
    theme_minimal() +
    theme(
      axis.text = element_text(size=12),
      strip.text = element_text(size=14),
      axis.title = element_text(size=14),
      legend.text = element_text(size=18),
      legend.title = element_text(size=18)
    ) + scale_color_manual(
      values = c(
        "Synchronous" = "black",
        "Asynchronous" = "tan2"
      )
    ) + geom_label(
      data = metrics_labels,
      aes(x = year_end_plot - 5, y = max_y *0.9, label = value, group = area),
      inherit.aes = FALSE,
      size = 5,
      hjust = 1,
      color = "burlywood",
      fill = "white", 
      label.size = 0.5,
      fontface = "bold"
    )
  
  }

plot_metric_by_parameter <- function(data, param, metric, bin_width, name_param) {
  
  param <- enquo(param)
  metric <- enquo(metric)
  
  summary_df <- data %>%
    mutate(param_val = !!param) %>%
    mutate(param_bin = cut(param_val, breaks = seq(floor(min(param_val)), ceiling(max(param_val)), by = bin_width))) %>%
    group_by(param_bin) %>%
    summarise(
      param_mid = mean(param_val),
      I_mean = mean(!!metric),
      I_d1 = quantile(!!metric, 0.1),
      I_d9 = quantile(!!metric, 0.9),
      I_q1 = quantile(!!metric, 0.25),
      I_q3 = quantile(!!metric, 0.75),
      I_min = min(!!metric),
      I_max = max(!!metric),
      .groups = "drop"
    )
  
  p <- ggplot(summary_df, aes(x = param_mid)) + 
    #geom_ribbon(aes(ymin = I_min, ymax = I_max), fill = "lightgrey", alpha = 0.4) +
    geom_ribbon(aes(ymin = I_d1, ymax = I_d9), fill = "lightblue", alpha = 0.7) +
    geom_ribbon(aes(ymin = I_q1, ymax = I_q3), fill = "steelblue", alpha = 0.4) +
    geom_line(aes(y = I_mean), color = "blue", size = 1) +
    labs(x = name_param) +
    theme_minimal() + theme(axis.title.y = element_blank()) +
    theme(axis.title.x = element_text(size = 10),
          axis.text.x = element_text(size = 7.5),
          axis.text.y = element_text(size = 7.5))
  
  return(p)
}
