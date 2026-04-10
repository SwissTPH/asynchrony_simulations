# ==============================================================================
#  SIS connected model (meta population) for malaria transmission
#
#  Model Implementation (deterministic version) with Odin library
#
#  created by : Bilal Benhana & Clara Champagne
#  originally created : 2024
#
#    _____         _           _______ _____  _    _ 
#   / ____|       (_)         |__   __|  __ \| |  | |
#  | (_____      ___ ___ ___     | |  | |__) | |__| |
#   \___ \ \ /\ / / / __/ __|    | |  |  ___/|  __  |
#   ____) \ V  V /| \__ \__ \    | |  | |    | |  | |
#  |_____/ \_/\_/ |_|___/___/    |_|  |_|    |_|  |_|
#
# ==============================================================================

#=================================== IMPORTS ===================================
library(odin)
library(tidyverse)

#===================== DETERMINISTIC MODEL IMPLEMENTATION ======================

connectedSIS_deterministic_generator <- odin::odin({
  
  ## Equations:
  deriv(x[]) <- r * ((1-x[i]) * (1-alpha_t[i]) * sum(pwlk[i, ]) - x[i])
  deriv(z[]) <-  r * (1-x[i]) * rho_t[i] * sum(pwlkN[i, ])
  
  ## Initial Conditions:
  initial(x[]) <- x0[i]
  initial(z[]) <- z0[i]
  
  ## Input:
  N[] <- user()
  W[,] <- user()
  alpha[,] <- user()
  rho[,] <- user()
  tt[] <- user()
  p[, ] <- user()
  x0[] <- user()
  z0[] <- user()
  R0[] <- user()
  r <- user()
  
  ## Output:
  output(s) <- s
  
  ## Variables definition:
  W_t[] <- interpolate(tt, W, "linear")
  alpha_t[] <- interpolate(tt, alpha, "linear")
  rho_t[] <- interpolate(tt, rho, "linear")
  s[] <- 1 - x[i]
  pN[, ] <- p[i,j] * N[i]
  pxN[ ,] <- pN[i,j] * x[i]
  K[] <- sum(pxN[, i]) / sum(pN[, i])
  pwlk[, ] <- p[i, j] * W_t[j] * R0[j] * K[j]
  pwlkN[, ] <- pwlk[i, j] * N[i]
  
  ## Dimensions:
  dim(N) <- user()
  n <- length(N)
  
  dim(W_t) <- n
  dim(alpha_t) <- n
  dim(rho_t) <- n
  dim(tt) <- user()
  dim(W) <- c(length(tt), length(W_t))
  dim(alpha) <- c(length(tt), length(alpha_t))
  dim(rho) <- c(length(tt), length(rho_t))
  dim(x) <- n
  dim(s) <- n
  dim(z) <- n
  dim(x0) <- n
  dim(z0) <- n
  dim(p) <- c(n, n)
  dim(pN) <- c(n, n)
  dim(pxN) <- c(n, n)
  dim(K) <- n
  dim(R0) <- n
  dim(pwlk) <- c(n, n)
  dim(pwlkN) <- c(n, n)
  
})

run_SISconnected_deterministic=function(parameters, n_days){
  
  t <- 0:n_days
  parameters$tt = 0:(n_days+1)
  n_foci=length(parameters$N)
  
  if(!is.matrix(parameters$W)){
    parameters$W= matrix(rep(parameters$W, n_days+2), nrow=n_days+2, ncol=n_foci, byrow = T)
  }
  
  if(!is.matrix(parameters$alpha)){
    parameters$alpha= matrix(rep(parameters$alpha, n_days+2), nrow=n_days+2, ncol=n_foci, byrow = T)
  }
  
  if(!is.matrix(parameters$rho)){
    parameters$rho= matrix(rep(parameters$rho, n_days+2), nrow=n_days+2, ncol=n_foci, byrow = T)
  }
  
  connectedSIS_model <- connectedSIS_deterministic_generator$new(user = parameters)
  
  y <- connectedSIS_model$run(t)
  
  
  ## Data preparation
  
  df <- as.data.frame(y)
  
  prevalence <- df %>%
    select(t, starts_with("x")) %>%
    gather(key = "variable", value = "value", -t)
  
  
  annual_incidence =  df %>%
    select(t, starts_with("z")) %>%
    gather(key = "variable", value = "value", -t)  %>%
    mutate(year=ceiling(t/365))%>%
    filter(t %%365 ==0) %>%
    group_by(variable)%>%
    mutate(value=c(0,diff(value)))
  annual_incidence$value[1] = annual_incidence$value[2]
  annual_incidence[annual_incidence$variable == "z[2]",]$value[1] <- annual_incidence[annual_incidence$variable == "z[2]",]$value[2]
  
  pop=data.frame(pop=parameters$N, variable=unique(annual_incidence$variable))
  
  annual_incidence%>% left_join(pop)%>%
    mutate(incidence=value*1000/pop)%>%
    rename("reportedCases"=value)
  
  return(list("simulation"=df, "annual_incidence"=annual_incidence, "prevalence"=prevalence))
}


create_timevarying_interventions=function(interventions, interv_times){
  if(interv_times[1]!=0){
    stop("please add 0 as first intervention time")
  }
  myvect=c()
  for(i in 1:(length(interv_times)-1)){
    myvect=c(myvect, rep(interventions[i,], (interv_times[i+1]-interv_times[i])))
  }
  output=matrix(myvect, 
                nrow=interv_times[length(interv_times)], ncol=ncol(interventions), byrow = T)
  return(output)
}



#===============================================================================
# 
# 
#                     ███████╗███╗░░██╗██████╗░
#                     ██╔════╝████╗░██║██╔══██╗
#                     █████╗░░██╔██╗██║██║░░██║
#                     ██╔══╝░░██║╚████║██║░░██║
#                     ███████╗██║░╚███║██████╔╝
#                     ╚══════╝╚═╝░░╚══╝╚═════╝░
#
#
#===============================================================================
