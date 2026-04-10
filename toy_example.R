# ==============================================================================
#  Toy Illustrative example of AIG using the SIS model
#
#                 Younes Iggidr
#                     2025
#
#    _____         _           _______ _____  _    _ 
#   / ____|       (_)         |__   __|  __ \| |  | |
#  | (_____      ___ ___ ___     | |  | |__) | |__| |
#   \___ \ \ /\ / / / __/ __|    | |  |  ___/|  __  |
#   ____) \ V  V /| \__ \__ \    | |  | |    | |  | |
#  |_____/ \_/\_/ |_|___/___/    |_|  |_|    |_|  |_|
#
# ==============================================================================


rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
options(scipen = 7)
source("utils.r")
n_foci=2
n_years=100
n_days=365*n_years


###############
## Parameters #
###############

## population
N=c(1000, 1000)

## connectivity matrix
p_1 = 0.9
p_2 = 0.99

## Interventions
start_year = 1
start_interv = start_year*365 
length_intervention = 5
nb_studied_cycles = 6

## vector control
W_1 = 0.5
W_2 = 0.5

## case management: effective coverage
alpha = rep(0, n_foci)

## effective reporting
rho = rep(1, n_foci)

## recovery rate
r=1/200

## R0
R0 = c(1.1, 1.1)

## initial prevalence in each foci
x0 = compute_equilibrium_prevalence(R0, N, r, p_1, p_2)
z0 = x0*r*365*1000

vars <- c("R0" = R0, "p_1" = p_1, 
          "p_2" = p_2, "W_1" = W_1, "W_2" = W_2, 
          "alpha" = alpha, "rho" = rho, "x0" = x0, "z0" = z0, "r" = r)

################
## Simulations #
################

simulations = simulate_SIR(n_days, start_interv, length_intervention, nb_studied_cycles, N, 
                           vars)
mysimul_A=simulations[[1]]
mysimul_B=simulations[[2]]
all_incidence=rbind(mysimul_A$annual_incidence %>% mutate(scenario="Synchronous"),
                    mysimul_B$annual_incidence %>% mutate(scenario="Asynchronous")
) %>%
  mutate(area=ifelse(variable=="z[1]", "Area1", "Area2"))

all_prevalence=rbind(mysimul_A$prevalence %>% mutate(scenario="Synchronous"),
                     mysimul_B$prevalence %>% mutate(scenario="Asynchronous"))

########################
## Metrics computation #
########################

x = compute_metrics(mysimul_B$annual_incidence$value, mysimul_A$annual_incidence$value, 
                    n_days, start_interv, length_intervention,
                    nb_studied_cycles)

#############################
## Prevalence Visualization #
#############################

visualize_prevalence(all_prevalence, start_interv, length_intervention)

############################
## Incidence Visualization #
############################


metrics_labels = data.frame(x[c("AIG_area1", "AIG_area2")])
colnames(metrics_labels) <- "value"
metrics_labels["area"] <- c("Area1","Area2")
metrics_labels <- metrics_labels %>% mutate(value = paste0("AIG in this area : ", value))

visualize_incidence(all_incidence, start_year, length_intervention, year_start_plot = 0, year_end_plot = 30, max_y = 170, gap = 5, metrics_labels = metrics_labels)

############################
## Time before elimination #
############################


elimination_time_synchronous_area1 = min(mysimul_A$annual_incidence[mysimul_A$annual_incidence$value < 10^(-5) & mysimul_A$annual_incidence$variable == "z[1]",]$t/365)
elimination_time_synchronous_area2 = min(mysimul_A$annual_incidence[mysimul_A$annual_incidence$value < 10^(-5) & mysimul_A$annual_incidence$variable == "z[2]",]$t/365)
elimination_time_synchronous = max(elimination_time_synchronous_area1, elimination_time_synchronous_area2)
elimination_time_asynchronous_area1 = min(mysimul_B$annual_incidence[mysimul_B$annual_incidence$value < 10^(-5) & mysimul_A$annual_incidence$variable == "z[1]",]$t/365)
elimination_time_asynchronous_area2 = min(mysimul_B$annual_incidence[mysimul_B$annual_incidence$value < 10^(-5) & mysimul_A$annual_incidence$variable == "z[2]",]$t/365)
elimination_time_asynchronous = max(elimination_time_asynchronous_area1, elimination_time_asynchronous_area2)


nb_case_before_elim_sync = sum(mysimul_A$annual_incidence$value[mysimul_A$annual_incidence$t/365 < elimination_time_synchronous])
nb_case_before_elim_sync_area1 = sum(mysimul_A$annual_incidence$value[mysimul_A$annual_incidence$t/365 < elimination_time_synchronous & mysimul_A$annual_incidence$variable == "z[1]"])
nb_case_before_elim_sync_area2 = sum(mysimul_A$annual_incidence$value[mysimul_A$annual_incidence$t/365 < elimination_time_synchronous & mysimul_A$annual_incidence$variable == "z[2]"])
nb_case_before_elim_async = sum(mysimul_B$annual_incidence$value[mysimul_B$annual_incidence$t/365 < elimination_time_asynchronous])
nb_case_before_elim_async_area1 = sum(mysimul_B$annual_incidence$value[mysimul_B$annual_incidence$t/365 < elimination_time_asynchronous & mysimul_A$annual_incidence$variable == "z[1]"])
nb_case_before_elim_async_area2 = sum(mysimul_B$annual_incidence$value[mysimul_B$annual_incidence$t/365 < elimination_time_asynchronous & mysimul_A$annual_incidence$variable == "z[2]"])

print(x)
print(elimination_time_asynchronous - elimination_time_synchronous)

print(elimination_time_synchronous)
print(elimination_time_synchronous_area1)
print(elimination_time_synchronous_area2)

print(elimination_time_asynchronous)
print(elimination_time_asynchronous_area1)
print(elimination_time_asynchronous_area2)


print(nb_case_before_elim_sync)
print(nb_case_before_elim_sync_area1)
print(nb_case_before_elim_sync_area2)
print(nb_case_before_elim_async)
print(nb_case_before_elim_async_area1)
print(nb_case_before_elim_async_area2)
