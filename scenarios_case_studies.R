rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
options(scipen = 7)
source("utils.r")
n_foci=2
n_days=100000


#############################
## Parameters of scenario A #
#############################

## population
N=c(1000, 1000)

## connectivity matrix
p_1 = 0.8
p_2 = 0.8

## Interventions
start_year = 2
start_interv = start_year*365 
length_intervention = 1.5
nb_studied_cycles = 6

## vector control
W_1 = 0.6
W_2 = 0.6

## case management: effective coverage
alpha = rep(0, n_foci)

## effective reporting
rho = rep(1, n_foci)

## recovery rate
r=1/200

## Initial prevalence
x0 = c(0.3, 0.3)

## R0
R0 = compute_R0(x0, N, r, p_1, p_2) #R0 = 1.43 for both

## Initial incidence
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


############
## Metrics #
############

A = compute_metrics(mysimul_B$annual_incidence$value, mysimul_A$annual_incidence$value, 
                    n_days, start_interv, length_intervention,
                    nb_studied_cycles)

############
## Figures #
############


metrics_labels = data.frame(A[c("AIG_area1", "AIG_area2")])
colnames(metrics_labels) <- "value"
metrics_labels["area"] <- c("Area1","Area2")
metrics_labels <- metrics_labels %>% mutate(value = paste0("AIG in this area : ", value))

inc_1 = visualize_incidence(all_incidence, start_year, length_intervention, year_start_plot = 0, year_end_plot = 25, max_y = 560, gap = 20, metrics_labels = metrics_labels)



#############################
## Parameters of scenario B #
#############################

## population
N=c(1000, 1000)

## connectivity matrix
p_1 = 0.9
p_2 = 0.99

## Interventions
start_year = 2
start_interv = start_year*365 
length_intervention = 3
nb_studied_cycles = 6

## vector control
W_1 = 0.5
W_2 = 0.7

## case management: effective coverage
alpha = rep(0, n_foci)

## effective reporting
rho = rep(1, n_foci)

## recovery rate
r=1/150

## Initial prevalence
x0 = c(0.2, 0.2)

## R0
R0 = compute_R0(x0, N, r, p_1, p_2) #R0 = 1.25 for both

## Initial incidence
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


############
## Metrics #
############


B = compute_metrics(mysimul_B$annual_incidence$value, mysimul_A$annual_incidence$value, 
                    n_days, start_interv, length_intervention,
                    nb_studied_cycles)

############
## Figures #
############


metrics_labels = data.frame(B[c("AIG_area1", "AIG_area2")])
colnames(metrics_labels) <- "value"
metrics_labels["area"] <- c("Area1","Area2")
metrics_labels <- metrics_labels %>% mutate(value = paste0("AIG in this area : ", value))

inc_2 = visualize_incidence(all_incidence, start_year, length_intervention, year_start_plot = 0, year_end_plot = 25, max_y = 500, gap = 20, metrics_labels = metrics_labels)


#############################
## Parameters of scenario C #
#############################

## population
N=c(1000, 1000)

## connectivity matrix
p_1 = 0.9
p_2 = 0.99

## Interventions
start_year = 2
start_interv = start_year*365 
length_intervention = 1.5
nb_studied_cycles = 6

## vector control
W_1 = 0.55
W_2 = 0.55

## case management: effective coverage
alpha = rep(0, n_foci)

## effective reporting
rho = rep(1, n_foci)

## recovery rate
r=1/60

## Initial prevalence
x0 = c(0.2, 0.2)

## R0
R0 = compute_R0(x0, N, r, p_1, p_2) #R0 = 1.25 for both

## Initial incidence
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


############
## Metrics #
############


C = compute_metrics(mysimul_B$annual_incidence$value, mysimul_A$annual_incidence$value, 
                    n_days, start_interv, length_intervention,
                    nb_studied_cycles)

############
## Figures #
############


metrics_labels = data.frame(C[c("AIG_area1", "AIG_area2")])
colnames(metrics_labels) <- "value"
metrics_labels["area"] <- c("Area1","Area2")
metrics_labels <- metrics_labels %>% mutate(value = paste0("AIG in this area : ", value))

inc_3 = visualize_incidence(all_incidence, start_year, length_intervention, year_start_plot = 0, year_end_plot = 25, max_y = 1250, gap = 40, metrics_labels = metrics_labels)


#############################
## Parameters of scenario D #
#############################

## population
N=c(1000, 1000)

## connectivity matrix
p_1 = 0.9/0.91
p_2 = 0.9/0.93

## Interventions
start_year = 1
start_interv = start_year*365 
length_intervention = 5
nb_studied_cycles = 6

## vector control
W_1 = 0.7
W_2 = 0.7

## case management: effective coverage
alpha = rep(0, n_foci)

## effective reporting
rho = rep(1, n_foci)

## recovery rate
r=1/60

## R0
R0 = c(1.05, 1.05)

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

############
## Metrics #
############


D = compute_metrics(mysimul_B$annual_incidence$value, mysimul_A$annual_incidence$value, 
                    n_days, start_interv, length_intervention,
                    nb_studied_cycles)

############
## Figures #
############


metrics_labels = data.frame(D[c("AIG_area1", "AIG_area2")])
colnames(metrics_labels) <- "value"
metrics_labels["area"] <- c("Area1","Area2")
metrics_labels <- metrics_labels %>% mutate(value = paste0("AIG in this area : ", value))

inc_4 = visualize_incidence(all_incidence, start_year, length_intervention, year_start_plot = 0, year_end_plot = 25, max_y = 300, gap = 10, metrics_labels = metrics_labels)



######################################################
## Patchwork of the scenarios' incidence in a column #
######################################################


inc_1 <- inc_1 + labs(x = NULL, y = "Annual incidence") 
inc_2 <- inc_2 + labs(x = NULL, y = "Annual incidence") + theme(strip.text = element_blank())
inc_3 <- inc_3 + labs(x = NULL, y = "Annual incidence") + theme(strip.text = element_blank())
inc_4 <- inc_4 + labs(x = "Years", y = "Annual incidence") + theme(strip.text = element_blank())

(inc_1 / inc_2 / inc_3 / inc_4) + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A', tag_suffix = ".") & 
  theme(legend.position = "bottom", plot.tag = element_text(face = "bold"))


###################
##  Print metrics #
###################

print(A)
print(B)
print(C)
print(D)
