# ==============================================================================
#  Uncertainty/Sensitivity Analysis and Decision Tree
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


rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source('utils.r')

##############
# LHS design #
##############

mysize=10000
X1 = data.frame(as.matrix( maximinLHS(n = mysize, k = 8)))
colnames(X1) = c("R0_1", "R0_2", "rinv", "time_intervention", "W_1", "W_2", "p_1", "p_2")
X1$rinv = floor(60 + X1$rinv * (201 - 60))
X1$time_intervention = floor(2*(0.5 + X1$time_intervention * 5))/2
X1$R0_1 = 0.9 + X1$R0_1 * (2.2 - 0.9)
X1$W_1 = 0.5 + X1$W_1 * (1-0.5)
X1$p_1 = 0.5 + X1$p_1 * 0.5
X1$R0_2 = 0.9 + X1$R0_2 * (2.2 - 0.9)
X1$W_2 = 0.5 + X1$W_2 * (1-0.5)
X1$p_2 = 0.5 + X1$p_2 * 0.5
X1$r = 1/X1$r


mysize=10000
X2 = data.frame(as.matrix( maximinLHS(n = mysize, k = 8)))
colnames(X2) = c("R0_1", "R0_2", "rinv", "time_intervention", "W_1", "W_2", "p_1", "p_2")
X2$rinv = floor(60 + X2$rinv * (201 - 60))
X2$time_intervention = floor(2*(0.5 + X2$time_intervention * 5))/2
X2$R0_1 = 0.9 + X2$R0_1 * (2.2 - 0.9)
X2$W_1 = 0.5 + X2$W_1 * (1-0.5)
X2$p_1 = 0.5 + X2$p_1 * 0.5
X2$R0_2 = 0.9 + X2$R0_2 * (2.2 - 0.9)
X2$W_2 = 0.5 + X2$W_2 * (1-0.5)
X2$p_2 = 0.5 + X2$p_2 * 0.5
X2$r = 1/X2$r

########################
# UNCERTAINTY AND PRCC #
########################

# Simulations database creation
df = metrics_computation(X1)

df$p_1 = 1-df$p_1
df$p_2 = 1-df$p_2
df$W_1 = 1-df$W_1
df$W_2 = 1-df$W_2

# Evolution of AIG by year as a function of a parameter

## R0_1
bin_width <- 0.1
p1 <- plot_metric_by_parameter(df, R0_1, AIG_year_area1, bin_width, "R0 in area 1")

## R0_2
bin_width <- 0.1
p2 <- plot_metric_by_parameter(df, R0_2, AIG_year_area1, bin_width, "R0 in area 2")

## p_1
bin_width <- 0.05
p3 <- plot_metric_by_parameter(df, p_1, AIG_year_area1, bin_width, expression(p["1,2"]))

## p_2
bin_width <- 0.05
p4 <- plot_metric_by_parameter(df, p_2, AIG_year_area1, bin_width, expression(p["2,1"]))

## RC_1
bin_width <- 0.05
p5 <- plot_metric_by_parameter(df, RC_1, AIG_year_area1, bin_width, expression(R[C]~" in area 1"), T)

## RC_2
bin_width <- 0.05
p6 <- plot_metric_by_parameter(df, RC_2, AIG_year_area1, bin_width, expression(R[C]~" in area 2"))

## W_1
bin_width <- 0.05
p7 <- plot_metric_by_parameter(df, W_1, AIG_year_area1, bin_width, "Intervention efficiency in area 1")

## W_2
bin_width <- 0.05
p8 <- plot_metric_by_parameter(df, W_2, AIG_year_area1, bin_width, "Intervention efficiency in area 2")

## rinv
bin_width <- 10
p9 <- plot_metric_by_parameter(df, rinv, AIG_year_area1, bin_width, "Duration of recovery")

## time_intervention
bin_width <- 0.5
p10 <- plot_metric_by_parameter(df, time_intervention, AIG_year_area1, bin_width, "Duration of the intervention")

## Grid plot
lineplots <- (
  (p1 | p2) / (p3 | p4) / 
    (p5 | p6) / (p7 | p8 ) / 
    (p9 | p10)
) + plot_layout(guides='collect')+ plot_annotation(
  title = 'C.', theme = theme(plot.title = element_text(face = "bold", size = 10)))


lineplots <- (
  (p1 | p2 | p3 | p4 | p5) / (p6 | p7 | p8 | p9 | p10)
) + plot_layout(guides='collect')+ plot_annotation(
  title = 'C.', theme = theme(plot.title = element_text(face = "bold", size = 10)))

lineplots <- (
  (p1 | p2 | p3 | p4) / (p5 | p6) / (p7 | p8 | p9 | p10)
) + plot_layout(guides='collect', axis_titles='collect_y')+ plot_annotation(
  title = 'C.', theme = theme(plot.title = element_text(face = "bold", size = 10)))

# Histogram of AIG per year distribution

first_centile <- quantile(df$AIG_year_area1, 0.01)
last_centile <- quantile(df$AIG_year_area1, 0.99)
quantile75 <- quantile(df$AIG_year_area1, 0.75)
quantile90 <- quantile(df$AIG_year_area1, 0.90)

df_filtered <- df[df$AIG_year_area1 >= first_centile & df$AIG_year_area1 <= last_centile,]

hist <- ggplot(df_filtered, aes(x = AIG_year_area1)) +
  geom_histogram(binwidth = 5, boundary = 0, fill = "skyblue", color = "white") +
  labs(x = "AIG per year in Area 1\n(cases per year per 1000 individuals)", 
       y = "Frequency",
       title = "A.") +
  theme(axis.title.y = element_blank(),
        plot.title = element_text(face = "bold", size = 10),
        axis.title.x = element_text(size = 10),
        axis.text.x = element_text(size = 7.5),
        axis.text.y = element_text(size = 7.5)) +
  scale_x_continuous(breaks = seq(0, 200, by = 10)) +
  geom_vline(aes(xintercept = quantile75))+
  annotate("text",x=quantile75+11, y=2000,label="q75")+
  geom_vline(aes(xintercept = quantile90))+
  annotate("text",x=quantile90+11, y=2000,label="q90")+
  theme_bw()
  

# Histogram low & high

first_decile <- quantile(df$AIG_year_area1, 0.1)
last_decile <- quantile(df$AIG_year_area1, 0.9)

df$label <- cut(df$AIG_year_area1,
                breaks = c(-Inf, first_decile, last_decile, Inf),
                labels = c("Low", "Medium", "High"),
                right = FALSE)

theme_no_y <- theme(
  axis.title.y = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank()
)

df_filtered <- df %>% filter(label != "Medium")

hist_rinv <- ggplot(df_filtered, aes(x = rinv, fill = label)) +
  geom_histogram(position = "identity", alpha = 0.8, boundary = 60, binwidth = 5, color = "black") +
  labs(title = "Histogram of the recovery time",
       fill = "Magnitude of AIG",
       x = "Recovery time (in days)") +
  scale_fill_manual(values = c("High" = "gold", "Low" = "green3")) +
  theme_no_y

hist_W1 <- ggplot(df_filtered, aes(x = W_1, fill = label)) +
  geom_histogram(position = "identity", alpha = 0.8, boundary = 0, binwidth = 0.02, color = "black") +
  labs(title = "Histogram of intervention efficiency in area 1",
       fill = "Magnitude of AIG", 
       x = "Intervention efficiency") +
  theme_no_y +
  scale_fill_manual(values = c("High" = "gold", "Low" = "green3")) +
  theme(legend.position = "none")

hist_W2 <- ggplot(df_filtered, aes(x = W_2, fill = label)) +
  geom_histogram(position = "identity", alpha = 0.8, boundary = 0, binwidth = 0.02, color = "black") +
  labs(title = "Histogram of intervention efficiency in area 2",
       x = "Intervention efficiency") +
  theme_no_y +
  scale_fill_manual(values = c("High" = "gold", "Low" = "green3")) +
  theme(legend.position = "none")

hist_RC_1 <- ggplot(df_filtered, aes(x = RC_1, fill = label)) +
  geom_histogram(position = "identity", alpha = 0.8, boundary = 1, binwidth = 0.1, color = "black") +
  labs(title = expression("Histogram of" ~ R[C] ~ "in area 1"),
       x = "Intervention efficiency") +
  theme_no_y +
  scale_fill_manual(values = c("High" = "gold", "Low" = "green3")) +
  theme(legend.position = "none")

hist_RC_2 <- ggplot(df_filtered, aes(x = RC_2, fill = label)) +
  geom_histogram(position = "identity", alpha = 0.8, boundary = 1, binwidth = 0.1, color = "black") +
  labs(title = expression("Histogram of" ~ R[C] ~ "in area 2"),
       x = "Intervention efficiency") +
  theme_no_y +
  scale_fill_manual(values = c("High" = "gold", "Low" = "green3")) +
  theme(legend.position = "none")

hist_rinv / (hist_W1 + hist_W2) / (hist_RC_1 + hist_RC_2) + plot_layout(guides = "collect") 


#########
# SOBOL #
#########

myvars <- c("rinv", "time_intervention", "R0_1", "R0_2", "W_1", "W_2", "p_1", "p_2")

#soboljansen_AIG = soboljansen(model = AIG1_computation, X1[myvars], X2[myvars], nboot = 100)

df_sobol_first_order = soboljansen_AIG$S
df_sobol_first_order$index = "first order"

df_sobol_total = soboljansen_AIG$T
df_sobol_total$index = "total"

sobol_AIG=rbind(df_sobol_first_order, df_sobol_total)

pl_first_total <- sobol_AIG %>%
  mutate(param = factor(X, levels = c("rinv", "time_intervention", "R0_1", 
                                      "R0_2", "W_1", "W_2", "p_1", "p_2"))) %>%
  ggplot()+
  geom_point(aes(x=original, y=param, col=index, shape=index), size=2)+
  geom_linerange(aes(xmin=min..c.i., xmax=max..c.i., y=param, col=index), size = 1)+
  labs(y='',x='Sobol indices', shape="", col="", title = "B.") +
  theme(plot.title = element_text(face = "bold", size = 10), axis.title.x = element_text(size = 10)) +
  scale_color_manual(values=c( "cyan3","darkblue"))+
  scale_y_discrete(labels = c("rinv"="Recovery time","time_intervention"="Length of intervention",
                              "R0_1"="R0 in area 1", "R0_2"="R0 in area 2",
                              "W_1"="Intervention efficiency in area 1", "W_2"="Intervention efficiency in area 2",
                              "p_1"=expression(p["1,2"]), "p_2"=expression(p["2,1"])))+
  theme_bw()+
  theme(legend.position = "bottom", legend.box="vertical", legend.text = element_text(size = 10), legend.title =  element_blank())+
  theme(strip.background =element_rect(fill="white", color="white"))+
  theme(legend.position=c(0.83,1),#c(0.83,1),
        legend.background=element_rect(colour = 1),
        legend.box = "horizontal")
pl_first_total
###################################
## Patchwork of article's charts  #
###################################

layout <- "
AABBBB
#CCCCC
"
pl_first_total

Figure3=(hist + pl_first_total) / free(wrap_elements(lineplots), side = "l") + plot_layout(heights = c(0.4, 1))
Figure3
ggsave(Figure3,filename ="Figure3.png" ,
       width=9.5, height=9)

lineplots


#################
# Decision Tree #
#################

df$label_ML <- ifelse(df$AIG_year_area1 >= last_decile, "High", "Low")
n_high <- sum(df$label_ML == "High")
n_low <- sum(df$label_ML == "Low")
n_total <- length(df$label_ML)


w <- ifelse(df$label_ML == "High", n_total/(2*n_high), n_total/(2*n_low))
model_tree <- rpart(label_ML ~ R0_1 + R0_2 + W_1 + W_2 + rinv + p_1 + p_2 + time_intervention, data = df, weights = w)

rename_vars <- c(
  "rinv" = "Duration of recovery",
  "W_1" = "Intervention efficiency in area 1",
  "W_2" = "Intervention efficiency in area 2",
  "time_intervention"  = "Duration of intervention",
  "R0_1" = "R0 in area 1",
  "<leaf>" = "<leaf>"
)

model_tree$frame$var <- rename_vars[ model_tree$frame$var ]
pal <- colorRampPalette(c("gold", "green3"))(4)
rpart.plot(model_tree, extra = 104, main="Decision tree", cex = 1, box.palette = pal)
