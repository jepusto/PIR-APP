rm(list = ls())
library(ggplot2)
setwd("Simulations/Simulation 4 - MTS penalty function")
load("Sim4MTScombined.Rdata")

# results_all$k_pretty <- as.factor(with(results_all, ifelse(k_priors == 0, "Unpenalized", paste("k prior = ", k_priors))))
# 
# results_all$k_pretty <- with(results_all, factor(k_pretty, levels = levels(k_pretty)[c(4,1:3)]))

simulation_plot <- function(data, statistic, graph = "bias",
                            breaks = c(-Inf, -1,-.50,-.15,-0.05,-0.02,0.02, 0.05,
                                       .15, .50 ,1, Inf),
                            labels = c("< -100%","-100% to -50%","-50% to -15%", 
                                       "-15% to -5%", "-5% to -2%","-2% to 2%",
                                       "2% to 5%", "5% to 15%", "15% to 50%",
                                       "50% to 100%", " > 100%")){
  
  data_sub <- subset(data, stat == statistic)
   
    data_sub$cuts <- cut(data_sub[,graph], breaks = breaks, labels = labels)
  
    qplot(phi, zeta, fill = cuts,
        geom = "tile",
        data = data_sub) +
    facet_grid(K_intervals + theta ~ k_priors, labeller = "label_both") +
    scale_y_continuous(limits = c(.005, .60), breaks=c(.05, .10, .20, .30, .40, .50)) + 
    scale_x_continuous(limits = c(0.005, .60), breaks=c(.05, .10, .20, .30, .40, .50)) +
    scale_color_brewer(type = "seq") + labs(fill = "Bias") + theme_bw() + 
    theme(axis.text.x = element_text(angle=45, hjust = 1))
}

simulation_plot(data = results, statistic = "logit phi")
simulation_plot(data = results, statistic = "logit phi", graph = "median_bias")
simulation_plot(data = results, statistic = "logit phi", graph = "var", breaks = c(-.01, .5, 1, 1.5, 2, 10, 100, 200, 300, 4000, labels = NULL)

simulation_plot(data = results, statistic = "log zeta")
simulation_plot(data = results, statistic = "log zeta", graph = "median_bias")

simulation_plot(data = results, statistic = "phi", labels = NULL)
simulation_plot(data = results, statistic = "phi", graph = "median_bias", labels = NULL)

simulation_plot(data = results, statistic = "zeta", labels = NULL)
simulation_plot(data = results, statistic = "zeta", graph = "median_bias", labels = NULL)

simulation_plot(data = results, statistic = "mu", labels = NULL)
simulation_plot(data = results, statistic = "mu", graph = "median_bias", labels = NULL)

simulation_plot(data = results, statistic = "lambda", labels = NULL)
simulation_plot(data = results, statistic = "lambda", graph = "median_bias", labels = NULL)

simulation_plot(data = results, statistic = "log mu")
simulation_plot(data = results, statistic = "log mu", graph = "median_bias")

simulation_plot(data = results, statistic = "log lambda")
simulation_plot(data = results, statistic = "log lambda", graph = "median_bias")