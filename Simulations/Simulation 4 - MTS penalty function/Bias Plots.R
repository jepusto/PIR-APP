rm(list = ls())
library(ggplot2)
setwd("Simulations/Simulation 4 - MTS penalty function")
load("sim4MTS_unp.Rdata")
load("sim4MTSpenalized.Rdata")
results_unp$theta <- 10
results_unp$k_priors <- 0
results_all <- rbind(results, results_unp)
results_unp$theta <- 20
results_all <- rbind(results_all, results_unp)

results_all$k_pretty <- as.factor(with(results_all, ifelse(k_priors == 0, "Unpenalized", paste("k prior = ", k_priors))))

results_all$k_pretty <- with(results_all, factor(k_pretty, levels = levels(k_pretty)[c(4,1:3)]))

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
    facet_grid(theta ~ k_pretty) +
    scale_y_continuous(limits = c(.005, .60), breaks=c(.05, .10, .20, .30, .40, .50)) + 
    scale_x_continuous(limits = c(0.005, .60), breaks=c(.05, .10, .20, .30, .40, .50)) +
    scale_color_brewer(type = "seq") + labs(fill = "Bias") + theme_bw() + 
    theme(axis.text.x = element_text(angle=45, hjust = 1))
}

simulation_plot(data = results_all, statistic = "logit phi")
simulation_plot(data = results_all, statistic = "logit phi", graph = "median_bias")

simulation_plot(data = results_all, statistic = "log zeta")
simulation_plot(data = results_all, statistic = "log zeta", graph = "median_bias")

simulation_plot(data = results_all, statistic = "phi", labels = NULL)
simulation_plot(data = results_all, statistic = "phi", graph = "median_bias", labels = NULL)

simulation_plot(data = results_all, statistic = "zeta", labels = NULL)
simulation_plot(data = results_all, statistic = "zeta", graph = "median_bias", labels = NULL)

simulation_plot(data = results_all, statistic = "mu", labels = NULL)
simulation_plot(data = results_all, statistic = "mu", graph = "median_bias", labels = NULL)

simulation_plot(data = results_all, statistic = "lambda", labels = NULL)
simulation_plot(data = results_all, statistic = "lambda", graph = "median_bias", labels = NULL)

simulation_plot(data = results_all, statistic = "log mu")
simulation_plot(data = results_all, statistic = "log mu", graph = "median_bias")

simulation_plot(data = results_all, statistic = "log lambda")
simulation_plot(data = results_all, statistic = "log lambda", graph = "median_bias")