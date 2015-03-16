library(ggplot2)
library(plyr)
rm(list = ls())
load("Simulations/Simulation 5 - PIR penalty Function/s5results_pen.Rdata")

K120results <- subset(s5results_pen, K_intervals == 120)

simulation_plot <- function(data, statistic, graph = "bias",
                            breaks = c(-Inf, -1,-.50,-.10,-0.02,0.02,
                                       .10, .50 ,1, Inf),
                            labels = c("< -100%","-100% to -50%","-50% to -10%", 
                                       "-10% to -2%","-2% to 2%",
                                       "2% to 10%", "10% to 50%",
                                       "50% to 100%", " > 100%")){
  
  data_sub <- subset(data, stat == statistic)
  
  data_sub$cuts <- cut(data_sub[,graph], breaks = breaks, labels = labels)
  
  print(qplot(phi, zeta, fill = cuts,
        geom = "tile",
        main = paste(graph, "small multiples for K =", data_sub$K_intervals[1]),
        data = data_sub) +
    facet_grid(theta ~ k_priors, labeller = "label_both") +
    scale_y_continuous(limits = c(.005, .60), breaks=c(.05, .10, .20, .30, .40, .50)) + 
    scale_x_continuous(limits = c(0.005, .9555), breaks=seq(.10, 1, .10)) +
    scale_color_brewer(type = "seq") + labs(fill = paste0("Bias in ", statistic)) + theme_bw() + 
    theme(axis.text.x = element_text(angle=45, hjust = 1)))
}

#log parameterizations
simulation_plot(K120results, statistic = "logit phi")
simulation_plot(K120results, statistic = "logit phi", graph = "median_bias")

simulation_plot(K120results, statistic = "log zeta")
simulation_plot(K120results, statistic = "log zeta", graph = "median_bias")

#natural parameterizations
simulation_plot(K120results, statistic = "phi", labels = c("< -1","-1 to -.50",
                                                           "-.5 to -.10", 
                                                           "-.10 to -.02",
                                                           "-.02 to .02",
                                                           ".02 to .10", 
                                                           ".10 to .50",
                                                           ".50 to 1", " > 1"))
simulation_plot(K120results, statistic = "phi", graph = "median_bias",
                labels = c("< -1","-1 to -.50",
                           "-.5 to -.10", 
                           "-.10 to -.02",
                           "-.02 to .02",
                           ".02 to .10", 
                           ".10 to .50",
                           ".50 to 1", " > 1"))

simulation_plot(K120results, statistic = "zeta", labels = c("< -1","-1 to -.50",
                                                                "-.5 to -.10", 
                                                                "-.10 to -.02",
                                                                "-.02 to .02",
                                                                ".02 to .10", 
                                                                ".10 to .50",
                                                                ".50 to 1", " > 1"))
simulation_plot(K120results, statistic = "zeta", graph = "median_bias",
                labels = c("< -1","-1 to -.50",
                           "-.5 to -.10", 
                           "-.10 to -.02",
                           "-.02 to .02",
                           ".02 to .10", 
                           ".10 to .50",
                           ".50 to 1", " > 1"))

# d_ply(s5results_pen, .(K_intervals), simulation_plot, statistic = "logit phi")
# d_ply(s5results_pen, .(K_intervals), simulation_plot, statistic = "logit phi", graph = "median_bias")
# 
# d_ply(s5results_pen, .(K_intervals), simulation_plot, statistic = "log zeta")
# d_ply(s5results_pen, .(K_intervals), simulation_plot, statistic = "log zeta", graph = "median_bias")



