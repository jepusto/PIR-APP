rm(list = ls())
library(ggplot2)
load("Simulations/PIR Simulations/PIR sample size sims BS.Rdata")
load("Simulations/PIR Simulations/PIR bootstrap results.Rdata")

qplot(factor(theta), bias, geom = "boxplot", data = subset(results, stat == "log zeta")) +
    facet_wrap(~K_intervals)

qplot(factor(theta), var, geom = "boxplot", data = subset(results, stat == "zeta")) +
  facet_wrap(~K_intervals)


qplot(factor(theta), bias, geom = "boxplot", data = subset(results, stat == "phi")) +
  facet_wrap(~K_intervals)
