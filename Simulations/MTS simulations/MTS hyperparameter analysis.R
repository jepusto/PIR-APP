rm(list=ls())
library(ggplot2)
library(plyr)
library(reshape2)

load("Simulations/MTS simulations/MTS hyperparameter sims.Rdata")

results <- within(results, {
  abs_bias <- abs(bias)
  rel_bias <- bias / ifelse(stat=="phi", phi, ifelse(stat=="zeta", zeta, NA))
  med_rel_bias <- median_bias / ifelse(stat=="phi", phi, ifelse(stat=="zeta", zeta, NA))
  rmse <- sqrt(bias^2 + var)
  alpha <- factor(k_priors)
})

names(results)[3] <- "K"

hyper_plot <- function(parm, measure, limits = NULL) {
  dat <- subset(results, stat == parm & (theta != Inf | alpha == 1) & (theta == Inf | alpha != 1)) 
  dat$y <- dat[,measure]
  ggplot(dat, aes(x = alpha, y = y)) + 
    geom_boxplot() + 
    coord_cartesian(ylim = limits) + 
    labs(y = measure) +
    facet_grid(K ~ theta, scales = "free", labeller = "label_both") + 
    theme_minimal()  
}

hyper_plot("phi", "bias")
hyper_plot("phi", "abs_bias")
hyper_plot("phi", "rmse")

hyper_plot("zeta", "bias", limits = c(-1,1))
hyper_plot("zeta", "rel_bias", limits = c(-1,1))
hyper_plot("zeta", "med_rel_bias", limits = c(-1,1))
hyper_plot("zeta", "abs_bias", limits = c(0,2))
hyper_plot("zeta", "rmse", limits = c(0,1))

hyper_plot("logit phi", "bias", limits = c(-1, 1))
hyper_plot("logit_phi", "abs_bias", limits = c(0, 1))
hyper_plot("logit_phi", "rmse", limits = c(0, 3))

hyper_plot("log zeta","bias", limits = c(-1, 1))
hyper_plot("log_zeta","abs_bias", limits = c(0, 1))
hyper_plot("log_zeta", "rmse", limits = c(0, 3))
