library(ggplot2)
library(plyr)
library(reshape2)

load("Simulations/MTS simulations/MTS hyperparameter sims.Rdata")

results <- within(results, {
  abs_bias <- abs(bias)
  rmse <- sqrt(bias^2 + var)
  alpha <- factor(k_priors)
})

names(results)[3] <- "K"

hyper_plot <- function(stat, measure, limits = NULL) {
  dat <- subset(results, stat == stat & (theta != Inf | alpha == 1) & (alpha != 1 | theta == Inf)) 
  dat$y <- dat[,measure]
  ggplot(dat, aes(x = alpha, y = y)) + 
    geom_boxplot() + 
    coord_cartesian(ylim = limits) + 
    labs(y = measure) +
    facet_grid(K ~ theta, scales = "free", labeller = "label_both") + 
    theme_minimal()  
}

hyper_plot("logit_phi", "bias", limits = c(-1, 1))
hyper_plot("log_zeta","bias", limits = c(-1, 1))
hyper_plot("logit_phi", "abs_bias", limits = c(0, 1))
hyper_plot("log_zeta","abs_bias", limits = c(0, 1))
hyper_plot("logit_phi", "rmse", limits = c(0, 3))
hyper_plot("log_zeta", "rmse", limits = c(0, 3))
