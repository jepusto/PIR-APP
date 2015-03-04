library(plyr)
library(reshape2)
setwd("Simulations/Simulation 4 - MTS penalty function")
rm(list = ls())
load("sim4MTS_unp2.Rdata")
load("sim4MTSpenalized2.Rdata")

results_unp2$theta <- 0
results_unp2$k_priors <- 0
results_all <- rbind(results2, results_unp2)

stats <- c("phi", "logit phi", "zeta", "log zeta")

results_sub <- subset(results_all, stat %in% stats)

mcerror <- ddply(results_sub, .(K_intervals, k_priors, theta, phi, zeta, stat), 
                 summarize, error = var/5000)

final_errors <- ddply(mcerror, .(K_intervals, stat, theta), summarize, max_error = max(error))

final_errors <- dcast(final_errors, stat +theta ~K_intervals, value.var = "max_error")


# median/max absolute bias
# median/max RMSE across phi/zeta, for a given level of K, by priors

rmse <- ddply(results_sub, .(K_intervals, k_priors, theta, phi, zeta, stat), summarize,
              RMSE = sqrt(var + bias^2))

rmse_summary <- ddply(rmse, .(K_intervals, k_priors, theta, stat), summarize,
                      max_rmse = max(RMSE), median_rmse = median(RMSE))
bias_summary <- ddply(results_sub, .(K_intervals, k_priors, theta, stat), summarize,
                      max_bias = max(abs(bias)), median_bias = median(abs(bias)))

summarized <- merge(bias_summary, rmse_summary, id.vars = c("K_intervals", "k_priors", "theta", "stat"))

max_bias <- dlply(summarized, .(K_intervals, stat), dcast, formula = k_priors + stat ~ theta, value.var = "max_bias")[1:8]
median_bias <- dlply(summarized, .(K_intervals, stat), dcast, formula = k_priors ~ theta, value.var = "median_bias")[1:8]
max_rmse <- dlply(summarized, .(K_intervals, stat), dcast, formula = k_priors ~ theta, value.var = "max_rmse")[1:8]
median_rmse <- dlply(summarized, .(K_intervals, stat), dcast, formula = k_priors ~ theta, value.var = "median_rmse")[1:8]
