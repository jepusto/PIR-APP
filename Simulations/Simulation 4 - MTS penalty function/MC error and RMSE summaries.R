library(plyr)
library(reshape2)
setwd("Simulations/Simulation 4 - MTS penalty function")
rm(list = ls())
load("sim4MTS_unp2.Rdata")
load("sim4MTSpenalized2.Rdata")

results_unp2$theta <- 0
results_unp2$k_priors <- 0
results_all <- rbind(results2, results_unp2)

mc_error <- function(data, statistic){

  data_sub <- subset(data, stat == statistic)
  
  mcerror <- ddply(data_sub, .(K_intervals, k_priors, theta, phi, zeta), summarize, 
                   error = var/5000)

  final_errors <- ddply(mcerror, .(K_intervals), summarize, max_error = max(error))
  final_errors
}

lphi <- cbind(stat = "logit phi", mc_error(data = results_all, statistic = "logit phi"))
phi <- cbind(stat = "phi", mc_error(data = results_all, statistic = "phi"))
lzeta <- cbind(stat = "log zeta", mc_error(data = results_all, statistic = "log zeta"))
zeta <- cbind(stat = "zeta", mc_error(data = results_all, statistic = "zeta"))

dcast(rbind(lphi, phi, lzeta, zeta), stat~K_intervals, value.var = "max_error")


# median/max absolute bias
# median/max RMSE across phi/zeta, for a given level of K, by priors

sim_summary <- function(data, statistic){
  
  data_sub <- subset(data, stat == statistic)
  
  rmse <- ddply(data_sub, .(K_intervals, k_priors, theta, phi, zeta), summarize,
        RMSE = sqrt(var + bias^2))
  
  rmse_summary <- ddply(rmse, .(K_intervals, k_priors, theta), summarize,
                        max_rmse = max(RMSE), median_rmse = median(RMSE))
  bias_summary <- ddply(data_sub, .(K_intervals, k_priors, theta), summarize,
                       max_bias = max(bias), median_bias = median(bias))
  
  summarized <- merge(bias_summary, rmse_summary, id.vars = c("K_intervals", "k_priors", "theta"))
  summarized
  
}

phi_sum <- sim_summary(data = results_all, stat = "phi")
zeta_sum <- sim_summary(data = results_all, stat = "zeta")

#phi max bias, median bias, max RMSE and median RMSE in order
dlply(phi_sum, .(K_intervals), dcast, formula = k_priors ~ theta, value.var = "max_bias")[1:2]
dlply(phi_sum, .(K_intervals), dcast, formula = k_priors ~ theta, value.var = "median_bias")[1:2]
dlply(phi_sum, .(K_intervals), dcast, formula = k_priors ~ theta, value.var = "max_rmse")[1:2]
dlply(phi_sum, .(K_intervals), dcast, formula = k_priors ~ theta, value.var = "median_rmse")[1:2]

#zeta max bias, median bias, max RMSE and median RMSE in order
dlply(zeta_sum, .(K_intervals), dcast, formula = k_priors ~ theta, value.var = "max_bias")[1:2]
dlply(zeta_sum, .(K_intervals), dcast, formula = k_priors ~ theta, value.var = "median_bias")[1:2]
dlply(zeta_sum, .(K_intervals), dcast, formula = k_priors ~ theta, value.var = "max_rmse")[1:2]
dlply(zeta_sum, .(K_intervals), dcast, formula = k_priors ~ theta, value.var = "median_rmse")[1:2]
