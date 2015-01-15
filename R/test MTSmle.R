library(ARPobservation)
library(plyr)
library(ggplot2)

rm(list=ls())
source("R/Brown Solomon Stevens.R")
source("R/PIR-APP.R")

iterations <- 10
phi <- 0.1
zeta <- 0.1
K <- 40

BS <- r_behavior_stream(n = iterations, mu = phi / zeta, lambda = (1 - phi) / zeta, 
                          F_event = F_exp(), F_interim = F_exp(), 
                          stream_length = K)
X <- momentary_time_recording(BS = BS, interval_length = 1, summarize = FALSE)

MLEs <- adply(X, 2, MTSmle, c = 1, transform = "exp")[,-1]
MLEs$zeta[is.na(MLEs$zeta)] <- max(MLEs$zeta, na.rm = TRUE)
missing_phi <- is.na(MLEs$phi)
MLEs$phi[missing_phi] <- colMeans(X)[missing_phi]
colMeans(MLEs)

MLEs_optim <- adply(X, 2, MTSmle, c = 1, 
               penalty_func = function(param) 0, 
               transform = "exp")[,-1]
colMeans(MLEs_optim)
cbind(MLEs, MLEs_optim)

PMLEs <- adply(X, 2, MTSmle, c = 1, 
                   penalty_func = Beta_Gamma(k_mu = 1.001, k_lambda = 1.001, theta_mu = 6, theta_lambda = 6, const = 2), 
                   transform = "exp")[,-1]
colMeans(PMLEs)

ests_wide <- cbind(MLEs, PMLEs)

MLEs$est <- "MLE"
PMLEs$est <- "PMLE"
ests_long <- rbind(MLEs, PMLEs)

ddply(ests_long, .(est), summarize, phi = median(phi), zeta = median(zeta))

qplot(phi, zeta, color = est, data = ests_long) + 
  geom_point(x = phi, y = zeta, color = "black", size = 4)
