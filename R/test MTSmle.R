library(ARPobservation)
library(plyr)
library(ggplot2)

rm(list=ls())
source("R/Brown Solomon Stevens.R")
source("R/PIR-APP.R")

iterations <- 5000
phi <- 0.4
zeta <- 0.4
K <- 100

BS <- r_behavior_stream(n = iterations, mu = phi / zeta, lambda = (1 - phi) / zeta, 
                          F_event = F_exp(), F_interim = F_exp(), 
                          stream_length = K)
X <- momentary_time_recording(BS = BS, interval_length = 1, summarize = FALSE)
BSS <- adply(X, 2, MTSmle, c = 1, transform = "exp")[,-1]
colMeans(BSS, na.rm = TRUE)

MLEs <- adply(X, 2, MTSmle, c = 1, 
               penalty_func = function(param, coding) 0, 
               transform = "exp")[,-1]
colMeans(MLEs)
apply(MLEs, 2, median)

PMLEs <- adply(X, 2, MTSmle, c = 1, 
                   penalty_func = Beta_Gamma(k_mu = 1.5, k_lambda = 1.5, theta_mu = 3, theta_lambda = 3, const = 2), 
                   transform = "exp")[,-1]
colMeans(PMLEs)
apply(PMLEs, 2, median)


ests_wide <- cbind(MLE = MLEs, PMLE = PMLEs)

ggplot(ests_wide, aes(MLE.phi, PMLE.phi)) + 
  geom_point(alpha = 0.1) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")

ggplot(ests_wide, aes(MLE.zeta, PMLE.zeta)) + 
  geom_point(alpha = 0.1) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")
last_plot() + coord_cartesian(x = c(0,1), y = c(0, 1))

MLEs$est <- "MLE"
PMLEs$est <- "PMLE"
ests_long <- rbind(MLEs, PMLEs)

ddply(ests_long, .(est), summarize, phi = mean(phi), zeta = mean(zeta))
medians <- ddply(ests_long, .(est), summarize, phi = median(phi), zeta = median(zeta))
medians

ggplot(ests_long, aes(phi, zeta, color = est)) + 
  geom_point(alpha = 0.1) + 
  geom_point(data = medians, size = 4) + 
  geom_point(x = phi, y = zeta, color = "black", size = 3) + 
  coord_cartesian(y = c(0, 1))
