library(xtable)
library(ARPobservation)
library(plyr)
library(reshape2)

rm(list=ls())
source("R/Brown Solomon Stevens.R")

iterations <- 100
phi <- 0.1
zeta <- 0.1
K <- 40

proportion_valid <- function(iterations, phi, zeta, K) {
  mu <- phi / zeta
  lambda <- (1 - phi) / zeta
  f <- function(mu, lambda, K) {
    BS <- r_behavior_stream(n = 1, mu = mu, lambda = lambda, 
                            F_event = F_exp(), F_interim = F_exp(), 
                            stream_length = K)
    X <- momentary_time_recording(BS = BS, interval_length = 1, summarize = FALSE)[,1]
    crossprod(cbind(1 - X[1:K], X[1:K]), cbind(1 - X[1 + 1:K], X[1 + 1:K]))
    MTSmle(X, c = 1)
  }
  MTS_mle <- data.frame(t(replicate(iterations, f(mu, lambda, K))))
  
  phi_finite <- is.finite(logit(MTS_mle$phi))
  zeta_finite <- is.finite(log(MTS_mle$zeta))
  #c(phi_valid = mean(phi_finite), zeta_valid = mean(zeta_finite), both_valid = mean(phi_finite * zeta_finite))
  mean(phi_finite * zeta_finite)
}


phi <- c(0.1, 0.2, 0.3, 0.4, 0.5)
zeta <- c(0.05, 0.1, 0.2, 0.25, 0.4, 0.5)
iterations <- 50

params <- expand.grid(phi = phi, zeta = zeta)
system.time(results <- maply(params, proportion_valid, K = 60, iterations = iterations))
