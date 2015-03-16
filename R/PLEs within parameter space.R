library(plyr)
library(ARPobservation)
rm(list=ls())
source("R/parallel setup.R")
source("R/PIR-APP.R")

MTS_proportion_valid <- function(iterations, phi, zeta, K, p_func, threshold) {
  
  mu <- phi / zeta
  lambda <- (1 - phi) / zeta
  
  f <- function(mu, lambda, K) {
    BS <- r_behavior_stream(n = 1, mu = mu, lambda = lambda, 
                            F_event = F_exp(), F_interim = F_exp(), 
                            stream_length = K)
    X <- momentary_time_recording(BS = BS, interval_length = 1, summarize = FALSE)[,1]
    MTSmle(X, c = 1, penalty_func = p_func)
  }
  
  MTS_mle <- data.frame(t(replicate(iterations, f(mu, lambda, K))))
  
  phi_finite <- abs(MTS_mle$logit_phi) < threshold
  zeta_finite <- abs(MTS_mle$log_zeta) < threshold
  mean(phi_finite * zeta_finite)
}

PIR_proportion_valid <- function(iterations, phi, zeta, K, p_func, threshold) {
  mu <- phi / zeta
  lambda <- (1 - phi) / zeta
  
  f <- function(mu, lambda, K) {
    BS <- r_behavior_stream(n = 1, mu = mu, lambda = lambda, 
                            F_event = F_exp(), F_interim = F_exp(), 
                            stream_length = K)
    U <- interval_recording(BS = BS, interval_length = 1, summarize = FALSE)[,1]
    PIRmle(U, c = 1, d = 0, coding = "PIR", penalty_func = p_func)
  }
  
  PIR_mle <- data.frame(t(replicate(iterations, f(mu, lambda, K))))
  
  phi_finite <- abs(PIR_mle$logit_phi) < threshold
  zeta_finite <- abs(PIR_mle$log_zeta) < threshold
  mean(phi_finite * zeta_finite)
}

source_func <- ls()

phi <- seq(0.1, 0.9, 0.1)
zeta <- c(0.02, 0.05, 0.10, 0.20, 0.25, 0.40, 0.50)
K <- 40
iterations <- 2000
threshold <- 8
p_func <- Beta_Gamma(k_mu = 1.1, k_lambda = 1.1, theta_mu = 20, theta_lambda = 20, const = 2)

params <- expand.grid(phi = phi, zeta = zeta)

cluster <- start_parallel(source_func)
clusterEvalQ(cluster, library(ARPobservation))
MTS_results <- maply(params, MTS_proportion_valid,
                     K = K, iterations = iterations, p_func = p_func, threshold = threshold,
                     .parallel = TRUE)
PIR_results <- maply(params, PIR_proportion_valid, 
                     K = K, iterations = iterations, p_func = p_func, threshold = threshold,
                     .parallel = TRUE)
stopCluster(cluster)
