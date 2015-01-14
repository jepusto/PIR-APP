library(ARPobservation)
library(plyr)

rm(list=ls())
source("R/parallel setup.R")
source("R/PIR-APP.R")

iterations <- 10
phi <- 0.1
zeta <- 0.02
K <- 60
threshold <- 8

mu <- phi / zeta
lambda <- (1 - phi) / zeta
BS <- r_behavior_stream(n = 1, mu = mu, lambda = lambda, 
                        F_event = F_exp(), F_interim = F_exp(), 
                        stream_length = K)
U <- interval_recording(BS = BS, interval_length = 1, summarize = FALSE)[,1]
PIRmle(U, c = 1, d = 0, coding = "PIR", transform = "exp")


PIR_proportion_valid <- function(iterations, phi, zeta, K, threshold) {
  mu <- phi / zeta
  lambda <- (1 - phi) / zeta
  
  f <- function(mu, lambda, K) {
    BS <- r_behavior_stream(n = 1, mu = mu, lambda = lambda, 
                            F_event = F_exp(), F_interim = F_exp(), 
                            stream_length = K)
    U <- interval_recording(BS = BS, interval_length = 1, summarize = FALSE)[,1]
    PIRmle(U, c = 1, d = 0, coding = "PIR", transform = "exp")
  }
  
  PIR_mle <- data.frame(t(replicate(iterations, f(mu, lambda, K))))
  
  phi_finite <- abs(logit(PIR_mle$phi)) < threshold
  zeta_finite <- abs(log(PIR_mle$zeta)) < threshold
  mean(phi_finite * zeta_finite)
}

source_func <- ls()


phi <- seq(0.1, 0.9, 0.1)
zeta <- c(0.02, 0.05, 0.10, 0.20, 0.25, 0.40, 0.50)
K <- 40
iterations <- 100
threshold <- 8
params <- expand.grid(phi = phi, zeta = zeta)

cluster <- start_parallel(source_func)
clusterEvalQ(cluster, library(ARPobservation))
system.time(PIR_results <- maply(params, PIR_proportion_valid, K = K, 
                     iterations = iterations, threshold = threshold,
                     .parallel = TRUE))
stopCluster(cluster)
