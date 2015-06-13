rm(list = ls())
library(plyr)
library(ARPobservation)
library(compiler)
library(Pusto)
source('r/PIR-APP.R')
enableJIT(3)

bootstrap_check <- function(phi, zeta, c, d, iterations, bootstraps, k_priors, theta){
  mu <- phi/zeta
  lambda <- (1-phi)/zeta
  
  U <- r_PIR(n = iterations, mu = mu, lambda = lambda, stream_length = 100, F_event = F_exp(),
             F_interim = F_exp(), interval_length = 1)

  CI <- adply(U, 1, PIRbootstrap, c = 1, d = 0, transform = "exp",
                     penalty_func = Beta_Gamma(k_mu = k_priors, 
                                               k_lambda = k_priors, 
                                               theta_mu = theta, 
                                               theta_lambda = theta,), iterations = bootstraps)
}

params <- data.frame(phi = c(.10, .45, .90), zeta = c(.5, .35, .2))

source_func <- ls()
cluster <- start_parallel(source_func, cores = 4)
clusterEvalQ(cluster, library(ARPobservation))
clusterEvalQ(cluster, library(plyr))
clusterEvalQ(cluster, library(compiler))
clusterEvalQ(cluster, enableJIT(3))

system.time(results <- mdply(params, bootstrap_check, c = 1, d = 0, k_priors = 1.5, theta = 10, iterations = 200, bootstraps = 2000, .parallel = TRUE))

stopCluster(cluster)

save(results, "Bootstrap Check.Rdata")