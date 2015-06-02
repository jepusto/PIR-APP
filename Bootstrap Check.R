library(plyr)
library(ARPobservation)
library(compiler)
library(Pusto)
source('r/PIR-APP.R')
enableJIT(3)

params$phi <- c(.10, .45, .90)
params$zeta <- c(.5, .35, .2)

params <- data.frame(phi = c(.10, .45, .90), zeta = c(.5, .35, .2))

covered <- function(U, zeta, k_priors, theta, bootstraps){
out <- PIRbootstrap(U = U, c = 1, d = 0, transform = "exp",
                    penalty_func = Beta_Gamma(k_mu = k_priors, 
                                              k_lambda = k_priors, 
                                              theta_mu = theta, 
                                              theta_lambda = theta,), iterations = bootstraps)

(out$CI_L[2] <= zeta & zeta <= out$CI_U[2])
}

bootstrap_check <- function(phi, zeta, c, d, iterations, bootstraps, k_priors, theta){
  mu <- phi/zeta
  lambda <- (1-phi)/zeta
  
  U <- r_PIR(n = iterations, mu = mu, lambda = lambda, stream_length = 100, F_event = F_exp(),
             F_interim = F_exp(), interval_length = 1)
  
  coverage <- apply(U, 1, covered, zeta = zeta, k_priors = 1.5, theta = 10, bootstraps = bootstraps)
  
  mean(coverage)
}

source_func <- ls()
system.time(results <- mdply(params, bootstrap_check, c = 1, d = 0, k_priors = 1.5, theta = 10, iterations = 10, bootstraps = 50))
