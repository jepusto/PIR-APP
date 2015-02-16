rm(list = ls())
library(plyr)
library(ARPobservation)
source("R/PIR-APP.R")

#-------------------------------------
#Data Generating Model
#-------------------------------------
#See ARPobservation package's r_MTS function

#-------------------------------------
#Performance statistics
#-------------------------------------

logit <- function(x) log(x) - log(1 - x)

expit <- function(x) 1 / (1 + exp(-x))

performance <- function(true_phi, true_zeta, ests){

  phi <- expit(ests[,1])
  zeta <- exp(ests[,2])
  
  means <- colMeans(ests, na.rm = TRUE)
  
  mu <- phi/zeta
  lambda <- (1-phi)/zeta
  true_mu <- true_phi/true_zeta
  true_lambda <- (1-true_phi)/true_zeta
  
  logmu <- log(mu)
  loglambda <- log(lambda)
  true_logmu <- log(true_mu)
  true_loglambda <- log(true_lambda)
  
  results <- matrix(NA, nrow = 8, ncol = 5, 
                    dimnames = list(c("logit_phi", "log_zeta", "phi", "zeta", "log_mu", "log_lambda", "mu", "lambda"), 
                                    c("bias", "median_bias", "var",  "phi_nas", "zeta_nas")))
  
  results[1,1] <- means[1] - logit(true_phi)
  results[1,2] <- median(ests[,1]) - logit(true_phi)
  results[1,3] <- var(ests[,1], na.rm = TRUE)
  results[,4]  <- mean(is.na(ests[,1]))
  results[,5]  <- mean(is.na(ests[,2]))
  
  results[2,1] <- means[2] - log(true_zeta)
  results[2,2] <- median(ests[,2], na.rm = TRUE) - logit(true_phi)
  results[2,3] <- var(ests[,2], na.rm = TRUE)

  
  results[3,1] <- mean(phi, na.rm = TRUE) - true_phi
  results[3,2] <- median(phi, na.rm = TRUE) - true_phi
  results[3,3] <- var(phi, na.rm = TRUE)

  results[4,1] <- mean(zeta, na.rm = TRUE) - true_zeta
  results[4,2] <- median(zeta, na.rm = TRUE) - true_zeta
  results[4,3] <- var(zeta, na.rm = TRUE)
 
  results[5,1] <- mean(logmu, na.rm = TRUE) - true_logmu
  results[5,2] <- median(logmu, na.rm = TRUE) - true_logmu
  results[5,3] <- var(logmu, na.rm = TRUE)

  
  results[6,1] <- mean(loglambda, na.rm = TRUE) - true_loglambda
  results[6,2] <- median(loglambda, na.rm = TRUE) - true_loglambda
  results[6,3] <- var(loglambda, na.rm = TRUE)
  
  results[7,1] <- mean(mu, na.rm = TRUE) - true_mu
  results[7,2] <- median(mu, na.rm = TRUE) - true_mu
  results[7,3] <- var(mu, na.rm = TRUE)
  
  results[8,1] <- mean(lambda, na.rm = TRUE) - true_lambda
  results[8,2] <- median(lambda, na.rm = TRUE) - true_lambda
  results[8,3] <- var(lambda, na.rm = TRUE)
  
  results
}


#-------------------------------------
#Simulation Driver
#-------------------------------------
# phi <- .20
# zeta <- .20
# K_intervals <- 90
# c <- 1
# k_priors <- 1.01
# theta <- 20
# iterations <- 4000

runSim <- function(phi, zeta, K_intervals, c, k_priors, theta, iterations, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  mu <- phi / zeta
  lambda <- (1 - phi) / zeta
  
  X <- r_MTS(n = iterations, mu = mu, lambda = lambda, stream_length = K_intervals * c, 
             F_event = F_exp(), F_interim = F_exp(), interval_length = c)
  
  ests <- adply(X, .margins = 1, .fun = MTSmle, c = 1, 
                penalty_func = Beta_Gamma(k_mu = k_priors, k_lambda = k_priors, 
                                          theta_mu = theta, theta_lambda = theta))[,-1]
  
  performance(true_phi = phi, true_zeta = zeta, ests = ests)
}

#-------------------------------------
# Design parameters
#-------------------------------------

K_intervals <- c(30, 90)
phi <- c(.05, .10, .20, .30, .40, .50)
zeta <- c(.05, .10, .20, .30, .40, .50)
k_priors <- c(1.01, 1.05, 1.1)
theta <- c(10, 20)
set.seed(8473695)


params <- expand.grid(phi = phi, zeta = zeta, K_intervals = K_intervals, k_priors = k_priors, theta = theta)
params$seed <- round(runif(nrow(params)) * 2^30)

system.time(results <- mdply(params, .fun = runSim, iterations = 100, c = 1))
results$stat <- c("logit_phi", "log_zeta", "phi", "zeta", "log_mu", "log_lambda", "mu", "lambda")