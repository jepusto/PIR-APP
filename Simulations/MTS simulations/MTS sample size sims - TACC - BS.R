rm(list = ls())
#-------------------------------------
#Data Generating Model
#-------------------------------------
#See ARPobservation package's r_MTS function

#-------------------------------------
#Performance statistics
#-------------------------------------

logit <- function(x) log(x) - log(1 - x)

logit_trunc <- function(x, trunc) {
  y <- pmin(pmax(x, trunc), 1 - trunc)
  log(y) - log(1 - y)
}

expit <- function(x) 1 / (1 + exp(-x))

performance <- function(true_phi, true_zeta, true_mu, true_lambda, ests, K) {
  
  all_trans <- cbind(ests[,1], ests[,2], ests[,3], 
                     logit(ests[,1]), log(ests[,2]), logit_trunc(ests[,3], trunc = 0.5 / K),
                     ests[,1]/ests[,2], (1-ests[,1])/ests[,2], 
                     log(ests[,1]/ests[,2]), log((1-ests[,1])/ests[2]))
  
  true_parms <- c(true_phi, true_zeta, true_phi,
                  logit(true_phi), log(true_zeta), logit(true_phi),
                  true_mu, true_lambda, log(true_mu), log(true_lambda))
  
  results <- data.frame(stat = c("phi", "zeta","mean", 
                                 "logit phi", "log zeta", "logit mean",
                                 "mu", "lambda", "log mu", "log lambda"))
  results$bias <- colMeans(all_trans, na.rm = TRUE) - true_parms
  results$median_bias <- apply(all_trans, 2, median) - true_parms
  results$var <- apply(all_trans, 2, var)
  results$missing <- apply(all_trans, 2, function(x) sum(is.na(x)))
  results$quant025 <- NA
  results$quant050 <- NA
  results$quant975 <- NA
  results[1,6:8] <- quantile(all_trans[,1], probs = c(0.025, 0.5, 0.975))
  results[2,6:8] <- quantile(all_trans[,2], probs = c(0.025, 0.5, 0.975))
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
                                          theta_mu = theta, theta_lambda = theta, const = 2), 
                                          transform = "exp")[,-1]
  MTSmean <- rowMeans(X[,-1])
    
  performance(true_phi = phi, true_zeta = zeta, 
              true_mu = mu, true_lambda = lambda, 
              ests = cbind(ests, mean = MTSmean), K = K_intervals)
}

#-------------------------------------
# Design parameters
#-------------------------------------
library(plyr)

K_intervals <- seq(10, 150, 10)
phi <- c(0, seq(0.01, 0.99, length.out = 50))
zeta <- c(0.00001, seq(0.005, 0.9, length.out = 50), 1:3)
k_priors <- c(1, 1.5)
theta <- c(10, Inf)
set.seed(20150316)

params <- expand.grid(phi = phi, zeta = zeta, K_intervals = K_intervals, k_priors = k_priors, theta = theta)
params <- subset(params, (k_priors == 1 & theta == Inf) | (k_priors == 1.5 & theta == 10))
params$seed <- round(runif(nrow(params)) * 2^30)
nrow(params)

#--------------------------------------------------------
# run simulations in parallel
#--------------------------------------------------------

source_func <- ls()

library(Rmpi)
library(snow)
library(foreach)
library(iterators)
library(doSNOW)
library(rlecuyer)
library(plyr)

# set up parallel processing
cluster <- getMPIcluster()
registerDoSNOW(cluster)

clusterEvalQ(cluster, source("PIR-APP.R"))
clusterEvalQ(cluster, library(plyr))
clusterEvalQ(cluster, library(ARPobservation))
clusterEvalQ(cluster, library(compiler))
clusterEvalQ(cluster, enableJIT(3))
clusterExport(cluster, source_func) 

system.time(results <- mdply(params, .fun = runSim, iterations = 10000, c = 1, .parallel = TRUE))
stopCluster(cluster)

save(results, file = "MTS sample size sims - BS.Rdata")

