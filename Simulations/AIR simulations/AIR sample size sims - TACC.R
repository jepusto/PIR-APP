rm(list = ls())
source("PIR-APP.R")

#-------------------------------------
#Data Generating Model
#-------------------------------------
#See ARPobservation package's r_MTS function

#-------------------------------------
#Performance statistics
#-------------------------------------

logit <- function(x) log(x) - log(1 - x)

expit <- function(x) 1 / (1 + exp(-x))

performance <- function(true_phi, true_zeta, true_mu, true_lambda, ests){
  
  all_trans <- cbind(ests[,1], ests[,2], 
                     logit(ests[,1]), log(ests[,2]), 
                     ests[,1]/ests[,2], (1-ests[,1])/ests[,2], 
                     log(ests[,1]/ests[,2]), log((1-ests[,1])/ests[2]))
  
  true_parms <- c(true_phi, true_zeta, 
                  logit(true_phi), log(true_zeta), 
                  true_mu, true_lambda, log(true_mu), log(true_lambda))
  
  results <- data.frame(stat = c("phi", "zeta", "logit phi", "log zeta", 
                                 "mu", "lambda", "log mu", "log lambda"))
  results$bias <- colMeans(all_trans, na.rm = TRUE) - true_parms
  results$median_bias <- apply(all_trans, 2, median) - true_parms
  results$var <- apply(all_trans, 2, var)
  results$missing <- apply(all_trans, 2, function(x) sum(is.na(x)))
  results
}


#-------------------------------------
#Simulation Driver
#-------------------------------------
# phi <- .20
# zeta <- .20
# K_intervals <- 40
# c <- 2
# d <- 0
# k_priors <- 1.5
# theta <- 10
# iterations <- 5000

runSim <- function(phi, zeta, K_intervals, c, k_priors, theta, iterations, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  mu <- phi / zeta
  lambda <- (1 - phi) / zeta
  
  XUW_array <- r_AIR(n = iterations, mu = mu, lambda = lambda, stream_length = K_intervals * c, 
             F_event = F_exp(), F_interim = F_exp(), interval_length = c)
  
  ests <- adply(XUW_array, .margins = 3, .fun = AIRmle, c = c, d = 0,
                penalty_func = Beta_Gamma(k_mu = k_priors, k_lambda = k_priors, 
                                          theta_mu = theta, theta_lambda = theta, const = 2), 
                                          transform = "exp")[,-1]
    
  performance(true_phi = phi, true_zeta = zeta, 
              true_mu = mu, true_lambda = lambda, 
              ests = ests)
}

#-------------------------------------
# Design parameters
#-------------------------------------
library(plyr)

K_intervals <- seq(5, 75, 5)
phi <- seq(.05, .50, .05)
zeta <- seq(.05, .50, .05)
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

clusterEvalQ(cluster, library(plyr))
clusterEvalQ(cluster, library(ARPobservation))
clusterEvalQ(cluster, library(compiler))
clusterEvalQ(cluster, enableJIT(3))
clusterExport(cluster, source_func) 

system.time(results <- mdply(params, .fun = runSim, iterations = 10000, c = 2, .parallel = TRUE))
stopCluster(cluster)

save(results, file = "AIR sample size sims.Rdata")

