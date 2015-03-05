rm(list = ls())
library(plyr)
library(ARPobservation)
source("PIR-APP.R")

#-------------------------------------
#Data Generating Model
#-------------------------------------
#See ARPobservation package's r_PIR function

#-------------------------------------
#Performance statistics
#-------------------------------------

logit <- function(x) log(x) - log(1 - x)

expit <- function(x) 1 / (1 + exp(-x))

performance <- function(true_phi, true_zeta, true_mu, true_lambda, ests){
  
  all_trans <- cbind(ests[,1], ests[,2], logit(ests[,1]), log(ests[,2]), 
                     ests[,1]/ests[,2], (1-ests[,1])/ests[,2], 
                     log(ests[,1]/ests[,2]), log((1-ests[,1])/ests[2]))
  
  true_parms <- c(true_phi, true_zeta, logit(true_phi), log(true_zeta), 
                  true_mu, true_lambda, log(true_mu), log(true_lambda))
  
  results <- data.frame(stat = c("phi", "zeta", "logit phi", "log zeta", "mu", "lambda", "log mu", "log lambda"))
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
# K_intervals <- 90
# c <- 1
# k_priors <- 1.01
# theta <- 20
# iterations <- 4000

runSim <- function(phi, zeta, K_intervals, c, rest_length, k_priors, theta, iterations, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  mu <- phi / zeta
  lambda <- (1 - phi) / zeta
  
  X <- r_PIR(n = iterations, mu = mu, lambda = lambda, 
             stream_length = K_intervals * c, F_event = F_exp(), 
             F_interim = F_exp(), interval_length = c, rest_length = 0)
  
  ests <- adply(X, .margins = 1, .fun = PIRmle, c = 1, d = rest_length, 
                penalty_func = Beta_Gamma(k_mu = k_priors, k_lambda = k_priors, 
                                          theta_mu = theta, theta_lambda = theta), 
                transform = "exp")[,-1]
  
  performance(true_phi = phi, true_zeta = zeta, true_mu = mu, true_lambda = lambda, ests = ests)
}

#-------------------------------------
# Design parameters
#-------------------------------------

K_intervals <- c(30, 60, 90, 120)
phi <- seq(.05, .95, .05)
zeta <- seq(.05, .50, .05)
k_priors <- c(1.001, 1.01, 1.05, 1.10)
theta <- c(5, 10, 20, 40)
set.seed(8473695)

print(dims <- c(length(K_intervals),length(phi),length(zeta),length(k_priors), length(theta)))
prod(dims)
params <- expand.grid(phi = phi, zeta = zeta, K_intervals = K_intervals, k_priors = k_priors, theta = theta)[sample(prod(dims)),]
params$seed <- round(runif(nrow(params)) * 2^30)


params <- expand.grid(phi = phi, zeta = zeta, K_intervals = K_intervals, k_priors = k_priors, theta = theta)
params$seed <- round(runif(nrow(params)) * 2^30)

#----------------------------------------------
# Run simulation in parallel on the TACC
#----------------------------------------------
source_func <- ls()

library(Rmpi)
library(snow)
library(foreach)
library(iterators)
library(doSNOW)
library(plyr)

# set up parallel processing
cluster <- getMPIcluster()
registerDoSNOW(cluster)

# export source functions
clusterExport(cluster, source_func)
clusterEvalQ(cluster, library(plyr))
clusterEvalQ(cluster, library(ARPobservation))
clusterEvalQ(cluster, library(compiler))
clusterEvalQ(cluster, enableJIT(3))

# execute simulations
system.time(s5results_pen <- mdply(params, .fun = runSim, c = 1, rest_length = 0, iterations = 5000, .parallel=TRUE))
stopCluster(cluster)

save(s5results_pen, file = "s5results_pen.Rdata")
