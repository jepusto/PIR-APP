rm(list = ls())
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
  
  performance(true_phi = phi, true_zeta = zeta, true_mu = mu, true_lambda = lambda, ests = ests)
}

#-------------------------------------
# Design parameters
#-------------------------------------
library(plyr)

K_intervals <- seq(10, 150, 10)
phi <- seq(.05, .50, .05)
zeta <- seq(.05, .50, .05)
k_priors <- c(1, 1.5)
theta <- c(10, Inf)
set.seed(20150316)

# check that penalty function zeros out when k_priors = 1, theta = Inf
params <- expand.grid(phi = phi, zeta = zeta, K_intervals = K_intervals, k_priors = k_priors, theta = theta)
params <- subset(params, (k_priors == 1 & theta == Inf) | (k_priors == 1.5 & theta == 10))
params$seed <- round(runif(nrow(params)) * 2^30)
nrow(params)

# # single core run
# system.time(results_single <- mdply(params, .fun = runSim, iterations = 5000, c = 1, .inform = TRUE))


# multi-core run

source_func <- ls()

start_parallel <- function(source_func) {
  require(parallel)
  require(foreach)
  require(iterators)
  require(doParallel)
  if (!is.na(pmatch("Windows", Sys.getenv("OS")))) {
    cluster <- makeCluster(detectCores() - 1, type = "SOCK")
    registerDoParallel(cluster)
    clusterEvalQ(cluster, library(plyr))
    clusterEvalQ(cluster, library(ARPobservation))
    clusterEvalQ(cluster, library(compiler))
    clusterEvalQ(cluster, enableJIT(3))
    clusterExport(cluster, source_func) 
  } else {
    registerDoParallel(cores=detectCores())
  }
  cluster
}

cluster <- start_parallel(source_func)
system.time(results <- mdply(params, .fun = runSim, iterations = 5000, c = 1, .parallel = TRUE))
stopCluster(cluster)

save(results, file = "Simulations/MTS simulations/MTS_sample_size_sim.Rdata")

# library(mailR)
# send.mail(from = "jepusto@gmail.com",
#           to = "pusto@austin.utexas.edu",
#           subject = "Simulation is done",
#           body = "All set!",
#           smtp = list(host.name = "smtp.gmail.com", 
#                       port = 465, 
#                       user.name = "jepusto", 
#                       passwd = "xiqzdlacycwuksdf", 
#                       ssl = TRUE),
#           authenticate = TRUE)
