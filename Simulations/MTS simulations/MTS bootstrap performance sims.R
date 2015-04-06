rm(list = ls())
library(ARPobservation)
library(plyr)
library(fields)
source("R/PIR-APP.R")

load("Simulations/MTS simulations/MTS sample size sims - BS.Rdata")
phiquant <- subset(results, stat == "phi")
zetaquant <- subset(results, stat == "zeta")
names(phiquant)[8:14] <- paste0("p", names(phiquant)[8:14])
names(zetaquant)[8:14] <- paste0("z", names(zetaquant)[8:14])
phiquant <- phiquant[, -7]
zetaquant <- zetaquant[, -7]
results <- merge(x = phiquant, y = zetaquant)
rm(list = c("phiquant", "zetaquant"))

#-------------------------------------
#Data Generating Model
#-------------------------------------
#See ARPobservation package's r_MTS function

#-------------------------------------
#Performance statistics
#-------------------------------------

BootstrapCIs <- function(ests, pctls) {
  
  pctls <- arrange(pctls, zeta, phi)
  phi_range <- unique(pctls$phi)
  zeta_range <- unique(pctls$zeta)
  
  # arrange percentiles into matrices
  zquant025 <- matrix(pctls$zquant025, length(phi_range), length(zeta_range))
  zquant975 <- matrix(pctls$zquant975, length(phi_range), length(zeta_range))
  pquant025 <- matrix(pctls$pquant025, length(phi_range), length(zeta_range))
  pquant975 <- matrix(pctls$pquant975, length(phi_range), length(zeta_range))
  
  # bilinear interpolation
  z025 <- interp.surface(list(x=phi_range, y = zeta_range, z = zquant025), ests)
  z975 <- interp.surface(list(x=phi_range, y = zeta_range, z = zquant975), ests)
  p025 <- interp.surface(list(x=phi_range, y = zeta_range, z = pquant025), ests)
  p975 <- interp.surface(list(x=phi_range, y = zeta_range, z = pquant975), ests)
  
  data.frame(p025, p975, z025, z975)
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
  
  bootCIs <- BootstrapCIs(ests = ests, pctls=results[results$K_intervals==K_intervals & results$k_priors==k_priors & results$theta == theta,])
  
  
  ispcovered <- with(bootCIs, (p025 < phi & phi < p975))
  iszcovered <- with(bootCIs, (z025 < zeta & zeta < z975))
  
  
  pNA <- mean(is.na(ispcovered))
  zNA <- mean(is.na(iszcovered))
  
  
  
  pOOR <- with(bootCIs, mean(p025 >= p975, na.rm=TRUE))
  zOOR <- with(bootCIs, mean(z025 >= z975, na.rm=TRUE))
  
  
  
  ispcovered[is.na(ispcovered)] <- FALSE
  iszcovered[is.na(iszcovered)] <- FALSE
  
  
  pcoverage <- mean(ispcovered)
  zcoverage <- mean(iszcovered)
  
  c(pcoverage = pcoverage, pNA = pNA, pOOR = pOOR, zcoverage = zcoverage, zNA = zNA, zOOR = zOOR)
}

#-------------------------------------
# Design parameters
#-------------------------------------

K_intervals <- seq(10, 150, 10)
phi <- seq(.05, .50, .05)
zeta <- seq(.05, .50, .05)
k_priors <- c(1, 1.5)
theta <- c(10, Inf)
set.seed(20150316)

params <- expand.grid(phi = phi, zeta = zeta, K_intervals = K_intervals, k_priors = k_priors, theta = theta)
params <- subset(params, (k_priors == 1 & theta == Inf) | (k_priors == 1.5 & theta == 10))
params$seed <- round(runif(nrow(params)) * 2^30)
nrow(params)

# # single core run
# system.time(results_single <- mdply(params, .fun = runSim, iterations = 10, c = 1, .inform = TRUE))


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
    clusterEvalQ(cluster, library(fields))
    clusterEvalQ(cluster, library(compiler))
    clusterEvalQ(cluster, enableJIT(3))
    clusterExport(cluster, source_func) 
  } else {
    registerDoParallel(cores=detectCores())
  }
  cluster
}

cluster <- start_parallel(source_func)
system.time(BSresults <- mdply(params, .fun = runSim, iterations = 50, c = 1, rest_length = 0, .parallel = TRUE))
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
