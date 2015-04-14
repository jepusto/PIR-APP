rm(list = ls())

load("PIR sample size sims BS.Rdata")
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
# K_intervals <- 30
# c <- 1
# rest_length <- 0
# k_priors <- 1
# theta <- Inf
# iterations <- 5000

runSim <- function(phi, zeta, K_intervals, c, rest_length, k_priors, theta, iterations, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  mu <- phi / zeta
  lambda <- (1 - phi) / zeta
  
  U <- r_PIR(n = iterations, mu = mu, lambda = lambda, 
             stream_length = K_intervals * c, F_event = F_exp(), 
             F_interim = F_exp(), interval_length = c, rest_length = 0)
  
  ests <- adply(U, .margins = 1, .fun = PIRmle, c = 1, d = rest_length, 
                penalty_func = Beta_Gamma(k_mu = k_priors, k_lambda = k_priors, 
                                          theta_mu = theta, theta_lambda = theta), 
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
phi <- seq(.05, .95, .05)
zeta <- seq(.05, .50, .05)
k_priors <- c(1, 1.5)
theta <- c(10, Inf)
set.seed(20150317)

params <- expand.grid(phi = phi, zeta = zeta, K_intervals = K_intervals, k_priors = k_priors, theta = theta)

params <- subset(params, (k_priors == 1 & theta == Inf) | (k_priors == 1.5 & theta == 10))
params <- params[sample(nrow(params)),]
params$seed <- round(runif(nrow(params)) * 2^30)
nrow(params)

#--------------------------------------------------------
# run simulations in parallel
#--------------------------------------------------------

source_func <- ls()
source_func
(obj_size <- sapply(source_func, object.size))
sum(obj_size)

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

clusterExport(cluster, source_func) 
clusterEvalQ(cluster, source("PIR-APP.R"))
clusterEvalQ(cluster, library(plyr))
clusterEvalQ(cluster, library(ARPobservation))
clusterEvalQ(cluster, library(fields))
clusterEvalQ(cluster, library(compiler))
clusterEvalQ(cluster, enableJIT(3))

system.time(BSresults <- mdply(params, .fun = runSim, iterations = 10000, c = 1, rest_length = 0, .parallel = TRUE))
stopCluster(cluster)

save(BSresults, file = "PIR sample size sims.Rdata")



