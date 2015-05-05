library(Pusto)
library(plyr)
library(ARPobservation)
rm(list=ls())
source("R/PIR-APP.R")
source_obj <- ls()
cluster <- start_parallel(source_obj)
clusterEvalQ(cluster, library(ARPobservation))

#------------------------------------------------------
# Momentary time sampling 
#------------------------------------------------------

phi <- 0.4
zeta <- 0.05
stream_length <- 500
c <- 1

MTS_info(phi, zeta, c_M)

MTS_dat <- r_MTS(n=10000, mu = phi / zeta, lambda = (1 - phi) / zeta, stream_length = stream_length, 
                 F_event = F_exp(), F_interim = F_exp(), interval_length = c, summarize = FALSE)

ests <- adply(MTS_dat, .margins = 1, .fun = MTSmle, c = c, transform = "exp", .parallel=TRUE)[,-1]
ests$rho <- with(ests, zeta / (phi * (1 - phi)))
ests$zeta <- NULL
(sim_cov <- cov(ests) * stream_length / c)
(info_cov <- chol2inv(chol(MTS_info(phi, zeta, c))))
max(abs(sim_cov - info_cov))

#------------------------------------------------------
# Augmented interval recording
#------------------------------------------------------

phi <- 0.4
zeta <- 0.6
zeta / (phi * (1 - phi))
stream_length <- 400
c <- 1.5
d <- 0.5

AIR_info(phi, zeta, c, d)

AIR_dat <- r_AIR(n=10000, mu = phi / zeta, lambda = (1 - phi) / zeta, stream_length = stream_length, 
                 F_event = F_exp(), F_interim = F_exp(), interval_length = c + d, rest_length = d)

ests <- adply(AIR_dat, .margins = 3, .fun = AIRmle, c = c, d = d, transform = "exp", .parallel = TRUE)[,-1]
ests$rho <- with(ests, zeta / (phi * (1 - phi)))
ests$zeta <- NULL
(sim_cov <- cov(ests) * stream_length / (c + d))
(info_cov <- chol2inv(chol(AIR_info(phi, zeta, c, d))))
max(abs(sim_cov - info_cov))


stopCluster(cluster)
