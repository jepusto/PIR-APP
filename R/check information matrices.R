library(plyr)
library(ARPobservation)
rm(list=ls())
source("R/PIR-APP.R")
source("R/APP-Information.R")

phi <- 0.4
zeta <- 0.05
stream_length <- 10000
c_M <- 1

MTS_info(phi, zeta, c_M)

MTS_dat <- r_MTS(n=10000, mu = phi / zeta, lambda = (1 - phi) / zeta, stream_length = stream_length, 
                 F_event = F_exp(), F_interim = F_exp(), interval_length = c_M, summarize = FALSE)

ests <- adply(MTS_dat, .margins = 1, .fun = MTSmle, c = c_M, transform = "exp")[,-1]
ests$rho <- with(ests, zeta / (phi * (1 - phi)))
ests$zeta <- NULL
(sim_cov <- cov(ests) * stream_length / c_M)
(info_cov <- chol2inv(chol(MTS_info(phi, zeta, c_M))))
max(abs(sim_cov - info_cov))
