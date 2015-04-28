library(plyr)
library(ARPobservation)
rm(list=ls())
source("R/PIR-APP.R")


BSS_vcov <- function(phi, zeta, c_M) {
  p0 <- p_0(c_M, phi, zeta)
  p1 <- p_1(c_M, phi, zeta)
  G_mat <- matrix(c((1 - p1) / (p0 + 1 - p1)^2, 1 / (c_M * (p1 - p0)), -p0 / (p0 + 1 - p1)^2, -1 / (c_M * (p1 - p0))), 2, 2)
  Sigma_mat <- diag(c(p0 * (1 - p0) / (1 - phi), p1 * (1 - p1) / phi))
  G_mat %*% Sigma_mat %*% t(G_mat)
}

MTS_info <- function(phi, zeta, c_M) {
  rho <- zeta / (phi * (1 - phi))
  p0 <- p_0(c_M, phi, zeta)
  p1 <- p_1(c_M, phi, zeta)
  I_phiphi <- (1 - exp(-rho * c_M)) * ((1 - phi) / (phi * (1 - p0)) + phi / ((1 - phi) * p1))
  I_phirho <- c_M * exp(-rho * c_M) * ((1 - phi) / (1 - p0) - phi / p1)
  I_rhorho <- phi * (1 - phi) * c_M^2 * exp(-2 * rho * c_M) * (1 + exp(-rho * c_M)) / ((1 - exp(-rho * c_M)) * (1 - p0) * p1)
  matrix(c(I_phiphi, I_phirho, I_phirho, I_rhorho), 2, 2)
}


phi <- 0.2
zeta <- 0.1
stream_length <- 10000
c_M <- 2

MTS_dat <- r_MTS(n=5000, mu = phi / zeta, lambda = (1 - phi) / zeta, stream_length = stream_length, 
                 F_event = F_exp(), F_interim = F_exp(), interval_length = c_M, summarize = FALSE)

ests <- adply(MTS_dat, .margins = 1, .fun = MTSmle, c = c_M, transform = "exp")[,-1]
ests$rho <- with(ests, zeta / (phi * (1 - phi)))
ests$zeta <- NULL
cov(ests) * stream_length / c_M
chol2inv(chol(MTS_info(phi, zeta, c_M)))
BSS_vcov(phi, zeta, c_M)
