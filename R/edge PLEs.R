rm(list=ls())
source("R/PIR-APP.R")

K <- 20
k_priors <- 1.5
theta <- 10

X_list <- list(X0 = rep(0, K + 1),
               X01 = c(rep(0, 1 + K / 2), rep(1, K / 2)),
               X1 = rep(1, K + 1))
U_list <- list(U0 = rep(0, K),
               U01 = c(rep(0, K / 2), rep(1, K / 2)),
               U1 = rep(1, K))

p_func <- Beta_Gamma(k_mu = k_priors, k_lambda = k_priors, theta_mu = theta, theta_lambda = theta, const = 2)

MTS_MLE <- sapply(X_list, MTSmle, c = 1, transform = "exp")
MTS_PLE <- sapply(X_list, MTSmle, c = 1, penalty_func = p_func, transform = "exp")
PIR_MLE <- sapply(U_list, PIRmle, c = 1, d = 0, transform = "exp")
PIR_PLE <- sapply(U_list, PIRmle, c = 1, d = 0, penalty_func = p_func, transform = "exp")
list(MTS_MLE = MTS_MLE, MTS_PLE = MTS_PLE, PIR_MLE = PIR_MLE, PIR_PLE = PIR_PLE)
  
x <- 0.5 * (sqrt(4 * (k_priors - 1)^2 + 4 * K * (2 * k_priors - 1) / (theta * k_priors)) - 2 * (k_priors - 1))
y <- 0.5 * (sqrt(4 * (k_priors - 1)^2 + 4 * K / theta) - 2 * (k_priors - 1))
(k_priors - 1) / (2 * k_priors + x - 1)
(k_priors - 1) / (2 * k_priors + y - 1)
x * (k_priors + x) / (2 * k_priors + x - 1) / K
x / K
y * (k_priors + y) / (2 * k_priors + y - 1) / K
y / K

