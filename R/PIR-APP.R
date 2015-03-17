#---------------------------
# transformation functions 
#---------------------------

# simple transform functions

logit <- function(x) log(x) - log(1 - x)
expit <- function(x) 1 / (1 + exp(-x))


# transform parameters from logit-phi, log-zeta to different scales

param_trans <- function(param, transform) {
  switch(transform,
         none = c(logit_phi = param[1], 
                  log_zeta = param[2]),
         exp = c(phi = expit(param[1]), 
                 zeta = exp(param[2])),
         means = c(mu = expit(param[1]) / exp(param[2]), 
                   lambda = (1 - expit(param[1])) / exp(param[2])),
         logmeans = c(log_mu = -log(1 + exp(-param[1])) - param[2], 
                      log_lambda = -log(exp(param[1]) + 1) - param[2]))
}


#---------------------------------------
# conditional probability functions
#---------------------------------------

p_0 <- function(x,phi,zeta) phi * (1 - exp(-1 * x * zeta / (phi * (1- phi))))
p_1 <- function(x,phi,zeta) phi + (1 - phi) * exp(-1 * x * zeta / (phi * (1- phi)))


#---------------------------------------------------
# penalty functions 
#---------------------------------------------------

# Note that priors are on (phi, 1 / zeta) when const = 1
# For inverse-gamma prior on zeta, use const = -1
# For priors on (mu, lambda), use const = 2

Beta_Gamma <- function(k_mu, k_lambda, theta_mu, theta_lambda, const = 1) {
  function(param, c, coding) {
    if (coding == "WIR") param[1] <- -param[1]
    -(k_mu - 1) * log(1 + exp(-param[1])) - 
      (k_lambda - 1) * log(exp(param[1]) + 1) - 
        (k_mu + k_lambda - const) * param[2] - 
          (expit(param[1]) / theta_mu + (1 - expit(param[1])) / theta_lambda) * exp(-param[2]) / c
  } 
}


# Normal-normal priors on log(mu), log(lambda)

Norm_mu_lambda <- function(g_mu, g_lambda, sigma_mu, sigma_lambda)
  function(param, c, coding) {
    if (coding == "WIR") param[1] <- -param[1]
    ((-log(1 + exp(-param[1])) - param[2] - g_mu - log(c))^2 / sigma_mu^2
     + (-log(exp(param[1]) + 1) - param[2] - g_lambda - log(c))^2 / sigma_lambda^2) / 2
  }


# Normal-normal priors on logit(phi), log(zeta)

Norm_phi_zeta <- function(g_phi, sigma_phi, g_zeta, sigma_zeta)
  function(param, c, coding) {
    if (coding == "WIR") param[1] <- -param[1]
    ((param[1] - g_phi)^2 / sigma_phi^2 + (param[2] - g_zeta - log(c))^2 / sigma_zeta^2) / 2
  }


#---------------------------------------
# log-likelihood function for MTS
#---------------------------------------

MTS_loglik <- function(param, Tmat, c) {
  phi <- expit(param[1])
  zeta <- exp(param[2])
  p0 <- p_0(c, phi, zeta)
  p1 <- p_1(c, phi, zeta)
  Tmat[1,1] * log(1 - p0) + Tmat[1,2] * log(p0) + Tmat[2,1] * log(1 - p1) + Tmat[2,2] * log(p1)
}

# penalized log-likelihood for MTS

MTS_loglik_pen <- function(param, Tmat, c, penalty_func)
  MTS_loglik(param, Tmat, c) + penalty_func(param, c, coding = "MTS")


#----------------------------------------------
# maximum likelihood estimation for MTS, 
# with optional penalty function
#----------------------------------------------

# transition matrix

Transitions <- function(X)  {
  N <- length(X)  
  crossprod(cbind(1 - X[1:(N-1)], X[1:(N-1)]), cbind(1 - X[2:N], X[2:N]))
}

# Brown, Solomon, & Stephens (1977) estimator

BSSest <- function(Tmat, c) {
  
  if (any(rowSums(Tmat) == 0)) return(c(NA, NA))
  
  P01 <- Tmat[1,2] / sum(Tmat[1,])
  P10 <- Tmat[2,1] / sum(Tmat[2,])  
  S <- P01 + P10
  
  if (S >= 1) {
    est <- c(logit(sum(Tmat[,2]) / sum(Tmat)), NA)
  } else {
    est <- c(logit(P01 / S), log(P01) + log(P10) + log(-log(1 - S)) - log(c) - 2 * log(S))
  }
   return(est)
}

# penalized maximum likelihood estimates for a single response string

MTSmle <- function(X, c, penalty_func = NULL,
                   phi_start = pmin(pmax(mean(X), 1 / length(X)), 1 - 1 / length(X)), zeta_start = .10, 
                   transform = "none") {
  
  if (sum(is.na(X)) > 0) return(c(phi = NA, zeta = NA))
  
  N <- length(X)
  
  # transition matrix
  Tmat <- crossprod(cbind(1 - X[1:(N-1)], X[1:(N-1)]), cbind(1 - X[2:N], X[2:N]))
  
  if (is.null(penalty_func)) {
    est <- BSSest(Tmat, c)
  } else {
    objective <- function(par) MTS_loglik_pen(par, Tmat = Tmat, c = c, penalty_func) 
    results <- optim(par = c(logit(phi_start), log(zeta_start)), fn = objective, 
                     control = list(fnscale = -1))
    est <- results$par
  }
  
  param_trans(est, transform = transform)
  
}

#---------------------------------------
# log-likelihood function for PIR
#---------------------------------------

# calculate conditional probabilities

PIRpsi <- function(phi, zeta, U, c, d) {
  
  K <- length(U)
  p1_cd <- p_1(c + d, phi, zeta)
  p0_cd <- p_0(c + d, phi, zeta)
  p0_d <- p_0(d, phi, zeta)
  e_neg_lambda <- exp(-1 * zeta * c / (1 - phi))
  
  psi <- vector(mode = "numeric", length = K)
  psi[1] <- phi
  for(i in 1:(K - 1)) {
    psi_lag <- psi[i]
    u <- U[i]
    psi[i + 1] <- ((psi_lag * p1_cd + (1 - psi_lag) * (p0_cd - p0_d * e_neg_lambda)) / (1 - (1-psi_lag) * e_neg_lambda))^u * p0_d^(1-u)
  }
  return(psi)
}


# log-likelihood for PIR

PIR_loglik <- function(param, U, c, d) {
  phi <- expit(param[1])
  zeta <- exp(param[2])
  psi <- PIRpsi(phi, zeta, U, c, d)
  loglik <- sum(U * log(1 - (1 - psi) * exp(-1 * zeta * c / (1 - phi))) + (1 - U) * (log(1 - psi) - zeta * c / (1-phi)))
  return(loglik)
}


# penalized log-likelihood for PIR

PIR_loglik_pen <- function(param, U, c, d, penalty_func, coding)
  PIR_loglik(param, U, c, d) + penalty_func(param, c, coding)


#----------------------------------------------
# maximum likelihood estimation for PIR, 
# with optional penalty function
#----------------------------------------------

# maximum likelihood estimates for a single response string

PIRmle <- function(U, c, d, coding = "PIR", penalty_func = NULL,
                   phi_start = max(mean(U) / 2, expit(-10)),
                   zeta_start = .10, transform = "none") {
  
  if(sum(is.na(U)) > 0) return(c(phi = NA, zeta = NA))
  
  if(coding == "WIR") {
    U <- 1 - U
    phi_start <- 1 - phi_start
  }
  
  if (is.null(penalty_func)) {
    objective <- function(par) PIR_loglik(par, U = U, c = c, d = d)
  } else {
    objective <- function(par) {
      PIR_loglik_pen(par, U = U, c = c, d = d, penalty_func, coding = coding)       
    }
  }
    
  results <- optim(par = c(logit(phi_start), log(zeta_start)), fn = objective, 
                   control = list(fnscale = -1))
  
  if(coding == "WIR") results$par[1] <- -results$par[1]
  
  param_trans(results$par, transform = transform)

}

#---------------------------------------
# log-likelihood function for AIR
#---------------------------------------

# transition probabilities for AIR data 

AIR_pi <- function(phi, zeta, c, d) {
  p0_d <- p_0(d, phi, zeta)
  p0_cd <- p_0(c + d, phi, zeta)
  exp_lambda <- exp(- zeta * c / (1 - phi))
  p1_d <- p_1(d, phi, zeta)
  p1_cd <- p_1(c + d, phi, zeta)
  exp_mu <- exp(- zeta * c / phi)
  c(pi_0000 = (1 - p0_d) * exp_lambda, pi_1000 = NA, 
              pi_0100 = p0_d * exp_lambda, pi_1100 = NA,
              pi_0010 = 1 - p0_cd - (1 - p0_d) * exp_lambda, pi_1010 = 1 - p1_cd - (1 - p1_d) * exp_mu, 
              pi_0110 = p0_cd - p0_d * exp_lambda, pi_1110 = p1_cd - p1_d * exp_mu, 
              pi_0001 = NA, pi_1001 = NA, 
              pi_0101 = NA, pi_1101 = NA,
              pi_0011 = NA, pi_1011 = (1 - p1_d) * exp_mu, 
              pi_0111 = NA, pi_1111 = p1_d * exp_mu)
}

# log-likelihood for AIR

AIR_loglik <- function(param, Tmat, c, d) {
  phi <- expit(param[1])
  zeta <- exp(param[2])
  trans_prob <- AIR_pi(phi, zeta, c, d)
  loglik <- sum(log(trans_prob) * Tmat, na.rm = TRUE)
  return(loglik)
}


# penalized log-likelihood for AIR

AIR_loglik_pen <- function(param, Tmat, c, d, penalty_func)
  AIR_loglik(param, Tmat, c, d) + penalty_func(param, c, coding = "AIR")


#----------------------------------------------
# maximum likelihood estimation for AIR, 
# with optional penalty function
#----------------------------------------------

# maximum likelihood estimates for a single response string

AIR_Tmat <- function(XUW) {
  K <- nrow(XUW) - 1
  MTS_lag <- factor(XUW[1:K, 1], levels = 0:1)
  table(MTS_lag, 
        factor(XUW[-1,1], levels = 0:1),
        factor(XUW[-1,2], levels = 0:1),
        factor(XUW[-1,3], levels = 0:1))
}

AIRmle <- function(XUW, c, d, penalty_func = NULL,
                   phi_start = min(max(mean(XUW[,1]), expit(-10)), 1 - expit(-10)),
                   zeta_start = .10, transform = "none") {
  
  if (sum(is.na(XUW[-1,])) + is.na(XUW[1,1]) > 0) return(c(phi = NA, zeta = NA))
  
  Tmat <- AIR_Tmat(XUW)
    
  if (is.null(penalty_func)) {
    objective <- function(par) AIR_loglik(par, Tmat = Tmat, c = c, d = d)
  } else {
    objective <- function(par) {
      AIR_loglik_pen(par, Tmat = Tmat, c = c, d = d, penalty_func)       
    }
  }
  
  results <- optim(par = c(logit(phi_start), log(zeta_start)), fn = objective, 
                   control = list(fnscale = -1))
    
  param_trans(results$par, transform = transform)
  
}

# mu <- 6
# lambda <- 9
# c <- 0.8
# d <- 0.2
# XUW <- r_AIR(n = 1, mu = mu, lambda = lambda, stream_length = 100000,
#              F_event = F_exp(), F_interim = F_exp(),
#              interval_length = c + d, rest_length = d)[,,1]
# 
# phi <- mu / (mu + lambda)
# zeta <- 1 / (mu + lambda)
# Tmat <- AIR_Tmat(XUW)
# Tmat[1,,,] <- Tmat[1,,,] / sum(Tmat[1,,,])
# Tmat[2,,,] <- Tmat[2,,,] / sum(Tmat[2,,,])
# trans_prob <- AIR_pi(phi, zeta, c, d)
# check_pi <- data.frame(pi = trans_prob, m = as.vector(Tmat))
# with(check_pi, all(m[is.na(pi)] == 0))
# subset(check_pi, !is.na(pi))
# 
# AIRmle(XUW, c, d, transform = "exp")
# c(phi, zeta)