library(ARPobservation)
library(plyr)

#---------------------------
# transformation functions 
#---------------------------

# simple transform functions

logit <- function(x) log(x) - log(1 - x)
expit <- function(x) 1 / (1 + exp(-x))


# transform functions to get phi and zeta

trans_phizeta <- function(pars) c(phi = expit(pars[1]), zeta = exp(pars[2]))


# transform functions to get mu and lambda

trans_mulambda <- function(pars) {
  phi <- expit(pars[1])
  zeta <- exp(pars[2])
  mu <- as.numeric(phi/zeta)
  lambda <- as.numeric((1-phi)/zeta)
  trans <- c(mu = mu, lamdba = lambda)
  trans
}


# transform functions to get log mu and log lambda

trans_logmulambda <- function(pars) {
  trans <- log(trans_mulambda(pars))
  names(trans) <- c("log_mu","log_lambda")
  trans
}


#---------------------------------------
# conditional probability functions
#---------------------------------------

p_0 <- function(x,phi,zeta) phi * (1- exp(-1 * x * zeta / (phi * (1- phi))))
p_1 <- function(x,phi,zeta) phi + (1 - phi) * exp(-1 * x * zeta / (phi * (1- phi)))



#---------------------------------------
# log-likelihood functions 
#---------------------------------------

# calculate interval probabilities

PIRpsi <- function(phi, zeta, U, c, d) {
  
  psi <- vector(mode = "numeric", length = length(U))
  psi[1] <- phi
  for(i in 2:length(U)){
    psi[i] <- ((psi[i-1] * p_1(c + d,phi,zeta) + (1 - psi[i-1]) * (p_0(c + d,phi,zeta) - p_0(d, phi, zeta) * exp(-1 * zeta * c / (1-phi)))) / (1 - (1-psi[i-1]) * exp(-1 * zeta * c/(1-phi))))^U[i-1] * p_0(d,phi,zeta)^(1-U[i-1]) 
  }
  return(psi)
}


# log-likelihood of a particular U

PIR_loglik <- function(phi, zeta, U, c, d) {
  psi <- PIRpsi(phi, zeta, U, c, d)
  loglik <- sum(U * log(1 - (1 - psi) * exp(-1 * zeta * c / (1-phi))) + 
                  (1 - U) * (log(1 - psi) - zeta * c / (1-phi)))
  return(loglik)
}

# penalized log-likelihood

PIR_pll <- function(phi, zeta, U, c, d, mu_hyp, lambda_hyp, theta_hyp) {
  log_lik <- PIR_loglik(phi = phi, zeta = zeta, U = U, c = c, d = d)
  sigma_sq <- log(theta_hyp + 1)
  penalty <- dnorm(log(phi / zeta), mean = log(mu_hyp), sd = sqrt(sigma_sq), log = TRUE) + 
    dnorm(log((1 - phi) / zeta), mean = log(lambda_hyp), sd = sqrt(sigma_sq), log = TRUE)
  c(log_lik = log_lik, penalty = penalty, pll = log_lik + penalty)
}



# function to be maximized - optim requires all maximized values to be contained in a single vector

ll2 <- function(pars, U, c, d, ...)  
  PIR_loglik(phi = expit(pars[1]), zeta = exp(pars[2]), U = U, c = c, d = d)

pll2 <- function(pars, U, c, d, mu_hyp, lambda_hyp, theta_hyp, ...) {
  PIR_pll(phi = expit(pars[1]), zeta = exp(pars[2]), U = U, c = c, d = d,
      mu_hyp = mu_hyp, lambda_hyp = lambda_hyp, theta_hyp = theta_hyp)["pll"]
} 


#---------------------------------------
# maximum likelihood estimation
#---------------------------------------

# get maxmimum likelihood estimates for a single response string

PIRmle <- function(U, c, d, coding = "PIR", objective = ll2,
                   phi_start = .50, zeta_start = .10, transform = "none", ...) {
  
  if(sum(is.na(U)) > 0) return(c(phi = NA, zeta = NA))
  
  if(coding == "WIR") U <- 1 - U
  
  results <- optim(par = c(logit(phi_start), log(zeta_start)), fn = objective, 
                   U = U, c = c, d = d, 
                   control = list(fnscale = -1), ...)
  
  if(coding == "WIR") results$par[1] <- -results$par[1]
  
  
  if (transform == "none") mle <- c(logit_phi = results$par[1], log_zeta = results$par[2])
  
  if (transform == "exp") mle <- trans_phizeta(results$par)
    
  if (transform == "means") mle <- trans_mulambda(results$par)
    
  if (transform == "logmeans") mle <- trans_logmulambda(results$par)
  
  return(mle)
  
}

#---------------------------------------
# bootstrap confidence intervals
#---------------------------------------

#generate a single sample of behavior for a given value of phi and zeta using ARPobservation

generateSample <- function(phi, zeta, c, d, intervals, coding) {  
  mu <- phi/zeta
  lambda <- (1-phi)/zeta
  intervalLength <- c + d
  
  BS <- r_behavior_stream(n = 1, mu = mu, lambda = lambda, 
                          F_event = F_exp(), F_interim = F_exp(), 
                          stream_length = intervalLength * intervals)
  partial <- (coding == "PIR")
  interval_recording(BS = BS, interval_length = intervalLength, rest_length = d, partial = partial, summarize = FALSE)
}

#generate "iterations" number of samples of behavior for a given phi/zeta

generateSamples <- function(phi, zeta, c, d, intervals, coding, iterations) {
  samples <- replicate(iterations,
                       generateSample(phi = phi, zeta = zeta, c = c, d = d, 
                                      intervals = intervals, coding = coding))[,1,]
  t(samples)
}

# bootstrap a confidence interval for a given phi/zeta

PIRbootstrap <- function(phi, zeta, c, d, intervals, 
                         coding = "PIR", objective = ll2,
                         phi_start = .50, zeta_start = .10, 
                         iterations = 2000, p = .05, ...) {
  
  U <- generateSamples(phi = phi, zeta = zeta, c = c, d = d, 
                       intervals = intervals, coding = coding, iterations = iterations)
  
  ests <- t(apply(X = U, MARGIN = 1, FUN = PIRmle, c = c, d = d, 
                  phi_start = phi_start, zeta_start = zeta_start, 
                  coding = coding, objective = objective, ...))
  
  phi_CI <- quantile(ests[,1], probs = c(p/2, 1-p/2))
  zeta_CI <- quantile(ests[,2], probs = c(p/2, 1-p/2))
  
  return(list(phi_CI = phi_CI, zeta_CI = zeta_CI))
}

# get both the maximum likelihood estimates and CIs for a given response string

PIRmleCI <- function(U, c, d, coding = "PIR", objective = ll2,
                     phi_start = .50, zeta_start = .10, 
                     transform = "none", iterations = 2000, p = .05, ...){
  
  ests <- PIRmle(U = U, c = c, d = d, coding = coding, objective = objective,
                 phi_start = phi_start, zeta_start = zeta_start, ...)
  
  CI <- PIRbootstrap(phi = expit(ests[1]), zeta = exp(ests[2]), 
                     c = c, d = d, intervals = length(U), 
                     coding = coding, objective = objective, 
                     iterations = iterations, p = p, ...)
  
  if(transform == "none"){ 
    return(list(phi = ests[1], phi_CI = CI$phi_CI, 
           zeta = ests[2], zeta_CI = CI$zeta_CI))
  }
  
  if(transform == "exp"){
    return(list(phi = expit(ests[1]), phi_CI = expit(CI$phi_CI),
                zeta = exp(ests[2]), zeta_CI = exp(CI$zeta_CI)))
  }

}

