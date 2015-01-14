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


#---------------------------------------------------
# penalty functions and penalized log-likelihood
#---------------------------------------------------

# Note that priors are on (phi, 1 / zeta) when const = 1
# For inverse-gamma prior on zeta, use const = -1
# For priors on (mu, lambda), use const = 2

Beta_Gamma <- function(k_mu, k_lambda, theta, const = 1) {
  function(param) {
    -(k_mu - 1) * log(1 + exp(-param[1])) - 
      (k_lambda - 1) * (param[1] + log(1 + exp(-param[1]))) - 
      (k_mu + k_lambda - const) * param[2] - exp(-param[2]) / theta
  }  
}

# Normal-normal priors on log(mu), log(lambda)

Normal_Normal <- function(mu_hyp, lambda_hyp, sigma_hyp)
  function(param) {
    ((log(param[1]) - log(param[2]) - log(mu_hyp))^2 
      + (log(1 - param[1]) - log(param[2]) - log(lambda_hyp))^2 ) / (2 * sigma_hyp^2)
  }
        
PIR_loglik_pen <- function(param, U, c, d, penalty_func)
  PIR_loglik(param, U, c, d) + penalty_func(param)

#-------------------------------------------------------------------
# maximum likelihood estimation, with optional penalty function
#-------------------------------------------------------------------

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
    objective <- function(par) PIR_loglik_pen(par, U = U, c = c, d = d, penalty_func) 
  }
    
  results <- optim(par = c(logit(phi_start), log(zeta_start)), fn = objective, 
                   control = list(fnscale = -1))
  
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
                         coding = "PIR", penalty_func = NULL,
                         iterations = 2000, p = .05) {
  
  U <- generateSamples(phi = phi, zeta = zeta, c = c, d = d, 
                       intervals = intervals, coding = coding, iterations = iterations)
  
  ests <- t(apply(X = U, MARGIN = 1, FUN = PIRmle, c = c, d = d, 
                  coding = coding, penalty_func = penalty_func))
  
  phi_CI <- quantile(ests[,1], probs = c(p/2, 1-p/2))
  zeta_CI <- quantile(ests[,2], probs = c(p/2, 1-p/2))
  
  return(list(phi_CI = phi_CI, zeta_CI = zeta_CI))
}

# get both the maximum likelihood estimates and CIs for a given response string

PIRmleCI <- function(U, c, d, coding = "PIR", penalty_func = NULL,
                     phi_start = max(mean(U) / 2, expit(-10)), zeta_start = .10, 
                     transform = "none", iterations = 2000, p = .05){
  
  ests <- PIRmle(U = U, c = c, d = d, coding = coding, penalty_func = penalty_func,
                 phi_start = phi_start, zeta_start = zeta_start)
  
  CI <- PIRbootstrap(phi = expit(ests[1]), zeta = exp(ests[2]), 
                     c = c, d = d, intervals = length(U), 
                     coding = coding, penalty_func = penalty_func, 
                     iterations = iterations, p = p)
  
  if(transform == "none"){ 
    return(list(phi = ests[1], phi_CI = CI$phi_CI, 
           zeta = ests[2], zeta_CI = CI$zeta_CI))
  }
  
  if(transform == "exp"){
    return(list(phi = expit(ests[1]), phi_CI = expit(CI$phi_CI),
                zeta = exp(ests[2]), zeta_CI = exp(CI$zeta_CI)))
  }

}

