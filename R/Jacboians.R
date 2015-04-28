#Functions that create Jacobian matrices from logit(phi) and log(zeta)

jacob_phizeta <- function(pars) {
  J <- matrix(data = c(0,0,0,0), nrow = 2, ncol = 2)
  J[1,1] <- exp(-pars[1])/(1 + exp(-pars[1]))^2
  J[2,2] <- exp(pars[2])
  return(J)
}


jacob_mulambda <- function(pars) {
  J <- matrix(data = c(0,0,0,0), nrow = 2, ncol = 2)
  J[1,1] <- exp(-pars[1]) / (1 + exp(-pars[1]))^2 / exp(pars[2])
  J[1,2] <- -exp(-pars[2]) / (1 + exp(-pars[1]))
  J[2,1] <- -exp(pars[1]) / (1 + exp(pars[1]))^2 / exp(pars[2])
  J[2,2] <- -exp(-pars[2]) / (1 + exp(pars[1]))
  
  return(J)
}


jacob_logmulambda <-function(pars) {
  J <- matrix(data = c(0,-1,0,-1), nrow = 2, ncol = 2)
  J[1,1] <- 1 / (1 + exp(pars[1]))
  J[1,2] <- -1 / (1 + exp(-pars[1]))
  return(J)
}
