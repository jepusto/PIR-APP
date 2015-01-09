logit <- function(x) log(x) - log(1 - x)
expit <- function(x) 1 / (1 + exp(-x))

MTSmle <- function(U, c) {

  if (sum(is.na(U)) > 0) return(c(phi = NA, zeta = NA))
  
  N <- length(U)

  # transition matrix
  K <- crossprod(cbind(1 - U[1:(N-1)], U[1:(N-1)]), cbind(1 - U[2:N], U[2:N]))
  
  if (any(rowSums(K) == 0)) return(c(phi = NA, zeta = NA))
  
  P01 <- K[1,2] / sum(K[1,])
  P10 <- K[2,1] / sum(K[2,])  
  S <- P01 + P10
  
  if (S >= 1) {
    phi <- sum(K[,2]) / (N - 1)
    zeta <- NA
  } else {
    phi <- P01 / S
    zeta <- - P01 * P10 * log(1 - S) / (c * S^2)
  }
  return(c(phi = phi, zeta = zeta))    
}

MTSmle_cov <- function(phi, zeta, c, intervals) {
  
  muI <- zeta / phi
  lambdaI <- zeta / (1 - phi)
  
  g <- exp((muI + lambdaI) * c)
  
  G_matrix <- matrix(nrow = 2, ncol = 2)
  G_matrix[1,1] <- muI + lambdaI * (g-1) / (c * (lambdaI + muI))
  G_matrix[1,2] <- -lambdaI * (1 - (g - 1) / (c * (lambdaI + muI)))
  G_matrix[2,1] <- -muI * (1 - (g - 1) / (c * (lambdaI + muI)))
  G_matrix[2,2] <- lambdaI + muI * (g - 1) / (c * (lambdaI + muI))
  
  G <- G_matrix * g  / (g - 1)
  
  S_matrix <- matrix(0, nrow = 2, ncol = 2)
  S_matrix[1,1] <- lambdaI * (g + lambdaI / muI)
  S_matrix[2,2] <- muI * (g + muI / lambdaI)
  
  S <- S_matrix * (g-1) / ((lambdaI + muI) * g^2)
  
  J <- matrix(nrow = 2, ncol = 2)
  J[1,1] <- 1 / lambdaI
  J[1,2] <- -1 / muI
  J[2,1] <- muI / (lambdaI * (muI + lambdaI))
  J[2,2] <- lambdaI / (muI * (muI + lambdaI))
  
  J %*% G %*% S %*% t(G) %*% t(J) / intervals
}

MTSmle_CI <- function(U, c, p = .05) {
  ests <- MTSmle(U = U, c = c)
  
  if (is.na(ests[2])) {
    return(list(phi = ests[1], phi_CI = c(NA,NA), 
                zeta = ests[2], zeta_CI = c(NA,NA)))
  }
  
  V_logparms <- MTSmle_cov(phi = ests[1], zeta = ests[2], c = c, intervals = length(U))
  sd_logparms <- sqrt(diag(V_logparms))
  
  zconf <- qnorm(c(p/2, 1-p/2))
  
  lphiCI <- logit(ests[1]) + zconf * sd_logparms[1]
  lzetaCI <- log(ests[2]) + zconf * sd_logparms[2]
  
  phi_CI <- expit(lphiCI)
  zeta_CI <- exp(lzetaCI)
  
  return(list(phi = ests[1], phi_CI = phi_CI, 
              zeta = ests[2], zeta_CI = zeta_CI))
}
