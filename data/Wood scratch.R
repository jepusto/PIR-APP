rm(list = ls())
source("R/PIR-APP.R")
wood <- read.csv("data/Wood et al data.csv", na.strings = ".", stringsAsFactors = FALSE)

PIR_wrap <- function(data_line, c = 1, d = 0, k_prior = 1, theta_prior = Inf){
  
  U <- data_line[4:70]
  U <- as.numeric(U[(U == 1|U == 0) & !is.na(U)])
  
  PIRmle(U = U, c = c, d = d, coding = data_line[3], 
         penalty_func = Beta_Gamma(k_mu = k_prior, k_lambda = k_prior, 
                                   theta_mu = theta_prior, 
                                   theta_lambda = theta_prior), transform = "exp")
}

results <- t(apply(wood, 1, PIR_wrap))

results2 <- t(apply(wood, 1, PIR_wrap, k_prior = 1.5, theta_prior = 10))

MLE <- cbind(wood, results)
PLE <- cbind(wood, results2)
MLE$estimator = "MLE"
PLE$estimator = "PLE"

wood_long <- rbind(MLE, PLE)

colnames(results) <- c("phi_mle", "zeta_mle")
colnames(results2) <- c("phi_ple", "zeta_ple")

wood_wide <- cbind(wood, results, results2)

head(wood_long)
head(wood_wide)

qplot(phi, zeta, color = factor(Measurement.Type), data = wood_long)

qplot(phi_mle, phi_ple, color = factor(Measurement.Type), data = wood_wide)
qplot(zeta_mle, zeta_ple, color = factor(Measurement.Type), data = wood_wide)
