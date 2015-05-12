
library(ARPobservation)
library(ggplot2)
library(plyr)
source("R/PIR-APP.R")
source("R/parallel setup.R")
wood <- read.csv("data/Wood et al data.csv", na.strings = ".", stringsAsFactors = FALSE)

PIR_wrap <- function(data_line, c = 1, d = 0, k_prior = 1, theta_prior = Inf){
  
  U <- data_line[4:70]
  U <- as.numeric(U[(U == 1|U == 0) & !is.na(U)])
  
  results <- PIRbootstrap(U = U, c = c, d = d, coding = data_line[3], 
         penalty_func = Beta_Gamma(k_mu = k_prior, k_lambda = k_prior, 
                                   theta_mu = theta_prior, 
                                   theta_lambda = theta_prior), transform = "exp", iterations = 100)
  
  summary <- mean(U)
  
  results$summary <- summary
  results
}


mle <- adply(wood, 1, PIR_wrap)

ple <- adply(wood, 1, PIR_wrap, k_prior = 1.5, theta_prior = 10)

mle$estimator <- "mle"
ple$estimator <- "ple"

wood <- rbind(mle, ple)


qplot(summary, est, ymin = CI_L, ymax = CI_U, color = factor(Measurement.Type), geom = "pointrange", data = subset(wood, parm == "phi" & estimator == "mle"))

qplot(summary, est, ymin = CI_L, ymax = CI_U, color = factor(Measurement.Type), geom = "pointrange", data = subset(wood, parm == "phi" & estimator == "ple"))


qplot(summary, est, ymin = CI_L, ymax = CI_U, color = factor(Measurement.Type), geom = "pointrange", data = subset(wood, parm == "phi"))+ facet_wrap(~estimator) 


qplot(summary, est, ymin = CI_L, ymax = CI_U, color = factor(Measurement.Type), geom = "pointrange", data = subset(wood, parm == "zeta" & estimator == "mle"))
qplot(summary, est, ymin = CI_L, ymax = CI_U, color = factor(Measurement.Type), geom = "pointrange", data = subset(wood, parm == "zeta" & estimator == "ple"))
