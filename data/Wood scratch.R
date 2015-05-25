library(ARPobservation)
library(ggplot2)
library(plyr)
library(reshape2)
source("R/PIR-APP.R")
source("R/parallel setup.R")
wood <- read.csv("data/Wood et al data.csv", na.strings = ".", stringsAsFactors = FALSE)

PIR_wrap <- function(data_line, c = 1, d = 0, k_prior = 1, theta_prior = Inf){
  
  U <- data_line[4:70]
  U <- as.numeric(U[(U == 1|U == 0) & !is.na(U)])
  
  results <- PIRbootstrap(U = U, c = c, d = d, coding = data_line[3], 
         penalty_func = Beta_Gamma(k_mu = k_prior, k_lambda = k_prior, 
                                   theta_mu = theta_prior, 
                                   theta_lambda = theta_prior), transform = "exp", iterations = 10)
  
  summary <- mean(U)
  
  results$summary <- summary
  results
}


source_func <- ls()
cluster <- start_parallel(source_func)
clusterEvalQ(cluster, library(ARPobservation))

system.time(ple <- adply(wood, 1, PIR_wrap, k_prior = 1.5, theta_prior = 10, .parallel = TRUE))

stopCluster(cluster)

ple <- ple[,-(4:70)]

cdr <- read.csv("data/transcribed.csv")[,-1]

ple <- merge(ple, cdr)

PIR <- subset(ple, Measurement.Type == "PIR")
WIR <- subset(ple, Measurement.Type == "WIR")

names(PIR)[5:9] <- paste0("PIR", names(PIR[5:9]))
names(WIR)[5:9] <- paste0("WIR", names(WIR[5:9]))

ple_wide <- merge(PIR[,-3],WIR[,-3])

pl3 <- dcast(ple2, Classroom.ID + Child + parm ~ Measurement.Type, value.var = c("parm", "est", "sd", "CI_L", "CI_U"))


head(ple)

qplot(summary, est, ymin = CI_L, ymax = CI_U, color = factor(Measurement.Type), geom = "pointrange", data = subset(wood, parm == "phi" & estimator == "mle"))

qplot(summary, est, ymin = CI_L, ymax = CI_U, color = factor(Measurement.Type), geom = "pointrange", data = subset(wood, parm == "phi" & estimator == "ple"))


qplot(summary, est, ymin = CI_L, ymax = CI_U, color = factor(Measurement.Type), geom = "pointrange", data = subset(wood, parm == "phi"))+ facet_wrap(~estimator) 


qplot(summary, est, ymin = CI_L, ymax = CI_U, color = factor(Measurement.Type), geom = "pointrange", data = subset(wood, parm == "zeta" & estimator == "mle"))
qplot(summary, est, ymin = CI_L, ymax = CI_U, color = factor(Measurement.Type), geom = "pointrange", data = subset(wood, parm == "zeta" & estimator == "ple"))
