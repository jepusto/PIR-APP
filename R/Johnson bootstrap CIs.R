library(plyr)
rm(list=ls())
source("R/parallel setup.R")
source("R/PIR-APP.R")

johnson <- subset(read.csv("data/Interval recording data from Johnson dissertation.csv"), Behavior == 1)
names(johnson)[4] <- "Clip"

johnson <- within(johnson, {
  Clip <- as.factor(Clip)
  levels(Method) <- c("MTS","PIR","WIR")
  levels(Clip) <- paste("Clip", levels(Clip))
  Prevalence <- ObsInts / TotalInts  
})

MTS <- subset(johnson, Method == "MTS")
PIR <- subset(johnson, Method == "PIR")
WIR <- subset(johnson, Method == "WIR")

MTmeans <- ddply(MTS, .variables= "Clip", summarize, Prevalence = mean(Prevalence))

source_func <- ls()

iterations <- 100
p <- 0.05
penalty_func <- Beta_Gamma(k_mu = 1.5, k_lambda = 1.5, theta_mu = 10, theta_lambda = 10, const = 2)

# set up parallel processing
cluster <- start_parallel(source_func)
clusterEvalQ(cluster, library(ARPobservation))

MTS_PLEs <- adply(as.matrix(MTS[,9:48]), 1, MTSbootstrap, 
                  c = 1, penalty_func = penalty_func, iterations = iterations, 
                  transform = "exp", .parallel = TRUE)
PIR_PLEs <- adply(as.matrix(PIR[,9:48]), 1, PIRbootstrap, coding = "PIR",
                  c = 1, d = 0, penalty_func = penalty_func, iterations = iterations, 
                  transform = "exp", .parallel = TRUE)
WIR_PLEs <- adply(as.matrix(WIR[,9:48]), 1, PIRbootstrap, coding = "WIR",
                  c = 1, d = 0, penalty_func = penalty_func, iterations = iterations, 
                  transform = "exp", .parallel = TRUE)

stopCluster(cluster)

johnson$X1 <- rownames(johnson)
PLEs <- merge(subset(johnson, select = c(X1, Clip, Method, Prevalence)), rbind(MTS_PLEs, PIR_PLEs, WIR_PLEs))
phi <- subset(PLEs, parm == "phi", select = -parm)
zeta <- subset(PLEs, parm == "zeta", select = -parm)

ggplot(phi, aes(x = Prevalence, y = est, ymin = CI_L, ymax = CI_U, color = Method, shape = Method)) + 
  geom_pointrange(position = position_jitter(w = 0.02)) + 
  geom_hline(data = MTmeans, aes(yintercept = Prevalence), linetype = "dashed") + 
  facet_wrap( ~ Clip, ncol = 3) + 
  coord_cartesian(xlim = c(-0.02,1.02), ylim = c(-0.02,1.02)) + 
  scale_x_continuous(breaks = seq(0.2, 1, 0.2)) +
  labs(x = "Proportion of intervals", y = "PLE prevalence estimate", 
       color = "Recording system", shape = "Recording system") + 
  theme_bw() 

phi_summary <- ddply(phi, .(Clip, Method), summarize, Summary_prop = mean(Prevalence), PLE = expit(mean(logit(est))))
dcast(melt(phi_summary, id.vars = c("Clip","Method")), Clip ~ variable  + Method)
