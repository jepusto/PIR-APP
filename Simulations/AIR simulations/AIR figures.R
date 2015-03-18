library(dplyr)
library(reshape2)
library(ggplot2)

rm(list=ls())

load("Simulations/AIR simulations/AIR sample size sims.Rdata")
names(results)[3] <- "K"

results <- within(results, {
  rmse <- sqrt(bias^2 + var)
  rmse_trunc <- pmin(rmse, 2)
  rrmse <- ifelse(stat == "zeta", rmse / zeta, NA)
  alpha <- factor(k_priors)
  estimator <- ifelse(k_priors==1 & theta==Inf, "MLE", "PLE")
  stat <- ordered(stat, c("phi","logit phi","zeta", "log zeta","mu","log mu","lambda","log lambda"))
})

with(results, table(stat, estimator))


#--------------------------------------
# Minimum sample size requirements
#--------------------------------------

SS <- within(subset(results, stat %in% c("phi","logit phi","zeta","log zeta")), {
  meets <- (stat == "phi" & abs(bias) < 0.03) | (stat == "logit phi" & abs(bias) < 0.05) |
    (stat == "zeta" & abs(bias) / zeta < 0.1) | (stat == "log zeta" & abs(bias) < 0.1)
  criterion <- factor(stat, labels = c("Absolute Bias(phi) < 0.03", "Absolute Bias(logit phi) < 0.05",
                                       "Absolute Relative Bias(zeta) < 0.10", "Absolute Bias(log zeta) < 0.10"))
})

SS_min <- filter(SS, meets | K == 75) %>%
  group_by(estimator, criterion, phi, zeta) %>%
  summarize(min_K = min(K))

SS_min <- within(SS_min, {
  min_K <- factor(min_K)
  levels(min_K)[15] <- "75+"
  min_K_fac <- min_K
  levels(min_K_fac) <- c(seq(5,15,5), rep(seq(20,75,10), each = 2))
})

table(SS_min$min_K_fac)

ggplot(subset(SS_min, estimator != "mean"), 
       aes(phi, zeta, fill = min_K_fac, label = min_K)) +
  geom_tile(alpha = 0.8) + 
  geom_text(size = 4) + 
  scale_fill_brewer(type = "seq", palette = 3) + 
  facet_grid(estimator ~ criterion) + 
  labs(fill = "Minimum sample size") + 
  theme_bw() + theme(legend.position = "none")

ggplot(subset(SS_min, estimator == "PLE"),
       aes(phi, zeta, fill = min_K_fac, label = min_K)) +
  geom_tile(alpha = 0.8) + 
  geom_text(size = 4) + 
  scale_fill_brewer(type = "seq", palette = 3) + 
  facet_wrap(~ criterion) + 
  labs(fill = "Minimum sample size") + 
  theme_bw() + theme(legend.position = "none")


#--------------------------------------
# Compare RMSE
#--------------------------------------

load("Simulations/AIR simulations/AIR sample size sims.Rdata")
results$proc <- "AIR"
results$K_intervals <- results$K_intervals * 2
results_AIR <- results
results$phi <- 1 - results$phi
results_AIR <- rbind(results_AIR, subset(results, phi > 0.5))
with(results_AIR, table(phi, zeta))

load("Simulations/MTS simulations/MTS sample size sims.Rdata")
results$proc <- "MTS"
results_MTS <- subset(results, ! stat %in% c("mean","logit mean"))
results$phi <- 1 - results$phi
results_MTS <- rbind(results_MTS, subset(results, phi > 0.5 & ! stat %in% c("mean","logit mean")))
with(results_MTS, table(phi, zeta))

load("Simulations/PIR simulations/PIR sample size sims.Rdata")
results$proc <- "PIR"
results_PIR <- subset(results, ! stat %in% c("mean","logit mean"))
with(results_PIR, table(phi, zeta))

results <- rbind(results_MTS, results_PIR, results_AIR)
names(results)[3] <- "K"
with(results, table(phi, zeta))

results <- within(results, {
  rmse <- sqrt(bias^2 + var)
  rmse_trunc <- pmin(rmse, 2)
  rrmse <- ifelse(stat == "zeta", rmse / zeta, NA)
  alpha <- factor(k_priors)
  estimator <- ifelse(k_priors==1 & theta==Inf, "MLE", "PLE")
  stat <- ordered(stat, c("phi","logit phi","zeta", "log zeta","mu","log mu","lambda","log lambda"))
})

with(results, table(stat, proc, estimator))


# RMSE of phi

ggplot(subset(results, K %in% c(30,60,90,120) & stat == "phi" & estimator == "PLE"),
       aes(x = proc, y = rmse_trunc, fill = proc)) + 
  geom_boxplot(position = "dodge") + 
  facet_grid(. ~ K, scales = "free", labeller = "label_both") + 
  labs(y = "Root mean squared error") + 
  theme_bw()

# Relative RMSE of zeta

ggplot(subset(results, K %in% c(30,60,90,120) & stat == "zeta" & estimator == "PLE"),
       aes(x = proc, y = rrmse, fill = proc)) + 
  geom_boxplot(position = "dodge") + 
  facet_grid(. ~ K, scales = "free") + 
  labs(y = "Relative root mean squared error") + 
  theme_bw()

# RMSE on logit/log scale

ggplot(subset(results, K %in% c(30,60,90,120) & stat %in% c("logit phi","log zeta") & estimator == "PLE"),
       aes(x = proc, y = rmse, fill = proc)) + 
  geom_boxplot(position = "dodge") + 
  facet_grid(stat ~ K, scales = "free") + 
  labs(y = "Root mean squared error") + 
  theme_bw()
