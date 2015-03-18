library(dplyr)
library(reshape2)
library(ggplot2)

rm(list=ls())

load("Simulations/MTS simulations/MTS sample size sims.Rdata")
names(results)[3] <- "K"

results <- subset(results, !(stat %in% c("mean","logit mean")) | (k_priors!=1 & theta!=Inf))

results <- within(results, {
  rmse <- sqrt(bias^2 + var)
  rmse_trunc <- pmin(rmse, 2)
  rrmse <- ifelse(stat == "zeta", rmse / zeta, NA)
  alpha <- factor(k_priors)
  estimator <- ifelse(k_priors==1 & theta==Inf, "MLE", ifelse(stat %in% c("mean","logit mean"), "mean", "PLE"))
  stat[stat=="logit mean"] <- "logit phi"
  stat[stat=="mean"] <- "phi"
  stat <- ordered(stat, c("phi","logit phi","zeta", "log zeta","mu","log mu","lambda","log lambda"))
})

with(results, table(stat, estimator))


# RMSE of phi

ggplot(subset(results, K %in% seq(30,150,30) & stat == "phi"),
       aes(x = estimator, y = rmse_trunc, fill = estimator)) + 
  geom_boxplot(position = "dodge") + 
  facet_grid(. ~ K, scales = "free", labeller = "label_both") + 
  labs(y = "Root mean squared error") + 
  theme_bw() + theme(legend.position = "none")

# Relative RMSE of zeta

ggplot(subset(results, K %in% seq(30,150,30) & stat == "zeta"),
       aes(x = estimator, y = rrmse, fill = estimator)) + 
  geom_boxplot(position = "dodge") + 
  coord_cartesian(ylim = c(0,3)) + 
  facet_grid(. ~ K, scales = "free") + 
  labs(y = "Relative root mean squared error") + 
  theme_bw() + theme(legend.position = "none")

# RMSE of logit-phi, log-zeta

ggplot(subset(results, K %in% seq(30,150,30) & stat %in% c("logit phi","log zeta")),
       aes(x = estimator, y = rmse, fill = estimator)) + 
  geom_boxplot(position = "dodge") + 
  coord_cartesian(ylim = c(0,2)) + 
  facet_grid(stat ~ K, scales = "free") + 
  labs(y = "Root mean squared error") + 
  theme_bw() + theme(legend.position = "none")

#--------------------------------------
# Minimum sample size requirements
#--------------------------------------

SS <- within(subset(results, stat %in% c("phi","logit phi","zeta","log zeta")), {
  meets <- (stat == "phi" & abs(bias) < 0.03) | (stat == "logit phi" & abs(bias) < 0.05) |
            (stat == "zeta" & abs(bias) / zeta < 0.1) | (stat == "log zeta" & abs(bias) < 0.1)
  criterion <- factor(stat, labels = c("Absolute Bias(phi) < 0.03", "Absolute Bias(logit phi) < 0.05",
                                       "Absolute Relative Bias(zeta) < 0.10", "Absolute Bias(log zeta) < 0.10"))
})

SS_min <- filter(SS, meets | K == 150) %>%
  group_by(estimator, criterion, phi, zeta) %>%
  summarize(min_K = min(K))

SS_min <- within(SS_min, {
  min_K <- factor(min_K)
  levels(min_K)[15] <- "150+"
  min_K_fac <- min_K
  levels(min_K_fac) <- c(seq(10,30,10), rep(seq(40,150,20), each = 2))
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

