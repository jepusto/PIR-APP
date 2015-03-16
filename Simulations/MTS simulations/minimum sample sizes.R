library(ggplot2)
library(reshape2)
rm(list=ls())
load("Simulations/MTS simulations/MTS_sample_size_sim.Rdata")

results$estimator <- with(results, ifelse(k_priors == 1 & theta == Inf, "MLE", "PLE"))


# relative bias in zeta

zeta <- subset(results, stat=="zeta")
zeta <- within(zeta, {
  rel_bias <- bias / zeta 
  RB_round <- pmax(pmin(rel_bias, 1), -0.5)
})

print(dcast(zeta, estimator + zeta ~ phi, value.var = "RB_round"), digits = 2)

ggplot(zeta, aes(x = phi, y = zeta, fill = RB_round)) + 
  geom_tile() +
  facet_grid(. ~ estimator)

# bias in phi

phi <- subset(results, stat == "phi")
phi$bias_r <- round(phi$bias, 3)
print(dcast(phi, estimator + zeta ~ phi, value.var = "bias_r"), digits = 2)
