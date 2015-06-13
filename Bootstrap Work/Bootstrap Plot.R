rm(list = ls())
library(ggplot2)
library(plyr)
load("Bootstrap Work/Bootstrap Check.RData")

results$covered <- with(results, ifelse(parm == "phi", phi >= CI_L & phi <= CI_U, 
                                        zeta >= CI_L & zeta <= CI_U))

ddply(results, .(phi, zeta, parm), summarize, sum(covered)/200)

results$panel <- with(results,paste("phi = ", phi, "& zeta = ", zeta))

phi <- subset(results, parm == "phi")
zeta <- subset(results, parm == "zeta")

zeta <- zeta[order(zeta$phi, zeta$est),]
zeta$case <- 1:200

true_zeta <- ddply(results, .(phi), summarize, phi = phi[1], zeta = zeta[1])

CIs <- ggplot(data = zeta, aes(x = est, y = case, xmin = CI_L, xmax = CI_U))

CIs + geom_point(color = "grey") + 
      geom_errorbarh(color = "grey") + 
      geom_vline(aes(xintercept = zeta), color = "red")+
      facet_wrap(~panel)

phi <- phi[order(phi$phi, phi$est),]
phi$case <- 1:200

true_zeta <- ddply(results, .(phi), summarize, phi = phi[1], zeta = zeta[1])

CIs <- ggplot(data = phi, aes(x = est, y = case, xmin = CI_L, xmax = CI_U))

CIs + geom_point(color = "grey") + 
  geom_errorbarh(color = "grey") + 
  geom_vline(aes(xintercept = phi), color = "red")+
  facet_wrap(~panel)



