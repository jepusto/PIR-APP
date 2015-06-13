rm(list = ls())
library(ggplot2)
library(plyr)
load("Bootstrap Check.RData")

zeta <- subset(results, parm == "zeta")

zeta <- zeta[order(zeta$phi, zeta$est),]
zeta$case <- 1:200

true_zeta <- ddply(results, .(phi), summarize, phi = phi[1], zeta = zeta[1])

CIs <- ggplot(data = zeta, aes(x = est, y = case, xmin = CI_L, xmax = CI_U))

CIs + geom_point(color = "grey") + 
      geom_errorbarh(color = "grey") + 
      geom_vline(aes(xintercept = zeta), color = "red")+
      facet_wrap(~phi)
