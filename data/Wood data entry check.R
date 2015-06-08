rm(list = ls())
Danny <- read.csv("data/Wood et al data.csv", na.strings = ".")
Chris <- read.csv("data/Runyon Wood Data 2 (6_5_2015).csv", na.string = ".")

levels(Chris$Measurement.Type) <- c("PIR", "WIR")

all.equal(Danny, Chris)

#The first error for "Child" is fine, this wasn't numbered and therefore we had to provide our own, disagreement is not an issue

#The second error is a disagreement about the type of missing that we aren't differentiating between anyway.