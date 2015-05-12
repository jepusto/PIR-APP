Danny <- read.csv("data/Wood et al data.csv", na.strings = ".")
Chris <- read.csv("data/Runyon Wood Data.csv", na.string = ".")

levels(Chris$Measurement.Type) <- c("PIR", "WIR")

all.equal(Danny, Chris)
