#!/bin/Rscript

require("ggplot2")

# svg("mygraphic.svg")

data <- read.csv("da-points.csv")
qplot(data$x, data$y, color=data$color) + geom_abline(intercept = 180, slope = 1) + geom_abline(intercept = 0, slope = 1) + geom_abline(intercept = -180, slope = 1)

# dev.off()
