# clear workspace
rm(list = ls())

# import necessary libraries
library(latex2exp)

### We generate the smoothed survival function plot for Binomial distribution with n=4 and p=0.5
# Number of trials
n <- 4

# Probability of success in each trial
p <- 0.5

# Values for which to calculate probabilities
x <- 0:n

# Calculate survival probabilities (1 - CDF)
survival_values <- 1 - pbinom(x, size = n, prob = p)

# Create a step function
survival_step <- stepfun(x, c(1, survival_values), f = 1, right = TRUE)

# output the plot with dpi 300 and jpeg format
jpeg("./output/smoothed.survival.png", width = 8, height = 6, units = "in", res = 300)
par(oma = c(0, 0, 0, 0))

# Plot the survival function
plot(survival_step,
      main = "",
      xlab = expression(italic(v)), ylab = "Survival Probability", verticals = FALSE, xaxt = "n",cex.lab=1.4
)

lines(c(-1, x), c(1, survival_values), lty = 2)

# add a line for y=0.5
abline(h = 0.5, lty = 3)

# add a point for the median with triangle pointing down
points(2, 1 - pbinom(1, size = n, prob = p), pch = 2, cex = 2)
points(1.5, 0.5, pch = 0, cex = 2)

# add a line for x=2 and x=1.5 from the x-axis to the curve
segments(2, -0.1, 2, 1 - pbinom(1, size = n, prob = p), lty = 3)
segments(1.5, -0.1, 1.5, 0.5, lty = 3)

axis(side = 1, at = c(1.5, 2), labels = c(TeX("$Q_{0.5}^{{m}}$"),TeX("$Q_{0.5}$")),cex.axis=1.3)
axis(side = 2, at = c(0.5), labels = c(0.5))

dev.off()
