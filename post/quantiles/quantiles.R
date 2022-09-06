# Project:   blogdown
# Objective: Describe quantiles theory and R function
# Author:    Edoardo Costantini
# Created:   2022-09-06
# Modified:  2022-09-06

# Set a seed

set.seed(20220906)

# Sample from a normal distribution

age <- rnorm(1e5, mean = 27, sd = 2)

# Define the 1st and 2nd quartile

quartiles <- quantile(age, probs = c(.25, .75))

# Plot density distribution

plot(density(age),
     main = "Quartiles for the probability distribution of age",
     xlab = NA)

# Costum x ticks
axis(side = 1, at = c(27, round(quartiles, 1)), labels = TRUE)

# Add points for quantiles
points(c(quartiles, median(age)),
        y = rep(0, length(quartiles)+1))
points(c(quartiles, median(age)),
       y = dnorm(c(quartiles, median(age)), mean = mean(age), sd = sd(age)))

# Add segments to devide plot
segments(x0 = quartiles[1], y0 = 0,
         x1 = quartiles[1], y1 = dnorm(quartiles[1],
                                       mean = mean(age), sd = sd(age)))
segments(x0 = median(age), y0 = 0,
         x1 = median(age), y1 = max(dnorm(age,
                                       mean = mean(age), sd = sd(age))))
segments(x0 = quartiles[2], y0 = 0,
         x1 = quartiles[2], y1 = dnorm(quartiles[2],
                                       mean = mean(age), sd = sd(age)))

# Add quartile labels
text(x = quartiles[1],
     y = -.005,
     "1st quartile")
text(x = quartiles[2],
     y = -.005,
     "3rd quartile")

# Add percentage under the curve labels
text(x = c(24, 30, 26.3, 27.7),
     y = c(.03, .03, .06, .06),
     "25 %")
