# Project:   blogdown
# Objective: Describe how boxplots work
# Author:    Edoardo Costantini
# Created:   2022-09-06
# Modified:  2022-09-06

# Generate some age variable for a university master programme
set.seed(20220906)
age <- round(rnorm(1e3, mean = 26, sd = 2), 0)

# Add some extreme outliers
age <- c(15, 20, 32, 38, age)

# Look at the boxplot
boxplot(age)

# Compute boxplot statistics manually ------------------------------------------

# Compute the median
med <- median(age)

# Compute 1st and 3rd quartiles
qnt <- quantile(age, probs = c(.25, .75))

# Compute Interquartile range
IQR <- diff(qnt)[[1]]

# Compute fences bounds
C <- 1.5 # range multiplier
fences <- c(lwr = qnt[[1]] - C * IQR, upr = qnt[[2]] + C * IQR)

# Put together the boxplot stats
bxstats <- sort(c(med = med, qnt, f = fences))

# Define names for boxplot statistics
bxstats_names <- c(ifelse(C == 1.5,
                          "inner fence \n lower bound",
                          "outer fence \n lower bound"),
                   "1st quartile",
                   "median",
                   "3rd quartile",
                   ifelse(C == 1.5,
                          "inner fence \n upper bound",
                          "outer fence \n upper bound"))
names(bxstats) <- bxstats_names

# Compute boxplot statistics with R function -----------------------------------

bxstats_auto <- boxplot.stats(age, coef = C)$stats

# Compare automatic and manual apporches

data.frame(manual = bxstats, auto = bxstats_auto)

# Visualize --------------------------------------------------------------------

# Allow two plots one next to the other
par(mfrow = c(1, 2))

# Plot C = 1.5
C <- 1.5
boxplot(age, range = C, main = paste0("Boxplot of age (C = ", C, ")"))

# Add arrows pointings to statistics
arrows(x0 = .565, y0 = boxplot.stats(age, coef = C)$stats,
       x1 = c(.875, rep(.765, 3), .875), y1 = boxplot.stats(age, coef = C)$stats)
# Add labels of statistics
text(x = rep(.55, 5),
     y = boxplot.stats(age, coef = C)$stats,
     labels = c("inner fence \n lower bound",
                "1st quartile",
                "median",
                "3rd quartile",
                "inner fence \n upper bound"),
     adj = 1)

# Plot C = 3
C <- 3
boxplot(age, range = C, main = paste0("Boxplot of age (C = ", C, ")"))

# Add arrows pointings to statistics
arrows(x0 = .565, y0 = boxplot.stats(age, coef = C)$stats,
       x1 = c(.875, rep(.765, 3), .875), y1 = boxplot.stats(age, coef = C)$stats)
# Add labels of statistics
text(x = rep(.55, 5),
     y = boxplot.stats(age, coef = C)$stats,
     labels = c("outer fence \n lower bound",
                "1st quartile",
                "median",
                "3rd quartile",
                "outer fence \n upper bound"),
     adj = 1)
