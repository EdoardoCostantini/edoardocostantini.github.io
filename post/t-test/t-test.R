# Project:   blogdown
# Objective: TODO
# Author:    Edoardo Costantini
# Created:   2022-09-06
# Modified:  2022-09-06


# Side of the test -------------------------------------------------------------

# Types of t-test --------------------------------------------------------------

# one sample

# two sample paired

# etc

# Group sizes ------------------------------------------------------------------

# Compare manual and automatic computation in R --------------------------------

###          ###
### Overview ###
###          ###

## You will practice inferential testing, some basic EDA techniques, missing
## data descriptives, and outlier analysis.

## You will need four datasets to complete the following tasks: "bfi_clean.rds",
## "tests.rds", "bfiOE.rds", and "airQual.rds". These datasets are saved in the
## "data" directory for this set of lab materials.


###                   ###
### Tasks / Questions ###
###                   ###

rm(list = ls(all = TRUE))

##--Preliminaries-------------------------------------------------------------##

## 1) If you have not already done so, use the "install.packages" function to
##    install the "mice" package.

install.packages("mice", repos = "http://cloud.r-project.org")

## 2) Use the "library" function to load the "mice" and "MASS" packages.

library(mice)
library(MASS)

## 3) Use the "paste0" function and the "readRDS" function to load the four
##    datasets into memory.

dataDir <- "../data/"

bfi1    <- readRDS(paste0(dataDir, "bfi_clean.rds"))
bfi2    <- readRDS(paste0(dataDir, "bfiOE.rds"))
tests   <- readRDS(paste0(dataDir, "tests.rds"))
airQual <- readRDS(paste0(dataDir, "airQual.rds"))


##--Testing-------------------------------------------------------------------##

### Use the "bfi_clean" data to complete the following analyses/tasks:

## 4a) Conduct a t-test to check for gender differences in mean levels of
##     "agree" (i.e., agreeableness). Do not assume equal variances.

out1.1 <- t.test(agree ~ gender, data = bfi1)
out1.1

## 4b) What is the value of the estimated mean difference in "agree"?

diff(out1.1$estimate)

### OR ###

diff(rev(out1.1$estimate))

## 4c) What is the value of the estimated t-statistic?

out1.1$statistic

## 4d) Is the estimated mean difference significant at the alpha = 0.05 level?

ifelse(out1.1$p.value < 0.05, "YES", "NO")

## 5a) Conduct a t-test to check for gender differences in mean levels of
##     "agree" (i.e., agreeableness). Assume equal variances.

out1.2 <- t.test(agree ~ gender, data = bfi1, var.equal = TRUE)
out1.2

## 5b) What is the value of the estimated mean "agree" for males?

out1.2$estimate[1]

## 5c) What is the value of the estimated t-statistic?

out1.2$statistic

## 5d) What is the 95% CI for the estimated mean difference?

out1.2$conf.int

## 5e) Is the t-statistic you found here different from the one you computed in
##     Q1? If so, why?

### YES. Although the estimated mean difference hasn't changed, the SE in Q1 was
### corrected to control for the unequal group variances, so the estimated
### t-statistic in Q1 was smaller.

### EXTRA: Compute everything manually
###        You can compute all of the values obtained with the t.test()
###        function manually.
###        This is extra content that might help you understand important
###        statistical comcepts. However, you do not need to know how to
###        this for the exam / quiz.
# First get all the information you need on the two samples (men and women)
mean_1 <- mean(bfi1$agree[bfi1$gender == "male"])   # mean of agree for men
s2_1   <- var(bfi1$agree[bfi1$gender == "male"])    # variance of agree for men
n_1    <- sum(bfi1$gender == "male")                # sample size of men

mean_2 <- mean(bfi1$agree[bfi1$gender == "female"]) # mean of agree for women
s2_2   <- var(bfi1$agree[bfi1$gender == "female"])  # variance of agree for women
n_2    <- sum(bfi1$gender == "female")              # sample size of women

est_h0 <- 0 # value of the estimate according to the null hypothesis

# Second, compute the estiamte of interest (difference between group means)
est <- mean_1 - mean_2

# Third, compute the pooled standard deviation
pool_sd <- sqrt( ( (n_1-1) * s2_1 + (n_2-1) * s2_2)/(n_1 + n_2 - 2))

# Fourth, compute the variability measure (standard error of the estimate)
variability <- pool_sd * sqrt(1/n_1 + 1/n_2)

# Fifth, compute the t-statistic
t_statistic <- (est - est_h0)/variability

# Sixth, compute the degrees of freedom
df <- n_1 + n_2 - 2

# Seventh, obtain the pvalue by checking the probability of observing the
# t_statistic you obsereved or something more extreme on a t-distribution
# with degrees of freedom df
pvalue <- pt(q = t_statistic, df = df) * 2 # * 2 beacuse it's a two sided test!

# Check the value is the same as the one you computed with the t.test() function
c(fast = out1.2$p.value, manual = pvalue)

## 6a) Test for a significant Pearson correlation between "agree" and "neuro"
##    (i.e., neuroticism).

out1.3 <- with(bfi1, cor.test(agree, neuro))
out1.3

## 6b) What is the value of the correlation coefficient?

out1.3$estimate

## 6c) Is this correlation significantly different from zero at the alpha = 0.05
##     level?

ifelse(out1.3$p.value < 0.05, "YES", "NO")

## EXTRA: Compute everything manually
r <- cor(bfi1$agree, bfi1$neuro)            # correlation coefficient
n <- nrow(bfi1)                             # sample size
t_statistic <- r * sqrt((n-2)/(1-r^2))      # t-statistic
df <- n - 2                                 # degrees of freedom
pvalue <- pt(q = t_statistic, df = df) * 2  # * 2 beacuse it's a two sided test!
c(fast = out1.3$p.value, manual = pvalue)   # check the results are the same

## 7a) Test the hypothesis that the correlation between "consc"
##     (i.e., conscientiousness) and "neuro" is less than zero.

out1.4 <- with(bfi1, cor.test(consc, neuro, alternative = "less"))
out1.4

## 7b) What is the value of the estimated correlation coefficient?

out1.4$estimate

## 7c) What is the 95% CI for the estimated correlation coefficient?

out1.4$conf.int

## 7d) Is this correlation significantly less than zero at the alpha = 0.05
##     level?

ifelse(out1.4$p.value < 0.05, "YES", "NO")

## EXTRA: Compute everything manually
r <- cor(bfi1$consc, bfi1$neuro)            # correlation coefficient
n <- nrow(bfi1)                             # sample size
t_statistic <- r * sqrt((n-2)/(1-r^2))      # t statistic
df <- n - 2                                 # degrees of freedom
pvalue <- pt(q = t_statistic, df = df)      # no * 2 this time! It's one sided.
c(fast = out1.4$p.value, manual = pvalue)   # check the results are the same
