# Project:   blogdown
# Objective: small script describing the key aspects of cross entropy
# Author:    Edoardo Costantini
# Created:   2022-04-22
# Modified:  2022-05-12

# Prepare environment ----------------------------------------------------------

# Packages
library(nnet)
library(MLmetrics)  # for LogLoss() function
library(FactoMineR) # for tab.disjonctif() function

# Default rounding for this sessino
options("digits" = 5)

# Fit mulinomial logistic model ------------------------------------------------

# Fit model
glm_mln <- multinom(Species ~ Sepal.Length, data = iris)

# Obtain p and p_har -----------------------------------------------------------

# store true labels in a matrix p
p <- FactoMineR::tab.disjonctif(iris$Species)

# check it
head(p)

# obtain predictions
p_hat <- predict(glm_mln, type = "probs")

# check it
head(p_hat)

# Compute CE with a loop -------------------------------------------------------

# Define parameters
N <- nrow(iris) # sample size
K <- nlevels(iris$Species) # number of classes

# Create storing object for CE
CE <- 0

# Compute CE with a loop
for (i in 1:N){
  for (k in 1:K){
    CE <- CE - p[i, k] * log(p_hat[i, k])
  }
}

# Print the value of CE
CE

# Compute CE using the matrices directly ---------------------------------------
ce <- -sum(diag(p %*% t(log(p_hat))))

# Print the value of ce
ce

# Binary cross entropy ---------------------------------------------------------

# Fit model
glm_log <- glm(am ~ hp + wt,
               family = binomial(link = 'logit'),
               data = mtcars)

# store true labels in a matrix p
p <- FactoMineR::tab.disjonctif(mtcars$am)

# obtain predicted probabilites in matrix form
pred_probs <- predict(glm_log, type = "response")
p_hat <- cbind(k_0 = 1 - pred_probs,
               k_1 = pred_probs)

# check the first few rows of p
head(p)

# check the first few rows of p_hat
head(p_hat)

# Compute CE using the matrices directly
ce <- -sum(diag(p %*% t(log(p_hat))))

# Print the value of ce
ce

# Express as average
ce / nrow(mtcars)

# Compute binary CE with MLmetrics implementation ------------------------------

# Obtain vector of true probabilities
p_vec <- mtcars$am

# Obtain vector of predicted probabilities
p_hat_vec <- predict(glm_log, type = "response")

# Compute and print binary CE with MLmetrics implementation
MLmetrics::LogLoss(y_pred = p_hat_vec,
                   y_true = p_vec)