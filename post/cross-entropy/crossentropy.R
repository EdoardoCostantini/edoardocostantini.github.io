# Project:   blogdown
# Objective: small script describing the key aspects of cross entropy
# Author:    Edoardo Costantini
# Created:   2022-04-22
# Modified:  2022-04-22

library(nnet)

# Data -------------------------------------------------------------------------

dat_bin <- ToothGrowth
dat_cat <- iris

# Binary cross entropy ---------------------------------------------------------

# Fit model
glm_log <- glm(supp ~ .,
               family = binomial(link = 'logit'),
               data = ToothGrowth)

# Predictions
preds_log <- predict(glm_log, type = "response")

# Compute entropy
p <- as.numeric(ToothGrowth$supp)-1 #
phat <- preds_log

x <- 0
for (i in 1:length(p)){
  x <- x + (p[i] * log(phat[i]))
}

- x

# Multi-categorical cross entropy ----------------------------------------------

# Fit model
glm_mln <- multinom(Species ~ ., data = iris)

# Predictions
preds_mln <- predict(glm_mln, type = "probs")

# Compute entropy
p <- FactoMineR::tab.disjonctif(iris$Species)
phat <- preds_mln

x <- 0
for (i in 1:nrow(p)){
  for (j in 1:ncol(p)){
    x <- x + (p[i, ] * log(phat[i, ]))
  }
}

- sum(x)

# Fast way to compute both -----------------------------------------------------

ce <- - sum(
  diag(
    p %*% t(log(phat))
  )
)