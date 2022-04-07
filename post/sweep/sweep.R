# Project:   knowledgeBase
# Objective: Understand the sweep operator
# Author:    Edoardo Costantini
# Created:   2021-11-19
# Modified:  2021-11-26
# Notion:    https://lavish-hollyhock-981.notion.site/Sweep-operator-a2147fb983594364a71b23398fbc0856
# Other ref:

# Packages and functions
source("sweepGoodnight.R") # main sweep function
library(ISR3)                  # for SWP functions
library(fastmatrix)            # alternative sweep

# ---------------------------------------------------------------------- #
# Little Rubin Dataset
# ---------------------------------------------------------------------- #

# Input data
  dat <- as.data.frame(
    matrix(
      data = c(7, 1, 11, 11, 7, 11, 3, 1, 2, 21, 1, 11, 10, 26,
               29, 56, 31, 52, 55, 71 ,31, 54, 47, 40, 66, 68,
               6, 15, 8, 8, 6, 9, 17, 22, 18, 4, 23, 9, 8,
               60, 52, 20, 47, 33, 22,6,44,22,26,34,12,12,
               78.5, 74.3, 104.3, 87.6, 95.9, 109.2, 102.7,
               72.5, 93.1, 115.9, 83.8, 113.3, 109.4),
      ncol = 5
    )
  )
  n <- nrow(dat)
  p <- ncol(dat)

# Create augmented crossproducts matrix T0 (matrix of sufficinet stats)
  # aka Uncorrected Sum of Squares and Crossproducts matrix (USSCP)
  dat_aug <- as.matrix(cbind(int = 1, dat))
  T0 <- crossprod(dat_aug)
  colSums(dat_aug)
  t(as.matrix(dat)) %*% as.matrix(dat)

# Create G-matrix from Little & Rubin p. 149
  G <- matrix(NA, ncol = p + 1, nrow = p + 1)
  G[, 1] <- colMeans(dat_aug)
  G[1, ] <- colMeans(dat_aug)
  G[-1,-1] <- crossprod(as.matrix(dat)) / n
  dimnames(G) <- dimnames(T0)

# Augmented covariance matrix from Little & Rubin p. 149 (eq 7.20)
  # Option 1: By sweeping G
  sweepGoodnight(G, 1)
  colMeans(dat) # means are there
  cov(dat) * (n - 1) / n # covariance matrix is there

  # Option 2: Make the augmented covariance matrix my self
  theta <- matrix(NA, ncol = p + 1, nrow = p + 1)
  theta[, 1] <- colMeans(dat_aug)
  theta[1, ] <- colMeans(dat_aug)
  theta[-1,-1] <- cov(dat) * (n-1) / n
  theta[1, 1] <- -1
  dimnames(theta) <- dimnames(T0)

# Fit some models
# Multivariate intercept only model
  # Define the dvs
  dvs <- c("V5", "V4")

  # Complicated but flexible way of writing the formula
  formula_lm <- paste0("cbind(",
                       paste0(dvs, collapse = ", "),
                       ")",
                       " ~ 1")

  # Fit the model with the MLM
  mlm0 <- lm(formula_lm, data = dat)

  # Do it with sweep shortcut
  coef(mlm0)
  sweepGoodnight(T0, 1)["int", dvs]
  sweepGoodnight(G, 1)["int", dvs]
  theta["int", dvs]

# Any multivariate model
# Sweep any variables as predictor of any variable
  # Define the dvs
  dvs <- c("V1", "V3", "V4")

  # Define the predictors
  preds <- c("V2", "V5")

  # Complicated but flexible way of writing the formula
  formula_lm <- paste0("cbind(",
                       paste0(dvs, collapse = ", "),
                       ")",
                       " ~ ",
                       preds = paste0(preds, collapse = " + "))

  # Fit the model with the MLM
  mlm0 <- lm(formula_lm, data = dat)
  coef(mlm0)

  # Do it with sweep shortcut
  sweep_over <- which(colnames(T0) %in% preds)
  sweepGoodnight(T0, c(1, sweep_over))[c("int", preds), dvs]
  sweepGoodnight(G, c(1, sweep_over))[c("int", preds), dvs]
  sweepGoodnight(theta, sweep_over)[c("int", preds), dvs]

  # Compare matrices
  T0 / n
  G
  ISR3::RSWP(theta, 1) # reverse sweep!

  ISR3::SWP(G, 1) # reverse sweep!
  theta

# ---------------------------------------------------------------------- #
# Draw data
# ---------------------------------------------------------------------- #

# Draw simple data
  n <- 50
  p <- 5
  X <- data.frame(MASS::mvrnorm(n, rep(0, p), diag(p)))
  X_aug <- cbind(int = 1, X)
  lm(X1 ~ ., X)

# Make crossproduct matrix
  T0 <- crossprod(as.matrix(X_aug))

# Create G-matrix as defined in Little Rubin p. 149
  G <- matrix(NA, ncol = p + 1, nrow = p + 1)
  G[, 1] <- colMeans(X_aug)
  G[1, ] <- colMeans(X_aug)
  G[-1,-1] <- crossprod(as.matrix(X)) / n
  dimnames(G) <- dimnames(T0)

# Augmented covariance matrix from Little & Rubin p. 149 (eq 7.20)
  theta <- matrix(NA, ncol = p + 1, nrow = p + 1)
  theta[, 1] <- colMeans(X_aug)
  theta[1, ] <- colMeans(X_aug)
  theta[-1,-1] <- cov(X) * (n-1) / n
  theta[1, 1] <- -1
  dimnames(theta) <- dimnames(T0)

# Multivaraite regressions
  # Define the dvs
  dvs <- c("X1", "X5")

  # Define the predictors
  preds <- c("X3", "X4")

  # LM
  # Complicated but flexible way of writing the formula
  formula_lm <- paste0("cbind(",
                       paste0(dvs, collapse = ", "),
                       ")",
                       " ~ ",
                       preds = paste0(preds, collapse = " + "))

  # Fit the model with the alternatives
  mlm0 <- lm(formula_lm, data = X)
  coef(mlm0)

  # Sweep T0
  sweep_over <- which(colnames(T0) %in% preds)
  sweepGoodnight(T0, c(1, sweep_over))[c("int", preds), dvs]

  # Sweep G
  sweepGoodnight(G, c(1, sweep_over))[c("int", preds), dvs]

  # Sweep theta (notice we do not sweep over 1)
  sweepGoodnight(theta, sweep_over)[c("int", preds), dvs]
