# Project:   blogdown
# Objective: Study Supervised PCR a la BairEtAl
# Author:    Edoardo Costantini
# Created:   2022-07-27
# Modified:  2022-07-27
# Notes:     https://tibshirani.su.domains/superpc/tutorial.html

# Load automatic package -------------------------------------------------------

library("superpc")

# Run example functions --------------------------------------------------------

set.seed(464)

x <- matrix(rnorm(1000*100), ncol = 100)
v1 <- svd(x[1:80,])$v[,1]
y <- 2 + 5 * v1 + .05 * rnorm(100)

xtest <- x
ytest <- 2 + 5 * v1 + .05 * rnorm(100)

featurenames <- paste("X", as.character(1:1000),sep="")

# Create train and test data
data.train <- list(x = x, y = y, featurenames = featurenames)
data.test <- list(x = xtest, y = ytest, featurenames = featurenames)

# train  the model. This step just computes the scores for each feature
train.obj <- superpc.train(data.train, type = "regression")

# cross-validate the model
cv.obj <- superpc.cv(train.obj, data.train)

#plot the cross-validation curves. From this plot we see that the 1st
# principal component is significant and the best threshold  is around 0.7

superpc.plotcv(cv.obj)

ls(cv.obj)

cv.obj$v.preval

fit.cts <- superpc.predict(train.obj,
                           data.train,
                           data.test,
                           threshold = 0.7,
                           n.components = 3,
                           prediction.type = "continuous")

# Manual approach --------------------------------------------------------------

x_m <- t(x)
  colnames(x_m) <- featurenames

x_m_test <- matrix(rnorm(25*1000), ncol = 1e3)
  colnames(x_m_test) <- featurenames

y_m <- y
theta = seq(0.01, .99, by = .01)
npcs <- 5
nfolds <- 10

  # Obtain R-squared for all simple linear regression models
  r2_vec <- apply(x_m, 2, function(j) {
    sqrt(summary(lm(y_m ~ j))$r.squared)
  })

  # DEfine predictor groups (pred groups) based on different theta
  pred_groups <- lapply(theta, function(m) {
    preds <- colnames(x_m)[r2_vec >= m]
    if (length(preds) >= 1) {
      preds
    } else {
      NULL
    }
  })

  # If theta used lead only to empty pred groups, say so
  if(all(sapply(pred_groups, is.null)) == TRUE){
    stop(
      paste0(
        "The threshold values used are too high. Try using a lower range."
      )
    )
  }

  # Drop empty pred_groups slots
  pred_groups <- pred_groups[!sapply(pred_groups, is.null)]

  # Drop possible duplicated pred_groups slots
  pred_groups <- unique(pred_groups)

  # Drop preds groups that are smaller than required npcs
  pred_groups <- pred_groups[sapply(pred_groups, length) >= npcs]

  # If there is no pred group with enough predictors for the required npcs, say so
  if(length(pred_groups) == 0){
    stop(
      paste0(
        "There is no threshold value that can select enough predictors to extract ",
        npcs, " PCs. Try using a smaller npcs or lower theta."
      )
    )
  }

  # Create a partition vector
  part <- sample(rep(1:nfolds, ceiling(nrow(x_m) / nfolds)))[1:nrow(x_m)]

  # Obtain Cross-validation error
  cve_obj <- lapply(pred_groups, function(set) {
    .spcrCVE(
      dv = y_m,
      pred = x_m[, set, drop = FALSE],
      K = nfolds,
      part = part,
      npcs = npcs
    )
  })

  # Extract CVEs
  cve <- sapply(cve_obj, "[[", 1)
  preds_active <- pred_groups[[which.min(cve)]]

  # Train PCR on dotxobs sample
  pcr_out <- pls::pcr(
    y_m ~ x_m[, preds_active, drop = FALSE],
    ncomp = npcs,
    scale = TRUE,
    center = TRUE,
    validation = "none"
  )

  # Define sigma
  RSS <- sqrt(sum(pcr_out$residuals^2))
  sigma <- RSS / (nrow(x_m) - npcs - 1)

  # Get prediction on (active) missing part
  yhat <- predict(
    object = pcr_out,
    newdata = x_m_test[, preds_active, drop = FALSE],
    ncomp = npcs,
    type = "response"
  )

  # Add noise for imputation uncertainty
  imputes <- yhat + rnorm(sum(wy)) * sigma

  .spcrCVE <- function(dv, pred, part, K = 10, npcs = 1) {
    # Input examples
    # dv   = mtcars[, 1]
    # pred = mtcars[, -1]
    # K    = 10
    # npcs = 5
    # part = sample(rep(1 : K, ceiling(nrow(mtcars) / K)))[1 : nrow(mtcars)]

    # Define a safe number of pcs
    q <- min(npcs, ncol(pred))

    # Create an empty storing object
    mse <- rep(NA, K)

    # Loop over K folds
    for (k in 1:K) {

      # Partition data:
      Xtr <- pred[part != k, , drop = FALSE]
      Xva <- pred[part == k, , drop = FALSE]
      ytr <- dv[part != k]
      yva <- dv[part == k]

      # Calibrate PCR on training datest
      pcr_out <- pls::pcr(
        ytr ~ Xtr,
        ncomp = q,
        scale = TRUE,
        center = TRUE,
        validation = "none"
      )

      # Get prediction on validation data set
      yva_hat <- predict(pcr_out, newdata = Xva, ncomp = q, type = "response")

      # Save MSE
      mse[k] <- MLmetrics::MSE(
        y_pred = yva_hat,
        y_true = yva
      )

    }

    # Return the CVE:
    cve <- sum(mse * (table(part) / length(part)))

    # Return
    return(list(
      cve = cve,
      npcs = q
    ))

  }