# Quick study of mice ampute

    library(mice)
    source("./content/post/ampute/simMissingness.R")

# Generate data to ampute ------------------------------------------------------

    p <- 4
    S <- matrix(.4, nrow = p, ncol = p)
    diag(S) <- 1
    X <- MASS::mvrnorm(
        n = 1e3, mu = rep(10, p),
        Sigma = S
    )

# Proportion of missingness

    # Ampute with ampute
    X_MA <- ampute(data = X, prop = .5, mech = "MCAR")

    # Ampute with stepwise univariate amputation (SUA)
    X_SUA <- apply(X, 2, function(j){
        j[sample(1:nrow(X), nrow(X)*.5)] <- NA
        j
    })

    # Variable wise prop of missingness
    colMeans(is.na(X_MA$amp))
    colMeans(is.na(X_SUA))

    # complete cases with missingness
    mean(rowSums(is.na(X_MA$amp)) == 0)
    mean(rowSums(is.na(X_SUA)) == 0)

# Can this be MAR? -------------------------------------------------------------

    X_MA <- ampute(data = X, prop = .5, mech = "MAR")
    X_MA$amp
    X_MA$patterns
    X_MA$weights

# Measures of effect size in MAR mechanisms ------------------------------------

    # Miss pats
    mp <- matrix(0, 1, p)
    mp[, 3:4] <- 1

    # Weight
    wt <- matrix(0, 1, p)
    wt[, 3:4] <- 1

    # Ampute with ampute
    X_MA <- ampute(
        data = X, 
        prop = .5, 
        mech = "MAR",
        patterns = mp,
        weights = wt
    )
    X_MA$amp
    X_MA$patterns
    X_MA$weights

    # Ampute with simMissingnes()
    nas <- simMissingness(
        pm = .3,
        data = X,
        preds = c(3, 4),
        type = "high",
        beta = NULL
    )

    # Make a data.frame for prediction
    dat <- data.frame(
        X_MA = is.na(X_MA$amp[, 1]),
        X_SUA = nas,
        X = X[, 3:4]
    )
    
    # Fit logistic models
    glm_MA <- glm(X_MA ~ X.1 + X.2, data = dat, family = "binomial")
    glm_SUA <- glm(X_SUA ~ X.1 + X.2, data = dat, family = "binomial")
    yhat_MA <- predict(glm_MA, dat, type = "response")
    yhat_SUA <- predict(glm_SUA, dat, type = "response")

    # R-square
    c(
        MA = (1 - glm_MA$deviance / glm_MA$null.deviance) * 100,
        sua = (1 - glm_SUA$deviance / glm_SUA$null.deviance) * 100
    )
    
    # AUC
    c(
        MA = pROC::auc(dat$X_MA, yhat_MA),
        sua = pROC::auc(dat$X_SUA, yhat_SUA)
    )