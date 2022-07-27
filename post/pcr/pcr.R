# Project:   blogdown
# Objective: Principal component regression
# Author:    Edoardo Costantini
# Created:   2022-07-01
# Modified:  2022-07-01

# Prepare environment ----------------------------------------------------------

    # Load packages
    library(pls)

# Preapre data -----------------------------------------------------------------

    # mtcars
    X <- as.matrix(mtcars[, -1])
    y <- mtcars[, 1]
    n <- nrow(X) # number of observations
    p <- ncol(X) # number of variables
    q <- 7

    # Process data -------------------------------------------------------------
    # X <- scale(X, center = TRUE, scale = TRUE)
    # y <- scale(y, center = FALSE, scale = FALSE)

    # Devide training and test data
    tr <- sample(x = 1:n, size = floor(2/3*n))
    Xtr <- scale(X[tr, ], center = TRUE, scale = TRUE)
    ytr <- scale(y[tr], center = TRUE, scale = FALSE)
    Xte <- X[-tr, ]
    yte <- y[-tr]

    # - Intercept checks -------------------------------------------------------

    lm_raw <- lm(y ~ ., data.frame(y = y[tr], Xtr))
    lm_cen <- lm(y ~ ., data.frame(y = ytr, Xtr))

    # Fitted values
    data.frame(
        raw = fitted(lm_raw),
        centered = fitted(lm_cen) + mean(y[tr]),
        diff = round(fitted(lm_raw) - (fitted(lm_cen) + mean(y[tr])), 3)
    )

    # Predict new values
    data.frame(
        raw = predict(lm_raw, newdata = data.frame(Xte)),
        centered = predict(lm_cen, newdata = data.frame(Xte)) + mean(y[tr]),
        diff = round(
            predict(lm_raw, newdata = data.frame(Xte)) -
                (predict(lm_cen, newdata = data.frame(Xte)) + mean(y[tr])),
            3
        )
    )

    # - Scaling checks ---------------------------------------------------------
    
    Xte_sc <- sapply(
        1:p,
        function(j) {
            (X[-tr, j] - mean(X[tr, j])) / sd(X[tr, j])
        }
    )
    colnames(Xte_sc) <- colnames(X)

    # Fit models on scaled and not scaled data
    lm_raw <- lm(y ~ ., data.frame(y = y[tr], X[tr, ]))
    lm_sca <- lm(y ~ ., data.frame(y = y[tr], Xtr))

    # Get predictions on scalde and not scaled data
    predict(lm_raw, newdata = data.frame(Xte))
    predict(lm_sca, newdata = data.frame(Xte_sc))

# Fit PCR models ---------------------------------------------------------------

    # - Fit with pls::pcr ------------------------------------------------------

    pcr_out <- pls::pcr(
        ytr ~ Xtr,
        ncomp = q,
        scale = FALSE,
        center = FALSE,
        validation = "none"
    )

    # - Fit PCR manually -------------------------------------------------------

    uSVt    <- svd(Xtr)
    u       <- uSVt$u[, 1:q]                               # orthonormal basis in the column space of X
    Sigma   <- diag(uSVt$d)[1:q, 1:q]                      # diagonal scaling matrix
    V       <- uSVt$v[, 1:q]                               # loadings
    Ts      <- Xtr %*% V                                   # component scores
    Ts_te   <- Xte %*% V                                   # component scores new data
    Bhat    <- V %*% solve(t(Ts) %*% Ts) %*% t(Ts) %*% ytr # regression coefficients
    Bhat_v2 <- V %*% solve(Sigma) %*% t(u) %*% ytr

# Check equivalence ------------------------------------------------------------

    ls(pcr_out)
    
    # - Scores -----------------------------------------------------------------

    round(pcr_out$scores - Ts, 5)

    # - Loadings ---------------------------------------------------------------

    pcr_out$projection - V

    # - B regression coefficients ----------------------------------------------

    data.frame(
        pls.pack = pcr_out$coefficients[, , q],
        manual1 = Bhat,
        manual2 = Bhat_v2,
        diff = round(
            pcr_out$coefficients[, , q] - Bhat, 5
        )
    )

    # - Fitted values ----------------------------------------------------------

    data.frame(
        pls.pack = pcr_out$fitted.values[, , q],
        manual.Ts = Ts %*% coef(lm(ytr ~ -1 + Ts)),
        manual.Xtr = Xtr %*% Bhat,
        diff = round(
            pcr_out$fitted.values[, , q] - Xtr %*% Bhat, 5
        )
    )

    # - Residuals --------------------------------------------------------------

    data.frame(
        pls.pack = pcr_out$residuals[, , q],
        manual.Ts = ytr - Ts %*% coef(lm(ytr ~ -1 + Ts)),
        manual.Xtr = ytr - Xtr %*% Bhat,
        diff = round(
            pcr_out$residuals[, , q] - (ytr - Xtr %*% Bhat), 5
        )
    )

    # - Project new data -------------------------------------------------------

    round(
        predict(pcr_out, newdata = Xte_sc, ncomp = 1:q, type = "scores") - Xte_sc %*% V,
        5
    )

    # - Prediction of new data -------------------------------------------------

    data.frame(
        pls.pack = predict(pcr_out, newdata = Xte, ncomp = q, type = "response"),
        manual.Xte = Xte %*% Bhat,
        manual.Ts_te = Ts_te %*% coef(lm(ytr ~ -1 + Ts)),
        diff = round(
            drop(predict(pcr_out, newdata = Xte, ncomp = q, type = "response")) - Xte %*% Bhat,
            5
        )
    )

    # - Prediction of new data on correct scale --------------------------------

    # Fit model on uncetered y as input
    pcr_out_rawy <- pls::pcr(
        y[tr] ~ X[tr, ],
        scale = TRUE,
        center = TRUE,
        ncomp = q,
        validation = "none"
    )

    # Predict new data with function
    pred_pls <- predict(pcr_out_rawy, newdata = X[-tr, ], ncomp = q, type = "response")

    # Predict on new data manual
    B0 <- drop(mean(y[tr]) - colMeans(Xtr) %*% Bhat)
    pred_manscale <- B0 + Xte_sc %*% Bhat

    # Comparison
    data.frame(
        pls.pack = pred_pls,
        manual.1step = Xte_sc %*% Bhat + mean(y[tr]),
        manual.2step = pred_manscale,
        diff = round(drop(pred_pls) - drop(pred_manscale), 5)
    )