# Fit a linear model -----------------------------------------------------------

    lm_fit <- lm(mpg ~ cyl + hp + wt, data = mtcars)

# Compute the residual standard error manually ---------------------------------

    # Define elements of the formula
    n <- nrow(mtcars) # sample size
    k <- 3 # number of parameters (regression coefficients)
    yhat <- fitted(lm_fit) # fitted y values
    y <- mtcars$mpg

    # Compute rse
    rse <- sqrt(sum((y - yhat)^2) / (n - k - 1))

    # Print rse
    rse

# residual standard error from lm output ---------------------------------------

    # Use the sigma function to extract it from an lm object
    sigma(lm_fit)

    # Compare with the manual computation
    sigma(lm_fit) - rse

# Check other computation of degrees of freedom

    lm_fit$df.residual
    summary(lm_fit)$df[2]
    df.residual(lm_fit)
    n - k - 1

# Center the dependent variable

    y <- mtcars[, 1] # mpg
    X <- as.matrix(mtcars[, -1])
    n <- nrow(X)
    p <- ncol(X)

    # Regular LM fit with intercept
    lm.fit <- lm(y ~ X)
    lm.fit$df.residual
    n - p - 1

    # Centered Y
    lm.fit.sc <- lm(scale(y, center = TRUE, scale = FALSE) ~ - 1 + X)
    lm.fit.sc$df.residual
    n - p

    # Compare gits
    summary(lm.fit)$r.squared
    summary(lm.fit.sc)$r.squared

    # Compute rse
    rse <- sqrt(sum((y - mean(y) - fitted(lm.fit.sc) )^2) / (n - p))
    rse

    # Use the sigma function to extract it from an lm object
    sigma(lm.fit.sc)
