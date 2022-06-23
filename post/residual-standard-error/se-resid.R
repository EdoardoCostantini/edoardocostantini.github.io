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

# Check other things

    df.residual(lm_fit)
    (n - k - 1)

# PLS degrees of freedom

library(plsdof)

data(Boston)
X <- as.matrix(Boston[, -14])
y <- as.vector(Boston[, 14])

dim()

my.pls1 <- pls.model(X, y, m = 5, compute.DoF = TRUE)
my.pls1

my.pls1$DoF