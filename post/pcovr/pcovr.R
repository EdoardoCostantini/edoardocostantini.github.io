# pcover estimation

# Fit regular pcover with every parameter fixed

library("PCovR")

# Load data --------------------------------------------------------------------
data(alexithymia)

colMeans(alexithymia$X)
colMeans(alexithymia$Y)

apply(alexithymia$X, 2, var)
apply(alexithymia$Y, 2, var)

# Subset data
X_raw <- alexithymia$X
y_raw <- alexithymia$Y[, 1, drop = FALSE]

n <- nrow(X_raw)
p <- ncol(X_raw)

# Scale data
X <- scale(X_raw)# * (n - 1) / n
y <- scale(y_raw)# * (n - 1) / n

# Define parameters
alpha <- .5
npcs <- 5

# Estimation -------------------------------------------------------------------

    # Estiamte with PCovR function
    out <- pcovr_est(
        X = X,
        Y = y,
        a = alpha,
        r = npcs # fixed number of components
    )

    # Estiamte manually Vervolet version
    Hx <- X %*% solve(S) %*% t(X)
    G_vv <- alpha * X %*% t(X) / sum(X^2) + (1 - alpha) * Hx %*% y %*% t(y) %*% Hx / sum(y^2)
    EG_vv <- eigen(G_vv) # eigen-decomposition of matrix
    T_vv <- EG_vv$vectors[, 1:npcs]

    # Estiamte manually de Jong version
    beta <- alpha * sum(y)^2 / (alpha * sum(y)^2 + (1 - alpha) * sum(X)^2)
    G_jv <- beta * X %*% t(X) + (1 - beta) * Y_hat %*% t(Y_hat)
    EG_jv <- eigen(G_jv) # eigen-decomposition of matrix
    T_jv <- EG_jv$vectors[, 1:npcs]

# Compare results --------------------------------------------------------------

    # T scores
    Ts <- list(
        PCovR = head(out$Te),
        Vervolet = head(T_vv),
        PCovR_man = head(X %*% out$W),
        deJong = head(T_jv)
    )

    # Weights
    W <- list(
        PCovR = out$W,
        Vervolet = solve(t(X) %*% X) %*% t(X) %*% T_vv,
        deJong = solve(t(X) %*% X) %*% t(X) %*% T_jv
    )

    # Px
    t(out$Te) %*% X
    t(out$W) %*% t(X) %*% X
    out$Px

    # Py
    cbind(
        Py = out$Py,
        TtY = t(out$Te) %*% y,
        WtXtY = t(out$W) %*% t(X) %*% y
    )

    # B
    cbind(
        B = out$B,
        WPY = out$W %*% out$Py,
        WWtXtY = out$W %*% t(out$W) %*% t(X) %*% y
    )

# Maximum likelihood tuning of alpha -------------------------------------------

    # Fit PCovR
    out <- pcovr_out <- pcovr(
        X = X_raw,
        Y = y_raw,
        rot = "none",
        R = npcs, # fixed number of components
        modsel = "seq" # fastest option
    )

    # Compute error ratio with function
    err <- ErrorRatio(
        X = X,
        Y = y,
        Rmin = npcs,
        Rmax = npcs
    )

    # Compute error ratio components
    lm_mod <- lm(y ~ -1 + X)
    ery <- 1 - summary(lm_mod)$r.squared

    Rmin <- 1
    sing <- svd(X)
    vec <- Rmin:Rmax
    vec <- c(vec[1] - 1, vec, vec[length(vec)] + 1)
    VAF <- c(0, cumsum(sing$d^2) / sum(sing$d^2))
    VAF <- VAF[vec + 1]
    scr <- array(NA, c(1, length(vec)))
    for (u in 2:(length(vec) - 1)) {
        scr[, u] <- (VAF[u] - VAF[u - 1]) / (VAF[u + 1] - VAF[u])
    }
    erx <- 1 - VAF[which.max(scr)]

    # Find alpha ML
    alpha_ML <- sum(X^2) / (sum(X^2) + sum(Y^2) * erx / ery)

    # Compare to one found by package
    out$a - alpha_ML