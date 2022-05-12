# Project:   blogdown
# Objective: Initial notes on EM for factor analysis
# Author:    Edoardo Costantini
# Created:   2022-05-10
# Modified:  2022-05-10

rm(list = ls())

# Packages ----------------------------------------------------------------

library(psych)
library(fastmatrix)

sweepGoodnight <- function (A, target){

  for(k in target){
    # Step 1: Let D = a_kk
    D <- A[k, k]

    # Step 2: Divide row k by D.
    A[k, ] <- A[k, ] / D

    # Step 3:
    # - For every other row i != k, let B = a_ik
    # - Subtract B \times row k from row i.
    # - set a_ik = -B/D.
    for(i in 1:nrow(A)){
      if(i != k){
        B <- A[i, k]
        A[i, ] <- A[i, ] - B * A[k, ]
        A[i, k] <- -1 * B / D
      }
    }
    # Step 4: Set a_kk = 1/D
    A[k, k] = 1/D
  }

  # Output
  return(A)
}

# Set up -----------------------------------------------------------------------

  # Baseline data
  v1 <- c(1,1,1,1,1,1,1,1,1,1,3,3,3,3,3,4,5,6)
  v2 <- c(1,2,1,1,1,1,2,1,2,1,3,4,3,3,3,4,6,5)
  v3 <- c(3,3,3,3,3,1,1,1,1,1,1,1,1,1,1,5,4,6)
  v4 <- c(3,3,4,3,3,1,1,2,1,1,1,1,2,1,1,5,6,4)
  v5 <- c(1,1,1,1,1,3,3,3,3,3,1,1,1,1,1,6,4,5)
  v6 <- c(1,1,1,2,1,3,3,3,4,3,1,1,1,2,1,6,5,4)
  Y <- cbind(v1,v2,v3,v4,v5,v6)
  n <- nrow(Y)
  p <- ncol(Y)
  q <- 3
  Ybar <- colMeans(Y)

# Psych version ----------------------------------------------------------------

  fa_out <- psych::fa(Y,
                      nfactors = 3,
                      covar = TRUE, # we want to work with the covariance matrix
                      rotate = "varimax",
                      fm = "ml")

# McNicholas -------------------------------------------------------------------

  # Constants
  Sxb <- cor(Y) # if we wanted to work with cor as default in psych::fa
  Sxb <- 1/(n-1) * t(Y - Ybar) %*% (Y - Ybar) # if data is centered, Cyy which is Sxb

  # Starting values
  Lambda <- matrix(.5, nrow = p, ncol = q)
  Psi <- diag(1, p)

  # EM
  for (i in 1:1e3){
    # E step
    Beta <- t(Lambda) %*% solve(Lambda %*% t(Lambda) + Psi)
    Theta <- diag(q) - Beta %*% Lambda + Beta %*% Sxb %*% t(Beta)

    # M step
    Lambda <- Sxb %*% t(Beta) %*% solve(Theta)
    Psi <- diag(diag(Sxb - Lambda %*% Beta %*% Sxb))
  }

  # Loadings
  varimax(Lambda)$loadings[]
  loadings(fa_out)[, c(2,3,1)]

  # Uniqueness Psi
  diag(Psi)
  fa_out$uniquenesses

  # Sigma = LL' + Psi
  Lambda %*% t(Lambda) + Psi
  loadings(fa_out)[, c(2,3,1)] %*% t(loadings(fa_out)[, c(2,3,1)]) + diag(fa_out$uniquenesses)

# RubinThayer1982 --------------------------------------------------------------

  # Constants
  a <- Ybar <- colMeans(Y)
  R <- diag(q)
  Ys <- scale(Y)
  Cyy <- 1/(n-1) * t(Ys) %*% (Ys) # if we want to work with correlation matrix
  Cyy <- 1/(n-1) * t(Y - Ybar) %*% (Y - Ybar) # if data is centered, Cyy which is Sxb

  # Initialize
  B <- B0 <- matrix(.5, nrow = q, ncol = p)
  Tau <- Tau0 <- diag(1, p)

  # EM
  for (i in 1:1e3){
    # E step estiamtes ---------------------------------------------------------
    # gamma = Beta
    gamma <- solve(Tau + t(B) %*% R %*% B) %*% t(B) %*% R
    # Delta = Theta - Beta %*% Sxb %*% t(Beta)
    Delta <- R - (R %*% B) %*% solve(Tau + t(B) %*% R %*% B) %*% t(B) %*% R

    # Cyz_hat
    Cyz_hat <- Cyy %*% gamma
    Czy_hat <- t(gamma) %*% Cyy

    # Czz_hat = Theta
    Czz_hat <- Delta + t(gamma) %*% Cyy %*% gamma

    # M step -------------------------------------------------------------------
    # B = Lambda
    B <- solve(Czz_hat) %*% t(Cyz_hat)
    # Taus = Psi
    taus <- diag(Cyy - Cyy %*% gamma %*% B)
    Tau <- diag(taus)
  }

  # Loadings
  lapply(list(
    RFA = loadings(fa_out)[, c(2,3,1)],
    MN2016 = varimax(Lambda)$loadings[],
    RT1982 = varimax(t(B))$loadings[, c(2, 1, 3)]
  ), round, 3)

  # Uniqueness
  sapply(list(
    RFA = fa_out$uniquenesses,
    MN2016 = diag(Psi),
    RT1982 = diag(Tau)
  ), round, 3)

# LiuRubin1998 -----------------------------------------------------------------

  # Constants
  a <- Ybar <- colMeans(Y)
  R <- diag(q)
  S_i <- 1:q
  Cyy <- 1/(n-1) * t(Y - Ybar) %*% (Y - Ybar) # if data is centered, Cyy which is Sxb
  # Ys <- scale(Y)
  # Cyy <- 1/(n-1) * t(Ys) %*% (Ys) # if we want to work with correlation matrix

  # Initialize
  B <- B0 <- matrix(.5, nrow = q, ncol = p)
  Sigma <- Sigma0 <- diag(1, p)

  # EM
  for (i in 1:1e3){

    # E step estiamtes ---------------------------------------------------------
    # Matrix
    imp_mat <- rbind(
      cbind(t(B) %*% R %*% B + Sigma, t(B) %*% R),
      cbind(R %*% B, R)
    )

    # Sweep matrix
    imp_mat_swept <- sweepGoodnight(imp_mat, target = 1:p)
    gamma <- imp_mat_swept[1:p, -c(1:p)]
    Delta <- imp_mat_swept[-c(1:p), -c(1:p)]

    # Cyz_hat
    Cyz_hat <- Cyy %*% gamma
    Czy_hat <- t(gamma) %*% Cyy

    # Czz_hat
    Czz_hat <- t(gamma) %*% Cyy %*% gamma + Delta

    # M step -------------------------------------------------------------------
    m_matrix <- rbind(
      cbind(Cyy, Cyz_hat),
      cbind(Czy_hat, Czz_hat)
    )
    m_matrix_swept <- sweepGoodnight(m_matrix, target = S_i + p)

    # Updates
    B <- m_matrix_swept[-c(1:p), 1:p]
    Sigma <- diag(diag(m_matrix_swept[1:p, 1:p]))
  }

  # Loadings
  lapply(list(
    RFA = loadings(fa_out)[, c(2,3,1)],
    MN2016 = varimax(Lambda)$loadings[],
    LR1982 = varimax(t(B))$loadings[, c(1, 3, 2)]
  ), round, 3)


# FA EM with missing values ----------------------------------------------------

  Y_miss <- mice::ampute(Y,
                         prop = .1,
                         patterns = matrix(c(1, 1, 0, rep(1, ncol(Y)-3),
                                             1, 0, 1, rep(1, ncol(Y)-3),
                                             1, 0, 0, rep(1, ncol(Y)-3)),
                                           ncol = ncol(Y),
                                           byrow = TRUE),
                         mech = "MAR",
                         type = "RIGHT"
  )

  # Define Missing data patterns
  patts <- mice::md.pattern(Y_miss$amp, plot = FALSE)
  R <- patts[-nrow(patts),-ncol(patts), drop = FALSE]
  R <- R[, colnames(Y), drop = FALSE]

  # Data dimensionality
  n <- nrow(Y)

  # Number of missing data patterns
  S <- nrow(R)
  # Columns observed for a given pattern
  # O <- apply(R, 1, function(x) {colnames(R)[x == 1]})
  O <- apply(R, 1, function(x) {x == 1}, simplify = FALSE)
  # Columns missings for a given pattern
  M <- apply(R, 1, function(x) {colnames(R)[x == 0]}) #
  # Define I matrices (which obs in which pattern)
  ry <- !is.na(Y_miss$amp)
  R_logi <- R == 1 # pattern config saved as True and False
  I <- vector("list", S)
  for (s in 1:S) {
    # s <- 1
    index <- NULL
    for (i in 1:n) {
      # i <- 1
      if(all.equal(ry[i, ], R_logi[s, ]) == TRUE) {
        index <- c(index, i)
      }
    }
    I[[s]] <- index
  }

  # Define sufficient statistics matrix (observed)
  Tobs_s <- vector("list", S)

  for (s in 1:S) {
    id_obs  <- I[[s]]
    dat_aug <- as.matrix(cbind(int = 1, Y[id_obs, , drop = FALSE]))
    Tobs_s[[s]] <- crossprod(dat_aug)

    # Fix NAs
    Tobs_s[[s]][is.na(Tobs_s[[s]])] <- 0
  }

  Tobs <- Reduce("+", Tobs_s)

  # Constants
  a <- Ybar <- colMeans(Y)
  R <- diag(q)
  S_i <- 1:q
  Cyy <- 1/(n-1) * t(Y - Ybar) %*% (Y - Ybar) # if data is centered, Cyy which is Sxb
  # Ys <- scale(Y)
  # Cyy <- 1/(n-1) * t(Ys) %*% (Ys) # if we want to work with correlation matrix

  # Initialize
  Phi <- diag(q)
  L <- L0 <- matrix(.5, nrow = p, ncol = q)
  Psi <- Psi0 <- diag(1, p)
  Sigma <- L %*% Phi %*% t(L) + Psi

  for (s in 1:S) {
    s <- 1
    id_obs  <- I[[s]]
    L[O[[s]], ]
    dat_aug <- as.matrix(cbind(int = 1, Y[id_obs, , drop = FALSE]))
    Tobs_s[[s]] <- crossprod(dat_aug)

    # Fix NAs
    Tobs_s[[s]][is.na(Tobs_s[[s]])] <- 0
  }


  # EM
  for (i in 1:1e3){

    # E step estiamtes ---------------------------------------------------------
    # Matrix
    imp_mat <- rbind(
      cbind(t(B) %*% R %*% B + Sigma, t(B) %*% R),
      cbind(R %*% B, R)
    )

    # Sweep matrix
    imp_mat_swept <- sweepGoodnight(imp_mat, target = 1:p)
    gamma <- imp_mat_swept[1:p, -c(1:p)]
    Delta <- imp_mat_swept[-c(1:p), -c(1:p)]

    # Cyz_hat
    Cyz_hat <- Cyy %*% gamma
    Czy_hat <- t(gamma) %*% Cyy

    # Czz_hat
    Czz_hat <- t(gamma) %*% Cyy %*% gamma + Delta

    # M step -------------------------------------------------------------------
    m_matrix <- rbind(
      cbind(Cyy, Cyz_hat),
      cbind(Czy_hat, Czz_hat)
    )
    m_matrix_swept <- sweepGoodnight(m_matrix, target = S_i + p)

    # Updates
    B <- m_matrix_swept[-c(1:p), 1:p]
    Sigma <- diag(diag(m_matrix_swept[1:p, 1:p]))
  }

  # Loadings
  lapply(list(
    RFA = loadings(fa_out)[, c(2,3,1)],
    MN2016 = varimax(Lambda)$loadings[],
    LR1982 = varimax(t(B))$loadings[, c(1, 3, 2)]
  ), round, 3)
