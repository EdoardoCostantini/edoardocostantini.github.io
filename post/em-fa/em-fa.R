# Project:   blogdown
# Objective: Initial notes on EM for factor analysis
# Author:    Edoardo Costantini
# Created:   2022-05-10
# Modified:  2022-05-10

rm(list = ls())

# Packages ---------------------------------------------------------------------

library(psych)
library(fastmatrix)

# Functions --------------------------------------------------------------------

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

augmentCov <- function(covmat, center,
                       dnames = list(c("int", colnames(covmat)),
                                       c("int", colnames(covmat)))
) {
  # Internals -------------------------------------------------------------

  # covmat = cov(mtcars)      # covariance matrix of a dataset
  # center = colMeans(mtcars) # vector of means of a dataset
  # dnames = list(c("int", colnames(covmat)),
  #                 c("int", colnames(covmat))) # list of row and column names

  # Body ------------------------------------------------------------------

  # Store the dimensionality of the data
  p <- ncol(covmat)

  # Define the structure of the augmented covariance matrix
  augCov <- matrix(rep(NA, (p + 1)^2 ),
                   ncol = (p + 1),
                   dimnames = dnames)

  # Assign the values to the correct slot
  augCov[1, 1] <- -1
  augCov[-1, 1] <- center
  augCov[1, -1] <- center
  augCov[-1,-1] <- covmat

  # Output
  return(augCov)
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


# FA EM with NAs: Naive use of augemnted matrix --------------------------------

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

  Y_m <- as.matrix(cbind(Y_miss$amp, matrix(rep(NA, n*q), ncol = q)))

  # Define Missing data patterns
  patts <- mice::md.pattern(Y_m, plot = FALSE)
  R <- patts[-nrow(patts),-ncol(patts), drop = FALSE]
  R <- R[, colnames(Y_m), drop = FALSE]

  # Data dimensionality
  n <- nrow(Y_m)

  # Number of missing data patterns
  S <- nrow(R)
  # Columns observed for a given pattern
  O <- apply(R, 1, function(x) {colnames(R)[x == 1]})
  # Columns missings for a given pattern
  M <- apply(R, 1, function(x) {colnames(R)[x == 0]}) #
  # Define I matrices (which obs in which pattern)
  ry <- !is.na(Y_m)
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
    dat_aug <- as.matrix(cbind(int = 1, Y_m[id_obs, , drop = FALSE]))
    Tobs_s[[s]] <- crossprod(dat_aug)

    # Fix NAs
    Tobs_s[[s]][is.na(Tobs_s[[s]])] <- 0
  }

  Tobs <- Reduce("+", Tobs_s)

  # Define intial theta0

  # Constants
  a <- Ybar <- colMeans(Y_m[1:p], na.rm = TRUE)
  R <- diag(q)
  S_i <- 1:q
  B <- B0 <- matrix(.5, nrow = q, ncol = p)
  Sigma <- Sigma0 <- diag(1, p)

  # Covariance matrix
  cov_mat <- rbind(
    cbind(t(B) %*% R %*% B + Sigma, t(B) %*% R),
    cbind(R %*% B, R)
  )

  # Augmented covariance matrix
  theta0 <- augmentCov(
    covmat = cov_mat,
    center = c(a, rep(0, q)),
    dnames = dimnames(Tobs)
  )

  ## EM algorithm
  # Starting values
  theta <- theta0

  iters <- 3

  # Iterations
  for (it in 1:iters) {
    # Reset T to info in the data
    # (will be updated based on new theta at every new iteration)
    #> E-step ####
    T <- Tobs
    # We only need to do this if there are missing data patterns other than
    # the fully observed pattern (s = 1)
    for (s in 1:S) {
      # Description: For every missing data patter, except the first one (complete data
      # missing data pattern)
      # s <- 1
      obs <- I[[s]]
      v_obs <- O[[s]]
      v_mis <- M[[s]]
      v_all <- colnames(Y_m)

      # Sweep theta over predictors for this missing data pattern
      # theta <- sweepGoodnight(theta, which(v_all %in% v_obs)+1)
      theta <- ISR3::SWP(theta, v_obs)
      
      # Define expectations (individual contributions)
      betas <- theta[c("int", v_obs), v_mis]
      cjs <- cbind(1, Y_m[obs, v_obs, drop = FALSE]) %*% betas

      # Update T matrix ##
      for (i in seq_along(obs)) {
        # i <- 1
        for (j in seq_along(v_mis)) {
          # j <- 1
          # Update for mean
          J <- which(v_all == v_mis[j])
          T[1, J + 1] <- T[1, J + 1] + cjs[i, j]
          T[J + 1, 1] <- T[1, J + 1]

          # Update for covariances w/ unobserved covariates for this id
          # (both j and k missing, includes covariances with itself k = j)
          for (k in seq_along(v_mis)) {
            # k <- 1
            K <- which(v_all == v_mis[k])
            if (K >= J) {
              T[K + 1, J + 1] <- T[K + 1, J + 1] + theta[K + 1, J + 1] + cjs[i, j] * cjs[i, k]
              T[J + 1, K + 1] <- T[K + 1, J + 1]
            }
          }
        }
      }
      # Make sure the means for the factors are 0s
      # T[1, -c(1:(p+1))] <- 0
      # T[-c(1:(p+1)), 1] <- 0
      theta <- ISR3::RSWP(theta, v_obs)
      # Note: this corresponds to the reverse sweep in the first
      # loop performed in the algorithm proposed by Schafer 1997.
      # It basically replaces the "if r_sj = 0 and theta_jj < 0".
      # For one E step, the covariance matrix used to compute individual
      # contirbutions in each missing data pattern is the same!
      
    }

    # > M-step ####
    theta <- ISR3::SWP((n^(-1) * T), 1)
  }

# FA EM with NAs: correct sufficient stats -------------------------------------
# Jamshidian1994 - An EM Algorithm for ML Factor Analysis with Missing Data
# Order data: o, m, l
patts <- mice::md.pattern(Y_m, plot = FALSE)
R <- patts[-nrow(patts), -ncol(patts), drop = FALSE]
Y_m <- Y_m[, colnames(R)]
p_o <- length(which(colSums(is.na(Y_m)) == 0))
p_m <- p - p_o

# Define starting values for everything
muo <- colMeans(Y_m[, 1:p_o])
Lo <- matrix(.5, nrow = p_o, ncol = q)
Psio <- diag(p_o)
Sigma00 <- Lo %*% Phi %*% t(Lo) + Psio
mum <- colMeans(Y_m[, (p_o+1):(p_o+p_m)], na.rm = TRUE)
Lm <- matrix(.5, nrow = p_m, ncol = q)
Psim <- diag(p_m)
Sigmamm <- Lm %*% Phi %*% t(Lm) + Psim
Phi <- diag(q)

# Compute expecationts
# Xbarstar
yms <- mum + Lm %*% Phi %*% t(Lo) %*% solve(Sigma00) %*% (t(Y_m[, 1:p_o]) - muo)
xbs <- 1/n * colSums(cbind(Y_m[, 1:p_o], t(yms)))

# Sstar
Eymym <- Sigmamm - Lm %*% Phi %*% t(Lo) %*% solve(Sigma00) %*% Lo %*% Phi %*% t(Lm) + yms %*% t(yms)
Exxt <- rbind(
  cbind(
    t(Y_m[, 1:p_o]) %*% Y_m[, 1:p_o], t(Y_m[, 1:p_o]) %*% t(yms)
  ),
  cbind(
    yms %*% Y_m[, 1:p_o], Eymym
  )
)
Ss <- 1 / n * Exxt

# fstar
fs <- Ef <- Phi %*% t(Lo) %*% solve(Sigma00) %*% (t(Y_m[, 1:p_o]) - muo)
fbs <- rowSums(Ef) / n

# Fstar
Eff <- Phi - Phi %*% t(Lo) %*% solve(Sigma00) %*% Lo %*% Phi + fs %*% t(fs)
Fs <- 1/n * Eff

# Vstar
Exf <- rbind(
  t(Y_m[, 1:p_o]) %*% t(fs),
  Lm %*% Phi - Lm %*% Phi %*% t(Lo) %*% solve(Sigma00) %*% Lo %*% Phi + yms %*% t(fs)
)
Vs <- 1 / n * Exf

# M steps
B <- Fs - fbs %*% t(fbs)
Lp <- Vs - xbs %*% t(fs)
G <- Ss - 2 * xbs %*% c(muo, mum) - 2 * Vs %*% 



t(t(Y_m[, 1:p_o]) - muo)
Y_m[18, 1:p_o] - muo

  # Constants
  a <- Ybar <- colMeans(Y)
  R <- diag(q)
  S_i <- 1:q
  # Cyy <- 1/(n-1) * t(Y - Ybar) %*% (Y - Ybar) # if data is centered, Cyy which is Sxb
  # Ys <- scale(Y)
  # Cyy <- 1/(n-1) * t(Ys) %*% (Ys) # if we want to work with correlation matrix

  # Initialize
  B <- B0 <- matrix(.5, nrow = q, ncol = p)
  Sigma <- Sigma0 <- diag(1, p)

  # E step estiamtes ---------------------------------------------------------
  # Covariance matrix
  cov_mat <- rbind(
    cbind(t(B) %*% R %*% B + Sigma, t(B) %*% R),
    cbind(R %*% B, R)
  )

  # Augmented covariance matrix
  Theta <- augmentCov(
    covmat = cov_mat,
    center = c(a, rep(0, q)),
    dnames = NULL
  )

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
