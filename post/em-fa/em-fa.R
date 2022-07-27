# Project:   blogdown
# Objective: Initial notes on EM for factor analysis
# Author:    Edoardo Costantini
# Created:   2022-05-10
# Modified:  2022-07-05

rm(list = ls())

# Packages ---------------------------------------------------------------------

library(psych)
library(fastmatrix)
library(ISR3)
library(dplyr)

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

ordLoad <- function(x) x[, order(x[1, ])]

# Set up -----------------------------------------------------------------------

  # Based on raw data ----------------------------------------------------------

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
  R <- diag(q)
  Ybar <- colMeans(Y)
  
  # Larger data
  bfi_complete <- bfi[rowSums(is.na(bfi)) == 0, 1:25]
  Y <- as.matrix(bfi_complete)
  n <- nrow(Y)
  p <- ncol(Y)
  q <- 5
  R <- diag(q)
  Ybar <- colMeans(Y)

  # Covariance based 
  Ycov <- cov(Y)                                  # Pysch and Stats
  Sxb <- 1 / (n - 1) * t(Y - Ybar) %*% (Y - Ybar) # McNicholas
  Cyy <- 1 / (n - 1) * t(Y - Ybar) %*% (Y - Ybar) # RT and LR

  # Correlation based
  Ycov  <- cor(Y)                     # Pysch and Stats
  Sxb   <- cor(Y)                     # McNicholas
  Ys <- scale(Y)
  Cyy <- 1 / (n - 1) * t(Ys) %*% (Ys) # RT and LR

  # Based on covariance directly -----------------------------------------------

  # Ycov <- Harman74.cor$cov   # Pysch and Stats
  # Sxb <- Y                # McNicholas
  # Cyy <- Y                # RT and LR
  # p <- ncol(Y)
  # q <- 4
  # R <- diag(q)
  
# Stats version ----------------------------------------------------------------
  
  stats_fa <- factanal(covmat = Ycov, factors = q, rotation = "varimax")
  loads_stats <- ordLoad(loadings(stats_fa)[])

# Psych version ----------------------------------------------------------------

  psych_fa <- psych::fa(Ycov,
    nfactors = q,
    covar = TRUE, # we want to work with the covariance matrix
    rotate = "varimax",
    fm = "ml"
  )
  loads_psych <- ordLoad(loadings(psych_fa)[])

# McNicholas -------------------------------------------------------------------

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
  loads_mcnic <- ordLoad(varimax(Lambda)$loadings[])
  uniqs_mcnic <- diag(Psi)

  # Sigma = LL' + Psi
  # Lambda %*% t(Lambda) + Psi
  # loadings(psych_fa) %*% t(loadings(psych_fa)) + diag(psych_fa$uniquenesses)

# RubinThayer1982 --------------------------------------------------------------

  # Initialize
  B <- B0 <- matrix(.5, nrow = q, ncol = p)
  Tau <- Tau0 <- diag(1, p)

  # EM
  for (i in 1:1e3){
    # E step estiamtes ---------------------------------------------------------
    delta <- solve(Tau + t(B) %*% R %*% B) %*% t(B) %*% R # delta = Beta
    Delta <- R - (R %*% B) %*% delta

    # Cyz_hat
    Cyz_hat <- Cyy %*% delta
    Czy_hat <- t(Cyz_hat)

    # Czz_hat = Theta
    Czz_hat <- Delta + t(delta) %*% Cyy %*% delta

    # M step -------------------------------------------------------------------
    B <- solve(Czz_hat) %*% t(Cyz_hat)      # B = Lambda
    taus <- diag(Cyy - Cyy %*% delta %*% B) # Taus = Psi
    Tau <- diag(taus)
  }

  # Store values
  loads_RT1982 <- ordLoad(varimax(t(B))$loadings[])
  uniqs_RT1982 <- diag(Tau)

# LiuRubin1998 -----------------------------------------------------------------

  # Constants
  a <- Ybar <- colMeans(Y)
  R <- diag(q)
  S_i <- 1:q

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
    B <- m_matrix_swept[-c(1:p), 1:p, drop = FALSE]
    Sigma <- diag(diag(m_matrix_swept[1:p, 1:p]))
  }

  # Store outputs
  loads_LR1982 <- ordLoad(varimax(t(B))$loadings[])
  uniqs_LR1982 <- diag(Sigma)

# Comparisons ------------------------------------------------------------------

  # Loadings Stats vs Psych
  round(loads_stats - loads_psych, 3) # Stats vs psych
  round(loads_psych - loads_mcnic, 3)
  round(loads_mcnic - loads_RT1982, 3)
  round(loads_RT1982 - loads_LR1982, 3)

  # Uniqueness
  round(
    data.frame(
      psych_vs_stats = psych_fa$uniquenesses - stats_fa$uniquenesses,
      psych_vs_mcnic = psych_fa$uniquenesses - uniqs_mcnic,
      mcnic_vs_RT1982 = uniqs_mcnic - uniqs_RT1982,
      RT1982_vs_LR1982 = uniqs_RT1982 - uniqs_LR1982
    ), 3
  )

# FA EM with NAs: Concatenate problems -----------------------------------------

  set.seed(1234)
  Y_miss <- mice::ampute(scale(Y),
                          prop = .3,
                          patterns = matrix(c(1, 1, 0, rep(1, ncol(Y)-3),
                                              1, 0, 1, rep(1, ncol(Y)-3),
                                              1, 0, 0, rep(1, ncol(Y)-3)),
                                            ncol = ncol(Y),
                                            byrow = TRUE),
                          mech = "MAR",
                          type = "RIGHT"
  )
  Y_miss <- as.matrix(Y_miss$amp)

  n <- nrow(Y_miss)
  p <- ncol(Y_miss)

  # Define Missing data patterns
  patts <- mice::md.pattern(Y_miss, plot = FALSE)
  R <- patts[-nrow(patts), -ncol(patts), drop = FALSE]
  R <- R[, colnames(Y_miss), drop = FALSE]

  # Number of missing data patterns
  S <- nrow(R)
  # Columns observed for a given pattern
  O <- apply(R, 1, function(x) {
    colnames(R)[x == 1]
  })
  # Columns missings for a given pattern
  M <- apply(R, 1, function(x) {
    colnames(R)[x == 0]
  }) #

  # Define I matrices (which obs in which pattern)
  ry <- !is.na(Y_miss)
  R_logi <- R == 1 # pattern config saved as True and False
  I <- vector("list", S)
  for (s in 1:S) {
    # s <- 1
    index <- NULL
    for (i in 1:n) {
      # i <- 1
      if (all.equal(ry[i, ], R_logi[s, ]) == TRUE) {
        index <- c(index, i)
      }
    }
    I[[s]] <- index
  }

  # Define sufficient statistics matrix (observed)
  Tobs_s <- vector("list", S)

  for (s in 1:S) {
    id_obs <- I[[s]]
    dat_aug <- as.matrix(cbind(int = 1, Y_miss[id_obs, , drop = FALSE]))
    Tobs_s[[s]] <- crossprod(dat_aug)

    # Fix NAs
    Tobs_s[[s]][is.na(Tobs_s[[s]])] <- 0
  }

  Tobs <- Reduce("+", Tobs_s)

  # Starting values
  # For NA part
  theta0 = matrix(rep(NA, (p+1)^2 ), ncol = (p+1),
                  dimnames = list(c("int", colnames(Y_miss)),
                                  c("int", colnames(Y_miss))
                                  ))
  theta0[, 1]   = c(-1, colMeans(Y_miss, na.rm = TRUE)) # T1 CC
  theta0[1, ] = c(-1, colMeans(Y_miss, na.rm = TRUE))
  theta0[-1, -1] = cov(Y_miss, use = "pairwise.complete.obs") * (n - 1) / n # T2 CC

  theta <- theta0

  # For FA part
  a <- a0 <- Ybar <- colMeans(Y_miss, na.rm = TRUE)
  R <- diag(q)
  # Ys <- scale(Y)
  # Cyy <- 1 / (n - 1) * t(Ys) %*% (Ys) # if we want to work with correlation matrix
  ccases <- rowSums(is.na(Y_miss)) == 0
  Cyy <- Cyy0 <- 1 / (sum(ccases) - 1) * t(Y_miss[ccases, ] - Ybar) %*% (Y_miss[ccases, ] - Ybar) # if data is centered, Cyy which is Sxb
  B <- B0 <- matrix(.5, nrow = q, ncol = p)
  Tau <- Tau0 <- diag(1, p)

  ## EM algorithm
  for (i in 1:1e3) {

    # E step estiamtes ---------------------------------------------------------
    # Get a new estimate of Cyy
    # Reset T to info in the data
    # (will be updated based on new theta at every new iteration)
    T <- Tobs
    if (S > 1) {
      # We only need to do this if there are missing data patterns other than
      # the fully observed pattern (s = 1)
      for (s in 2:S) {
        # Description: For every missing data patter, except the first one (complete data
        # missing data pattern)
        # s <- 2
        obs <- I[[s]]
        v_obs <- O[[s]]
        v_mis <- M[[s]]
        v_all <- colnames(Y_miss)

        # Sweep theta over predictors for this missing data pattern
        theta <- ISR3::SWP(theta, v_obs)

        # Define expectations (individual contributions)
        betas <- theta[c("int", v_obs), v_mis]
        cjs <- cbind(1, Y_miss[obs, v_obs, drop = FALSE]) %*% betas

        # Update T matrix ##
        for (i in seq_along(obs)) {
          for (j in seq_along(v_mis)) {
            # j <- 1
            # Update for mean
            J <- which(v_all == v_mis[j])
            T[1, J + 1] <- T[1, J + 1] + cjs[i, j]
            T[J + 1, 1] <- T[1, J + 1]

            # Update for covariances w/ observed covariates for this id
            # (for Ks observed for this id)
            for (k in seq_along(v_obs)) {
              # k <- 1
              K <- which(v_all == v_obs[k])
              T[K + 1, J + 1] <- T[K + 1, J + 1] + cjs[i, j] * Y_miss[obs[i], K]
              T[J + 1, K + 1] <- T[K + 1, J + 1]
            }

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
        theta <- ISR3::RSWP(theta, v_obs)
        # Note: this corresponds to the reverse sweep in the first
        # loop performed in the algorithm proposed by Schafer 1997.
        # It basically replaces the "if r_sj = 0 and theta_jj < 0".
        # For one E step, the covariance matrix used to compute individual
        # contirbutions in each missing data pattern is the same!
      }
    }

    # Use new Cyy
    Cyy <- theta[colnames(Y_miss), colnames(Y_miss)]

    delta <- solve(Tau + t(B) %*% R %*% B) %*% t(B) %*% R # delta = Beta
    Delta <- R - (R %*% B) %*% delta

    # Cyz_hat
    Cyz_hat <- Cyy %*% delta
    Czy_hat <- t(Cyz_hat)

    # Czz_hat = Theta
    Czz_hat <- Delta + t(delta) %*% Cyy %*% delta

    # M step -------------------------------------------------------------------
    # from regular EM for missing data
    theta <- SWP((n^(-1) * T), 1)

    # from regular EM for factor analysis
    B <- solve(Czz_hat) %*% t(Cyz_hat) # B = Lambda
    taus <- diag(Cyy - Cyy %*% delta %*% B) # Taus = Psi
    Tau <- diag(taus)
  }

  # Store values
  loads_FAMISS <- ordLoad(varimax(t(B))$loadings[])
  uniqs_FAMISS <- diag(Tau)

  # Loadings Stats vs Psych
  round(loads_RT1982 - loads_FAMISS, 3)
  round(loads_stats - loads_psych, 3) # Stats vs psych
  round(loads_psych - loads_mcnic, 3)
  round(loads_mcnic - loads_RT1982, 3)

  # Uniqueness
  round(
    data.frame(
      psych_vs_stats = psych_fa$uniquenesses - stats_fa$uniquenesses,
      psych_vs_mcnic = psych_fa$uniquenesses - uniqs_mcnic,
      mcnic_vs_RT1982 = uniqs_mcnic - uniqs_RT1982,
      RT1982_vs_LR1982 = uniqs_RT1982 - uniqs_FAMISS
    ), 3
  )


# FA EM with NAs: Naive use of augmented matrix --------------------------------

  set.seed(1234)
  # bfi_complete <- bfi[rowSums(is.na(bfi)) == 0, 1:25]
  # Y <- bfi_complete
  # q <- 5
  # n <- nrow(Y)
  # p <- ncol(Y)
  Y_miss <- mice::ampute(Y,
                         prop = .3,
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

  # Define initial theta0

  # Constants
  a <- Ybar <- colMeans(Y_m[, 1:p], na.rm = TRUE)
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

  iters <- 1e2

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

          # Update for covariances w/ observed covariates for this id
          # (for Ks observed for this id)
          for (k in seq_along(v_obs)) {
            # k <- 1
            K <- which(v_all == v_obs[k])
            T[K + 1, J + 1] <- T[K + 1, J + 1] + cjs[i, j] * Y[obs[i], K]
            T[J + 1, K + 1] <- T[K + 1, J + 1]
          }

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
      T[1, -c(1:(p+1))] <- 0
      T[-c(1:(p+1)), 1] <- 0
      T[-c(1:(p+1)), -c(1:(p+1))] <- 0

      theta <- ISR3::RSWP(theta, v_obs)
      theta[1, -c(1:(p+1))] <- 0
      theta[-c(1:(p+1)), 1] <- 0
      theta[-c(1:(p+1)), -c(1:(p+1))] <- R
      # Note: this corresponds to the reverse sweep in the first
      # loop performed in the algorithm proposed by Schafer 1997.
      # It basically replaces the "if r_sj = 0 and theta_jj < 0".
      # For one E step, the covariance matrix used to compute individual
      # contirbutions in each missing data pattern is the same!
      
    }

    # > M-step ####
    theta <- ISR3::SWP((n^(-1) * T), 1)
    theta[1, -c(1:(p+1))] <- 0
    theta[-c(1:(p+1)), 1] <- 0
    theta[-c(1:(p+1)), -c(1:(p+1))] <- R
  }

(Betas <- theta[2:(p+1), -c(1:(p+1))])

# FA EM with NAs: correct sufficient stats -------------------------------------

# Jamshidian1994 - An EM Algorithm for ML Factor Analysis with Missing Data

  x <- t(Y_m)
  xi <- x[, 1]

# Order data: o, m, l
  patts <- mice::md.pattern(Y_m, plot = FALSE)
  R <- patts[-nrow(patts), -ncol(patts), drop = FALSE]
  p_o <- length(which(colSums(is.na(Y_m)) == 0))
  p_m <- p - p_o
  x <- t(Y_m[, colnames(R)])

# Define starting values for everything
  muo     <- rowMeans(x[1:p_o, ])
  mum     <- rowMeans(x[(p_o+1):(p_o+p_m), ], na.rm = TRUE)
  Lo      <- matrix(.5, nrow = p_o, ncol = q)
  Lm      <- matrix(.5, nrow = p_m, ncol = q)
  yo      <- x[1:p_o, ]
  Psio    <- diag(p_o)
  Psim    <- diag(p_m)
  Phi     <- diag(q)
  Sigma00 <- Lo %*% Phi %*% t(Lo) + Psio
  Sigmamm <- Lm %*% Phi %*% t(Lm) + Psim

# E step Compute expecationts

  # Xbarstar
  yms <- mum + Lm %*% Phi %*% t(Lo) %*% solve(Sigma00) %*% (yo - muo)
  xbs <- 1/n * rowSums(rbind(yo, yms))

  # Sstar
  Eymym <- Sigmamm - Lm %*% Phi %*% t(Lo) %*% solve(Sigma00) %*% Lo %*% Phi %*% t(Lm) + yms %*% t(yms)
  Exxt <- rbind(
    cbind(
      yo %*% t(yo), yo %*% t(yms)
    ),
    cbind(
      yms %*% t(yo), Eymym
    )
  )
  Ss <- 1 / n * Exxt

  # fstar (fixed to 0)
  fs <- Ef <- Phi %*% t(Lo) %*% solve(Sigma00) %*% (yo - muo)
  fbs <- rowSums(Ef) / n
  fbs <- rep(0, q)

  # Fstar
  Eff <- Phi - Phi %*% t(Lo) %*% solve(Sigma00) %*% Lo %*% Phi + fs %*% t(fs)
  Fs <- 1/n * Eff

  # Vstar
  Exf <- rbind(
    yo %*% t(fs),
    Lm %*% Phi - Lm %*% Phi %*% t(Lo) %*% solve(Sigma00) %*% Lo %*% Phi + yms %*% t(fs)
  )
  Vs <- 1 / n * Exf

# M steps
  B <- Fs - fbs %*% t(fbs)
  Lp <- Vs - rbind(yo, yms) %*% t(fs) # <- this does not add up!
  G <- Ss - 2 * xbs %*% c(muo, mum) - 2 * Vs %*%

# Other experiments ------------------------------------------------------------

  # Constants
  a <- Ybar <- colMeans(Y)
  R <- diag(q)
  S_i <- 1:q

  # Initialize
  B <- B0 <- matrix(.5, nrow = q, ncol = p)
  Sigma <- Sigma0 <- diag(1, p)

  # EM
  for (i in 1:1e3) {

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
    B <- m_matrix_swept[-c(1:p), 1:p, drop = FALSE]
    Sigma <- diag(diag(m_matrix_swept[1:p, 1:p]))
  }

  # Store outputs
  loads_LR1982 <- ordLoad(varimax(t(B))$loadings[])
  uniqs_LR1982 <- diag(Sigma)