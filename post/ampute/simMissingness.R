### Title:    Simulate MAR Missingness via Logistic Regression
### Author:   Kyle M. Lang
### Created:  2019-11-06
### Modified: 2019-11-11

###--------------------------------------------------------------------------###

## Objective function for the logistic regression offset:
fOffset <- function(offset, eta, pm, type) {
    f <- switch(type,
                center = offset - abs(eta),
                tails  = abs(eta) - offset,
                offset + eta
                )
    
    (mean(plogis(f)) - pm)^2
}

###--------------------------------------------------------------------------###

## Optimize the logistic regression offset for a given value of the linear
## predictor (eta) to get a desired percent missing (pm):
optOffset <- function(pm, eta, type, tol = c(0.1, 0.001), maxIter = 10) {
    for(k in 1 : maxIter) {
        ## Define the search range:
        int <- k * range(eta)
        
        ## Optimize the objective over 'int':
        out <- optimize(f        = fOffset,
                        interval = int,
                        eta      = eta,
                        pm       = pm,
                        type     = type)
        
        ## Are we far enough from the boundary?
        dist   <- out$minimum - int
        check1 <- all(abs(dist) > tol[1] * diff(int))
        
        ## Are we within tolerance?
        check2 <- out$objective < tol[2]
        
        if(check1 & check2) break
    }
    ## Did we fail?
    if(!check1 | ! check2) stop("I could not optimize this function.")
    
    out
}

###--------------------------------------------------------------------------###
                                 
## Simulate a nonresponse vector:
simMissingness <- function(pm,
                           data,
                           preds = colnames(data),
                           type  = "high",
                           beta  = NULL)
{
    ## Define a trivial slope vector, if necessary:
    if(is.null(beta)) beta <- rep(1.0, length(preds))
    
    ## Compute (and standardize) the linear predictor:
    eta <- scale(
        as.matrix(data[ , preds]) %*% matrix(beta)
    )
    
    ## Optimize the offset:
    offset <- optOffset(pm = pm, eta = eta, type = type)$minimum
    
    ## Compute the probabilities of nonresponse:
    probs <- plogis(
        switch(type,
               high   = offset + eta,
               low    = offset - eta,
               center = offset - abs(eta),
               tails  = abs(eta) - offset
               )
    )
    
    ## Return a logical nonresponse vector:
    as.logical(rbinom(n = length(eta), size = 1, prob = probs))
}

###--------------------------------------------------------------------------###
