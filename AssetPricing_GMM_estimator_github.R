rm(list=ls())
library(sandwich)

# Two pass regression estimates for a simple cross sectional asset pricing model
# As in Cochrane 2005: 12.2 Cross-Sectional Regressions, equation 12.11
two_pass_estimates <- function(r, f, constant=FALSE){
  # r:            t x n matrix of asset returns (t: sample size; n: number of assets/portfolios)
  # f:            t x k matrix of factors (t: sample size; k: number of factors)
  # constant:     TRUE or FALSE; refers to a intercept in the second stage / cross sectional regression
  
  k = dim(f)[2] # number of factors
  
  if(is.null(k)){
    k=1
  }
  
  r <- as.matrix(r)
  f <- as.matrix(f)
  
  # TimeSeries regression, equationwise OLS of portfolio returns on factor values and a constant 
  TSmodel <- lm(r ~ 1 + f) 
  
  betas <- t(TSmodel$coefficients[-1,]) 
  r_means <- colMeans(r) # compute the mean return of the portfolios/assets 
  
  if(k==1){
    betas = t(betas)
    }
  
  # CrossSectional Regression
  if(constant==FALSE){
    CSmodel <- lm(r_means ~ 0 + betas) # without a constant
  } else if(constant==TRUE){
    CSmodel <- lm(r_means ~ 1 + betas) # with a constant  
  }else{
    print("ERROR: variable 'constant' must be TRUE or FALSE")
  }
  return(list(TS = TSmodel, CS = CSmodel))
}

# GMM estimation of the 2 Pass estimation of a Asset Pricing model with a constant in the cross sectional regression if constant=TRUE
# Mostly following: Burnside, Craig. 2011. "The Cross Section of Foreign Currency Risk Premia and Consumption Growth Risk: Comment." American Economic Review, 101 (7): 3456-76.
GMM_se <- function(r, f, constant=FALSE, adjust=FALSE){
  # r:            t x n matrix of asset returns (t: sample size; n: number of assets/portfolios), if columns have names they will appear in the result
  # f:            t x k matrix of factors (t: sample size; k: number of factors), if columns have names they will appear in the result
  # constant:     TRUE or FALSE; refers to a intercept in the second stage / cross sectional regression
  # adjust:       TRUE or FALSE; Should a finite sample adjustment be made when computing the Newey-West covariance matrix? Note: the adjustment = t/(t-#momentconditions), #momentconditions = n*(k+2)
  
  k <- dim(f)[2] # number of factors 
  
  if(is.null(k)){
    k=1
  }
  
  n <- dim(r)[2] # number of portfolios/assets
  t <- dim(r)[1] # sample size 
  x <- as.matrix(cbind(r,f)) # matrix of our data: columns 1 to n are portfolio return, columns n+1 to n+k are factor values
  
  # preestimation of parameters
  # the only use of preestimating the parameters is to safe computational time during the optimization routine
  # this is since the GMM parameter estimates recover the OLS regression estimates
  pre_estim <- two_pass_estimates(r, f, constant) 
  tet <- c(pre_estim$TS$coefficients,pre_estim$CS$coefficients) 
  #tet <- rep(0, (k+1)*n+k) # in case of no constant
  #tet <- rep(0, (k+1)*n+k+1) # in case of constant
  
  # the n(k + 2) x T vector of moment conditions
  u_t <- function(tet, x){
        
    f_tilde <- t(cbind(1, x[,(n+1):(n+k)])) # K+1 x t matrix of constant and factor values
    r <- x[,1:n] # define r to be a t x n matrix of portfolio returns
    beta_tilde <- matrix(tet[1:(n*(k+1))], k+1, n) # matrix of first pass betas and intercepts
    betaOLS <- t(beta_tilde[-1,]) # the betas without the constants
    
    if(k==1){
      betaOLS  <- t(betaOLS)
    }
    
    lambdas <- tet[(n*(k+1)+1):(length(tet))] # the lambdas and (the gamma in case of constant=TRUE)
    
    # first block of moment condition for the first observation
    u <- f_tilde[,1] %*% (r[1,1] - t(f_tilde[,1]) %*% beta_tilde[,1])
    
    # the remaining n-1 blocks for the first set of moment conditions for the first observation
    for(i in 2:n){
      u <- rbind(u, f_tilde[,1] %*% (r[1,i] - t(f_tilde[,1]) %*% beta_tilde[,i]))
    }
    # do the same for each of the remaining observations in time and bind columnwise
    for(m in 2:t){
      u1 <-  f_tilde[,m] %*% (r[m,1] - t(f_tilde[,m]) %*% beta_tilde[,1])
      
      for(i in 2:n){
        u1 <- rbind(u1, f_tilde[,m] %*% (r[m,i] - t(f_tilde[,m]) %*% beta_tilde[,i]))
      }
      u <- cbind(u, u1)
    }
    
    # the last set of moment equations (responsible for the cross sectional regression)
    if(constant==FALSE){
      u <- rbind(u, t(r) - matrix(betaOLS %*% lambdas, n, t))
    }else if(constant==TRUE){
      u <- rbind(u, t(r) - lambdas[1] - matrix(betaOLS %*% lambdas[-1], n, t))
    }
    return(u)
  }
  
  beta_tilde <- matrix(tet[1:(n*(k+1))], k+1, n) # (k+1) x n matrix of betas and intercepts
  betaOLS <- t(beta_tilde[-1,]) # n x k matrix of betas
  
  if(k==1){
    betaOLS  <- t(betaOLS)
  }
  
  # building the aT weighting matrix which produces the two pass estimates
  if(constant==FALSE){
    aT <- matrix(0, n*(k+1)+k, n*(k+2))
    aT[1:(n*(k+1)), 1:(n*(k+1))] <- diag(1, n*(k+1), n*(k+1))
    aT[(n*(k+1)+1):dim(aT)[1], (n*(k+1)+1):dim(aT)[2]] <- t(betaOLS)
  }else if(constant==TRUE){
    aT <- matrix(0, (n+1)*(k+1), n*(k+2)) # (n + 1)(k + 1) x [n(k + 2) x 1] matrix
    aT[1:(n*(k+1)), 1:(n*(k+1))] <- diag(1, n*(k+1), n*(k+1))
    X <- cbind(1, betaOLS)
    aT[(n*(k+1)+1):dim(aT)[1], (n*(k+1)+1):dim(aT)[2]] <- t(X)
  }
  
  # function to be minimized
  obj <- function(tet){
    gT <- aT%*%rowMeans(u_t(tet,x))
    return(t(gT)%*%gT)
  }
  # Estimating the parameters! This step is not really necessary since we are only interested in the standard errors, since the 
  # GMM estimator reproduces the exact coefficients from the two pass OLS regressions
  theta <- optim(tet, obj, method="BFGS")$par
  
  # Building the gradient of the moment condition function
  DgT <- function(theta, x){
  
    f_tilde <- t(cbind(1, x[,(n+1):(n+k)]))
    M_f_tilde <- (f_tilde %*% t(f_tilde))/t
    
    if(constant==FALSE){
      # initialize matrix
      dg <- matrix(0, n*(k+2), n*(k+1)+k)
      
      # upper left corner of gradient matrix
      dg[1:(n*(k+1)), 1:(n*(k+1))] <- -diag(1, n, n) %x% M_f_tilde
      lambdas <- theta[(n*(k+1)+1):(length(theta))]
      
      # lower left corner of gradient matrix
      dg[(n*(k+1)+1):dim(dg)[1], 1:(n*(k+1))] <- -diag(1, n, n) %x% t(c(0, lambdas))
      
      betas_a <- matrix(theta[1:(n*(k+1))], k+1, n) # matrix of betas and first pass intercepts
      betaOLS <- t(betas_a[-1,]) # the betas without the intercepts
      
      if(k==1){betaOLS  <- t(betaOLS)}
        dg[(n*(k+1)+1):dim(dg)[1], (n*(k+1)+1):dim(dg)[2]] <- -betaOLS
        
    }else if(constant==TRUE){
      # initilize matrix
      dg <- matrix(0, n*(k+2), n*(k+1)+k+1)
      
      # upper left corner of gradient matrix
      dg[1:(n*(k+1)), 1:(n*(k+1))] <- -diag(1, n, n) %x% M_f_tilde
      lambdas <- theta[(n*(k+1)+2):(length(theta))]
      
      # lower left corner of gradient matrix
      dg[(n*(k+1)+1):dim(dg)[1], 1:(n*(k+1))] <- -diag(1, n, n) %x% t(c(0, lambdas))
      
      betas_a <- matrix(theta[1:(n*(k+1))], k+1, n) # matrix of betas and first pass intercepts
      betaOLS <- t(betas_a[-1,]) # the betas without the intercepts
      
      if(k==1){
        betaOLS  <- t(betaOLS)
      }
      
      X <- cbind(1, betaOLS)
        
      dg[(n*(k+1)+1):dim(dg)[1], (n*(k+1)+1):dim(dg)[2]] <- -X
    }
    return(dg)
  }
  
  # computing the covariance matrix, here using the Newey-West Bartlett HAC estimator
  S <- lrvar(t(u_t(theta, x)), type = "Newey-West", prewhite = FALSE, adjust = adjust)
  
  ad <- solve(aT%*%DgT(theta,x))
  V <- ad%*%aT%*%S%*%t(aT)%*%t(ad)
  se = sqrt(diag(V))
  
  # returning results
  coef1pass <- matrix(theta[1:(n*(k+1))], k+1, n)
  colnames(coef1pass) <- colnames(r)
  if(k>1){rownames(coef1pass) <- c("intercept", colnames(f))}
  
  coef2pass <- theta[(n*(k+1)+1):length(theta)]
  
  se1pass <- matrix(se[1:(n*(k+1))], k+1, n)
  colnames(se1pass) <- colnames(r)
  if(k>1){
    rownames(se1pass) <- c("intercept", colnames(f))
  }
  
  se2pass <- se[(n*(k+1)+1):length(se)]
  
  if(constant==TRUE){
    names(se2pass) <- c("intercept", colnames(f))
    names(coef2pass) <- c("intercept", colnames(f))
  }else if(constant==FALSE){
    names(se2pass) <- colnames(f)
    names(coef2pass) <- colnames(f)
  }

  results <- list(time_series_reg = list(coefficients = coef1pass, standard_errors = se1pass),  cross_sectional_reg=list(coefficients = coef2pass, standard_errors = se2pass), raw=list(coefficients = theta, standard_errors = se))
  return(results)
}


###############################################################################
###############################################################################

# Example:

Data <- read.csv(paste0(wdir, "/PfReturns.csv"), header =TRUE, sep=";", dec=",", check.names=TRUE)
r <- Data[,2:6] #  t x n matrix matrix of asset returns, Note: named columns will appear in the result
f <- as.matrix(Data[,23:24]) # t x k matrix of factors
constant <- TRUE
adjust <-  TRUE # If the time-series sample is too small compared to the number of assets and factors the adjustment factor becomes negative resulting in failure in computing the covariance matrix
# adjust:       TRUE or FALSE; Should a finite sample adjustment be made when computing the Newey-West covariance matrix? Note: the adjustment = t/(t-#momentconditions), #momentconditions = n*(k+2)

GMM_results <- GMM_se(r, f, constant, adjust)

# first pass / Time-series Results
print(GMM_results$time_series_reg)
# second pass / Cross sectional results
print(GMM_results$cross_sectional_reg)


# Standard OLS two stage regression results
OLS_results <- two_pass_estimates(r, f, constant) 

# first pass / Time-series Results
print(summary(OLS_results$TS))
# second pass / Cross sectional results
print(summary(OLS_results$CS))
