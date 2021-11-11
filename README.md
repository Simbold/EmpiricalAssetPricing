# EmpiricalAssetPricing

GMM estimation of a standard cross sectional asset pricing model:
- This repository cotains one .R file which can be used to estimate a standard cross sectional asset pricing model via GMM
- The file contains two function:
  1. OLS_two_pass_estimates() which calculates the coefficient and standard errors using OLS
  2. GMM_two_pass_estimates() which calculates the standard errors using the GMM method (Note: the coefficient are the same from both methods)
- For further info see Cochrane 2005: 12.2 Cross-Sectional Regressions 
- The estimator alows for an intercept in cross-sectional regression (parameter: 'constant') 
- The covariance matrix of the moment conditions is computed using the Newey-West HAC estimator
- The parameter 'adjust' refers to a finite sample adjust for the Newey-West estimator the adjustment = t/(t-#momentconditions), #momentconditions = n*(k+2), where k=#factors, n=#assets, t=#sampleSize
