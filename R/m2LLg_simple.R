#-------------------------------------------------------------------------------
#
#          m2LL_simple
#
#-------------------------------------------------------------------------------

# minus 2 times REML log likelihood for partitioned data
m2LLg_simple <- function(theta, z, X, xycoords, corModel)
{
	# number of groups
  ngrps = length(z)
  # emply list to hold inverse of covariance matrices
  Sigma_i_inv = vector("list", ngrps)
  # empty matrices and scalars to hold sums
  p = dim(X[[1]])[2]
  Sxx = matrix(0, nrow = p, ncol = p)
  Sxy = matrix(0, nrow = p, ncol = 1)
  logDetV = 0
  logDetXViX = 0
  for(i in 1:ngrps) {
		# distance matrix for ith group, scaled by range parameter
	  dismat <- as.matrix(dist(xycoords[[i]]))/exp(theta[2])	
	  # if effective range is more than 4 times maximum distance, return 
	  # large value
		if(theta[2]/3 > 4*max(dismat)) return(1e+30)	  
	  # covariance matrix, scaled by partial sill
	  covMat <- exp(theta[1])*corModel(dismat) 
	  # add nugget effect and a small amount in case
	  # nugget effect is zero
    diag(covMat) = diag(covMat) + exp(theta[3]) + 1e-10
    # get inverse of covariance matrix for ith group and store it
    Sigma_i_inv[[i]] = solve(covMat)
    # intermediate matrix products
    ViX = Sigma_i_inv[[i]] %*% X[[i]]
    XViX = t(X[[i]]) %*% ViX
    # accumulated the sum of XViX
    Sxx = Sxx + XViX
    # accumulate the sum of X %*% Vi %*% y
    Sxy = Sxy + t(ViX) %*% z[[i]] 
    # accumulate determinants blockwise to get overall determinant
    logDetV = logDetV + as.numeric(determinant(covMat, 
      logarithm = TRUE)$modulus)
    logDetXViX = logDetXViX + as.numeric(determinant(XViX, 
      logarithm = TRUE)$modulus)
  }
	
	# global estimate of fixed effects after summing partition parts
  betaHat = solve(Sxx, Sxy)
  # get sum of squared residuals scaled by inverse covariance matrix
  rVir = 0
  for(i in 1:ngrps) 
    rVir = rVir + t(z[[i]] - X[[i]] %*% betaHat) %*% 
      Sigma_i_inv[[i]] %*% 
      (z[[i]] - X[[i]] %*% betaHat)

	# minus 2 time loglikelihood for REML (ignoring parts that don't contain
	# theta)
	minus2LL <- rVir + logDetV + logDetXViX
	as.numeric(minus2LL)
}

