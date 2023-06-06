#-------------------------------------------------------------------------------
#
#          m2LL
#
#-------------------------------------------------------------------------------

#' minus 2 times loglikelihood for case 12
#'
#' fits a geostatistical linear model
#'
#' @param theta vector of estimated covariance parameters
#' @param z vector of data
#' @param X design matrix for fixed effects
#' @param Z list of design matrices for each random effect 
#' @param xcoords vector with the x-coordinates
#' @param ycoords vector with the y-coordinates
#' @param estMeth estimation method.  Default is "REML" for restricted maximum likelihood.  Other options are "ML" for maximum likelihood
#'
#' @return minus 2 times the loglikelihood
#'
#' @author Jay Ver Hoef
m2LLgprof <- function(theta, z, X, xycoords, estMeth)
{
  if(any(abs(theta) > 10)) return(1e+30)
  p = length(X[[1]][1,])
  ngrps = length(z)
  qrlist = vector("list", ngrps)
  Sxx = matrix(0, nrow = p, ncol = p)
  Sxy = matrix(0, nrow = p, ncol = 1)
  logDetV = 0
  logDetXViX = 0
  df = 0
  for(i in 1:ngrps) {
		ni = length(z[[i]])
	  dismat <- as.matrix(dist(xycoords[[i]]))/exp(theta[2])	
	  covMat <- exp(theta[1])*corModelExponential(dismat)/(1 + exp(theta[1])) 
    diag(covMat) = rep(1, times = ni) # this adds the nugget effect
    qrlist[[i]] = qr(covMat, LAPACK = TRUE)
    ViX = solve(qrlist[[i]], X[[i]])
    XViX = crossprod(X[[i]],ViX)
    Sxx = Sxx + XViX
    Sxy = Sxy + t(crossprod(z[[i]],ViX)) 
    logDetV = logDetV + sum(log(abs(diag(qr.R(qrlist[[i]])))))
    logDetXViX = logDetXViX + as.numeric(determinant(XViX, 
      logarithm = TRUE)$modulus)
#    df = df + ni - p
		df = df + ni
  }
	
  betaHat = solve(Sxx, Sxy)
  rVir = 0
  for(i in 1:ngrps) 
    rVir = rVir + sum((z[[i]] - 
      X[[i]] %*% betaHat) *
      solve(qrlist[[i]], (z[[i]] - 
      X[[i]] %*% betaHat)))


  df = df - p
	minus2LL <- logDetV + df*log(rVir)
	if(estMeth == "REML") minus2LL <- minus2LL + logDetXViX
	as.numeric(minus2LL)
}

