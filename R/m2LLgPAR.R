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
m2LLgPAR <- function(theta, z, X, xcoords, ycoords, 
	estMeth, grpindx)
{
	# theta = splmbd20out$theta
	# z = splmbd20out$z
	# X = splmbd20out$X
	# xcoords = splmbd20out$xcoords
	# ycoords = splmbd20out$ycoords
	# estMeth = 'REML'
	# grpindx = splmbd20out$grpindx

  if(theta[2] > 2.3) return(1e+30)
  p = length(X[1,])
  n = length(X[,1])
  Vilist = foreach(i=1:max(grpindx)) %dopar% {
	  dismat <- distGeoAni(xcoords[grpindx == i], ycoords[grpindx == i], 
      xcoords[grpindx == i], ycoords[grpindx == i], 
      rotate = 90, range = exp(theta[2]), minorp = 1)	
	  covMat <- exp(theta[1])*corModelExponential(dismat) +
      (exp(theta[3]) + 1e-6)*diag(dim(dismat)[1])
    QR = qr(covMat, LAPACK = TRUE)
		ViX = solve(QR, X[grpindx == i,])
		XViX = crossprod(X[grpindx == i,], ViX)
    list(Viz = solve(QR, z[grpindx == i,]), ViX = ViX, Sxx = XViX,
    Sxy = t(crossprod(z[grpindx == i], ViX)), 
    logDetV = sum(log(abs(diag(qr.R(QR))))),
    logDetXViX = as.numeric(determinant(XViX, 
      logarithm = TRUE)$modulus))
  }
  Sxx = Reduce('+',lapply(Vilist, function(x){x$Sxx}))
  Sxy = Reduce('+',lapply(Vilist, function(x){x$Sxy}))
  logDetV = Reduce('+',lapply(Vilist, function(x){x$logDetV}))
  logDetXViX = Reduce('+',lapply(Vilist, function(x){x$logDetXViX}))
  betaHat = solve(Sxx) %*% Sxy
  rVirlist = foreach(i=1:max(grpindx)) %dopar% {
		crossprod(z[grpindx == i] - X[grpindx == i,] %*% betaHat,
      Vilist[[i]]$Viz - Vilist[[i]]$ViX %*% betaHat)
  }
  rVir = Reduce('+',rVirlist)
	minus2LL <- logDetV + rVir + n*log(2*pi)
	if(estMeth == "REML") minus2LL <- minus2LL + logDetXViX - p*log(2*pi)
	return( as.numeric(minus2LL) )
}

