#-------------------------------------------------------------------------------
#
#           fefe
#
#-------------------------------------------------------------------------------

#' Fixed effects estimation for geostatistical linear model
#'
#' Fixed effects estimation for geostatistical linear model
#'
#' @param cope_out Object of class 'cope.'  Output from cope() function.
#' @param data Data set used with cope() function
#' @param spatial_model  Spatial autocorrelation models for errors.  The list of spatial autocorrelation 
#' models is "exponential","expRadon2","expRadon4","gaussian","stable",
#' "rationalQuad","cauchyGrav","cauchyMag","cauchy","circular","spherical",
#' "cubic","penta","cardinalSine","besselK","besselJ"  Default is "exponential".
#' Any names in the list not given above will be searched among the columns in the data set and
#' used as a factor variable for levels of a traditional random effect.
#' @param theta input the covariance parameters directly
#' @param z a vector of observed data
#' @param X design matrix for observed data
#' @param xcoords x-coordinates for observed data.  Same dimension as z.
#' @param ycoords y-coordinates for observed data.  Same dimension as z.
#' @param extrap extra parameter for those autocorrelation models that use it.
#' Currently, these are models "stable", "besselK", "besselJ".
#' @param use parallel processing?  TRUE or FALSE.
#' @param  subsampindx Column containing indexes for data partitioning.
#'
#' @return a list of class "splmm".  The functions "summary" and "print" are used to obtain and print a summary. "anova" returns just the analysis of variance table...
#'
#' @author Jay Ver Hoef
#' @export
fefe2 <- function(cope_out = NULL, data, spatial_model = corModelExponential, 
	theta = NULL, z = NULL, X = NULL, xcoords = NULL, ycoords = NULL, 
	extrap = NULL, par = FALSE, subsampindx = NULL) 
{

	theta = cope_out$theta
	xycoords = cope_out$xylist
	ViXlist = cope_out$ViXlist
	z = cope_out$zlist
	xycoords = cope_out$xylist
	qrlist = cope_out$qrlist
	X = cope_out$Xlist
	betaHat = cope_out$betaHat
	Sxx = cope_out$Sxx
 	p = dim(X)[2]
	ngrps = length(ViXlist)

			bhatList = vector("list", ngrps)
			covbList = vector("list", ngrps)
			Sxx = matrix(0, nrow = p, ncol = p)
			Sxy = matrix(0, nrow = p, ncol = 1)

			for(i in 1:ngrps) {
				dismat <- as.matrix(dist(xycoords[[i]]))/theta[2]	
				covMat <- theta[1]*corModelExponential(dismat)
				diag(covMat) = diag(covMat) + theta[3] + 1e-10*theta[1]
				QR = qr(covMat, LAPACK = TRUE)
				ViX = solve(QR, X[[i]])
				XViX = crossprod(X[[i]],ViX)
				XViz = crossprod(z[[i]],ViX)
				covbList[[i]] = solve(XViX)
				bhatList[[i]] = covbList[[i]] %*% t(XViz)
				Sxx = Sxx + XViX
				Sxy = Sxy + t(XViz)
			}
			bhat = solve(Sxx) %*% Sxy


	outpt <- list(
		bhat = bhat,
		betaHat = betaHat,
		bhatList = bhatList,
		covbList = covbList,
    Sxx = Sxx,
    Sxy = Sxy
	)
#	class(outpt) <- "fefe"
	outpt
}

