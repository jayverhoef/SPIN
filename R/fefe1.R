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
fefe1 <- function(cope_out = NULL) 
{
	theta = cope_out$theta
	xylist = cope_out$xylist
	ViXlist = cope_out$ViXlist
	betaHat = cope_out$betaHat
	Sxx = cope_out$Sxx
 	p = dim(ViXlist[[1]])[2]
	ngrps = length(ViXlist)

	Wxx = matrix(0, nrow = p, ncol = p)
	for(i in 1:(ngrps - 1)) {
		for(j in (i + 1):ngrps) {
				dismat = as.matrix(pdist(xylist[[i]],xylist[[j]]))/theta[2]
				covMatg <- theta[1]*corModelExponential(dismat) 
				Wxx = Wxx + t(ViXlist[[i]]) %*% covMatg %*% ViXlist[[j]]
			}
		}

	covb = solve(Sxx) + 2*solve(Sxx) %*% Wxx %*% solve(Sxx)

	outpt <- list(
		bhat = cope_out$betaHat,
		covb = covb,
    Wxx = Wxx
	)
#	class(outpt) <- "fefe"
	outpt
}


