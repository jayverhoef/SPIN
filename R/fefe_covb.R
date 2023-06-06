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
fefe_covb <- function(cope_out, compute_covb = TRUE,
  compute_covbalt1 = FALSE, compute_covbalt2 = FALSE) 
{
	covbList = NULL
	bhatList = NULL
  covbalt1 = NULL
  covbalt2 = NULL
  covb = NULL
  theta = cope_out$theta
  zlist = cope_out$zlist
  Xlist = cope_out$Xlist 
  xylist = cope_out$xylist 
  qrlist = cope_out$qrlist
  ViXlist = cope_out$ViXlist 
  Xlist = cope_out$Xlist
  bhat = cope_out$betaHat
  Sxx = cope_out$Sxx
  Sxxi = (theta[1] + theta[3])*solve(Sxx)
  ngrps = length(zlist)
	if(compute_covb == TRUE) {
    Wxx = Sxx*0
    for(i in 1:(ngrps -1)) {
      for(j in (i + 1):ngrps) {
        XViSigijVjX = t(ViXlist[[i]]) %*%
          (theta[1]*corModelExponential(as.matrix(pdist:::pdist(xylist[[i]],xylist[[j]]))/
            theta[2])) %*%
          ViXlist[[j]]/((theta[1] + theta[3])^2)
        Wxx = Wxx + XViSigijVjX + t(XViSigijVjX)
      }
    }
    covb = Sxxi + Sxxi %*% Wxx %*% Sxxi
  }
  if(compute_covbalt1 == TRUE | compute_covbalt2 == TRUE) {
    for(i in 1:ngrps) {
      XViX = t(Xlist[[i]]) %*% ViXlist[[i]]
      XViz = t(ViXlist[[i]]) %*% zlist[[i]]
      covbList[[i]] = (theta[1] + theta[3])*solve(XViX)
      bhatList[[i]] = covbList[[i]] %*% XViz/(theta[1] + theta[3])
    }
  }
  if(compute_covbalt1 == TRUE) {
    bhatbar = Reduce('+', bhatList)[,1]/ngrps
    covbalt1 = Sxx*0
    for(i in 1:ngrps) {
#      covbalt1 = covbalt1 + outer(as.vector((bhatList[[i]][,1] - bhatbar)),
#        as.vector((bhatList[[i]][,1] - bhatbar)))
      covbalt1 = covbalt1 + outer(as.vector((bhatList[[i]][,1] - bhat)),
        as.vector((bhatList[[i]][,1] - bhat)))
    }
    covbalt1= covbalt1/((ngrps - 1)*ngrps) # - 
#      outer(as.vector(bhatbar) - as.vector(bhat),
#      as.vector(bhatbar) - as.vector(bhat))
  }
  if(compute_covbalt2 == TRUE) covbalt2 = Reduce('+', covbList)/
        ngrps^2

	outpt <- list(
		bhat = bhat,
		covb = covb,
    covbalt1 = covbalt1,
    covbalt2 = covbalt2,
		bhatList = bhatList,
    covbList = covbList)
#	class(outpt) <- "fefe"
	outpt
}

