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
fefe <- function(cope_out = NULL, data, spatial_model = corModelExponential, 
	theta = NULL, z = NULL, X = NULL, xcoords = NULL, ycoords = NULL, 
	extrap = NULL, par = FALSE, subsampindx = NULL, compute_covbalt1 = FALSE,
  compute_covbalt2 = FALSE) 
{
	grpindx = subsampindx
	Vilist = NULL
	bhatList = NULL
	Sxx =  NULL
	Sxy =  NULL
	Wxx = NULL
	ViX = NULL
	Viz = NULL
  covbalt1 = NULL
  covbalt2 = NULL
	p = dim(X)[2]

	if(!is.null(cope_out)) {
		if(class(data) == 'SpatialPointsDataFrame') {
			DF = data@data
		}

		trms <- terms(cope_out$formula, data = pDF)
		respCol <- cope_out$respCol
		theta = cope_out$theta
		xcoords = cope_out$coords[,1]
		ycoords = cope_out$coords[,2]
		z = cope_out$z
		X = cope_out$X
  }
	if(is.null(subsampindx)) {
		dismat <- distGeoAni(xcoords, ycoords, xcoords, ycoords, 
			rotate = 90, range = theta[2], minorp = 1)	
		if(length(theta) == 4) extrap = theta[4]
		covMat <- theta[1]*spatial_model(dismat) +
				theta[3]*diag(dim(dismat)[1])
		QR = qr(covMat)
		ViX <- solve(covMat, X)
		Viz = solve(covMat,z)
		XViX <- crossprod(X, ViX)
		covb <- solve(XViX)
		bhat <- covb %*% crossprod(ViX, z)
	}
	if(!is.null(subsampindx)) {
		Wxxloop = function(k, ijpair, theta, X, Vilist, 
			xcoords, ycoords, grpindx)
		{ 
			Wxx.5 = t(X[grpindx == ijpair[k,1],]) %*% Vilist[[ijpair[k,1]]] %*% (theta[1]*
				corModelExponential(
					distGeoAni(
						xcoords[grpindx == ijpair[k,1]], ycoords[grpindx == ijpair[k,1]],
						xcoords[grpindx == ijpair[k,2]], ycoords[grpindx == ijpair[k,2]], 
						rotate = 90, range = theta[2], minorp = 1
					)
				)) %*% Vilist[[ijpair[k,2]]]%*% X[grpindx == ijpair[k,2],]
      Wxx.5 + t(Wxx.5)
		}
		if(par == TRUE) {
			Vilist = vector("list", max(grpindx)) 
			Sxx = matrix(0, nrow = p, ncol = p)
			Sxy = matrix(0, nrow = p, ncol = 1)
			for(i in 1:max(grpindx)) {
				dismat <- distGeoAni(xcoords[grpindx == i],
					ycoords[grpindx == i], xcoords[grpindx == i], 
					ycoords[grpindx == i], rotate = 90, 
					range = theta[2], minorp = 1)	
				covMatg <- theta[1]*corModelExponential(dismat) +
					theta[3]*diag(dim(dismat)[1])
				Vilist[[i]] = solve(covMatg)
				Sxx = Sxx + t(X[grpindx == i,]) %*% Vilist[[i]] %*% 
					X[grpindx == i,]
				Sxy = Sxy + t(X[grpindx == i,]) %*% Vilist[[i]] %*% 
					z[grpindx == i]
			}
			bhat = solve(Sxx) %*% Sxy

			ijpair = NULL
			for(i in 1:(max(grpindx)-1)) {
				ijpair = rbind(ijpair, 
					cbind(rep(i, times = max(grpindx) - i), (i+1):max(grpindx))
				)
			}
			Wxxlist = foreach(k=1:length(ijpair[,1])) %dopar% {
				Wxxloop(k, ijpair, theta, X, Vilist, xcoords, ycoords, grpindx)
			}
			Wxx = Reduce('+', Wxxlist)
			covb = solve(Sxx) + 2*solve(Sxx) %*% Wxx %*% solve(Sxx)
		}
		if(par == FALSE) {
			Vilist = vector("list", max(grpindx))
			bhatList = vector("list", max(grpindx))
      covbList = vector("list", max(grpindx))
			Sxx = matrix(0, nrow = p, ncol = p)
			Sxy = matrix(0, nrow = p, ncol = 1)
			for(i in 1:max(grpindx)) {
				dismat <- distGeoAni(xcoords[grpindx == i],
					ycoords[grpindx == i], xcoords[grpindx == i], 
					ycoords[grpindx == i], rotate = 90, 
					range = theta[2], minorp = 1)	
				covMatg <- theta[1]*corModelExponential(dismat) +
					theta[3]*diag(dim(dismat)[1])
				Vilist[[i]] = solve(covMatg)
				XViX = t(X[grpindx == i,]) %*% Vilist[[i]] %*% 
					X[grpindx == i,]
				XViz = t(X[grpindx == i,]) %*% Vilist[[i]] %*% 
					z[grpindx == i]
        covbList[[i]] = solve(XViX)
				bhatList[[i]] = cbind(covbList[[i]] %*% XViz, diag(solve(XViX)))
				Sxx = Sxx + XViX
				Sxy = Sxy + XViz
			}
			bhat = solve(Sxx) %*% Sxy
      if(compute_covbalt1 == TRUE) {
        bhatbar = Reduce('+', bhatList)[,1]/max(grpindx)
        covbalt1 = Sxx*0
        for(i in 1:max(grpindx)) {
          covbalt1 = covbalt1 + outer(as.vector((bhatList[[i]][,1] - bhatbar)),
            as.vector((bhatList[[i]][,1] - bhatbar)))
        }
        covbalt1= covbalt1/((max(grpindx) - 1)*max(grpindx)) - 
          outer(as.vector(bhatbar) - as.vector(bhat),
            as.vector(bhatbar) - as.vector(bhat))
      }
      if(compute_covbalt2 == TRUE) covbalt2 = Reduce('+', covbList)/
        max(grpindx)^2

			ijpair = NULL
			for(i in 1:(max(grpindx)-1)) 
				ijpair = rbind(ijpair, 
					cbind(rep(i, times = max(grpindx) - i), (i+1):max(grpindx))
				)
			Wxx = 0
			for(k in 1:length(ijpair[,1])) {
				Wxx = Wxx + Wxxloop(k, ijpair, theta, X, Vilist, 
					xcoords, ycoords, grpindx)
			}
			covb = solve(Sxx) + solve(Sxx) %*% Wxx %*% solve(Sxx)
		}
	}

	outpt <- list(
		bhat = bhat,
		covb = covb,
    covbalt1 = covbalt1,
    covbalt2 = covbalt2,
		Vilist = Vilist,
		bhatList = bhatList,
    covbList = covbList,
    Sxx = Sxx,
    Sxy = Sxy,
    Wxx = Wxx,
    ViX = ViX,
    Viz = Viz
	)
#	class(outpt) <- "fefe"
	outpt
}

