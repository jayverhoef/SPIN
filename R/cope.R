#-------------------------------------------------------------------------------
#
#           cope
#
#-------------------------------------------------------------------------------

#' Covariance parameter estimation for a geostatistical linear model
#'
#' Covariance parameter estimation for a geostatistical linear model
#'
#' @param formula an R linear model formula
#' @param data an data object with spatial coordinates.  If a plain data.frame, then xcoordscol and ycoordscol must be specified.  If an sp object, xcoordscol and ycoordscol should be NULL (the default).
#' @param x_column name, in quotes, of the column containing the x-coordinate (if it was not possible to obtain the coordinates from the data class). Default is NULL. 
#' @param y_column name, in quotes, of the column containing the x-coordinate (if it was not possible to obtain the coordinates from the data class). Default is NULL. 
#' @param estMeth estimation method.  Default is "REML" for restricted maximum likelihood.  Other options are "ML" for maximum likelihood
#' @param spatial_model spatial autocorrelation models 
#' for random errors.  The list of spatial autocorrelation 
#' models is "exponential","expRadon2","expRadon4","gaussian","stable",
#' "rationalQuad","cauchyGrav","cauchyMag","cauchy","circular","spherical",
#' "cubic","penta","cardinalSine","besselK","besselJ"  Default is "exponential".
#' @param random_formula a list of variance components, 
#' Any names in the list not given above will be searched among the columns in the data set and
#' used as a factor variable for levels of a traditional random effect.
#' @param useAnisotropy include anistropy in parameter estimation?  Default is "FALSE"
#' @param subsample_col A column of factors indicating grouping for use in subsampling for large data sets.  Default is "NULL," which creates a single grouping.

#'
#' @return a list of class "slm_cope".  The functions "summary" and "print" are used to obtain and print a summary. "anova" returns just the analysis of variance table...
#'
#' @author Jay Ver Hoef
#' @export
cope <- function(formula, data, x_column, y_column,
	spatial_model = 'exponential',
	random_formula = NULL, use.anistropy = FALSE,
  thetaini = c(2,1,0.02), estMeth = "REML", par = FALSE,
  subSampCol = NULL, optMeth = 'Nelder-Mead', profile = TRUE) 
{

	# ----------------------------------------------------------------------------
	# prepare data
	# ----------------------------------------------------------------------------

	#formula = as.formula(z ~ x)
	#data = d1
	#xcoordscol = 'xcoords'
	#ycoordscol = 'ycoords'
	#ngroups = 500
	trms <- terms(formula, data = data)
	respCol <- as.character(as.list(attr(trms,"variables")[-1]))[1]
	covList <- attr(trms,"term.labels")
	#design matrix
	X <- model.matrix(formula, data)
	z <- as.matrix(data[, respCol])
	#vectors of spatial coordinates
	xcoords <- data[,x_column]
	ycoords <- data[,y_column]

	n <- length(z)
	p <- sum(svd(X)$d>1e-10)
	if(is.null(subSampCol)) {grpindx = rep(1, times = n)
	} else { grpindx = data[,subSampCol] }

  xylist = vector(mode = "list", length = max(grpindx)) 
  zlist = vector(mode = "list", length = max(grpindx)) 
  Xlist = vector(mode = "list", length = max(grpindx)) 
  for(i in 1:max(grpindx)) {
    xylist[[i]] = cbind(xcoords[grpindx == i], ycoords[grpindx == i])
    zlist[[i]] = z[grpindx == i]
    Xlist[[i]] = X[grpindx == i,]
	}
	#undebug(m2LLgPAR)
	if(par == FALSE & profile == FALSE) {
	starttime = Sys.time()
	parmest <-optim(thetaini, m2LLg, z = zlist, X = Xlist, 
			xycoords = xylist, 
			estMeth = 'REML', method = optMeth)
	stoptime = Sys.time()
	theta <- exp(parmest$par)
	}
	if(par == FALSE & profile == TRUE) {
		starttime = Sys.time()
		if(length(thetaini) != 2) return('thetaini must be of length 2')
		parmest <-optim(thetaini, m2LLgprof, z = zlist, X = Xlist, 
				xycoords = xylist, 
				estMeth = 'REML', method = optMeth)
		stoptime = Sys.time()
		theta = parmest$par
		X = Xlist
		xycoords = xylist
		z = zlist
		p = length(X[[1]][1,])
		ngrps = length(z)
		qrlist = vector("list", ngrps)
		ViXlist = vector("list", ngrps)
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
			ViXlist[[i]] = ViX
			XViX = crossprod(X[[i]],ViX)
			Sxx = Sxx + XViX
			Sxy = Sxy + t(crossprod(z[[i]],ViX)) 
			logDetV = logDetV + sum(log(abs(diag(qr.R(qrlist[[i]])))))
			logDetXViX = logDetXViX + as.numeric(determinant(XViX, 
				logarithm = TRUE)$modulus)
			df = df + ni - p
		}
		
		betaHat = solve(Sxx, Sxy)
		rVir = 0
		for(i in 1:ngrps) 
			rVir = rVir + sum((z[[i]] - 
				X[[i]] %*% betaHat) *
				solve(qrlist[[i]], (z[[i]] - 
				X[[i]] %*% betaHat)))
		theta <- c(rVir/(df)*exp(theta[1])/(1 + exp(theta[1])),
			exp(theta[2]),
			rVir/(df)/(1 + exp(theta[1])))
	}
  if(par == TRUE) {
	parmest <-optim(log(thetaini), m2LLgPAR, z = z, X = X, 
			xcoords = xcoords, ycoords = ycoords, 
			estMeth = 'REML', grpindx = grpindx)
	}
	m2LL <- parmest$value
	list(theta = theta, m2LL = m2LL, zlist = zlist, Xlist = Xlist, 
		xylist = xylist, qrlist = qrlist,	ViXlist = ViXlist, 
		parmest = parmest, betaHat = betaHat, Sxx = Sxx)
}
