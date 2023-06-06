# exponential autocorrelation model
corModelExponential <- function(distance.matrix)
{
	exp(-distance.matrix) 
}

# spherical autocorrelation model
corModelSpherical <- function(distance.matrix)
{
	CovMat <- (1 - 1.5*distance.matrix + 0.5*distance.matrix^3)
	CovMat[distance.matrix > 1] <- 0
	CovMat
}

#-------------------------------------------------------------------------------
#
#          cope_simple
#
#-------------------------------------------------------------------------------

cope_simple <- function(formula, data, x_column, y_column,
	corModel = corModelExponential,
  thetaini = c(2,1,0.02), subSampCol = NULL) 
{

	# ----------------------------------------------------------------------------
	# prepare data
	# ----------------------------------------------------------------------------

	# extract information from formula and data to create design matrix and
	# response variable
	trms <- terms(formula, data = data)
	respCol <- as.character(as.list(attr(trms,"variables")[-1]))[1]
	covList <- attr(trms,"term.labels")
	# design matrix
	X <- model.matrix(formula, data)
	# response variable
	z <- as.matrix(data[, respCol])
	# vectors of spatial coordinates
	xcoords <- data[,x_column]
	ycoords <- data[,y_column]

	# total number of rows in data set
	n <- length(z)
	# number of linearly independent columns in X
	p <- sum(svd(X)$d>1e-10)
	# grouping variable for data partitioning/covariance matrix blocking
	grpindx = data[,subSampCol]

	# emply lists for partitioned objects
  xylist = vector(mode = "list", length = max(grpindx)) 
  zlist = vector(mode = "list", length = max(grpindx)) 
  Xlist = vector(mode = "list", length = max(grpindx)) 
  # put useful items in the lists
  for(i in 1:max(grpindx)) {
		# spatial coordinates by group
    xylist[[i]] = cbind(xcoords[grpindx == i], ycoords[grpindx == i])
    # response variable by group
    zlist[[i]] = z[grpindx == i]
    # design matrix by group
    Xlist[[i]] = X[grpindx == i,]
	}

	# ----------------------------------------------------------------------------
	# estimate covariance parameters
	# ----------------------------------------------------------------------------

	# optimize minus 2 times the loglikelihood
	# all parameters will be optimized on the log scale
	parmest = optim(log(thetaini), m2LLg_simple, z = zlist, X = Xlist, 
			xycoords = xylist, corModel = corModel)
	theta <- exp(parmest$par)
	theta = c(theta[1], theta[3], theta[2])
	names(theta) = c('de', 'ie', 'range')
	theta
}
