#-------------------------------------------------------------------------------
#
#           geostatSim
#
#-------------------------------------------------------------------------------

#' Simulate geostatistical data on set of given locations
#'
#' simulates geostatistical data on set of given locations
#'
#' @param loc.data data.frame with x- and y-coordinates of locations for simulated data 
#' @param xcol name of the column in loc.data with x-coordinates, default is "x"
#' @param ycol name of the column loc.data with y-coordinates, default is "y"
#' @param parsil partial sill of autocorrelation model, default = 1
#' @param range range of autocorrelation model, default = 1
#' @param nugget range of autocorrelation model, default = 0
#' @param minorp proportion of range in x direction to that of y direction for unrotated anisotropic model, default = 1
#' @param rotate rotation of anisotropic axes, default = 90
#' @param CorModel autocorrelation model, default = "Exponential".  Other possibilities are "Spherical".
#'
#' @return data.frame of three columns, the original loc.data appended with a 3rd column of simulated geostatistical data
#'
#' @author Jay Ver Hoef
#' @export
geostatSim <- function(loc.data, xcol = "x", ycol = "y", 
	parsil = 1, range = 1, nugget = 0,
	minorp = 1, rotate = 90, extrap = NULL,
	CorModel = "Exponential")
{
	xcoord <- loc.data[,xcol]
	ycoord <- loc.data[,ycol]
	n <- length(xcoord)
	dismat <- distGeoAni(xcoord, ycoord, xcoord, ycoord, rotate, range, minorp)
	# compute correlation matrix for scaled distance matrix
		if(CorModel == "Exponential") CovMat <- corModelExponential(dismat)
#		if(CorModel == "ExpRadon2") CovMat <- CorModel.ExpRadon2(dismat)
#		if(CorModel == "ExpRadon4") CovMat <- CorModel.ExpRadon4(dismat)
#		if(CorModel == "Gaussian") CovMat <- CorModel.Gaussian(dismat)
#		if(CorModel == "Stable") CovMat <- CorModel.Stable(dismat, extrap)
#		if(CorModel == "RationalQuad") CovMat <- CorModel.RationalQuad(dismat)
#		if(CorModel == "CauchyGrav") CovMat <- CorModel.CauchyGrav(dismat)
#		if(CorModel == "CauchyMag") CovMat <- CorModel.CauchyMag(dismat)
#		if(CorModel == "Cauchy") CovMat <- CorModel.Cauchy(dismat, extrap)
#		if(CorModel == "Circular") CovMat <- CorModel.Circular(dismat)
		if(CorModel == "Spherical") CovMat <- corModelSpherical(dismat)
#		if(CorModel == "Cubic") CovMat <- CorModel.Cubic(dismat)
#		if(CorModel == "Penta") CovMat <- CorModel.Penta(dismat)
#		if(CorModel == "CardinalSine") CovMat <- CorModel.CardinalSine(dismat)
#		if(CorModel == "BesselK") CovMat <- CorModel.BesselK(dismat, extrap)
#		if(CorModel == "BesselJ") CovMat <- CorModel.BesselJ(dismat, extrap)
	CovMat <- parsil*CovMat + diag(nugget, nrow = n, ncol = n)
	data.frame(loc.data, z = t(chol(CovMat))%*%rnorm(n))
}

