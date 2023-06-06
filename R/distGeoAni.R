#-------------------------------------------------------------------------------
#
#           distGeoAni
#
#-------------------------------------------------------------------------------

#' Compute anistropy corrected distance between two sets of data
#'
#' computes anistropy corrected distance between two sets of data
#'
#' @param xrow vector with x-coordinates that will form rows of distance matrix 
#' @param yrow vector with y-coordinates that will form rows of distance matrix, must be of same length as xrow
#' @param xcol vector with x-coordinates that will form columns of distance matrix 
#' @param ycol vector with y-coordinates that will form columns of distance matrix, must be of same length as xcol
#' @param rotate rotation of anisotropic axes, default = 0
#' @param range range of autocorrelation model, default = 1
#' @param minorp proportion of range in x direction to that of y direction for unrotated anisotropic model, default = 1
#'
#' @return matrix of distances
#'
#' @author Jay Ver Hoef
distGeoAni <- function(xrow, yrow, xcol, ycol, rotate = 0, range = 1, minorp = 1)
{
	# total number of observations for each set of coordinates  
		n.rows <- length(xrow)
		n.cols <- length(xcol)
	# expand all x-coordinates
		sxr <- matrix(xrow, ncol = 1) %*% 
			matrix(rep(1,times = n.cols), nrow = 1)
		sxc <- matrix(rep(1,times = n.rows), ncol = 1) %*% 
			matrix(xcol, nrow = 1)
		syr <- matrix(yrow,ncol = 1) %*% 
			matrix(rep(1,times = n.cols), nrow = 1)
		syc <- matrix(rep(1,times = n.rows), ncol = 1) %*% 
			matrix(ycol, nrow = 1)
	# find difference in coordinates between all pairwise locations
		sxdif <- sxr - sxc
		sydif <- syr - syc
	# rotate coordinates
		newx <- cos(rotate*.0174533)*sxdif - sin(rotate*.0174533)*sydif
		newy <- sin(rotate*.0174533)*sxdif + cos(rotate*.0174533)*sydif
	# scale coordinates by minor and major axes */
		newx <- newx/(range*minorp)
		newy <- newy/range
	# compute distance for the scaled and rotated coordinates */
		sqrt(newx^2 + newy^2)
}

