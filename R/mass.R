#-------------------------------------------------------------------------------
#
#           mass
#
#-------------------------------------------------------------------------------

#' Makes subsamples for faster model fitting
#'
#' Makes subsamples for faster model fitting
#'
#' @param x vector of x-coordinates
#' @param y vector of y-coordinates
#' @param ngroups number of groups to make
#' @param gsize group size, only used for 'subsemble' method
#' @param ssmeth  subsampling method.  Valid methods are 'random' for 
#' randomly assigning locations to groups, 'compKmean' for using K-means
#' on coordinates to create compact spatial groups, 'zimmer' for 
#' using K-means to create compact spatial groups, but then randomly
#' re-assigning 10 percent in each group to get more distant distances,
#' 'subsemble' for ngroups of size gsize as in Barbian and Assuncao
#' (2017).
#'
#' @return a matrix where the first column is the row index of the x
#' (and y) vector, and the second column is group membership.
#'
#' @author Jay Ver Hoef
#' @export
mass = function(x, y, ngroups, gsize = NULL, ssmeth)
{
	used = 0
	n = length(x)
  if(ssmeth == 'random') {
		gsize = trunc(length(x)/ngroups)
		nxtra = length(x) - gsize*ngroups
		grpindx = kronecker((1:(ngroups-nxtra)),rep(1, times = gsize))
		if(nxtra > 0) grpindx = c(grpindx,
			kronecker(((ngroups - nxtra + 1):ngroups), 
			rep(1, times = gsize + 1)))
		grpindx = grpindx[order(runif(length(grpindx)))]
		grpindx = cbind(1:n,grpindx)
		colnames(grpindx) = c('row.indx', 'grpindx')
		used = 1
	}
	if(ssmeth == 'compKmean') {
		grpindx <- kmeans(as.data.frame(cbind(x,y)), ngroups,
			iter.max = 500)$cluster
		grpindx = cbind(1:n,grpindx)
		colnames(grpindx) = c('row.indx', 'grpindx')
		used = 1
	}
	if(ssmeth == 'zimmer') {
		grpindx <- kmeans(as.data.frame(cbind(x,y)), ngroups,
			iter.max = 500)$cluster
		gsamp = sample(1:length(grpindx), round(length(grpindx)/10))
			grpindx[gsamp] = 
			grpindx[gsamp][order(runif(round(length(grpindx)/10)))]
		grpindx = cbind(1:n,grpindx)
		colnames(grpindx) = c('row.indx', 'grpindx')
		used = 1
	}
	if(ssmeth == 'subsemble') {
		if(is.null(gsize)) stop(
			'Need group size (gsize argument) for subsemble method')
		xy = cbind(x,y)
		row.names(xy) = 1:n
		xytemp = xy
		grpindx = NULL
		for(i in 1:ngroups) {
			grpi = sample(1:(dim(xytemp)[1]),1)
			dxy = xytemp[grpi,, drop = FALSE]
			pxy = xytemp[!((1:dim(xytemp)[1]) %in% grpi),]
			nearxy = knn(data = pxy, query = dxy, k = gsize - 1)
			indxi = c(rownames(dxy), rownames(pxy)[nearxy$nn.idx])
			xytemp = xytemp[!(rownames(xytemp) %in% indxi),]
			grpindx = rbind(grpindx,
				cbind(as.integer(as.character(indxi)), 
				rep(i, times = length(indxi))))
		}
		colnames(grpindx) = c('row.indx', 'grpindx')
		used = 1
	}
	if(used == 0) stop(
		'No valid ssmeth was given (random, compKmean, zimmer, subsemble)')
  grpindx
}
