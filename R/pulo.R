#-------------------------------------------------------------------------------
#
#               pulo
#
#-------------------------------------------------------------------------------

#' Predictions for partitioned spatial linear models using local estimate of fixed effects
#'
#' Predictions for partitioned spatial linear models using local estimate of fixed effects
#'
#' @param theta covariance parameters to be used in estimating fixed effects
#' @param z a vector of observed data
#' @param X design matrix for observed data
#' @param Xp design matrix for prediction data
#' @param xcoords x-coordinates for observed data.  Same dimension as z.
#' @param ycoords y-coordinates for observed data.  Same dimension as z.
#' @param xcoordsp x-coordinates for prediction data.  Same number of rows as Xp.
#' @param ycoordsp y-coordinates for prediction data.  Same number of rows as Xp.
#' @param predmeth prediction method.  'claskrig' uses all data. 'nearnei' uses nearest neighbor method.
#' @param nNN number of nearest neighbors to use when predmeth='nearnei'.
#' @param use parallel processing?  TRUE or FALSE.
#'
#' @details These are the universal kriging equations found in Cressie (1993, pg. 154-155).
#'
#' @return \code{predict.splm} produces a data.frame with 2 columns for predictions, and prediction standard errors for missing response variable, one with the response variable name appended by \code{.pred} for predictions and by \code{.predSE} for the prediction standard errors.
#'
#' @references \cite{Cressie, N.A.C. (1993) Statistics for Spatial Data. Wiley.}
#'
#' @author Jay Ver Hoef
#' @export

pulo <- function(theta, z, X, Xp, xcoords, ycoords, xcoordsp, ycoordsp,
  predmeth = 'nearnei', nNN = 50, par = FALSE)
{
	
	if(predmeth == 'claskrig') {
		dismat <- distGeoAni(xcoords, ycoords, xcoords, ycoords, 
			rotate = 90, range = theta[2], minorp = 1)	
		covMat <- theta[1]*corModelExponential(dismat) +
				theta[3]*diag(dim(dismat)[1])
		dismat <- distGeoAni(xcoords, ycoords, xcoordsp, 
				ycoordsp, rotate = 90, 
				range = theta[2], minorp = 1)
		Vpred <- theta[1]*corModelExponential(dismat)
		Vi = solve(covMat)
		ViX <- Vi %*% X
		XViX <- crossprod(X, ViX)
		covb <- solve(XViX)
		bhat <- covb %*% crossprod(ViX, z)
		sill <- theta[1] + theta[3]
		preds <- matrix(NA, nrow = length(xcoordsp), ncol = 2)
		preds[,1] <- apply(as.vector((Vi %*% z)) * Vpred, 2, sum) +
			Xp %*% bhat - t(Vpred) %*% Vi %*% X %*% bhat	
		preds[,2] <- sqrt(rep(sill, times = nrow(Xp)) - 
			apply((Vi %*% Vpred) * Vpred, 2, sum) +
			apply((covb %*% t(Xp)) * t(Xp), 2, sum) -
			2*apply((covb %*% t(Xp)) * (t(X) %*% Vi %*% Vpred), 2, sum) +
			apply((covb %*% t(X) %*% Vi %*% Vpred) * (t(X) %*% Vi %*% Vpred), 2, sum))
		colnames(preds) <- c(paste('z', ".pred", sep = ""), 
			paste('z', ".predSE", sep = ""))

	}
	if(predmeth == 'nearnei') {
		dxy = as.matrix(cbind(xcoords,ycoords))
		pxy = as.matrix(cbind(xcoordsp,ycoordsp))
		nearxy = knn(data = dxy, query = pxy, k = nNN)
		if(par == FALSE) {
		  preds <- matrix(NA, nrow = length(xcoordsp), ncol = 2)
			for(i in 1:length(xcoordsp)) {
				dismat <- distGeoAni(xcoords[nearxy$nn.idx[i,]], 
					ycoords[nearxy$nn.idx[i,]], xcoordsp[i], 
					ycoordsp[i], rotate = 90, 
					range = theta[2], minorp = 1)
				Vpred <- theta[1]*corModelExponential(dismat)
				dismat <- distGeoAni(xcoords[nearxy$nn.idx[i,]], 
					ycoords[nearxy$nn.idx[i,]], xcoords[nearxy$nn.idx[i,]], 
					ycoords[nearxy$nn.idx[i,]], rotate = 90, 
					range = theta[2], minorp = 1)
				V <- theta[1]*corModelExponential(dismat)
				diag(V) = diag(V) + theta[3]					
				Vi <- solve(V)
				Xi <- X[nearxy$nn.idx[i,],]
				zi <- z[nearxy$nn.idx[i,]]
				ViX <- Vi %*% Xi
				XViX <- crossprod(Xi, ViX)
				covb <- solve(XViX)
				bhat <- covb %*% crossprod(ViX, zi)
				sill <- theta[1] + theta[3]
				preds[i,1] = sum(as.vector((Vi %*% zi)) * Vpred) +
					Xp[i,] %*% bhat - t(Vpred) %*% Vi %*% Xi %*% bhat
				preds[i,2] = sqrt(sill -
					sum((Vi %*% Vpred) * Vpred) +
					sum((covb %*% Xp[i,]) * t(Xp[i,, drop = F])) -
					2*sum((covb %*% Xp[i,]) * (t(Xi) %*% Vi %*% Vpred)) +
					sum((covb %*% t(Xi) %*% Vi %*% Vpred) * 
					(t(Xi) %*% Vi %*% Vpred)))
			}
		}
		if(par == TRUE) {
			NNloop = function(i, theta, X, Xp, xcoords, ycoords,
				xcoordsp, ycoordsp, nearxy) 
			{
		  predi <- matrix(NA, nrow = 1, ncol = 2)
				dismat <- distGeoAni(xcoords[nearxy$nn.idx[i,]], 
					ycoords[nearxy$nn.idx[i,]], xcoordsp[i], 
					ycoordsp[i], rotate = 90, 
					range = theta[2], minorp = 1)
				Vpred <- theta[1]*corModelExponential(dismat)
				dismat <- distGeoAni(xcoords[nearxy$nn.idx[i,]], 
					ycoords[nearxy$nn.idx[i,]], xcoords[nearxy$nn.idx[i,]], 
					ycoords[nearxy$nn.idx[i,]], rotate = 90, 
					range = theta[2], minorp = 1)
				V <- theta[1]*corModelExponential(dismat)
				diag(V) = diag(V) + theta[3]					
				Vi <- solve(V)
				Xi <- X[nearxy$nn.idx[i,],]
				zi <- z[nearxy$nn.idx[i,]]
				ViX <- Vi %*% Xi
				XViX <- crossprod(Xi, ViX)
				covb <- solve(XViX)
				bhat <- covb %*% crossprod(ViX, zi)
				sill <- theta[1] + theta[3]
				predi[1,1] = sum(as.vector((Vi %*% zi)) * Vpred) +
					Xp[i,] %*% bhat - t(Vpred) %*% Vi %*% Xi %*% bhat
				predi[1,2] = sqrt(sill -
					sum((Vi %*% Vpred) * Vpred) +
					sum((covb %*% Xp[i,]) * t(Xp[i,, drop = F])) -
					2*sum((covb %*% Xp[i,]) * (t(Xi) %*% Vi %*% Vpred)) +
					sum((covb %*% t(Xi) %*% Vi %*% Vpred) * 
					(t(Xi) %*% Vi %*% Vpred)))
				predi
			}
			NNout = foreach(j=1:length(xcoordsp)) %dopar% {
				NNloop(j, theta, X, Xp, xcoords, ycoords, xcoordsp, 
					ycoordsp, nearxy)
			}
			preds = t(matrix(unlist(NNout), nrow = 2))
		}
	}
	preds
}


