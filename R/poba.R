#-------------------------------------------------------------------------------
#
#           poba
#
#-------------------------------------------------------------------------------

#' Prediction over block aggregrations for geostatistical models
#'
#' Prediction over block aggregrations for geostatistical models
#'
#' @param theta covariance parameters to be used in estimating fixed effects
#' @param z a vector of observed data
#' @param X design matrix for observed data
#' @param Xp design matrix for prediction data
#' @param xcoords x-coordinates for observed data.  Same dimension as z.
#' @param ycoords y-coordinates for observed data.  Same dimension as z.
#' @param xcoordsp x-coordinates for prediction data.  Same number of rows as Xp.
#' @param ycoordsp y-coordinates for prediction data.  Same number of rows as Xp.
#' @param predmeth prediction method.  'alldata' uses all data. 'NNdata' uses nearest neighbor method.
#' @param nNN number of nearest neighbors to use when predmeth='NNdata'.
#' @param use parallel processing?  TRUE or FALSE.
#'
#' @return a list of two scalars; one is the block prediction, the other the standard error for that prediction.
#'
#' @author Jay Ver Hoef
#' @export
poba <- function(theta, z, X, Xp, xcoords, ycoords, xcoordsp, ycoordsp,
  predmeth = 'alldata', nNN = 50, par = FALSE)
{
	N = length(xcoordsp)
	n = length(xcoords)
	if(predmeth == 'alldata') {
		dismat <- distGeoAni(xcoords, ycoords, xcoords, ycoords, 
			rotate = 90, range = theta[2], minorp = 1)	
		Voo <- theta[1]*corModelExponential(dismat) +
				theta[3]*diag(n)
		dismat <- distGeoAni(xcoords, ycoords, xcoordsp, 
				ycoordsp, rotate = 90, 
				range = theta[2], minorp = 1)
		Vop <- theta[1]*corModelExponential(dismat)
		dismat <- distGeoAni(xcoordsp, ycoordsp, xcoordsp, 
				ycoordsp, rotate = 90, 
				range = theta[2], minorp = 1)
		Vpp <- theta[1]*corModelExponential(dismat) +
			theta[3]*diag(N)
		Vi = solve(Voo)
		ViX <- Vi %*% X
		XViX <- crossprod(X, ViX)
		covb <- solve(XViX)
		Wall = Xp %*% covb %*% t(ViX) +
			t(Vop) %*% Vi %*%(diag(n) - X %*% covb %*% t(ViX))
		v_star = apply(Wall,2,mean)
		v = rep(1/N, times = N)
		Yp_bar_hat = v_star %*% z
		Yp_bar_hat_se = sqrt(t(v_star) %*% Voo %*% v_star - 
			2*t(v_star) %*% Vop %*% v + t(v) %*% Vpp %*% v)
	} 
	if(predmeth == 'NNdata') {
		dxy = as.matrix(cbind(xcoords,ycoords))
		pxy = as.matrix(cbind(xcoordsp,ycoordsp))
		nearxy = knn(data = dxy, query = pxy, k = nNN)
		v_star = rep(0, times = n)
		Sig_pp_v = rep(0, times = N)
		for(i in 1:N) {
			dismat = as.matrix(pdist(dxy[nearxy$nn.idx[i,],],pxy[i,]))/theta[2]
			Vpred <- theta[1]*corModelExponential(dismat)
			dismat = as.matrix(dist(dxy[nearxy$nn.idx[i,],]))/theta[2]
			V <- theta[1]*corModelExponential(dismat)
			diag(V) = diag(V) + theta[3]			
			Vi <- solve(V)
			Xi <- X[nearxy$nn.idx[i,],]
			zi <- z[nearxy$nn.idx[i,]]
			ViX <- Vi %*% Xi
			XViX <- crossprod(Xi, ViX)
			covb <- solve(XViX)
			bhatnoz = covb %*% t(ViX)
			v_star[nearxy$nn.idx[i,]] = v_star[nearxy$nn.idx[i,]] + 
				Xp[i,] %*% bhatnoz + t(Vpred) %*% Vi - 
				(t(Vpred) %*% ViX) %*% bhatnoz
			dismat = as.matrix(pdist(pxy[i,],pxy))/theta[2]
			Sig_pp_v[i] = theta[1]*mean(corModelExponential(dismat)) + 
			  #one term in the vector will have a nugget variance
				theta[3]/length(xcoordsp)		
		}
		v_star = v_star/N
		
		Sig_op_v = rep(0, times = n)
		Sig_oo_vstar = rep(0, times = n)
		for(i in 1:n) {
			dismat = as.matrix(pdist(dxy[i,],pxy))/theta[2]
			Sig_op_v[i] = theta[1]*mean(corModelExponential(dismat)) 
			dismat = as.matrix(pdist(dxy[i,],dxy))/theta[2]
      Sig_oo_vstar[i] = theta[1]*sum(corModelExponential(dismat)*v_star) +
				#ith term in the vector will have an additional nugget variance
				v_star[i]*theta[3]
		}	
		v = rep(1/N, times = N)
		Yp_bar_hat = v_star %*% z
		Yp_bar_hat_se = sqrt(sum(Sig_oo_vstar*v_star) - 
			2*sum(Sig_op_v*v_star) + 
			mean(Sig_pp_v))
	}
	list(block_pred = Yp_bar_hat, block_pred_se = Yp_bar_hat_se)
}

