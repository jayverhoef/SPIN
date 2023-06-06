#Rscript sim_scen1_500pts.R

library('spNNGP')
library('spPart')
library('viridis')
library('classInt')
library('nabor')
library('slm')

#library('doParallel')
#registerDoParallel(cores=7)

coeftable = function(fefeout)
{
	EVplusVE = NULL
	analytic = cbind(fefeout$bhat, sqrt(diag(fefeout$covb)))
	if(!is.null(fefeout$bhatList))
	EVplusVE = cbind(apply(simplify2array(fefeout$bhatList), 1:2, mean)[,1],
	sqrt((apply(simplify2array(fefeout$bhatList), 1:2, mean)[,2]/length(fefeout$bhatList) + 
		apply(simplify2array(fefeout$bhatList), 1:2, var)[,1]/length(fefeout$bhatList))/2))
	list(analytic = analytic, EVplusVE = EVplusVE)
}

set.seed(7011)
nsim = 1000
scen1_1000pt_NNGP_est1 = matrix(NA, nrow = nsim, ncol = 3)
scen1_1000pt_NNGP_est1se = matrix(NA, nrow = nsim, ncol = 3)
scen1_1000pt_NNGP_est2 = matrix(NA, nrow = nsim, ncol = 3)
scen1_1000pt_NNGP_est2se = matrix(NA, nrow = nsim, ncol = 3)
scen1_1000pt_NNGP_RMSPE = matrix(NA, nrow = nsim, ncol = 3)
scen1_1000pt_NNGP_PI90 = matrix(NA, nrow = nsim, ncol = 3)
scen1_1000pt_NNGP_time = matrix(NA, nrow = nsim, ncol = 3)
# scen1_1000pt_NNGP_blk_pred = matrix(NA, nrow = nsim, ncol = 3)
# scen1_1000pt_NNGP_blk_se = matrix(NA, nrow = nsim, ncol = 3)
iter = 1
for(iter in 1:nsim) {
#-----------------------------------------------------------------------
#                       SIMULATE
#-----------------------------------------------------------------------
  cat("\n", "Iteration: ", iter)
  cat("\n", "Simulate")
  startiter = Sys.time()
	n = 1000
	x = runif(n)
	y = runif(n)
	xysyst = pointSimSyst(nrow = 40, ncol = 40, lower.x.lim = 0, upper.x.lim = 1, 
			lower.y.lim = 0, upper.y.lim = 1)
	x = c(x, xysyst[,1])
	y = c(y, xysyst[,2])

	xyz1 = geostatSim(loc.data = data.frame(x = x, y = y), xcol = "x", 
		ycol = "y", parsil = 10, range = .5, nugget = .01, 
		CorModel = 'Spherical')
	z1std = (xyz1[,3] - mean(xyz1[,3]))/sqrt(var(xyz1[,3]))
	xyz2 = geostatSim(loc.data = data.frame(x = x, y = y), xcol = "x", 
		ycol = "y", parsil = 10, range = .5, nugget = .01, 
		CorModel = 'Spherical')
	z2std = (xyz2[,3] - mean(xyz2[,3]))/sqrt(var(xyz2[,3]))

	X1 = as.matrix(cbind(rep(1, times = n + 1600), rnorm(n + 1600), z2std))
	Z1 =  1*X1[,1]+ 1*X1[,2] + 1*X1[,3] + z1std + .1*rnorm(n + 1600)

	xp = x[(n+1):(n + 1600)]
	x = x[1:n]
	yp = y[(n+1):(n + 1600)]
	y = y[1:n]
	X1p = X1[(n+1):(n + 1600),]
	X1 = X1[1:n,]
	Z1p = Z1[(n+1):(n + 1600)]
	Z1 = Z1[1:n]


	gcomp = mass(x, y, ngroups = round(1000/50), ssmeth = 'compKmean')
		
	d2 = data.frame(x = x, y = y, X1 = X1[,2], X2 =  X1[,3], z = Z1, 
  gcomp = gcomp[,2])

	dp = data.frame(x = xp, y = yp, X1 = X1p[,2], X2 = X1p[,3], z = Z1p)

#-----------------------------------------------------------------------
#                      NNGP
#-----------------------------------------------------------------------
cat("\n")
  phi.range = exp((-12:12)/4)
  phi.range = phi.range - min(phi.range) + .001
  alpha.range = exp((-12:12)/4)
  alpha.range = alpha.range - min(alpha.range) + .001
  g = length(phi.range)

  kronecker(rep(1, times = g), phi.range)
  theta.alpha <- cbind(
    kronecker(rep(1, times = g), phi.range),
    kronecker(alpha.range,rep(1, times = g))
  )
  colnames(theta.alpha) <- c("phi", "alpha")
  sigma.sq.IG <- c(2, 10)

	starttime = Sys.time()
  NNGPout_rmspe <- spConjNNGP(Z1 ~ X1-1, coords=cbind(x,y), n.neighbors = 10,
                  X.0 = X1p, coords.0 = cbind(xp,yp),
                  k.fold = 5, score.rule = "rmspe",
                  n.omp.threads = 1,
                  theta.alpha = theta.alpha, sigma.sq.IG = sigma.sq.IG, 
                  cov.model = 'exponential')
	stoptime = Sys.time()
	NNGP_time = difftime(stoptime, starttime, units="secs")
  

#-----------------------------------------------------------------------
#                    Full Covariance Matrix Methods
#-----------------------------------------------------------------------

	starttime = Sys.time()
  cat("\n", "Full Cov")
	cope_all_slm = slm:::cope(z ~ X1 + X2, data = d2, x_column = 'x', y_column = 'y')
  fefe_all_slm = slm:::fefe(cope_all_slm)
  betaHat_all = fefe_all_slm$bhat
  covb_all = fefe_all_slm$covb
	cope_all_pulo_krig = spPart:::pulo(cope_all_slm$theta, Z1, X1, X1p, x, y, xp, yp,
		predmeth = 'claskrig', par = FALSE)
	stoptime = Sys.time()
	FullCov_time = difftime(stoptime, starttime, units="secs")

#-----------------------------------------------------------------------
#                  Partitioned Matrix using 42
#-----------------------------------------------------------------------

  cat("\n", "Partitioned Cov")
  
	starttime = Sys.time()
	cope_comp = spPart:::cope(z ~ X1 + X2, data = d2, x_column = 'x', y_column = 'y', 
		subSampCol = 'gcomp', thetaini = c(3,3))
	cope_comp_fefe_comp = spPart:::fefe(theta = cope_comp$theta, z = Z1, 
		X = X1, subsampindx = gcomp[,2], xcoords = x, ycoords = y)
  betaHat_comp = cope_comp_fefe_comp$bhat
  covb_comp = cope_comp_fefe_comp$covb
	cope_comp_fefe_comp_puloNN = spPart:::pulo1(theta = cope_comp$theta, 
    betahat = betaHat_comp, covb = covb_comp, z = Z1, X = X1, Xp = X1p, 
    xcoords = x, ycoords = y, xcoordsp = xp, ycoordsp = yp, 
    nNN = 50, predmeth = 'nearnei', par = FALSE)
	stoptime = Sys.time()
	PartCov_time = difftime(stoptime, starttime, units="secs")

#-----------------------------------------------------------------------
#                  Store Results
#-----------------------------------------------------------------------


	scen1_1000pt_NNGP_est1[iter,1] = betaHat_all[2]
 	scen1_1000pt_NNGP_est1[iter,2] = betaHat_comp[2]
  scen1_1000pt_NNGP_est1[iter,3] = NNGPout_rmspe$beta.hat[2]
	scen1_1000pt_NNGP_est1se[iter,1] = sqrt(covb_all[2,2])
	scen1_1000pt_NNGP_est1se[iter,2] = sqrt(covb_comp[2,2])
 	scen1_1000pt_NNGP_est1se[iter,3] = sqrt(diag(NNGPout_rmspe$beta.var)[2])
	scen1_1000pt_NNGP_est2[iter,1] = betaHat_all[3]
 	scen1_1000pt_NNGP_est2[iter,2] = betaHat_comp[3]
  scen1_1000pt_NNGP_est2[iter,3] = NNGPout_rmspe$beta.hat[3]
	scen1_1000pt_NNGP_est2se[iter,1] = sqrt(covb_all[3,3])
	scen1_1000pt_NNGP_est2se[iter,2] = sqrt(covb_comp[3,3])
 	scen1_1000pt_NNGP_est2se[iter,3] = sqrt(diag(NNGPout_rmspe$beta.var)[3])
	scen1_1000pt_NNGP_RMSPE[iter,1] = mean((Z1p - cope_all_pulo_krig[,1])^2)
	scen1_1000pt_NNGP_RMSPE[iter,2] = mean((Z1p - cope_comp_fefe_comp_puloNN[,1])^2)
  scen1_1000pt_NNGP_RMSPE[iter,3] = mean((Z1p - NNGPout_rmspe$y.0.hat)^2) 
	scen1_1000pt_NNGP_PI90[iter,1] = mean(cope_all_pulo_krig[,1] - 
		qnorm(.95)*cope_all_pulo_krig[,2] < Z1p & Z1p < cope_all_pulo_krig[,1] + 
		qnorm(.95)*cope_all_pulo_krig[,2])
	scen1_1000pt_NNGP_PI90[iter,2] = mean(cope_comp_fefe_comp_puloNN[,1] - 
		qnorm(.95)*cope_comp_fefe_comp_puloNN[,2] < Z1p & Z1p < cope_comp_fefe_comp_puloNN[,1] + 
		qnorm(.95)*cope_comp_fefe_comp_puloNN[,2])
	scen1_1000pt_NNGP_PI90[iter,3] = mean(NNGPout_rmspe$y.0.hat - 
		qnorm(.95)*sqrt(NNGPout_rmspe$y.0.hat.var) < Z1p & Z1p < NNGPout_rmspe$y.0.hat + 
		qnorm(.95)*sqrt(NNGPout_rmspe$y.0.hat.var))
  scen1_1000pt_NNGP_time[iter,3] = NNGP_time
  scen1_1000pt_NNGP_time[iter,1] = FullCov_time
  scen1_1000pt_NNGP_time[iter,2] = PartCov_time
  

	
	stopiter = Sys.time()
	itertime = difftime(stopiter, startiter, units="mins")
  cat("\n", "iteration time: ", itertime, "\n")
}

basepath = '/media/jay/data/desktop_data/2019_papers/fastREML/spPart_package/spPart/data/'
save(scen1_1000pt_NNGP_est1, file = paste0(basepath,'scen1_1000pt_NNGP_est1.rda'))
save(scen1_1000pt_NNGP_est1se, file = paste0(basepath,'scen1_1000pt_NNGP_est1se.rda'))
save(scen1_1000pt_NNGP_est2, file = paste0(basepath,'scen1_1000pt_NNGP_est2.rda'))
save(scen1_1000pt_NNGP_est2se, file = paste0(basepath,'scen1_1000pt_NNGP_est2se.rda'))
save(scen1_1000pt_NNGP_RMSPE, file = paste0(basepath,'scen1_1000pt_NNGP_RMSPE.rda'))
save(scen1_1000pt_NNGP_PI90, file = paste0(basepath,'scen1_1000pt_NNGP_PI90.rda'))
save(scen1_1000pt_NNGP_time, file = paste0(basepath,'scen1_1000pt_NNGP_time.rda'))


