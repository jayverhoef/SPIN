#Rscript sim_scen1_500pts.R

basepath = '/home/xverhoef/data/2019_papers/fastREML/testingScripts/'
source(paste0(basepath,'distGeoAni.R'))
source(paste0(basepath,'geostatSim.R'))
source(paste0(basepath,'sim2DSinSurf.R'))
source(paste0(basepath,'pointSimSyst.R'))
source(paste0(basepath,'corModels.R'))
source(paste0(basepath,'mass.R'))
source(paste0(basepath,'cope.R'))
source(paste0(basepath,'fefe.R'))
source(paste0(basepath,'pulo.R'))
source(paste0(basepath,'poba.R'))
source(paste0(basepath,'m2LLg.R'))
source(paste0(basepath,'m2LLgPAR.R'))

library('spPart')
library('viridis')
library('classInt')
library('nabor')

library('doParallel')
registerDoParallel(cores=7)

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

set.seed(2100)
nsim = 2000
scen1_500pts_est = matrix(NA, nrow = nsim, ncol = 9)
scen1_500pts_estse = matrix(NA, nrow = nsim, ncol = 9)
scen1_500pts_RMSPE = matrix(NA, nrow = nsim, ncol = 7)
scen1_500pts_PI90 = matrix(NA, nrow = nsim, ncol = 7)
scen1_500pts_blk_pred = matrix(NA, nrow = nsim, ncol = 8)
scen1_500pts_blk_se = matrix(NA, nrow = nsim, ncol = 7)

for(iter in 1:nsim) {
#-----------------------------------------------------------------------
#                       SIMULATE
#-----------------------------------------------------------------------
  cat("\n", "Simulation: ", iter)
  cat("\n", "Simulate")
  startiter = Sys.time()
	n = 500
	x = runif(n)
	y = runif(n)
	xysyst = pointSimSyst(nrow = 40, ncol = 40, lower.x.lim = 0, upper.x.lim = 1, 
			lower.y.lim = 0, upper.y.lim = 1)
	x = c(x, xysyst[,1])
	y = c(y, xysyst[,2])

	xyz1 = geostatSim(loc.data = data.frame(x = x, y = y), xcol = "x", 
		ycol = "y", parsil = 10, range = .5, nugget = .01, 
		CorModel = 'Spherical')
	zstd = (xyz1[,3] - mean(xyz1[,3]))/sqrt(var(xyz1[,3]))
	xyz2 = geostatSim(loc.data = data.frame(x = x, y = y), xcol = "x", 
		ycol = "y", parsil = 10, range = .5, nugget = .01, 
		CorModel = 'Spherical')
	bowl = -((x-0.5)^2 + (y - 0.5)^2)
	bowl = (bowl - mean(bowl))/sqrt(var(bowl))

	X1 = as.matrix(cbind(rep(1, times = n + 1600), rnorm(n + 1600)))
	X2 = (xyz2[,3] - mean(xyz2[,3]))/sqrt(var(xyz2[,3]))
	Z1 =  1*X1[,1]+ 1*X1[,2] + zstd
	Z2 =  1*X1[,1]+ 1*X1[,2] + zstd + bowl
	Z3 =  1*X1[,1]+ 1*X2 + zstd
	Z4 =  1*X1[,1]+ 1*X2 + zstd + bowl

	xp = x[(n+1):(n + 1600)]
	x = x[1:n]
	yp = y[(n+1):(n + 1600)]
	y = y[1:n]
	X1p = X1[(n+1):(n + 1600),]
	X1 = X1[1:n,]
	X2p = X2[(n+1):(n + 1600)]
	X2 = X2[1:n]
	Z1p = Z1[(n+1):(n + 1600)]
	Z1 = Z1[1:n]
	Z2p = Z2[(n+1):(n + 1600)]
	Z2 = Z2[1:n]
	Z3p = Z3[(n+1):(n + 1600)]
	Z3 = Z3[1:n]
	Z4p = Z4[(n+1):(n + 1600)]
	Z4 = Z4[1:n]

	# gi1 = mass(x, y, ngroups = 5, gsize = 60, ssmeth = 'subsemble')

	gi2 = mass(x, y, ngroups = 10, ssmeth = 'random')

	gi3 = mass(x, y, ngroups = 10, ssmeth = 'compKmean')
		
	gi4 = mass(x, y, ngroups = 10, ssmeth = 'zimmer')

	#d1 = data.frame(x = x[gi1[,1]], y = y[gi1[,1]], X1 = X1[gi1[,1],2],
	#	X2 = X2[gi1[,1]], z1 = Z1[gi1[,1]], z2 = Z2[gi1[,1]], z3 = Z3[gi1[,1]],
	#	z4 = Z4[gi1[,1]])

	d2 = data.frame(x = x, y = y, X1 = X1[,2], X2 = X2, z1 = Z1, z2 = Z2, 
		z3 = Z3, z4 = Z4, grpindx.rand = gi2[,2], grpindx.comp = gi3[,2],
		grpindx.zimm = gi4[,2], grpindx.one = rep(1, times = n))

	dp = data.frame(x = xp, y = yp, X1 = X1p[,2], X2 = X2p, z1 = Z1p, 
		z2 = Z2p, z3 = Z3p, z4 = Z4p)

#-----------------------------------------------------------------------
#                       COPE
#-----------------------------------------------------------------------
  cat("\n", "cope")
	cope_all = cope(z1 ~ X1, data = d2, x_column = 'x', y_column = 'y', 
		subSampCol = 'grpindx.one')

	cope_rand = cope(z1 ~ X1, data = d2, x_column = 'x', y_column = 'y', 
		subSampCol = 'grpindx.rand')

	cope_comp = cope(z1 ~ X1, data = d2, x_column = 'x', y_column = 'y', 
		subSampCol = 'grpindx.comp')

	cope_zimm = cope(z1 ~ X1, data = d2, x_column = 'x', y_column = 'y', 
		subSampCol = 'grpindx.zimm')

#-----------------------------------------------------------------------
#                      FEFE
#-----------------------------------------------------------------------

  cat("\n", "fefe")

	cope_all_fefe_all = fefe(theta = cope_all$theta, z = Z1, 
		X = X1, xcoords = x, ycoords = y)

	cope_rand_fefe_all = fefe(theta = cope_rand$theta, z = Z1, 
		X = X1, xcoords = x, ycoords = y)

	cope_rand_fefe_rand = fefe(theta = cope_rand$theta, z = Z1, 
		X = X1, subsampindx = gi2[,2], xcoords = x, ycoords = y)

	cope_comp_fefe_all = fefe(theta = cope_comp$theta, z = Z1, 
		X = X1, xcoords = x, ycoords = y)

	cope_comp_fefe_rand = fefe(theta = cope_comp$theta, z = Z1, 
		X = X1, subsampindx = gi2[,2], xcoords = x, ycoords = y)

	cope_comp_fefe_comp = fefe(theta = cope_comp$theta, z = Z1, 
		X = X1, subsampindx = gi3[,2], xcoords = x, ycoords = y)

	cope_zimm_fefe_all = fefe(theta = cope_zimm$theta, z = Z1, 
		X = X1, xcoords = x, ycoords = y)

	cope_zimm_fefe_rand = fefe(theta = cope_zimm$theta, z = Z1, 
		X = X1, subsampindx = gi2[,2], xcoords = x, ycoords = y)

	cope_zimm_fefe_zimm = fefe(theta = cope_zimm$theta, z = Z1, 
		X = X1, subsampindx = gi4[,2], xcoords = x, ycoords = y)

	scen1_500pts_est[iter,1] = coeftable(cope_all_fefe_all)$analytic[2,1]
	scen1_500pts_estse[iter,1] = coeftable(cope_all_fefe_all)$analytic[2,2]
	scen1_500pts_est[iter,2] = coeftable(cope_rand_fefe_all)$analytic[2,1]
	scen1_500pts_estse[iter,2] = coeftable(cope_rand_fefe_all)$analytic[2,2]
	scen1_500pts_est[iter,3] = coeftable(cope_rand_fefe_rand)$analytic[2,1]
	scen1_500pts_estse[iter,3] = coeftable(cope_rand_fefe_rand)$analytic[2,2]
	scen1_500pts_est[iter,4] = coeftable(cope_rand_fefe_all)$analytic[2,1]
	scen1_500pts_estse[iter,4] = coeftable(cope_rand_fefe_all)$analytic[2,2]
	scen1_500pts_est[iter,5] = coeftable(cope_comp_fefe_rand)$analytic[2,1]
	scen1_500pts_estse[iter,5] = coeftable(cope_comp_fefe_rand)$analytic[2,2]
	scen1_500pts_est[iter,6] = coeftable(cope_comp_fefe_comp)$analytic[2,1]
	scen1_500pts_estse[iter,6] = coeftable(cope_comp_fefe_comp)$analytic[2,2]
	scen1_500pts_est[iter,7] = coeftable(cope_zimm_fefe_all)$analytic[2,1]
	scen1_500pts_estse[iter,7] = coeftable(cope_zimm_fefe_all)$analytic[2,2]
	scen1_500pts_est[iter,8] = coeftable(cope_zimm_fefe_rand)$analytic[2,1]
	scen1_500pts_estse[iter,8] = coeftable(cope_zimm_fefe_rand)$analytic[2,2]
	scen1_500pts_est[iter,9] = coeftable(cope_zimm_fefe_zimm)$analytic[2,1]
	scen1_500pts_estse[iter,9] = coeftable(cope_zimm_fefe_zimm)$analytic[2,2]

#-----------------------------------------------------------------------
#                      PULO
#-----------------------------------------------------------------------
  cat("\n", "pulo")
  
	cope_all_pulo_krig = pulo(cope_all$theta, Z1, X1, X1p, x, y, xp, yp,
		predmeth = 'claskrig', par = FALSE)

	cope_rand_pulo_krig = pulo(cope_rand$theta, Z1, X1, X1p, x, y, xp, yp,
		predmeth = 'claskrig', par = FALSE)

	cope_rand_pulo_nene = pulo(cope_rand$theta, Z1, X1, X1p, x, y, xp, yp,
		predmeth = 'nearnei', par = FALSE)

	cope_comp_pulo_krig = pulo(cope_comp$theta, Z1, X1, X1p, x, y, xp, yp,
		predmeth = 'claskrig', par = FALSE)

	cope_comp_pulo_nene = pulo(cope_comp$theta, Z1, X1, X1p, x, y, xp, yp,
		predmeth = 'nearnei', par = FALSE)

	cope_zimm_pulo_krig = pulo(cope_zimm$theta, Z1, X1, X1p, x, y, xp, yp,
		predmeth = 'claskrig', par = FALSE)
		
	cope_zimm_pulo_nene = pulo(cope_zimm$theta, Z1, X1, X1p, x, y, xp, yp,
		predmeth = 'nearnei', par = FALSE)

		
	scen1_500pts_RMSPE[iter,1] = mean((Z1p - cope_all_pulo_krig[,1])^2)
	scen1_500pts_PI90[iter,1] = mean(cope_all_pulo_krig[,1] - 
		qnorm(.95)*cope_all_pulo_krig[,2] < Z1p & Z1p < cope_all_pulo_krig[,1] + 
		qnorm(.95)*cope_all_pulo_krig[,2])
	scen1_500pts_RMSPE[iter,2] = mean((Z1p - cope_rand_pulo_krig[,1])^2)
	scen1_500pts_PI90[iter,2] = mean(cope_rand_pulo_krig[,1] - 
		qnorm(.95)*cope_rand_pulo_krig[,2] < Z1p & Z1p < cope_rand_pulo_krig[,1] + 
		qnorm(.95)*cope_rand_pulo_krig[,2])
	scen1_500pts_RMSPE[iter,3] = mean((Z1p - cope_rand_pulo_nene[,1])^2)
	scen1_500pts_PI90[iter,3] = mean(cope_rand_pulo_nene[,1] - 
		qnorm(.95)*cope_rand_pulo_nene[,2] < Z1p & Z1p < cope_rand_pulo_nene[,1] + 
		qnorm(.95)*cope_rand_pulo_nene[,2])
	scen1_500pts_RMSPE[iter,4] = mean((Z1p - cope_comp_pulo_krig[,1])^2)
	scen1_500pts_PI90[iter,4] = mean(cope_comp_pulo_krig[,1] - 
		qnorm(.95)*cope_comp_pulo_krig[,2] < Z1p & Z1p < cope_comp_pulo_krig[,1] + 
		qnorm(.95)*cope_comp_pulo_krig[,2])
	scen1_500pts_RMSPE[iter,5] = mean((Z1p - cope_comp_pulo_nene[,1])^2)
	scen1_500pts_PI90[iter,5] = mean(cope_comp_pulo_nene[,1] - 
		qnorm(.95)*cope_comp_pulo_nene[,2] < Z1p & Z1p < cope_comp_pulo_nene[,1] + 
		qnorm(.95)*cope_comp_pulo_nene[,2])
	scen1_500pts_RMSPE[iter,6] = mean((Z1p - cope_zimm_pulo_krig[,1])^2)
	scen1_500pts_PI90[iter,6] = mean(cope_zimm_pulo_krig[,1] - 
		qnorm(.95)*cope_zimm_pulo_krig[,2] < Z1p & Z1p < cope_zimm_pulo_krig[,1] + 
		qnorm(.95)*cope_zimm_pulo_krig[,2])
	scen1_500pts_RMSPE[iter,7] = mean((Z1p - cope_zimm_pulo_nene[,1])^2)
	scen1_500pts_PI90[iter,7] = mean(cope_zimm_pulo_nene[,1] - 
		qnorm(.95)*cope_zimm_pulo_nene[,2] < Z1p & Z1p < cope_zimm_pulo_nene[,1] + 
		qnorm(.95)*cope_zimm_pulo_nene[,2])

#-----------------------------------------------------------------------
#                      POBA
#-----------------------------------------------------------------------

  cat("\n", "poba")
  
	copa_all_poba_all = poba(cope_all$theta, Z1, X1, X1p, x, y, xp, yp,
		predmeth = 'alldata')

	copa_rand_poba_all = poba(cope_rand$theta, Z1, X1, X1p, x, y, xp, yp,
		predmeth = 'alldata')
		
	copa_rand_poba_NN = poba(cope_rand$theta, Z1, X1, X1p, x, y, xp, yp,
		predmeth = 'NNdata')

	copa_comp_poba_all = poba(cope_comp$theta, Z1, X1, X1p, x, y, xp, yp,
		predmeth = 'alldata')
		
	copa_comp_poba_NN = poba(cope_comp$theta, Z1, X1, X1p, x, y, xp, yp,
		predmeth = 'NNdata')

	copa_zimm_poba_all = poba(cope_zimm$theta, Z1, X1, X1p, x, y, xp, yp,
		predmeth = 'alldata')
		
	copa_zimm_poba_NN = poba(cope_zimm$theta, Z1, X1, X1p, x, y, xp, yp,
		predmeth = 'NNdata')
	
	scen1_500pts_blk_pred[iter,1] = copa_all_poba_all$block_pred
	scen1_500pts_blk_se[iter,1] = copa_all_poba_all$block_pred_se
	scen1_500pts_blk_pred[iter,2] = copa_rand_poba_all$block_pred
	scen1_500pts_blk_se[iter,2] = copa_rand_poba_all$block_pred_se
	scen1_500pts_blk_pred[iter,3] = copa_rand_poba_NN$block_pred
	scen1_500pts_blk_se[iter,3] = copa_rand_poba_NN$block_pred_se
	scen1_500pts_blk_pred[iter,4] = copa_comp_poba_all$block_pred
	scen1_500pts_blk_se[iter,4] = copa_comp_poba_all$block_pred_se
	scen1_500pts_blk_pred[iter,5] = copa_comp_poba_NN$block_pred
	scen1_500pts_blk_se[iter,5] = copa_comp_poba_NN$block_pred_se
	scen1_500pts_blk_pred[iter,6] = copa_zimm_poba_all$block_pred
	scen1_500pts_blk_se[iter,6] = copa_zimm_poba_all$block_pred_se
	scen1_500pts_blk_pred[iter,7] = copa_zimm_poba_NN$block_pred
	scen1_500pts_blk_se[iter,7] = copa_zimm_poba_NN$block_pred_se
	scen1_500pts_blk_pred[iter,8] = mean(Z1p)
	
	stopiter = Sys.time()
	itertime = difftime(stopiter, startiter, units="mins")
  cat("\n", "iteration time: ", itertime)
}
cat("\n")

basepath = '/home/xverhoef/data/2019_papers/fastREML/testingScripts/temp/'

save(scen1_500pts_est, file = paste0(basepath,'scen1_500pts_est.rda'))
save(scen1_500pts_estse, file = paste0(basepath,'scen1_500pts_estse.rda'))
save(scen1_500pts_RMSPE, file = paste0(basepath,'scen1_500pts_RMSPE.rda'))
save(scen1_500pts_PI90, file = paste0(basepath,'scen1_500pts_PI90.rda'))
save(scen1_500pts_blk_pred, file = paste0(basepath,'scen1_500pts_blk_pred.rda'))
save(scen1_500pts_blk_se, file = paste0(basepath,'scen1_500pts_blk_se.rda'))

