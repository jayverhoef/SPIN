#Rscript sim_scen1_500pts.R
#q(save = 'no')
#R
library('spPart')
library('nabor')
library('spNNGP')
library('pdist')

ScriptPath = paste0('/mnt/ExtraDrive1/Work/desktop_data/2019_papers',
	'/fastREML/spPart_package/spPart/R/')
source(paste0(ScriptPath, 'mass.R'))
source(paste0(ScriptPath, 'pointSimSyst.R'))

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

nsim = 1000

effBlockType_est1 = matrix(NA, nrow = nsim, ncol = 9)
effBlockType_est1se = matrix(NA, nrow = nsim, ncol = 9)
effBlockType_est2 = matrix(NA, nrow = nsim, ncol = 9)
effBlockType_est2se = matrix(NA, nrow = nsim, ncol = 9)
effBlockType_RMSPE = matrix(NA, nrow = nsim, ncol = 9)
effBlockType_PI90 = matrix(NA, nrow = nsim, ncol = 9)

set.seed(1091)
for(iter in 1:nsim) {
#-----------------------------------------------------------------------
#                       SIMULATE
#-----------------------------------------------------------------------
  cat("\n", "Simulation: ", iter)
  cat("\n", "Simulate")
  startiter = Sys.time()

  n = 1000 + round((runif(1)^2)*9000)

	x = runif(n)
	y = runif(n)
	xysyst = pointSimSyst(nrow = 40, ncol = 40, lower.x.lim = 0, upper.x.lim = 1, 
			lower.y.lim = 0, upper.y.lim = 1)
	x = c(x, xysyst[,1])
	y = c(y, xysyst[,2])

if(n <= 2000) {
  range = runif(1)*2 + .01
	xyz1 = geostatSim(loc.data = data.frame(x = x, y = y), xcol = "x", 
		ycol = "y", parsil = 10, range = range, nugget = .0001, 
		CorModel = 'Spherical')
	zstd1 = (xyz1[,3] - mean(xyz1[,3]))/sqrt(var(xyz1[,3]))
	xyz2 = geostatSim(loc.data = data.frame(x = x, y = y), xcol = "x", 
		ycol = "y", parsil = 10, range = range, nugget = .0001, 
		CorModel = 'Spherical')
	zstd2 = (xyz2[,3] - mean(xyz2[,3]))/sqrt(var(xyz2[,3]))
}
if(n > 2000) {
	xyz1 = sim2DSinSurf(x, y)
	zstd1 = (xyz1[,3] - mean(xyz1[,3]))/sqrt(var(xyz1[,3]))
	xyz2 = sim2DSinSurf(x, y)
	zstd2 = (xyz2[,3] - mean(xyz2[,3]))/sqrt(var(xyz2[,3]))  
}


  # design matrix with X1 spatially independent and X1 spatially patterned
	X = as.matrix(cbind(rep(1, times = n + 1600), rnorm(n + 1600), zstd2))
  colnames(X) = c('int','X1','X2')

  # make response variable with random proportion of nugget
  prop = runif(1)
	Z =  1*X[,1]+ 1*X[,2] + 1*X[,3] + 
    10*prop*zstd1 + 10*(1 - prop)*rnorm(n + 1600, 0, 1)

	xp = x[(n+1):(n + 1600)]
	x = x[1:n]
	yp = y[(n+1):(n + 1600)]
	y = y[1:n]
	X1p = X[(n+1):(n + 1600),2]
	X1 = X[1:n,2]
	X2p = X[(n+1):(n + 1600),3]
	X2 = X[1:n,3]
	Zp = Z[(n+1):(n + 1600)]
	Z = Z[1:n]
  Xp = X[(n+1):(n + 1600),]
  X = X[1:n,]

  # make group sizes random between 25 and 225
  grpsize = round(25 + runif(1)*200)
	gir = mass(x, y, ngroups = round(n/grpsize), ssmeth = 'random')
	gic = mass(x, y, ngroups = round(n/grpsize), ssmeth = 'compKmean')
	giz = mass(x, y, ngroups = round(n/grpsize), ssmeth = 'zimmer')

	#d1 = data.frame(x = x[gi1[,1]], y = y[gi1[,1]], X1 = X1[gi1[,1],2],
	#	X2 = X2[gi1[,1]], z1 = Z1[gi1[,1]], z2 = Z2[gi1[,1]], z3 = Z3[gi1[,1]],
	#	z4 = Z4[gi1[,1]])

	d2 = data.frame(x = x, y = y, X1 = X1, X2 = X2, z = Z, 
    grpindx.rand = gir[,2], grpindx.comp = gic[,2],
		grpindx.zimm = giz[,2], grpindx.one = rep(1, times = n))

	dp = data.frame(x = xp, y = yp, X1 = X1p, X2 = X2p, z = Zp)

#-----------------------------------------------------------------------
#                       COPE
#-----------------------------------------------------------------------
  cat("\n", "cope for sample size", n)

	starttime = Sys.time()
	cope_rand = cope(z ~ X1 + X2, data = d2, x_column = 'x', y_column = 'y', 
		subSampCol = 'grpindx.rand', thetaini = c(3,3))
	stoptime = Sys.time()
	difftime(stoptime, starttime, units="secs")

	starttime = Sys.time()
	cope_comp = cope(z ~ X1 + X2, data = d2, x_column = 'x', y_column = 'y', 
		subSampCol = 'grpindx.comp', thetaini = c(3,3))
	stoptime = Sys.time()
	difftime(stoptime, starttime, units="secs")

	starttime = Sys.time()
	cope_zimm = cope(z ~ X1 + X2, data = d2, x_column = 'x', y_column = 'y', 
		subSampCol = 'grpindx.zimm', thetaini = c(3,3))
	stoptime = Sys.time()
	difftime(stoptime, starttime, units="secs")

#-----------------------------------------------------------------------
#                      FEFE
#-----------------------------------------------------------------------

  cat("\n", "fefe")
	starttime = Sys.time()
	cope_rand_fefe_rand = fefe(theta = cope_rand$theta, z = Z, 
		X = X, subsampindx = gir[,2], xcoords = x, ycoords = y)
	stoptime = Sys.time()
	difftime(stoptime, starttime, units="secs")
	effBlockType_est1[iter,1] = coeftable(cope_rand_fefe_rand)$analytic[2,1]
	effBlockType_est1se[iter,1] = coeftable(cope_rand_fefe_rand)$analytic[2,2]
	effBlockType_est2[iter,1] = coeftable(cope_rand_fefe_rand)$analytic[3,1]
	effBlockType_est2se[iter,1] = coeftable(cope_rand_fefe_rand)$analytic[3,2]

	starttime = Sys.time()
	cope_rand_fefe_comp = fefe(theta = cope_rand$theta, z = Z, 
		X = X, subsampindx = gic[,2], xcoords = x, ycoords = y)
	stoptime = Sys.time()
	difftime(stoptime, starttime, units="secs")
	effBlockType_est1[iter,2] = coeftable(cope_rand_fefe_comp)$analytic[2,1]
	effBlockType_est1se[iter,2] = coeftable(cope_rand_fefe_comp)$analytic[2,2]
	effBlockType_est2[iter,2] = coeftable(cope_rand_fefe_comp)$analytic[3,1]
	effBlockType_est2se[iter,2] = coeftable(cope_rand_fefe_comp)$analytic[3,2]
	
	starttime = Sys.time()
	cope_rand_fefe_zimm = fefe(theta = cope_rand$theta, z = Z, 
		X = X, subsampindx = giz[,2], xcoords = x, ycoords = y)
	stoptime = Sys.time()
	difftime(stoptime, starttime, units="secs")
	effBlockType_est1[iter,3] = coeftable(cope_rand_fefe_zimm)$analytic[2,1]
	effBlockType_est1se[iter,3] = coeftable(cope_rand_fefe_zimm)$analytic[2,2]
	effBlockType_est2[iter,3] = coeftable(cope_rand_fefe_zimm)$analytic[3,1]
	effBlockType_est2se[iter,3] = coeftable(cope_rand_fefe_zimm)$analytic[3,2]

	starttime = Sys.time()
	cope_comp_fefe_rand = fefe(theta = cope_comp$theta, z = Z, 
		X = X, subsampindx = gir[,2], xcoords = x, ycoords = y)
	stoptime = Sys.time()
	difftime(stoptime, starttime, units="secs")
	effBlockType_est1[iter,4] = coeftable(cope_comp_fefe_rand)$analytic[2,1]
	effBlockType_est1se[iter,4] = coeftable(cope_comp_fefe_rand)$analytic[2,2]
	effBlockType_est2[iter,4] = coeftable(cope_comp_fefe_rand)$analytic[3,1]
	effBlockType_est2se[iter,4] = coeftable(cope_comp_fefe_rand)$analytic[3,2]

	starttime = Sys.time()
	cope_comp_fefe_comp = fefe(theta = cope_comp$theta, z = Z, 
		X = X, subsampindx = gic[,2], xcoords = x, ycoords = y)
	stoptime = Sys.time()
	difftime(stoptime, starttime, units="secs")
	effBlockType_est1[iter,5] = coeftable(cope_comp_fefe_comp)$analytic[2,1]
	effBlockType_est1se[iter,5] = coeftable(cope_comp_fefe_comp)$analytic[2,2]
	effBlockType_est2[iter,5] = coeftable(cope_comp_fefe_comp)$analytic[3,1]
	effBlockType_est2se[iter,5] = coeftable(cope_comp_fefe_comp)$analytic[3,2]

	starttime = Sys.time()
	cope_comp_fefe_zimm = fefe(theta = cope_comp$theta, z = Z, 
		X = X, subsampindx = giz[,2], xcoords = x, ycoords = y)
	stoptime = Sys.time()
	difftime(stoptime, starttime, units="secs")
	effBlockType_est1[iter,6] = coeftable(cope_comp_fefe_zimm)$analytic[2,1]
	effBlockType_est1se[iter,6] = coeftable(cope_comp_fefe_zimm)$analytic[2,2]
	effBlockType_est2[iter,6] = coeftable(cope_comp_fefe_zimm)$analytic[3,1]
	effBlockType_est2se[iter,6] = coeftable(cope_comp_fefe_zimm)$analytic[3,2]

	starttime = Sys.time()
	cope_zimm_fefe_rand = fefe(theta = cope_zimm$theta, z = Z, 
		X = X, subsampindx = gir[,2], xcoords = x, ycoords = y)
	stoptime = Sys.time()
	difftime(stoptime, starttime, units="secs")
	effBlockType_est1[iter,7] = coeftable(cope_zimm_fefe_rand)$analytic[2,1]
	effBlockType_est1se[iter,7] = coeftable(cope_zimm_fefe_rand)$analytic[2,2]
	effBlockType_est2[iter,7] = coeftable(cope_zimm_fefe_rand)$analytic[3,1]
	effBlockType_est2se[iter,7] = coeftable(cope_zimm_fefe_rand)$analytic[3,2]

	starttime = Sys.time()
	cope_zimm_fefe_comp = fefe(theta = cope_zimm$theta, z = Z, 
		X = X, subsampindx = gic[,2], xcoords = x, ycoords = y)
	stoptime = Sys.time()
	difftime(stoptime, starttime, units="secs")
	effBlockType_est1[iter,8] = coeftable(cope_zimm_fefe_comp)$analytic[2,1]
	effBlockType_est1se[iter,8] = coeftable(cope_zimm_fefe_comp)$analytic[2,2]
	effBlockType_est2[iter,8] = coeftable(cope_zimm_fefe_comp)$analytic[3,1]
	effBlockType_est2se[iter,8] = coeftable(cope_zimm_fefe_comp)$analytic[3,2]

	starttime = Sys.time()
	cope_zimm_fefe_zimm = fefe(theta = cope_zimm$theta, z = Z, 
		X = X, subsampindx = giz[,2], xcoords = x, ycoords = y)
	stoptime = Sys.time()
	difftime(stoptime, starttime, units="secs")
	effBlockType_est1[iter,9] = coeftable(cope_zimm_fefe_zimm)$analytic[2,1]
	effBlockType_est1se[iter,9] = coeftable(cope_zimm_fefe_zimm)$analytic[2,2]
	effBlockType_est2[iter,9] = coeftable(cope_zimm_fefe_zimm)$analytic[3,1]
	effBlockType_est2se[iter,9] = coeftable(cope_zimm_fefe_zimm)$analytic[3,2]

#-----------------------------------------------------------------------
#                      PULO
#-----------------------------------------------------------------------
  cat("\n", "pulo")
	starttime = Sys.time()
	cope_rand_fefe_rand_puloNN = pulo1(theta = cope_rand$theta, 
    betahat = cope_rand_fefe_rand$bhat, covb = cope_rand_fefe_rand$covb,
    z = Z, X = X, Xp = Xp, xcoords = x, ycoords = y, 
    xcoordsp = xp, ycoordsp = yp, predmeth = 'nearnei', par = FALSE)
	stoptime = Sys.time()
	difftime(stoptime, starttime, units="secs")
	effBlockType_RMSPE[iter,1] = mean((Zp - cope_rand_fefe_rand_puloNN[,1])^2)
	effBlockType_PI90[iter,1] = mean(cope_rand_fefe_rand_puloNN[,1] - 
		qnorm(.95)*cope_rand_fefe_rand_puloNN[,2] < Zp & Zp < cope_rand_fefe_rand_puloNN[,1] + 
		qnorm(.95)*cope_rand_fefe_rand_puloNN[,2])

	starttime = Sys.time()
	cope_rand_fefe_comp_puloNN = pulo1(theta = cope_rand$theta, 
    betahat = cope_rand_fefe_comp$bhat, covb = cope_rand_fefe_comp$covb,
    z = Z, X = X, Xp = Xp, xcoords = x, ycoords = y, 
    xcoordsp = xp, ycoordsp = yp, predmeth = 'nearnei', par = FALSE)
	stoptime = Sys.time()
	difftime(stoptime, starttime, units="secs")
	effBlockType_RMSPE[iter,2] = mean((Zp - cope_rand_fefe_comp_puloNN[,1])^2)
	effBlockType_PI90[iter,2] = mean(cope_rand_fefe_comp_puloNN[,1] - 
		qnorm(.95)*cope_rand_fefe_comp_puloNN[,2] < Zp & Zp < cope_rand_fefe_comp_puloNN[,1] + 
		qnorm(.95)*cope_rand_fefe_comp_puloNN[,2])
	 
	starttime = Sys.time()
	cope_rand_fefe_zimm_puloNN = pulo1(theta = cope_rand$theta, 
    betahat = cope_rand_fefe_zimm$bhat, covb = cope_rand_fefe_zimm$covb,
    z = Z, X = X, Xp = Xp, xcoords = x, ycoords = y, 
    xcoordsp = xp, ycoordsp = yp, predmeth = 'nearnei', par = FALSE)
	stoptime = Sys.time()
	difftime(stoptime, starttime, units="secs")
	effBlockType_RMSPE[iter,3] = mean((Zp - cope_rand_fefe_zimm_puloNN[,1])^2)
	effBlockType_PI90[iter,3] = mean(cope_rand_fefe_zimm_puloNN[,1] - 
		qnorm(.95)*cope_rand_fefe_zimm_puloNN[,2] < Zp & Zp < cope_rand_fefe_zimm_puloNN[,1] + 
		qnorm(.95)*cope_rand_fefe_zimm_puloNN[,2])

	starttime = Sys.time()
	cope_comp_fefe_rand_puloNN = pulo1(theta = cope_comp$theta, 
    betahat = cope_comp_fefe_rand$bhat, covb = cope_comp_fefe_rand$covb,
    z = Z, X = X, Xp = Xp, xcoords = x, ycoords = y, 
    xcoordsp = xp, ycoordsp = yp, predmeth = 'nearnei', par = FALSE)
	stoptime = Sys.time()
	difftime(stoptime, starttime, units="secs")
	effBlockType_RMSPE[iter,4] = mean((Zp - cope_comp_fefe_rand_puloNN[,1])^2)
	effBlockType_PI90[iter,4] = mean(cope_comp_fefe_rand_puloNN[,1] - 
		qnorm(.95)*cope_comp_fefe_rand_puloNN[,2] < Zp & Zp < cope_comp_fefe_rand_puloNN[,1] + 
		qnorm(.95)*cope_comp_fefe_rand_puloNN[,2])

	starttime = Sys.time()
	cope_comp_fefe_comp_puloNN = pulo1(theta = cope_comp$theta, 
    betahat = cope_comp_fefe_comp$bhat, covb = cope_comp_fefe_comp$covb,
    z = Z, X = X, Xp = Xp, xcoords = x, ycoords = y, 
    xcoordsp = xp, ycoordsp = yp, predmeth = 'nearnei', par = FALSE)
	stoptime = Sys.time()
	difftime(stoptime, starttime, units="secs")
	effBlockType_RMSPE[iter,5] = mean((Zp - cope_comp_fefe_comp_puloNN[,1])^2)
	effBlockType_PI90[iter,5] = mean(cope_comp_fefe_comp_puloNN[,1] - 
		qnorm(.95)*cope_comp_fefe_comp_puloNN[,2] < Zp & Zp < cope_comp_fefe_comp_puloNN[,1] + 
		qnorm(.95)*cope_comp_fefe_comp_puloNN[,2])

	starttime = Sys.time()
	cope_comp_fefe_zimm_puloNN = pulo1(theta = cope_comp$theta, 
    betahat = cope_comp_fefe_zimm$bhat, covb = cope_comp_fefe_zimm$covb,
    z = Z, X = X, Xp = Xp, xcoords = x, ycoords = y, 
    xcoordsp = xp, ycoordsp = yp, predmeth = 'nearnei', par = FALSE)
	stoptime = Sys.time()
	difftime(stoptime, starttime, units="secs")
	effBlockType_RMSPE[iter,6] = mean((Zp - cope_comp_fefe_zimm_puloNN[,1])^2)
	effBlockType_PI90[iter,6] = mean(cope_comp_fefe_zimm_puloNN[,1] - 
		qnorm(.95)*cope_comp_fefe_zimm_puloNN[,2] < Zp & Zp < cope_comp_fefe_zimm_puloNN[,1] + 
		qnorm(.95)*cope_comp_fefe_zimm_puloNN[,2])
    
	starttime = Sys.time()
	cope_zimm_fefe_rand_puloNN = pulo1(theta = cope_zimm$theta, 
    betahat = cope_zimm_fefe_rand$bhat, covb = cope_zimm_fefe_rand$covb,
    z = Z, X = X, Xp = Xp, xcoords = x, ycoords = y, 
    xcoordsp = xp, ycoordsp = yp, predmeth = 'nearnei', par = FALSE)
	stoptime = Sys.time()
	difftime(stoptime, starttime, units="secs")
	effBlockType_RMSPE[iter,7] = mean((Zp - cope_zimm_fefe_rand_puloNN[,1])^2)
	effBlockType_PI90[iter,7] = mean(cope_zimm_fefe_rand_puloNN[,1] - 
		qnorm(.95)*cope_zimm_fefe_rand_puloNN[,2] < Zp & Zp < cope_zimm_fefe_rand_puloNN[,1] + 
		qnorm(.95)*cope_zimm_fefe_rand_puloNN[,2])

	starttime = Sys.time()
	cope_zimm_fefe_comp_puloNN = pulo1(theta = cope_zimm$theta, 
    betahat = cope_zimm_fefe_comp$bhat, covb = cope_zimm_fefe_comp$covb,
    z = Z, X = X, Xp = Xp, xcoords = x, ycoords = y, 
    xcoordsp = xp, ycoordsp = yp, predmeth = 'nearnei', par = FALSE)
	stoptime = Sys.time()
	difftime(stoptime, starttime, units="secs")
	effBlockType_RMSPE[iter,8] = mean((Zp - cope_zimm_fefe_comp_puloNN[,1])^2)
	effBlockType_PI90[iter,8] = mean(cope_zimm_fefe_comp_puloNN[,1] - 
		qnorm(.95)*cope_zimm_fefe_comp_puloNN[,2] < Zp & Zp < cope_zimm_fefe_comp_puloNN[,1] + 
		qnorm(.95)*cope_zimm_fefe_comp_puloNN[,2])

	starttime = Sys.time()
	cope_zimm_fefe_zimm_puloNN = pulo1(theta = cope_zimm$theta, 
    betahat = cope_zimm_fefe_zimm$bhat, covb = cope_zimm_fefe_zimm$covb,
    z = Z, X = X, Xp = Xp, xcoords = x, ycoords = y, 
    xcoordsp = xp, ycoordsp = yp, predmeth = 'nearnei', par = FALSE)
	stoptime = Sys.time()
	difftime(stoptime, starttime, units="secs")
	effBlockType_RMSPE[iter,9] = mean((Zp - cope_zimm_fefe_zimm_puloNN[,1])^2)
	effBlockType_PI90[iter,9] = mean(cope_zimm_fefe_zimm_puloNN[,1] - 
		qnorm(.95)*cope_zimm_fefe_zimm_puloNN[,2] < Zp & Zp < cope_zimm_fefe_zimm_puloNN[,1] + 
		qnorm(.95)*cope_zimm_fefe_zimm_puloNN[,2])
    
		
	stopiter = Sys.time()
	itertime = difftime(stopiter, startiter, units="mins")
  cat("\n", "iteration time: ", itertime)
}

basepath = '/media/jay/data/2019_papers/fastREML/spPart_package/spPart/data/'
save(effBlockType_est1, file = paste0(basepath,'effBlockType_est1.rda'))
save(effBlockType_est1se, file = paste0(basepath,'effBlockType_est1se.rda'))
save(effBlockType_est2, file = paste0(basepath,'effBlockType_est2.rda'))
save(effBlockType_est2se, file = paste0(basepath,'effBlockType_est2se.rda'))
save(effBlockType_RMSPE, file = paste0(basepath,'effBlockType_RMSPE.rda'))
save(effBlockType_PI90, file = paste0(basepath,'effBlockType_PI90.rda'))
