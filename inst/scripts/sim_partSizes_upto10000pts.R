#Rscript sim_scen1_500pts.R
#q(save = 'no')
#R
library('spPart')
library('nabor')
library('spNNGP')
library('pdist')


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

nsim = 400

partSizes_est1 = matrix(NA, nrow = nsim, ncol = 16)
partSizes_est1se = matrix(NA, nrow = nsim, ncol = 16)
partSizes_est2 = matrix(NA, nrow = nsim, ncol = 16)
partSizes_est2se = matrix(NA, nrow = nsim, ncol = 16)
partSizes_RMSPE = matrix(NA, nrow = nsim, ncol = 16)
partSizes_PI90 = matrix(NA, nrow = nsim, ncol = 16)
partSizes_copetime = matrix(NA, nrow = nsim, ncol = 4)
partSizes_fefetime = matrix(NA, nrow = nsim, ncol = 4)


set.seed(1092)
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

  # make group sizes of 25, 50, 100, and 200 using compact
	gi25 = mass(x, y, ngroups = round(n/25), ssmeth = 'random')
	gi50 = mass(x, y, ngroups = round(n/50), ssmeth = 'compKmean')
	gi100 = mass(x, y, ngroups = round(n/100), ssmeth = 'zimmer')
	gi200 = mass(x, y, ngroups = round(n/200), ssmeth = 'zimmer')

	#d1 = data.frame(x = x[gi1[,1]], y = y[gi1[,1]], X1 = X1[gi1[,1],2],
	#	X2 = X2[gi1[,1]], z1 = Z1[gi1[,1]], z2 = Z2[gi1[,1]], z3 = Z3[gi1[,1]],
	#	z4 = Z4[gi1[,1]])

	d2 = data.frame(x = x, y = y, X1 = X1, X2 = X2, z = Z, 
    grpindx.25 = gi25[,2], grpindx.50 = gi50[,2],
		grpindx.100 = gi100[,2], grpindx.200 = gi200[,2])

	dp = data.frame(x = xp, y = yp, X1 = X1p, X2 = X2p, z = Zp)

#-----------------------------------------------------------------------
#                       COPE
#-----------------------------------------------------------------------
  cat("\n", "cope for sample size", n)

	starttime = Sys.time()
	cope_25 = cope(z ~ X1 + X2, data = d2, x_column = 'x', y_column = 'y', 
		subSampCol = 'grpindx.25', thetaini = c(3,3))
	stoptime = Sys.time()
	partSizes_copetime[iter,1] = difftime(stoptime,
    starttime, units="secs")

	starttime = Sys.time()
	cope_50 = cope(z ~ X1 + X2, data = d2, x_column = 'x', y_column = 'y', 
		subSampCol = 'grpindx.50', thetaini = c(3,3))
	stoptime = Sys.time()
	partSizes_copetime[iter,2] = difftime(stoptime,
    starttime, units="secs")

	starttime = Sys.time()
	cope_100 = cope(z ~ X1 + X2, data = d2, x_column = 'x', y_column = 'y', 
		subSampCol = 'grpindx.100', thetaini = c(3,3))
	stoptime = Sys.time()
	partSizes_copetime[iter,3] = difftime(stoptime,
    starttime, units="secs")

	starttime = Sys.time()
	cope_200 = cope(z ~ X1 + X2, data = d2, x_column = 'x', y_column = 'y', 
		subSampCol = 'grpindx.200', thetaini = c(3,3))
	stoptime = Sys.time()
	partSizes_copetime[iter,4] = difftime(stoptime,
    starttime, units="secs")

#-----------------------------------------------------------------------
#                      FEFE
#-----------------------------------------------------------------------

  cat("\n", "fefe")
  
	starttime = Sys.time()
	cope25_fefe25 = fefe(theta = cope_25$theta, z = Z, 
		X = X, subsampindx = gi25[,2], xcoords = x, ycoords = y)
	stoptime = Sys.time()
	partSizes_est1[iter,1] = coeftable(cope25_fefe25)$analytic[2,1]
	partSizes_est1se[iter,1] = coeftable(cope25_fefe25)$analytic[2,2]
	partSizes_est2[iter,1] = coeftable(cope25_fefe25)$analytic[3,1]
	partSizes_est2se[iter,1] = coeftable(cope25_fefe25)$analytic[3,2]
  partSizes_fefetime[iter,1] = difftime(stoptime,
    starttime, units="secs")

	starttime = Sys.time()
	cope25_fefe50 = fefe(theta = cope_25$theta, z = Z, 
		X = X, subsampindx = gi50[,2], xcoords = x, ycoords = y)
	stoptime = Sys.time()
	partSizes_est1[iter,2] = coeftable(cope25_fefe50)$analytic[2,1]
	partSizes_est1se[iter,2] = coeftable(cope25_fefe50)$analytic[2,2]
	partSizes_est2[iter,2] = coeftable(cope25_fefe50)$analytic[3,1]
	partSizes_est2se[iter,2] = coeftable(cope25_fefe50)$analytic[3,2]
  partSizes_fefetime[iter,2] = difftime(stoptime,
    starttime, units="secs")

	starttime = Sys.time()
	cope25_fefe100 = fefe(theta = cope_25$theta, z = Z, 
		X = X, subsampindx = gi100[,2], xcoords = x, ycoords = y)
	stoptime = Sys.time()
	partSizes_est1[iter,3] = coeftable(cope25_fefe100)$analytic[2,1]
	partSizes_est1se[iter,3] = coeftable(cope25_fefe100)$analytic[2,2]
	partSizes_est2[iter,3] = coeftable(cope25_fefe100)$analytic[3,1]
	partSizes_est2se[iter,3] = coeftable(cope25_fefe100)$analytic[3,2]
  partSizes_fefetime[iter,3] = difftime(stoptime,
    starttime, units="secs")

	starttime = Sys.time()
	cope25_fefe200 = fefe(theta = cope_25$theta, z = Z, 
		X = X, subsampindx = gi200[,2], xcoords = x, ycoords = y)
	stoptime = Sys.time()
	partSizes_est1[iter,4] = coeftable(cope25_fefe200)$analytic[2,1]
	partSizes_est1se[iter,4] = coeftable(cope25_fefe200)$analytic[2,2]
	partSizes_est2[iter,4] = coeftable(cope25_fefe200)$analytic[3,1]
	partSizes_est2se[iter,4] = coeftable(cope25_fefe200)$analytic[3,2]
  partSizes_fefetime[iter,4] = difftime(stoptime,
    starttime, units="secs")
	
	starttime = Sys.time()
	cope50_fefe25 = fefe(theta = cope_50$theta, z = Z, 
		X = X, subsampindx = gi25[,2], xcoords = x, ycoords = y)
	stoptime = Sys.time()
	partSizes_est1[iter,5] = coeftable(cope50_fefe25)$analytic[2,1]
	partSizes_est1se[iter,5] = coeftable(cope50_fefe25)$analytic[2,2]
	partSizes_est2[iter,5] = coeftable(cope50_fefe25)$analytic[3,1]
	partSizes_est2se[iter,5] = coeftable(cope50_fefe25)$analytic[3,2]

	starttime = Sys.time()
	cope50_fefe50 = fefe(theta = cope_50$theta, z = Z, 
		X = X, subsampindx = gi50[,2], xcoords = x, ycoords = y)
	stoptime = Sys.time()
	partSizes_est1[iter,6] = coeftable(cope50_fefe50)$analytic[2,1]
	partSizes_est1se[iter,6] = coeftable(cope50_fefe50)$analytic[2,2]
	partSizes_est2[iter,6] = coeftable(cope50_fefe50)$analytic[3,1]
	partSizes_est2se[iter,6] = coeftable(cope50_fefe50)$analytic[3,2]

	starttime = Sys.time()
	cope50_fefe100 = fefe(theta = cope_50$theta, z = Z, 
		X = X, subsampindx = gi100[,2], xcoords = x, ycoords = y)
	stoptime = Sys.time()
	partSizes_est1[iter,7] = coeftable(cope50_fefe100)$analytic[2,1]
	partSizes_est1se[iter,7] = coeftable(cope50_fefe100)$analytic[2,2]
	partSizes_est2[iter,7] = coeftable(cope50_fefe100)$analytic[3,1]
	partSizes_est2se[iter,7] = coeftable(cope50_fefe100)$analytic[3,2]

	starttime = Sys.time()
	cope50_fefe200 = fefe(theta = cope_50$theta, z = Z, 
		X = X, subsampindx = gi200[,2], xcoords = x, ycoords = y)
	stoptime = Sys.time()
	partSizes_est1[iter,8] = coeftable(cope50_fefe200)$analytic[2,1]
	partSizes_est1se[iter,8] = coeftable(cope50_fefe200)$analytic[2,2]
	partSizes_est2[iter,8] = coeftable(cope50_fefe200)$analytic[3,1]
	partSizes_est2se[iter,8] = coeftable(cope50_fefe200)$analytic[3,2]

	starttime = Sys.time()
	cope100_fefe25 = fefe(theta = cope_100$theta, z = Z, 
		X = X, subsampindx = gi25[,2], xcoords = x, ycoords = y)
	stoptime = Sys.time()
	partSizes_est1[iter,9] = coeftable(cope100_fefe25)$analytic[2,1]
	partSizes_est1se[iter,9] = coeftable(cope100_fefe25)$analytic[2,2]
	partSizes_est2[iter,9] = coeftable(cope100_fefe25)$analytic[3,1]
	partSizes_est2se[iter,9] = coeftable(cope100_fefe25)$analytic[3,2]

	starttime = Sys.time()
	cope100_fefe50 = fefe(theta = cope_100$theta, z = Z, 
		X = X, subsampindx = gi50[,2], xcoords = x, ycoords = y)
	stoptime = Sys.time()
	partSizes_est1[iter,10] = coeftable(cope100_fefe50)$analytic[2,1]
	partSizes_est1se[iter,10] = coeftable(cope100_fefe50)$analytic[2,2]
	partSizes_est2[iter,10] = coeftable(cope100_fefe50)$analytic[3,1]
	partSizes_est2se[iter,10] = coeftable(cope100_fefe50)$analytic[3,2]

	starttime = Sys.time()
	cope100_fefe100 = fefe(theta = cope_100$theta, z = Z, 
		X = X, subsampindx = gi100[,2], xcoords = x, ycoords = y)
	stoptime = Sys.time()
	partSizes_est1[iter,11] = coeftable(cope100_fefe100)$analytic[2,1]
	partSizes_est1se[iter,11] = coeftable(cope100_fefe100)$analytic[2,2]
	partSizes_est2[iter,11] = coeftable(cope100_fefe100)$analytic[3,1]
	partSizes_est2se[iter,11] = coeftable(cope100_fefe100)$analytic[3,2]

	starttime = Sys.time()
	cope100_fefe200 = fefe(theta = cope_100$theta, z = Z, 
		X = X, subsampindx = gi200[,2], xcoords = x, ycoords = y)
	stoptime = Sys.time()
	partSizes_est1[iter,12] = coeftable(cope100_fefe200)$analytic[2,1]
	partSizes_est1se[iter,12] = coeftable(cope100_fefe200)$analytic[2,2]
	partSizes_est2[iter,12] = coeftable(cope100_fefe200)$analytic[3,1]
	partSizes_est2se[iter,12] = coeftable(cope100_fefe200)$analytic[3,2]

	starttime = Sys.time()
	cope200_fefe25 = fefe(theta = cope_200$theta, z = Z, 
		X = X, subsampindx = gi25[,2], xcoords = x, ycoords = y)
	stoptime = Sys.time()
	partSizes_est1[iter,13] = coeftable(cope200_fefe25)$analytic[2,1]
	partSizes_est1se[iter,13] = coeftable(cope200_fefe25)$analytic[2,2]
	partSizes_est2[iter,13] = coeftable(cope200_fefe25)$analytic[3,1]
	partSizes_est2se[iter,13] = coeftable(cope200_fefe25)$analytic[3,2]

	starttime = Sys.time()
	cope200_fefe50 = fefe(theta = cope_200$theta, z = Z, 
		X = X, subsampindx = gi50[,2], xcoords = x, ycoords = y)
	stoptime = Sys.time()
	partSizes_est1[iter,14] = coeftable(cope200_fefe50)$analytic[2,1]
	partSizes_est1se[iter,14] = coeftable(cope200_fefe50)$analytic[2,2]
	partSizes_est2[iter,14] = coeftable(cope200_fefe50)$analytic[3,1]
	partSizes_est2se[iter,14] = coeftable(cope200_fefe50)$analytic[3,2]

	starttime = Sys.time()
	cope200_fefe100 = fefe(theta = cope_200$theta, z = Z, 
		X = X, subsampindx = gi100[,2], xcoords = x, ycoords = y)
	stoptime = Sys.time()
	partSizes_est1[iter,15] = coeftable(cope200_fefe100)$analytic[2,1]
	partSizes_est1se[iter,15] = coeftable(cope200_fefe100)$analytic[2,2]
	partSizes_est2[iter,15] = coeftable(cope200_fefe100)$analytic[3,1]
	partSizes_est2se[iter,15] = coeftable(cope200_fefe100)$analytic[3,2]

	starttime = Sys.time()
	cope200_fefe200 = fefe(theta = cope_200$theta, z = Z, 
		X = X, subsampindx = gi200[,2], xcoords = x, ycoords = y)
	stoptime = Sys.time()
	partSizes_est1[iter,16] = coeftable(cope200_fefe200)$analytic[2,1]
	partSizes_est1se[iter,16] = coeftable(cope200_fefe200)$analytic[2,2]
	partSizes_est2[iter,16] = coeftable(cope200_fefe200)$analytic[3,1]
	partSizes_est2se[iter,16] = coeftable(cope200_fefe200)$analytic[3,2]

#-----------------------------------------------------------------------
#                      PULO
#-----------------------------------------------------------------------

  cat("\n", "pulo")
  
	starttime = Sys.time()
	cope25_fefe25_puloNN = pulo1(theta = cope_25$theta, 
    betahat = cope25_fefe25$bhat, covb = cope25_fefe25$covb,
    z = Z, X = X, Xp = Xp, xcoords = x, ycoords = y, 
    xcoordsp = xp, ycoordsp = yp, predmeth = 'nearnei', par = FALSE)
	stoptime = Sys.time()
	difftime(stoptime, starttime, units="secs")
	partSizes_RMSPE[iter,1] = mean((Zp - cope25_fefe25_puloNN[,1])^2)
	partSizes_PI90[iter,1] = mean(cope25_fefe25_puloNN[,1] -
		qnorm(.95)*cope25_fefe25_puloNN[,2] < Zp & Zp < cope25_fefe25_puloNN[,1] + 
		qnorm(.95)*cope25_fefe25_puloNN[,2])

	starttime = Sys.time()
	cope25_fefe50_puloNN = pulo1(theta = cope_25$theta, 
    betahat = cope25_fefe50$bhat, covb = cope25_fefe50$covb,
    z = Z, X = X, Xp = Xp, xcoords = x, ycoords = y, 
    xcoordsp = xp, ycoordsp = yp, predmeth = 'nearnei', par = FALSE)
	stoptime = Sys.time()
	difftime(stoptime, starttime, units="secs")
	partSizes_RMSPE[iter,2] = mean((Zp - cope25_fefe50_puloNN[,1])^2)
	partSizes_PI90[iter,2] = mean(cope25_fefe50_puloNN[,1] -
		qnorm(.95)*cope25_fefe50_puloNN[,2] < Zp & Zp < cope25_fefe50_puloNN[,1] + 
		qnorm(.95)*cope25_fefe50_puloNN[,2])

	starttime = Sys.time()
	cope25_fefe100_puloNN = pulo1(theta = cope_25$theta, 
    betahat = cope25_fefe100$bhat, covb = cope25_fefe100$covb,
    z = Z, X = X, Xp = Xp, xcoords = x, ycoords = y, 
    xcoordsp = xp, ycoordsp = yp, predmeth = 'nearnei', par = FALSE)
	stoptime = Sys.time()
	difftime(stoptime, starttime, units="secs")
	partSizes_RMSPE[iter,3] = mean((Zp - cope25_fefe100_puloNN[,1])^2)
	partSizes_PI90[iter,3] = mean(cope25_fefe100_puloNN[,1] -
		qnorm(.95)*cope25_fefe100_puloNN[,2] < Zp & Zp < cope25_fefe100_puloNN[,1] + 
		qnorm(.95)*cope25_fefe100_puloNN[,2])

	starttime = Sys.time()
	cope25_fefe200_puloNN = pulo1(theta = cope_25$theta, 
    betahat = cope25_fefe200$bhat, covb = cope25_fefe200$covb,
    z = Z, X = X, Xp = Xp, xcoords = x, ycoords = y, 
    xcoordsp = xp, ycoordsp = yp, predmeth = 'nearnei', par = FALSE)
	stoptime = Sys.time()
	difftime(stoptime, starttime, units="secs")
	partSizes_RMSPE[iter,4] = mean((Zp - cope25_fefe200_puloNN[,1])^2)
	partSizes_PI90[iter,4] = mean(cope25_fefe200_puloNN[,1] -
		qnorm(.95)*cope25_fefe200_puloNN[,2] < Zp & Zp < cope25_fefe200_puloNN[,1] + 
		qnorm(.95)*cope25_fefe200_puloNN[,2])

	starttime = Sys.time()
	cope50_fefe25_puloNN = pulo1(theta = cope_50$theta, 
    betahat = cope50_fefe25$bhat, covb = cope50_fefe25$covb,
    z = Z, X = X, Xp = Xp, xcoords = x, ycoords = y, 
    xcoordsp = xp, ycoordsp = yp, predmeth = 'nearnei', par = FALSE)
	stoptime = Sys.time()
	difftime(stoptime, starttime, units="secs")
	partSizes_RMSPE[iter,5] = mean((Zp - cope50_fefe25_puloNN[,1])^2)
	partSizes_PI90[iter,5] = mean(cope50_fefe25_puloNN[,1] -
		qnorm(.95)*cope50_fefe25_puloNN[,2] < Zp & Zp < cope50_fefe25_puloNN[,1] + 
		qnorm(.95)*cope50_fefe25_puloNN[,2])

	starttime = Sys.time()
	cope50_fefe50_puloNN = pulo1(theta = cope_50$theta, 
    betahat = cope50_fefe50$bhat, covb = cope50_fefe50$covb,
    z = Z, X = X, Xp = Xp, xcoords = x, ycoords = y, 
    xcoordsp = xp, ycoordsp = yp, predmeth = 'nearnei', par = FALSE)
	stoptime = Sys.time()
	difftime(stoptime, starttime, units="secs")
	partSizes_RMSPE[iter,6] = mean((Zp - cope50_fefe50_puloNN[,1])^2)
	partSizes_PI90[iter,6] = mean(cope50_fefe50_puloNN[,1] -
		qnorm(.95)*cope50_fefe50_puloNN[,2] < Zp & Zp < cope50_fefe50_puloNN[,1] + 
		qnorm(.95)*cope50_fefe50_puloNN[,2])

	starttime = Sys.time()
	cope50_fefe100_puloNN = pulo1(theta = cope_50$theta, 
    betahat = cope50_fefe100$bhat, covb = cope50_fefe100$covb,
    z = Z, X = X, Xp = Xp, xcoords = x, ycoords = y, 
    xcoordsp = xp, ycoordsp = yp, predmeth = 'nearnei', par = FALSE)
	stoptime = Sys.time()
	difftime(stoptime, starttime, units="secs")
	partSizes_RMSPE[iter,7] = mean((Zp - cope50_fefe100_puloNN[,1])^2)
	partSizes_PI90[iter,7] = mean(cope50_fefe100_puloNN[,1] -
		qnorm(.95)*cope50_fefe100_puloNN[,2] < Zp & Zp < cope50_fefe100_puloNN[,1] + 
		qnorm(.95)*cope50_fefe100_puloNN[,2])

	starttime = Sys.time()
	cope50_fefe200_puloNN = pulo1(theta = cope_50$theta, 
    betahat = cope50_fefe200$bhat, covb = cope50_fefe200$covb,
    z = Z, X = X, Xp = Xp, xcoords = x, ycoords = y, 
    xcoordsp = xp, ycoordsp = yp, predmeth = 'nearnei', par = FALSE)
	stoptime = Sys.time()
	difftime(stoptime, starttime, units="secs")
	partSizes_RMSPE[iter,8] = mean((Zp - cope50_fefe200_puloNN[,1])^2)
	partSizes_PI90[iter,8] = mean(cope50_fefe200_puloNN[,1] -
		qnorm(.95)*cope50_fefe200_puloNN[,2] < Zp & Zp < cope50_fefe200_puloNN[,1] + 
		qnorm(.95)*cope50_fefe200_puloNN[,2])

	starttime = Sys.time()
	cope100_fefe25_puloNN = pulo1(theta = cope_100$theta, 
    betahat = cope100_fefe25$bhat, covb = cope100_fefe25$covb,
    z = Z, X = X, Xp = Xp, xcoords = x, ycoords = y, 
    xcoordsp = xp, ycoordsp = yp, predmeth = 'nearnei', par = FALSE)
	stoptime = Sys.time()
	difftime(stoptime, starttime, units="secs")
	partSizes_RMSPE[iter,9] = mean((Zp - cope100_fefe25_puloNN[,1])^2)
	partSizes_PI90[iter,9] = mean(cope100_fefe25_puloNN[,1] -
		qnorm(.95)*cope100_fefe25_puloNN[,2] < Zp & Zp < cope100_fefe25_puloNN[,1] + 
		qnorm(.95)*cope100_fefe25_puloNN[,2])

	starttime = Sys.time()
	cope100_fefe50_puloNN = pulo1(theta = cope_100$theta, 
    betahat = cope100_fefe50$bhat, covb = cope100_fefe50$covb,
    z = Z, X = X, Xp = Xp, xcoords = x, ycoords = y, 
    xcoordsp = xp, ycoordsp = yp, predmeth = 'nearnei', par = FALSE)
	stoptime = Sys.time()
	difftime(stoptime, starttime, units="secs")
	partSizes_RMSPE[iter,10] = mean((Zp - cope100_fefe50_puloNN[,1])^2)
	partSizes_PI90[iter,10] = mean(cope100_fefe50_puloNN[,1] -
		qnorm(.95)*cope100_fefe50_puloNN[,2] < Zp & Zp < cope100_fefe50_puloNN[,1] + 
		qnorm(.95)*cope100_fefe50_puloNN[,2])

	starttime = Sys.time()
	cope100_fefe100_puloNN = pulo1(theta = cope_100$theta, 
    betahat = cope100_fefe100$bhat, covb = cope100_fefe100$covb,
    z = Z, X = X, Xp = Xp, xcoords = x, ycoords = y, 
    xcoordsp = xp, ycoordsp = yp, predmeth = 'nearnei', par = FALSE)
	stoptime = Sys.time()
	difftime(stoptime, starttime, units="secs")
	partSizes_RMSPE[iter,11] = mean((Zp - cope100_fefe100_puloNN[,1])^2)
	partSizes_PI90[iter,11] = mean(cope100_fefe100_puloNN[,1] -
		qnorm(.95)*cope100_fefe100_puloNN[,2] < Zp & Zp < cope100_fefe100_puloNN[,1] + 
		qnorm(.95)*cope100_fefe100_puloNN[,2])

	starttime = Sys.time()
	cope100_fefe200_puloNN = pulo1(theta = cope_100$theta, 
    betahat = cope100_fefe200$bhat, covb = cope100_fefe200$covb,
    z = Z, X = X, Xp = Xp, xcoords = x, ycoords = y, 
    xcoordsp = xp, ycoordsp = yp, predmeth = 'nearnei', par = FALSE)
	stoptime = Sys.time()
	difftime(stoptime, starttime, units="secs")
	partSizes_RMSPE[iter,12] = mean((Zp - cope100_fefe200_puloNN[,1])^2)
	partSizes_PI90[iter,12] = mean(cope100_fefe200_puloNN[,1] -
		qnorm(.95)*cope100_fefe200_puloNN[,2] < Zp & Zp < cope100_fefe200_puloNN[,1] + 
		qnorm(.95)*cope100_fefe200_puloNN[,2])
    
	starttime = Sys.time()
	cope200_fefe25_puloNN = pulo1(theta = cope_200$theta, 
    betahat = cope200_fefe25$bhat, covb = cope200_fefe25$covb,
    z = Z, X = X, Xp = Xp, xcoords = x, ycoords = y, 
    xcoordsp = xp, ycoordsp = yp, predmeth = 'nearnei', par = FALSE)
	stoptime = Sys.time()
	difftime(stoptime, starttime, units="secs")
	partSizes_RMSPE[iter,13] = mean((Zp - cope200_fefe25_puloNN[,1])^2)
	partSizes_PI90[iter,13] = mean(cope200_fefe25_puloNN[,1] -
		qnorm(.95)*cope200_fefe25_puloNN[,2] < Zp & Zp < cope200_fefe25_puloNN[,1] + 
		qnorm(.95)*cope200_fefe25_puloNN[,2])

	starttime = Sys.time()
	cope200_fefe50_puloNN = pulo1(theta = cope_200$theta, 
    betahat = cope200_fefe50$bhat, covb = cope200_fefe50$covb,
    z = Z, X = X, Xp = Xp, xcoords = x, ycoords = y, 
    xcoordsp = xp, ycoordsp = yp, predmeth = 'nearnei', par = FALSE)
	stoptime = Sys.time()
	difftime(stoptime, starttime, units="secs")
	partSizes_RMSPE[iter,14] = mean((Zp - cope200_fefe50_puloNN[,1])^2)
	partSizes_PI90[iter,14] = mean(cope200_fefe50_puloNN[,1] -
		qnorm(.95)*cope200_fefe50_puloNN[,2] < Zp & Zp < cope200_fefe50_puloNN[,1] + 
		qnorm(.95)*cope200_fefe50_puloNN[,2])

	starttime = Sys.time()
	cope200_fefe100_puloNN = pulo1(theta = cope_200$theta, 
    betahat = cope200_fefe100$bhat, covb = cope200_fefe100$covb,
    z = Z, X = X, Xp = Xp, xcoords = x, ycoords = y, 
    xcoordsp = xp, ycoordsp = yp, predmeth = 'nearnei', par = FALSE)
	stoptime = Sys.time()
	difftime(stoptime, starttime, units="secs")
	partSizes_RMSPE[iter,15] = mean((Zp - cope200_fefe100_puloNN[,1])^2)
	partSizes_PI90[iter,15] = mean(cope200_fefe100_puloNN[,1] -
		qnorm(.95)*cope200_fefe100_puloNN[,2] < Zp & Zp < cope200_fefe100_puloNN[,1] + 
		qnorm(.95)*cope200_fefe100_puloNN[,2])

	starttime = Sys.time()
	cope200_fefe200_puloNN = pulo1(theta = cope_200$theta, 
    betahat = cope200_fefe200$bhat, covb = cope200_fefe200$covb,
    z = Z, X = X, Xp = Xp, xcoords = x, ycoords = y, 
    xcoordsp = xp, ycoordsp = yp, predmeth = 'nearnei', par = FALSE)
	stoptime = Sys.time()
	difftime(stoptime, starttime, units="secs")
	partSizes_RMSPE[iter,16] = mean((Zp - cope200_fefe200_puloNN[,1])^2)
	partSizes_PI90[iter,16] = mean(cope200_fefe200_puloNN[,1] -
		qnorm(.95)*cope200_fefe200_puloNN[,2] < Zp & Zp < cope200_fefe200_puloNN[,1] + 
		qnorm(.95)*cope200_fefe200_puloNN[,2])
		
	stopiter = Sys.time()
	itertime = difftime(stopiter, startiter, units="mins")
  cat("\n", "iteration time: ", itertime)
}

basepath = '/media/jay/data/2019_papers/fastREML/spPart_package/spPart/data/'
save(partSizes_est1, file = paste0(basepath,'partSizes_est1.rda'))
save(partSizes_est1se, file = paste0(basepath,'partSizes_est1se.rda'))
save(partSizes_est2, file = paste0(basepath,'partSizes_est2.rda'))
save(partSizes_est2se, file = paste0(basepath,'partSizes_est2se.rda'))
save(partSizes_RMSPE, file = paste0(basepath,'partSizes_RMSPE.rda'))
save(partSizes_PI90, file = paste0(basepath,'partSizes_PI90.rda'))
save(partSizes_copetime, file = paste0(basepath,'partSizes_copetime.rda'))
save(partSizes_fefetime, file = paste0(basepath,'partSizes_fefetime.rda'))
