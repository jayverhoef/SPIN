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

nsim = 100

fefealt_est1 = matrix(NA, nrow = nsim, ncol = 17)
fefealt_est2 = matrix(NA, nrow = nsim, ncol = 17)

set.seed(4009)
for(iter in 1:nsim) {
#-----------------------------------------------------------------------
#                       SIMULATE
#-----------------------------------------------------------------------
  cat("\n", "Simulation: ", iter)
  cat("\n", "Simulate")
  startiter = Sys.time()

  n = 20000 + round((runif(1)^2)*80000)

	x = runif(n)
	y = runif(n)
	xysyst = pointSimSyst(nrow = 40, ncol = 40, lower.x.lim = 0, upper.x.lim = 1, 
			lower.y.lim = 0, upper.y.lim = 1)
	x = c(x, xysyst[,1])
	y = c(y, xysyst[,2])


	xyz1 = sim2DSinSurf(x, y)
	zstd1 = (xyz1[,3] - mean(xyz1[,3]))/sqrt(var(xyz1[,3]))
	xyz2 = sim2DSinSurf(x, y)
	zstd2 = (xyz2[,3] - mean(xyz2[,3]))/sqrt(var(xyz2[,3]))  

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
  cat("\n", "cope25 for sample size", n)
	cope_25 = cope(z ~ X1 + X2, data = d2, x_column = 'x', y_column = 'y', 
		subSampCol = 'grpindx.25', thetaini = c(3,3))

  cat("\n", "cope50 for sample size", n)
	cope_50 = cope(z ~ X1 + X2, data = d2, x_column = 'x', y_column = 'y', 
		subSampCol = 'grpindx.50', thetaini = c(3,3))

  cat("\n", "cope100 for sample size", n)
	cope_100 = cope(z ~ X1 + X2, data = d2, x_column = 'x', y_column = 'y', 
		subSampCol = 'grpindx.100', thetaini = c(3,3))

  cat("\n", "cope200 for sample size", n)
	cope_200 = cope(z ~ X1 + X2, data = d2, x_column = 'x', y_column = 'y', 
		subSampCol = 'grpindx.200', thetaini = c(3,3))

#-----------------------------------------------------------------------
#                      FEFE
#-----------------------------------------------------------------------
  
  cat("\n", "fefe25 for sample size", n)
	fefe_covbalt25 = fefe_covb(cope_out = cope_25, compute_covb = TRUE,
    compute_covbalt1 = TRUE, compute_covbalt2 = TRUE)

  cat("\n", "fefe50 for sample size", n)
  #undebug(spPart:::fefe_covb)
	fefe_covbalt50 = fefe_covb(cope_out = cope_50, compute_covb = TRUE,
    compute_covbalt1 = TRUE, compute_covbalt2 = TRUE)

  cat("\n", "fefe100 for sample size", n)
	fefe_covbalt100 = fefe_covb(cope_out = cope_100, compute_covb = TRUE,
    compute_covbalt1 = TRUE, compute_covbalt2 = TRUE)

  cat("\n", "fefe200 for sample size", n)
	fefe_covbalt200 = fefe_covb(cope_out = cope_200, compute_covb = TRUE,
    compute_covbalt1 = TRUE, compute_covbalt2 = TRUE)

	fefealt_est1[iter,1] = fefe_covbalt25$bhat[2]
	fefealt_est1[iter,2] = sqrt(fefe_covbalt25$covb[2,2])
	fefealt_est1[iter,3] = sqrt(fefe_covbalt25$covbalt1[2,2])
	fefealt_est1[iter,4] = sqrt(fefe_covbalt25$covbalt2[2,2])
	fefealt_est1[iter,5] = fefe_covbalt50$bhat[2]
	fefealt_est1[iter,6] = sqrt(fefe_covbalt50$covb[2,2])
	fefealt_est1[iter,7] = sqrt(fefe_covbalt50$covbalt1[2,2])
	fefealt_est1[iter,8] = sqrt(fefe_covbalt50$covbalt2[2,2])
	fefealt_est1[iter,9] = fefe_covbalt100$bhat[2]
	fefealt_est1[iter,10] = sqrt(fefe_covbalt100$covb[2,2])
	fefealt_est1[iter,11] = sqrt(fefe_covbalt100$covbalt1[2,2])
	fefealt_est1[iter,12] = sqrt(fefe_covbalt100$covbalt2[2,2])
	fefealt_est1[iter,13] = fefe_covbalt200$bhat[2]
	fefealt_est1[iter,14] = sqrt(fefe_covbalt200$covb[2,2])
	fefealt_est1[iter,15] = sqrt(fefe_covbalt200$covbalt1[2,2])
	fefealt_est1[iter,16] = sqrt(fefe_covbalt200$covbalt2[2,2])
	fefealt_est1[iter,17] = n

	fefealt_est2[iter,1] = fefe_covbalt25$bhat[3]
	fefealt_est2[iter,2] = sqrt(fefe_covbalt25$covb[3,3])
	fefealt_est2[iter,3] = sqrt(fefe_covbalt25$covbalt1[3,3])
	fefealt_est2[iter,4] = sqrt(fefe_covbalt25$covbalt2[3,3])
	fefealt_est2[iter,5] = fefe_covbalt50$bhat[3]
	fefealt_est2[iter,6] = sqrt(fefe_covbalt50$covb[3,3])
	fefealt_est2[iter,7] = sqrt(fefe_covbalt50$covbalt1[3,3])
	fefealt_est2[iter,8] = sqrt(fefe_covbalt50$covbalt2[3,3])
	fefealt_est2[iter,9] = fefe_covbalt100$bhat[3]
	fefealt_est2[iter,10] = sqrt(fefe_covbalt100$covb[3,3])
	fefealt_est2[iter,11] = sqrt(fefe_covbalt100$covbalt1[3,3])
	fefealt_est2[iter,12] = sqrt(fefe_covbalt100$covbalt2[3,3])
	fefealt_est2[iter,13] = fefe_covbalt200$bhat[3]
	fefealt_est2[iter,14] = sqrt(fefe_covbalt200$covb[3,3])
	fefealt_est2[iter,15] = sqrt(fefe_covbalt200$covbalt1[3,3])
	fefealt_est2[iter,16] = sqrt(fefe_covbalt200$covbalt2[3,3])
	fefealt_est2[iter,17] = n
  cat("\n")
}

fefealt_est1_9 = fefealt_est1
fefealt_est2_9 = fefealt_est2

basepath = '/media/jay/data/desktop_data/2019_papers/fastREML/spPart_package/spPart/data/'
save(fefealt_est1_9, file = paste0(basepath,'fefealt_est1_9.rda'))
save(fefealt_est2_9, file = paste0(basepath,'fefealt_est2_9.rda'))
