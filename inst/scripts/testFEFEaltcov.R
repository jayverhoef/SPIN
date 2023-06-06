#Rscript sim_scen1_500pts.R
q(save = 'no')
R
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


  n = 10000

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

  # make group sizes 
  grpsize = 50
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

	starttime = Sys.time()
	cope_comp = cope(z ~ X1 + X2, data = d2, x_column = 'x', y_column = 'y', 
		subSampCol = 'grpindx.comp', thetaini = c(3,3))
	stoptime = Sys.time()
	difftime(stoptime, starttime, units="secs")

#-----------------------------------------------------------------------
#                      FEFE
#-----------------------------------------------------------------------

 
	starttime = Sys.time()
  #undebug(fefe)
	cope_comp_fefe_rand = fefe(theta = cope_comp$theta, z = Z, 
		X = X, subsampindx = gic[,2], xcoords = x, ycoords = y,
    compute_covbalt1 = TRUE, compute_covbalt2 = TRUE)
	stoptime = Sys.time()
	difftime(stoptime, starttime, units="secs")
  
covb = cope_comp_fefe_rand$covb
covbalt1 = cope_comp_fefe_rand$covbalt1
covbalt2 = cope_comp_fefe_rand$covbalt2
covbalt2[1,1] = covbalt2[1,1]*n/grpsize

round(covb, digits = 5)
round(covbalt1, digits = 5)
round(covbalt2, digits = 5)

round(cbind(sqrt(diag(covb)),
  sqrt(diag(covbalt1)),
  sqrt(diag(covbalt2))),digits=5)
  
