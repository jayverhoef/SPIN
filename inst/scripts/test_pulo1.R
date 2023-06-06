#Rscript sim_scen1_500pts.R
#q(save = 'no')
#R

library('spNNGP')
library('spPart')
library('viridis')
library('classInt')
library('nabor')

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

set.seed(2001)
nsim = 20

scen1_RMSPE = matrix(NA, nrow = nsim, ncol = 3)
scen1_PI90 = matrix(NA, nrow = nsim, ncol = 3)

iter = 1
for(iter in 1:nsim) {
#-----------------------------------------------------------------------
#                       SIMULATE
#-----------------------------------------------------------------------
  cat("\n", "Simulation: ", iter)
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

	gi2 = mass(x, y, ngroups = 20, ssmeth = 'random')

	gi3 = mass(x, y, ngroups = 20, ssmeth = 'compKmean')
		
	gi4 = mass(x, y, ngroups = 20, ssmeth = 'zimmer')

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

	cope_comp = cope(z1 ~ X1, data = d2, x_column = 'x', y_column = 'y', 
		subSampCol = 'grpindx.comp', thetaini = c(3,3))


#-----------------------------------------------------------------------
#                      FEFE
#-----------------------------------------------------------------------

  cat("\n", "fefe")


	cope_comp_fefe_comp = fefe(theta = cope_comp$theta, z = Z1, 
		X = X1, subsampindx = gi3[,2], xcoords = x, ycoords = y)


#-----------------------------------------------------------------------
#                      PULO
#-----------------------------------------------------------------------
  cat("\n", "pulo")
  
	cope_comp_pulo_nene = pulo(cope_comp$theta, Z1, X1, X1p, x, y, xp, yp,
		predmeth = 'nearnei', par = FALSE)

	cope_comp_pulo1_nene = pulo1(cope_comp$theta, betahat = cope_comp_fefe_comp$bhat, 
		covb = cope_comp_fefe_comp$covb, Z1, X1, X1p, x, y, xp, yp,
		predmeth = 'nearnei', par = FALSE)

#-----------------------------------------------------------------------
#                      NNGP
#-----------------------------------------------------------------------

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

NNGPout_rmspe <- spConjNNGP(Z1 ~ X1-1, coords=cbind(x,y), n.neighbors = 29,
                  X.0 = X1p, coords.0 = cbind(xp,yp),
                  k.fold = 5, score.rule = "rmspe",
                  n.omp.threads = 1,
                  theta.alpha = theta.alpha, sigma.sq.IG = sigma.sq.IG, 
                  cov.model = 'exponential')

		
	scen1_RMSPE[iter,1] = mean((Z1p - cope_comp_pulo_nene[,1])^2)
	scen1_PI90[iter,1] = mean(cope_comp_pulo_nene[,1] - 
		qnorm(.95)*cope_comp_pulo_nene[,2] < Z1p & Z1p < cope_comp_pulo_nene[,1] + 
		qnorm(.95)*cope_comp_pulo_nene[,2])
	scen1_RMSPE[iter,2] = mean((Z1p - cope_comp_pulo1_nene[,1])^2)
	scen1_PI90[iter,2] = mean(cope_comp_pulo1_nene[,1] - 
		qnorm(.95)*cope_comp_pulo1_nene[,2] < Z1p & Z1p < cope_comp_pulo1_nene[,1] + 
		qnorm(.95)*cope_comp_pulo1_nene[,2])
	scen1_RMSPE[iter,3] = mean((Z1p - NNGPout_rmspe$y.0.hat)^2)
	scen1_PI90[iter,3] = mean(NNGPout_rmspe$y.0.hat - 
		qnorm(.95)*sqrt(NNGPout_rmspe$y.0.hat.var) < Z1p & Z1p < NNGPout_rmspe$y.0.hat + 
		qnorm(.95)*sqrt(NNGPout_rmspe$y.0.hat.var))

	stopiter = Sys.time()
	itertime = difftime(stopiter, startiter, units="mins")
  cat("\n", "iteration time: ", itertime)
}
cat("\n")

apply(scen1_RMSPE,2,mean)
apply(scen1_PI90,2,mean)

basepath = '/home/xverhoef/data/2019_papers/fastREML/spPart_package/spPart/data/'

save(scen1_RMSPE, file = paste0(basepath,'scen1_RMSPE.rda'))
save(scen1_PI90, file = paste0(basepath,'scen1_PI90.rda'))

#load(paste0(basepath,'scen1_1000pt_NNGP_est.rda'))

