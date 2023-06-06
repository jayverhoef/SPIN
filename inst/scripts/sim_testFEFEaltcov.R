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

nsim = 2000

covbType_est1 = matrix(NA, nrow = nsim, ncol = 1)
covbType_est1se = matrix(NA, nrow = nsim, ncol = 3)
covbType_est2 = matrix(NA, nrow = nsim, ncol = 1)
covbType_est2se = matrix(NA, nrow = nsim, ncol = 3)

set.seed(2095)
for(iter in 1:nsim) {
#-----------------------------------------------------------------------
#                       SIMULATE
#-----------------------------------------------------------------------
  cat("\n", "Simulation: ", iter)
  cat("\n", "Simulate")
  startiter = Sys.time()

	x = runif(n)
	y = runif(n)
	xyz1 = sim2DSinSurf(x, y)
	zstd1 = (xyz1[,3] - mean(xyz1[,3]))/sqrt(var(xyz1[,3]))
	xyz2 = sim2DSinSurf(x, y)
	zstd2 = (xyz2[,3] - mean(xyz2[,3]))/sqrt(var(xyz2[,3]))  


  # design matrix with X1 spatially independent and X1 spatially patterned
	X = as.matrix(cbind(rep(1, times = n), rnorm(n), zstd2))
  colnames(X) = c('int','X1','X2')

  # make response variable with random proportion of nugget
  prop = runif(1)
	Z =  1*X[,1]+ 1*X[,2] + 1*X[,3] + 
    10*prop*zstd1 + 10*(1 - prop)*rnorm(n, 0, 1)

	X1 = X[,2]
	X2 = X[,3]



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

#-----------------------------------------------------------------------
#                       COPE
#-----------------------------------------------------------------------

  cat("\n", "cope for sample size", n)

	starttime = Sys.time()
	cope_comp = cope(z ~ X1 + X2, data = d2, x_column = 'x', y_column = 'y', 
		subSampCol = 'grpindx.comp', thetaini = c(3,3))
	stoptime = Sys.time()
	difftime(stoptime, starttime, units="secs")

#-----------------------------------------------------------------------
#                      FEFE
#-----------------------------------------------------------------------
 cat("\n", "fefe")

 
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


ses = cbind(sqrt(diag(covb)),
  sqrt(diag(covbalt1)),
  sqrt(diag(covbalt2)))
 
covbType_est1[iter] = cope_comp_fefe_rand$bhat[2]
covbType_est1se[iter,] = ses[2,]
covbType_est2[iter] = cope_comp_fefe_rand$bhat[3]
covbType_est2se[iter,] = ses[3,]

}

basepath = '/media/jay/data/2019_papers/fastREML/spPart_package/spPart/data/'
save(covbType_est1, file = paste0(basepath,'covbType_est1.rda'))
save(covbType_est1se, file = paste0(basepath,'covbType_est1se.rda'))
save(covbType_est2, file = paste0(basepath,'covbType_est2.rda'))
save(covbType_est2se, file = paste0(basepath,'covbType_est2se.rda'))

mean(covbType_est1)
covbT1se = covbType_est1se
covbT2se = covbType_est2se
covbT1se[is.na(covbT1se[,2]),2] =
  covbT1se[is.na(covbT1se[,2]),3] 
covbT2se[is.na(covbT2se[,2]),2] =
  covbT2se[is.na(covbT2se[,2]),3] 

mean(covbType_est1 - qnorm(.95)*covbT1se[,1] < 1 & 1 <
  covbType_est1 + qnorm(.95)*covbT1se[,1])
mean(covbType_est1 - qnorm(.95)*covbT1se[,2] < 1 & 1 <
  covbType_est1 + qnorm(.95)*covbT1se[,2])
mean(covbType_est1 - qnorm(.95)*covbT1se[,3] < 1 & 1 <
  covbType_est1 + qnorm(.95)*covbT1se[,3])
  
mean(covbType_est2 - qnorm(.95)*covbT2se[,1] < 1 & 1 <
  covbType_est2 + qnorm(.95)*covbT2se[,1])
mean(covbType_est2 - qnorm(.95)*covbT2se[,2] < 1 & 1 <
  covbType_est2 + qnorm(.95)*covbT2se[,2])
mean(covbType_est2 - qnorm(.95)*covbT2se[,3] < 1 & 1 <
  covbType_est2 + qnorm(.95)*covbT2se[,3])
