#Rscript sim_scen1_500pts.R
q(save = 'no')
R
library('spPart')
library('nabor')
library('spNNGP')
library('pdist')

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

set.seed(7014)
nsim = 1000
test_pulo_poba_nNNsize = matrix(NA, nrow = nsim, ncol = 36)

iter = 1
for(iter in 1:nsim) {
#-----------------------------------------------------------------------
#                       SIMULATE
#-----------------------------------------------------------------------

  cat("\n", "Iteration: ", iter)

  U = runif(1)
  if(U < .5)
	n = 1000 + round(runif(1)*1000)
  if(U >= .5)
  n = 2000 + round(runif(1)*8000)
  cat("\n", "Simulate sample size: ", n)

	nrc = 40
	N = nrc^2
	x = runif(n)
	y = runif(n)
	xysyst = pointSimSyst(nrow = nrc, ncol = nrc, lower.x.lim = 0, upper.x.lim = 1, 
			lower.y.lim = 0, upper.y.lim = 1)
	x = c(x, xysyst[,1])
	y = c(y, xysyst[,2])

  if(U < .5) {
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

  if(U >= .5) {
    xyz1 = sim2DSinSurf(x, y)
    zstd1 = (xyz1[,3] - mean(xyz1[,3]))/sqrt(var(xyz1[,3]))
    xyz2 = sim2DSinSurf(x, y)
    zstd2 = (xyz2[,3] - mean(xyz2[,3]))/sqrt(var(xyz2[,3]))  
  }

  # design matrix with X1 spatially independent and X1 spatially patterned
	X = as.matrix(cbind(rep(1, times = n + nrc^2), rnorm(n + nrc^2), zstd2))
  colnames(X) = c('int','X1','X2')

  # make response variable with random proportion of nugget
  prop = runif(1)
	Z = 1*X[,1]+ 1*X[,2] + 1*X[,3] + 
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

	gi50 = mass(x[1:n], y[1:n], ngroups = n/50, ssmeth = 'compKmean')[,2]

  d1 = data.frame(z = Z, X1 = X1, X2 = X2, xcoords = x,
    ycoords = y, grpindx.comp = gi50)
	
#-----------------------------------------------------------------------
#                       COPE
#-----------------------------------------------------------------------

  cat("\n", "cope")
	
	start_time = Sys.time()
	#undebug(cope)
	cope_comp = cope(z ~ X1 + X2, data = d1, x_column = 'xcoords', y_column = 'ycoords', 
		subSampCol = 'grpindx.comp', thetaini = c(3,3))
	stop_time = Sys.time()
	cope_time = difftime(stop_time, start_time, units = 'mins')
	
#-----------------------------------------------------------------------
#                      PULO
#-----------------------------------------------------------------------

  cat("\n", "pulo25")

  # local estimates of betahat
	start_time = Sys.time()
	cope_c50_pulo_25 = pulo(cope_comp$theta, Z, X, Xp, x, y, xp, yp,
		predmeth = 'nearnei', par = FALSE, nNN = 25)
	stop_time = Sys.time()
	pulo_25_time = difftime(stop_time, start_time, units = 'mins')
  test_pulo_poba_nNNsize[iter,1] = mean((Zp - cope_c50_pulo_25[,1])^2)
	test_pulo_poba_nNNsize[iter,2] = mean(cope_c50_pulo_25[,1] -
		qnorm(.95)*cope_c50_pulo_25[,2] < Zp & Zp < cope_c50_pulo_25[,1] + 
		qnorm(.95)*cope_c50_pulo_25[,2])

	
  # global estimates of betahat (requires covb for global estimate)
	start_time = Sys.time()
	fefe_comp = fefe(theta = cope_comp$theta, z = Z, 
		X = X, subsampindx = gi50, xcoords = x, ycoords = y)
	cope_c50_pulo1_25 = pulo1(cope_comp$theta, 
    betahat = fefe_comp$bhat, covb = fefe_comp$covb, Z, X, Xp, x, y, xp, yp,
		predmeth = 'nearnei', nNN=25, par = FALSE)
	stop_time = Sys.time()
	pulo1_25_time = difftime(stop_time, start_time, units = 'mins')
  test_pulo_poba_nNNsize[iter,3] = mean((Zp - cope_c50_pulo1_25[,1])^2)
	test_pulo_poba_nNNsize[iter,4] = mean(cope_c50_pulo1_25[,1] -
		qnorm(.95)*cope_c50_pulo1_25[,2] < Zp & Zp < cope_c50_pulo1_25[,1] + 
		qnorm(.95)*cope_c50_pulo1_25[,2])

  cat("\n", "pulo50")

  # local estimates of betahat
	start_time = Sys.time()
	cope_c50_pulo_50 = pulo(cope_comp$theta, Z, X, Xp, x, y, xp, yp,
		predmeth = 'nearnei', par = FALSE, nNN = 50)
	stop_time = Sys.time()
	pulo_50_time = difftime(stop_time, start_time, units = 'mins')
  test_pulo_poba_nNNsize[iter,5] = mean((Zp - cope_c50_pulo_50[,1])^2)
	test_pulo_poba_nNNsize[iter,6] = mean(cope_c50_pulo_50[,1] -
		qnorm(.95)*cope_c50_pulo_50[,2] < Zp & Zp < cope_c50_pulo_50[,1] + 
		qnorm(.95)*cope_c50_pulo_50[,2])
	
  # global estimates of betahat (requires covb for global estimate)
	start_time = Sys.time()
	fefe_comp = fefe(theta = cope_comp$theta, z = Z, 
		X = X, subsampindx = gi50, xcoords = x, ycoords = y)
	cope_c50_pulo1_50 = pulo1(cope_comp$theta, 
    betahat = fefe_comp$bhat, covb = fefe_comp$covb, Z, X, Xp, x, y, xp, yp,
		predmeth = 'nearnei', nNN=50, par = FALSE)
	stop_time = Sys.time()
	pulo1_50_time = difftime(stop_time, start_time, units = 'mins')
  test_pulo_poba_nNNsize[iter,7] = mean((Zp - cope_c50_pulo1_50[,1])^2)
	test_pulo_poba_nNNsize[iter,8] = mean(cope_c50_pulo1_50[,1] -
		qnorm(.95)*cope_c50_pulo1_50[,2] < Zp & Zp < cope_c50_pulo1_50[,1] + 
		qnorm(.95)*cope_c50_pulo1_50[,2])

  cat("\n", "pulo100")

  # local estimates of betahat
	start_time = Sys.time()
	cope_c50_pulo_100 = pulo(cope_comp$theta, Z, X, Xp, x, y, xp, yp,
		predmeth = 'nearnei', par = FALSE, nNN = 100)
	stop_time = Sys.time()
	pulo_100_time = difftime(stop_time, start_time, units = 'mins')
  test_pulo_poba_nNNsize[iter,9] = mean((Zp - cope_c50_pulo_100[,1])^2)
	test_pulo_poba_nNNsize[iter,10] = mean(cope_c50_pulo_100[,1] -
		qnorm(.95)*cope_c50_pulo_100[,2] < Zp & Zp < cope_c50_pulo_100[,1] + 
		qnorm(.95)*cope_c50_pulo_100[,2])
	
  # global estimates of betahat (requires covb for global estimate)
	start_time = Sys.time()
	fefe_comp = fefe(theta = cope_comp$theta, z = Z, 
		X = X, subsampindx = gi50, xcoords = x, ycoords = y)
	cope_c50_pulo1_100 = pulo1(cope_comp$theta, 
    betahat = fefe_comp$bhat, covb = fefe_comp$covb, Z, X, Xp, x, y, xp, yp,
		predmeth = 'nearnei', nNN=100, par = FALSE)
	stop_time = Sys.time()
	pulo1_100_time = difftime(stop_time, start_time, units = 'mins')
  test_pulo_poba_nNNsize[iter,11] = mean((Zp - cope_c50_pulo1_100[,1])^2)
	test_pulo_poba_nNNsize[iter,12] = mean(cope_c50_pulo1_100[,1] -
		qnorm(.95)*cope_c50_pulo1_100[,2] < Zp & Zp < cope_c50_pulo1_100[,1] + 
		qnorm(.95)*cope_c50_pulo1_100[,2])

  cat("\n", "pulo200")

  # local estimates of betahat
	start_time = Sys.time()
	cope_c50_pulo_200 = pulo(cope_comp$theta, Z, X, Xp, x, y, xp, yp,
		predmeth = 'nearnei', par = FALSE, nNN = 200)
	stop_time = Sys.time()
	pulo_200_time = difftime(stop_time, start_time, units = 'mins')
  test_pulo_poba_nNNsize[iter,13] = mean((Zp - cope_c50_pulo_200[,1])^2)
	test_pulo_poba_nNNsize[iter,14] = mean(cope_c50_pulo_200[,1] -
		qnorm(.95)*cope_c50_pulo_200[,2] < Zp & Zp < cope_c50_pulo_200[,1] + 
		qnorm(.95)*cope_c50_pulo_200[,2])
	
  # global estimates of betahat (requires covb for global estimate)
	start_time = Sys.time()
	fefe_comp = fefe(theta = cope_comp$theta, z = Z, 
		X = X, subsampindx = gi50, xcoords = x, ycoords = y)
	cope_c50_pulo1_200 = pulo1(cope_comp$theta, 
    betahat = fefe_comp$bhat, covb = fefe_comp$covb, Z, X, Xp, x, y, xp, yp,
		predmeth = 'nearnei', nNN=200, par = FALSE)
	stop_time = Sys.time()
	pulo1_200_time = difftime(stop_time, start_time, units = 'mins')
  test_pulo_poba_nNNsize[iter,15] = mean((Zp - cope_c50_pulo1_200[,1])^2)
	test_pulo_poba_nNNsize[iter,16] = mean(cope_c50_pulo1_200[,1] -
		qnorm(.95)*cope_c50_pulo1_200[,2] < Zp & Zp < cope_c50_pulo1_200[,1] + 
		qnorm(.95)*cope_c50_pulo1_200[,2])

#-----------------------------------------------------------------------
#                      POBA
#-----------------------------------------------------------------------

  cat("\n", "poba25")
  
	start_time = Sys.time()
	cope_c50_poba_25 = poba(cope_comp$theta, Z, X, Xp, x, y, xp, yp,
		predmeth = 'NNdata', nNN = 25)
	stop_time = Sys.time()
	poba_25_time = difftime(stop_time, start_time, units = 'mins')
  test_pulo_poba_nNNsize[iter,17] = (mean(Zp) - cope_c50_poba_25$block_pred)^2
	test_pulo_poba_nNNsize[iter,18] = (cope_c50_poba_25$block_pred -
		qnorm(.95)*cope_c50_poba_25$block_pred_se < mean(Zp) & 
    mean(Zp) < cope_c50_poba_25$block_pred + 
    qnorm(.95)*cope_c50_poba_25$block_pred_se)*1
 
  cat("\n", "poba50")

	start_time = Sys.time()
	cope_c50_poba_50 = poba(cope_comp$theta, Z, X, Xp, x, y, xp, yp,
		predmeth = 'NNdata', nNN = 50)
	stop_time = Sys.time()
	poba_50_time = difftime(stop_time, start_time, units = 'mins')
  test_pulo_poba_nNNsize[iter,19] = (mean(Zp) - cope_c50_poba_50$block_pred)^2
	test_pulo_poba_nNNsize[iter,20] = (cope_c50_poba_50$block_pred -
		qnorm(.95)*cope_c50_poba_50$block_pred_se < mean(Zp) & 
    mean(Zp) < cope_c50_poba_50$block_pred + 
    qnorm(.95)*cope_c50_poba_50$block_pred_se)*1

  cat("\n", "poba100")

	start_time = Sys.time()
	cope_c50_poba_100 = poba(cope_comp$theta, Z, X, Xp, x, y, xp, yp,
		predmeth = 'NNdata', nNN = 100)
	stop_time = Sys.time()
	poba_100_time = difftime(stop_time, start_time, units = 'mins')
  test_pulo_poba_nNNsize[iter,21] = (mean(Zp) - cope_c50_poba_100$block_pred)^2
	test_pulo_poba_nNNsize[iter,22] = (cope_c50_poba_100$block_pred -
		qnorm(.95)*cope_c50_poba_100$block_pred_se < mean(Zp) & 
    mean(Zp) < cope_c50_poba_100$block_pred + 
    qnorm(.95)*cope_c50_poba_100$block_pred_se)*1

  cat("\n", "poba200")

	start_time = Sys.time()
	cope_c50_poba_200 = poba(cope_comp$theta, Z, X, Xp, x, y, xp, yp,
		predmeth = 'NNdata', nNN = 200)
	stop_time = Sys.time()
	poba_200_time = difftime(stop_time, start_time, units = 'mins')
  test_pulo_poba_nNNsize[iter,23] = (mean(Zp) - cope_c50_poba_200$block_pred)^2
	test_pulo_poba_nNNsize[iter,24] = (cope_c50_poba_200$block_pred -
		qnorm(.95)*cope_c50_poba_200$block_pred_se < mean(Zp) & 
    mean(Zp) < cope_c50_poba_200$block_pred + 
    qnorm(.95)*cope_c50_poba_200$block_pred_se)*1
    
  test_pulo_poba_nNNsize[iter,25] = pulo_25_time
  test_pulo_poba_nNNsize[iter,26] = pulo1_25_time
  test_pulo_poba_nNNsize[iter,27] = pulo_50_time
  test_pulo_poba_nNNsize[iter,28] = pulo1_50_time
  test_pulo_poba_nNNsize[iter,29] = pulo_100_time
  test_pulo_poba_nNNsize[iter,30] = pulo1_100_time
  test_pulo_poba_nNNsize[iter,31] = pulo_200_time
  test_pulo_poba_nNNsize[iter,32] = pulo1_200_time
  test_pulo_poba_nNNsize[iter,33] = poba_25_time
  test_pulo_poba_nNNsize[iter,34] = poba_50_time
  test_pulo_poba_nNNsize[iter,35] = poba_100_time
  test_pulo_poba_nNNsize[iter,36] = poba_200_time
}

basepath = '/media/jay/data/desktop_data/2019_papers/fastREML/spPart_package/spPart/data/'
save(test_pulo_poba_nNNsize, file = paste0(basepath,'test_pulo_poba_nNNsize.rda'))

means = apply(test_pulo_poba_nNNsize,2,mean)
pulo_poba_nNNsize_table = rbind(
  means[c(1,3,2,4,17,18,25,26,33)],
  means[c(5,7,6,8,19,20,27,28,34)],
  means[c(9,11,10,12,21,22,29,30,35)],
  means[c(13,15,14,16,23,24,31,32,36)]
)
pulo_poba_nNNsize_table[,c(7,8,9)] = pulo_poba_nNNsize_table[,c(7,8,9)]*60
library(xtable)
  print(
    xtable(pulo_poba_nNNsize_table, 
      align = c('l',rep('l', times = length(pulo_poba_nNNsize_table[1,]))),
      digits = c(0,1,1,3,3,4,3,1,1,1),
    ),
    size = 'footnotesize',
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
  )
