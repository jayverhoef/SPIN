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

iter = 1
set.seed(2003)
#-----------------------------------------------------------------------
#                       SIMULATE
#-----------------------------------------------------------------------

	cat("\n", "Rep number ", rep)
	n = iter*5000
  cat("\n", "Simulate sample size ", n)

	nrc = 100
	N = nrc^2
	x = runif(n)
	y = runif(n)
	xysyst = pointSimSyst(nrow = nrc, ncol = nrc, lower.x.lim = 0, upper.x.lim = 1, 
			lower.y.lim = 0, upper.y.lim = 1)
	x = c(x, xysyst[,1])
	y = c(y, xysyst[,2])

	xyz1 = sim2DSinSurf(x,y)
	zstd = (xyz1[,3] - mean(xyz1[,3]))/sqrt(var(xyz1[,3]))
	prop = runif(1)
	zstd = prop*zstd + (1 - prop)*rnorm(n + N, 0, 1)
	

#	X1 = as.matrix(cbind(rep(1, times = n + N), rnorm(n + N)))
	xyX1 = sim2DSinSurf(x,y)
	X1std = (xyX1[,3] - mean(xyX1[,3]))/sqrt(var(xyX1[,3]))
	X1 = as.matrix(cbind(rep(1, times = n + N), X1std))
	Z1 =  1*X1[,1]+ 1*X1[,2] + zstd

	xp = x[(n+1):(n + N)]
	x = x[1:n]
	yp = y[(n+1):(n + N)]
	y = y[1:n]
	X1p = X1[(n+1):(n + N),]
	X1 = X1[1:n,]
	Z1p = Z1[(n+1):(n + N)]
	Z1 = Z1[1:n]

	gi50 = mass(x[1:n], y[1:n], ngroups = n/50, ssmeth = 'compKmean')[,2]

  d1 = data.frame(z1 = Z1[1:n], x = X1[1:n,2], xcoords = x[1:n],
    ycoords = y[1:n], grpindx.comp = gi50)
	
#-----------------------------------------------------------------------
#                       COPE
#-----------------------------------------------------------------------

  cat("\n", "cope")
	
	start_time = Sys.time()
	#undebug(cope)
	cope_comp = cope(z1 ~ x, data = d1, x_column = 'xcoords', y_column = 'ycoords', 
		subSampCol = 'grpindx.comp', thetaini = c(3,3))
	stop_time = Sys.time()
	cope_time = difftime(stop_time, start_time, units = 'mins')
	
#-----------------------------------------------------------------------
#                      PULO
#-----------------------------------------------------------------------

  cat("\n", "pulo")

  # local estimates of betahat
	start_time = Sys.time()
	cope_c50_pulo_25 = pulo(cope_comp$theta, Z1, X1, X1p, x, y, xp, yp,
		predmeth = 'nearnei', par = FALSE, nNN = 25)
	stop_time = Sys.time()
	pulo_25_time = difftime(stop_time, start_time, units = 'mins')
	head(cope_c50_pulo_25)

	
  # global estimates of betahat (requires covb for global estimate)
	start_time = Sys.time()
	fefe_comp = fefe(theta = cope_comp$theta, z = Z1, 
		X = X1, subsampindx = gi3, xcoords = x, ycoords = y)
	cope_c50_pulo1_25 = pulo1(cope_comp$theta, 
    betahat = fefe_comp$bhat, covb = fefe_comp$covb, Z1, X1, X1p, x, y, xp, yp,
		predmeth = 'nearnei', nNN=25, par = FALSE)
	stop_time = Sys.time()
	pulo1_time = difftime(stop_time, start_time, units = 'mins')
	head(cope_c50_pulo1_25)


  cat("\n", "poba")
  
	start_time = Sys.time()
	cope_c50_poba_25 = poba(cope_comp$theta, Z1, X1, X1p, x, y, xp, yp,
		predmeth = 'NNdata', nNN = 25)
	stop_time = Sys.time()
	poba_time = difftime(stop_time, start_time, units = 'mins')
 
	start_time = Sys.time()
	cope_c50_poba_50 = poba(cope_comp$theta, Z1, X1, X1p, x, y, xp, yp,
		predmeth = 'NNdata', nNN = 50)
	stop_time = Sys.time()
	poba_time = difftime(stop_time, start_time, units = 'mins')

poba_time

#sqrt(copa_all_poba_all$junk1 - copa_all_poba_all$junk2 + copa_all_poba_all$junk3)
#sqrt(copa_all_poba1_all$junk1 - copa_all_poba1_all$junk2 + copa_all_poba1_all$junk3)
