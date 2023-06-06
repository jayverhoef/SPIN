#Rscript sim_scen1_500pts.R
q(save = 'no')
R

library('spPart')
library('nabor')
library('spNNGP')


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

nrep = 5
set.seed(7011)

BigSinSurf_time = NULL
iterList = c(500,1000,1500,2000,2500,3000,3500,4000,4500,5000,10000,15000,
  20000,30000,50000,75000,100000)
#iterList = c(500,1000,3000,20000)
overall_start = Sys.time()
for(rep in 1:nrep) {
for(iter in iterList) {
#-----------------------------------------------------------------------
#                       SIMULATE
#-----------------------------------------------------------------------
  cat("\n", "Replicate: ", rep)
	n = iter
  cat("\n", "Simulate sample size ", n)
  startiter = Sys.time()
	nrc = round(sqrt(n))
	N = nrc^2
	x = runif(n)
	y = runif(n)
	xysyst = pointSimSyst(nrow = nrc, ncol = nrc, lower.x.lim = 0, upper.x.lim = 1, 
			lower.y.lim = 0, upper.y.lim = 1)
	x = c(x, xysyst[,1])
	y = c(y, xysyst[,2])

	xyz1 = sim2DSinSurf(x,y) 
	z1std = (xyz1[,3] - mean(xyz1[,3]))/sqrt(var(xyz1[,3]))
	xyz2 = sim2DSinSurf(x,y) 
	z2std = (xyz2[,3] - mean(xyz2[,3]))/sqrt(var(xyz2[,3]))

	X1 = as.matrix(cbind(rep(1, times = n + N), rnorm(n + N), z2std))
	Z1 =  1*X1[,1]+ 1*X1[,2] + 1*X1[,3] + z1std + .1*rnorm(length(z1std),0,1)

	xp = x[(n+1):(n + N)]
	x = x[1:n]
	yp = y[(n+1):(n + N)]
	y = y[1:n]
	X1p = X1[(n+1):(n + N),]
	X1 = X1[1:n,]
	Z1p = Z1[(n+1):(n + N)]
	Z1 = Z1[1:n]

	starttime = Sys.time()
	gi3 = mass(x[1:n], y[1:n], ngroups = n/50, ssmeth = 'compKmean')[,2]
	stoptime = Sys.time()
	mapi_time = difftime(stoptime, starttime, units="secs")

  d1 = data.frame(z1 = Z1[1:n], x1 = X1[1:n,2], x2 = X1[1:n,3], 
  xcoords = x[1:n], ycoords = y[1:n], grpindx.comp = gi3)
	
#	library(viridis)
#	library(classInt)
#  f10 = classIntervals(dp$z1)
#  plot(dp$x, dp$y, pch = 19,
#		col = findColours(f10, viridis(length(f10[[2]]))), cex = .3)

#-----------------------------------------------------------------------
#                       COPE
#-----------------------------------------------------------------------
  cat("\n", "cope all")
  cope_all_time = NA
#	cope_all = cope(z1 ~ X1, data = d2, x_column = 'x', y_column = 'y', 
#		subSampCol = 'grpindx.one')

  if(n < 2501) {
    start_time = Sys.time()
    cope_all = cope(z1 ~ x1 + x2, data = d1, x_column = 'xcoords', y_column = 'ycoords', thetaini = c(3,3))
    stop_time = Sys.time()
    cope_all_time = difftime(stop_time, start_time, units="secs")
  }
	
  cat("\n", "cope compact")
	start_time = Sys.time()
	cope_comp = cope(z1 ~ x1 + x2, data = d1, x_column = 'xcoords', y_column = 'ycoords', subSampCol = 'grpindx.comp', thetaini = c(3,3))
	stop_time = Sys.time()
	cope_comp_time = difftime(stop_time, start_time, units="secs")



#-----------------------------------------------------------------------
#                      FEFE
#-----------------------------------------------------------------------

  cat("\n", "fefe")

	start_time = Sys.time()
  #undebug(fefe_covb)
	fefe_covb = fefe_covb(cope_out = cope_comp)
	stop_time = Sys.time()
	fefe_covb_time = difftime(stop_time, start_time, units="secs")
	coeftable(fefe_covb)

	start_time = Sys.time()
	fefe_covbalt = fefe_covb(cope_out = cope_comp, compute_covb = FALSE,
    compute_covbalt1 = TRUE, compute_covbalt2 = TRUE)
	stop_time = Sys.time()
	fefe_covbalt_time = difftime(stop_time, start_time, units="secs")
	


#-----------------------------------------------------------------------
#                      PULO
#-----------------------------------------------------------------------
  cat("\n", "pulo")
  
	start_time = Sys.time()
	cope200_fefe25_puloNN = pulo1(theta = cope_comp$theta, 
    betahat = fefe_covb$bhat, covb = fefe_covb$covb,
    z = Z1, X = X1, Xp = X1p, xcoords = x, ycoords = y, 
    xcoordsp = xp, ycoordsp = yp, predmeth = 'nearnei', par = FALSE)
	stop_time = Sys.time()
	pulo_time = difftime(stop_time, start_time, units="secs")
  
#-----------------------------------------------------------------------
#                      POBA
#-----------------------------------------------------------------------

#  cat("\n", "poba")
  		
#	start_time = Sys.time()
#	copa_comp_poba_NN = poba(cope_comp$theta, Z1, X1, X1p, x, y, xp, yp,
#		predmeth = 'NNdata')
#	stop_time = Sys.time()
#	poba_time = stop_time - start_time
	

################################################################################
################################################################################
#         spNNGP
################################################################################
################################################################################


phi.range = exp((-12:12)/4)
phi.range = phi.range - min(phi.range) + .001
alpha.range = exp((-12:12)/4)
alpha.range = alpha.range - min(alpha.range) + .001
g = length(phi.range)

theta.alpha <- cbind(
	kronecker(rep(1, times = g), phi.range),
	kronecker(alpha.range,rep(1, times = g))
)
colnames(theta.alpha) <- c("phi", "alpha")
sigma.sq.IG <- c(2, 10)

cat("\n", "NNGP","\n")

t
NNGP_time = difftime(stop_time, start_time, units="secs")

BigSinSurf_time = rbind(BigSinSurf_time,c(n, iter, rep,
  cope_all_time,
  mapi_time,
	cope_comp_time,
	fefe_covb_time,
  fefe_covbalt_time,
	pulo_time,
	NNGP_time
)
)
cat("\n")
}
}
colnames(BigSinSurf_time) = c('n', 'iter', 'rep', 'copeall', 'mapi','copecomp', 'fefe', 'fefealt', 'pulo','NNGP')

overall_end = Sys.time()

difftime(overall_end, overall_start, units="mins")

basepath = '/media/jay/data/desktop_data/2019_papers/fastREML/spPart_package/spPart/data/'

save(BigSinSurf_time, file = paste0(basepath,'BigSinSurf_time.rda'))

BSSall = aggregate(BigSinSurf_time, by = list(BigSinSurf_time[,'n']), FUN = mean)
#BSSall[5,'copeall'] = 1800
par(mar = c(5,5,1,1))
plot(BSSall[,'n'],BSSall[,'copeall'], type = 'l',
  ylab = 'Time (Secs)', xlab = 'Sample Size',
  cex.lab = 2, cex.axis = 1.5, lwd = 6, ylim = c(0,1800))
lines(BSSall[,'n'],(BSSall[,'mapi'] + BSSall[,'copecomp'] + 
  BSSall[,'fefe'] + BSSall[,'pulo']), 
  lwd = 6, col = 'green3')
lines(BSSall[,'n'],(BSSall[,'mapi'] + BSSall[,'copecomp'] + 
  BSSall[,'fefealt'] + BSSall[,'pulo']), 
  lwd = 6, lty = 2, col = 'green3')
lines(BSSall[,'n'],BSSall[,'NNGP'], 
  lwd = 6, col = 'red3')
legend(10000,1800, legend = c('Full Cov Matrix', 'Part. Cov Matrix', 'Part. Matrix AltCovb', 'NNGPconj'), lty = c(1,1,2,1), lwd = c(4,4,4,4),
  col = c('black','green3','green3','red3'), cex = 2.5)
