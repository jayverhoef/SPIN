#Rscript sim_scen1_500pts.R

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

storage = NULL
for(rep in 1:2) {
for(iter in 1:10) {

	cat("\n", "Rep number ", rep)
	n = iter*5000
  cat("\n", "Simulate sample size ", n)

	nrc = 200
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
#	xyX1 = sim2DSinSurf(x,y)
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

#	gi3 = mass(x[1:n], y[1:n], ngroups = n/50, ssmeth = 'compKmean')[,2]
	gi3 = mass(x[1:n], y[1:n], ngroups = n/50, ssmeth = 'random')[,2]

  d1 = data.frame(z1 = Z1[1:n], x = X1[1:n,2], xcoords = x[1:n],
    ycoords = y[1:n], grpindx.comp = gi3)
	
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
	
  cat("\n", "fefe1")
 
	start_time = Sys.time()
	cope_comp_fefe_comp = fefe(theta = cope_comp$theta, z = Z1, 
		X = X1, subsampindx = gi3, xcoords = x, ycoords = y)
	stop_time = Sys.time()
	fefe_time = difftime(stop_time, start_time, units = 'mins')

	start_time = Sys.time()
	#undebug(fefe1)
  cope_fefe1_comp = fefe1(cope_comp)
	stop_time = Sys.time()
	fefe1_time = difftime(stop_time, start_time, units = 'mins')
  round(sqrt(diag(cope_fefe1_comp$covb)),8)
  
	start_time = Sys.time()
  #undebug(fefe2)
	cope_fefe2_comp = fefe2(cope_comp)
	stop_time = Sys.time()
	fefe2_time = difftime(stop_time, start_time, units = 'mins')

p = dim(cope_fefe2_comp$bhatList[[1]])[1]
ngrps = length(cope_fefe2_comp$bhatList)
VE = matrix(0, nrow = p, ncol = p)
EV = VE
empbhat = rep(0, times = p)
for(i in 1:ngrps)
	empbhat = empbhat  + cope_fefe2_comp$bhatList[[i]]
empbhat = empbhat/ngrps
bias2 = (empbhat - cope_fefe2_comp$bhat) %*% t(empbhat - cope_fefe2_comp$bhat)
for(i in 1:ngrps) {
	VE = VE + (cope_fefe2_comp$bhatList[[i]] - empbhat) %*%
		t(cope_fefe2_comp$bhatList[[i]] - empbhat)
	EV = EV + cope_fefe2_comp$covbList[[i]]
}	
VE = VE/ngrps/(ngrps-1)
EV = EV/ngrps^2

cbind(empbhat,
 sqrt(diag(VE - bias2)))
cbind(empbhat,
 sqrt(diag(EV)))
 
coeftable(cope_comp_fefe_comp)
cbind(empbhat,
 sqrt(diag((VE - bias2)/2 + EV/2)))

sqrt(diag(EV*ngrps))[1]

prop

	storage = rbind(storage, 
	 c(n, cope_time, cope_comp$parmest$counts['function'], fefe1_time))
}
}

storage

plot(storage[40:49,c(1,4)], pch = 19, cex = 2, col = 'red')
points(storage[40:49,1:2], pch = 19, cex = 2, col = 'blue')

plot(storage[,c(1,3)])

