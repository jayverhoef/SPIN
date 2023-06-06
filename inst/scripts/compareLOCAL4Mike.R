#Rscript sim_scen1_500pts.R
#q(save = 'no')
#R
library('spmodel')
library('viridis')
library('classInt')

ScriptPath = paste0('/mnt/ExtraDrive1/Work/desktop_data/2019_papers',
	'/fastREML/spPart_package/spPart/R/')
source(paste0(ScriptPath, 'cope_simple.R'))
source(paste0(ScriptPath, 'm2LLg_simple.R'))

# a function to color the simulated data to visuall inspect how autocorrelated
# they are
colorpoints = function(x, y, z, nclass, cex = 1)
{
	cip = classIntervals(z, n = nclass, style = 'fisher')
	palp = viridis(nclass)
	cip_colors = findColours(cip, palp)
	old.par = par(mar = c(5,5,1,1))
	plot(x, y, col = cip_colors, pch = 19, cex = cex,
		xlab = 'x', ylab = 'y', cex.lab = 2, cex.axis = 1.5)
	par(old.par)
}

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#              FIRST SIMULATION, SINGLE GROUP
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
	
	# set seed for reproducibility
	set.seed(1001)
	#set sample size
	n = 200
	# simulate random spatial x-coordinates
	xcoord = runif(n)
	# simulate random spatial x-coordinates
	ycoord = runif(n)
	# make a data.frame from coordinates
	xyDF = data.frame(x = xcoord, y = ycoord)
	
	# chose a range parameter and simulate pure autocorrelation with variance 1
  range = 1
  # specify covariance values
	spcov_params_val <- spcov_params("spherical", 
		de = 1, ie = .00001, range = range)
	# simulate the data for response variable
	z = sprnorm(spcov_params_val, data = xyDF, 
		xcoord = x, ycoord = y)
	xyz_4y = data.frame(xyDF, z = z)
	# simulate an autocorrelated covariate as well
	z = sprnorm(spcov_params_val, data = xyDF, 
		xcoord = x, ycoord = y)
	xyz_4x = data.frame(xyDF, z = z)

  # design matrix with X1 spatially independent and X1 spatially patterned
	X = as.matrix(cbind(rep(1, times = n), rnorm(n), xyz_4x$z))
  colnames(X) = c('int','X1','X2')

  # make response variable with a specified standard deviation for
  # autocorrelated error and another for independent error
  de = 10
  ie = 1
	Z =  1*X[,1]+ 1*X[,2] + 1*X[,3] + 
    de*xyz_4y$z + ie*rnorm(n, 0, 1)

  # make group sizes random with 25 groups
  ngroups = 1
	gsize = trunc(length(xcoord)/ngroups)
	nxtra = length(xcoord) - gsize*ngroups
	grpindx = kronecker((1:(ngroups-nxtra)),rep(1, times = gsize))
	if(nxtra > 0) grpindx = c(grpindx,
		kronecker(((ngroups - nxtra + 1):ngroups), 
		rep(1, times = gsize + 1)))
	grpindx = grpindx[order(runif(length(grpindx)))]

	# create a data.frame of the data
	d2 = data.frame(x = xcoord, y = ycoord, 
		X1 = X[,2], X2 = X[,3], z = Z, grpindx = grpindx)

	# plots the simulated spatial data
	colorpoints(xcoord, ycoord, xyz_4y$z, nclass = 8)
	colorpoints(xcoord, ycoord, X[,3], nclass = 8)
	colorpoints(xcoord, ycoord, Z, nclass = 8)

#-----------------------------------------------------------------------
#                ESTIMATE COVARIANCE PARAMETERS
#-----------------------------------------------------------------------

	source(paste0(ScriptPath, 'cope_simple.R'))
	source(paste0(ScriptPath, 'm2LLg_simple.R'))
	#undebug(cope_simple)
	#undebug(m2LLg_simple)
	# Jay's covariance estimation function is called cope_simple()
	# use true covariance values as intial values
	cope_out = cope_simple(z ~ X1 + X2, data = d2, x_column = 'x', y_column = 'y', 
		subSampCol = 'grpindx', thetaini = c(100, 1, 1))
	# set spmodel initial values as those estimated by cope_simple()
	spcov_ini = spcov_initial('exponential', 
		range = cope_out['range'], ie = cope_out['ie'], 
		 de = cope_out['de']) 
	# fit regular spmodel
	spmodel_out_nolocal = splm(z ~ X1 + X2, data = d2, xcoord = 'x', ycoord = 'y',
		spcov_initial = spcov_initial, 
#		local = list(index = d2$grpindx),
		control = list(reltol = 1e-12), estmethod = 'reml')
	# fit the big data model using local option in spmodel
	spmodel_out = splm(z ~ X1 + X2, data = d2, xcoord = 'x', ycoord = 'y',
		spcov_initial = spcov_initial, 
		local = list(index = d2$grpindx),
		control = list(reltol = 1e-12), estmethod = 'reml')
	# compare covariance parameter estimation using Jay's code vs. spmodel
	summary(spmodel_out_nolocal)
	summary(spmodel_out)
	cope_out
	
###
###  THEY ALL MATCH, AT LEAST PRETTY CLOSELY
###

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#              REPEAT SIMULATION, SINGLE GROUP
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
	
	# set seed for reproducibility
	set.seed(1002)
	#set sample size
	n = 200
	# simulate random spatial x-coordinates
	xcoord = runif(n)
	# simulate random spatial x-coordinates
	ycoord = runif(n)
	# make a data.frame from coordinates
	xyDF = data.frame(x = xcoord, y = ycoord)
	
	# chose a range parameter and simulate pure autocorrelation with variance 1
  range = 1
  # specify covariance values
	spcov_params_val <- spcov_params("spherical", 
		de = 1, ie = .00001, range = range)
	# simulate the data for response variable
	z = sprnorm(spcov_params_val, data = xyDF, 
		xcoord = x, ycoord = y)
	xyz_4y = data.frame(xyDF, z = z)
	# simulate an autocorrelated covariate as well
	z = sprnorm(spcov_params_val, data = xyDF, 
		xcoord = x, ycoord = y)
	xyz_4x = data.frame(xyDF, z = z)

  # design matrix with X1 spatially independent and X1 spatially patterned
	X = as.matrix(cbind(rep(1, times = n), rnorm(n), xyz_4x$z))
  colnames(X) = c('int','X1','X2')

  # make response variable with a specified standard deviation for
  # autocorrelated error and another for independent error
  de = 10
  ie = 1
	Z =  1*X[,1]+ 1*X[,2] + 1*X[,3] + 
    de*xyz_4y$z + ie*rnorm(n, 0, 1)

  # make group sizes random with 25 groups
  ngroups = 1
	gsize = trunc(length(xcoord)/ngroups)
	nxtra = length(xcoord) - gsize*ngroups
	grpindx = kronecker((1:(ngroups-nxtra)),rep(1, times = gsize))
	if(nxtra > 0) grpindx = c(grpindx,
		kronecker(((ngroups - nxtra + 1):ngroups), 
		rep(1, times = gsize + 1)))
	grpindx = grpindx[order(runif(length(grpindx)))]

	# create a data.frame of the data
	d2 = data.frame(x = xcoord, y = ycoord, 
		X1 = X[,2], X2 = X[,3], z = Z, grpindx = grpindx)

	# plots the simulated spatial data
	colorpoints(xcoord, ycoord, xyz_4y$z, nclass = 8)
	colorpoints(xcoord, ycoord, X[,3], nclass = 8)
	colorpoints(xcoord, ycoord, Z, nclass = 8)

#-----------------------------------------------------------------------
#                ESTIMATE COVARIANCE PARAMETERS
#-----------------------------------------------------------------------

	source(paste0(ScriptPath, 'cope_simple.R'))
	source(paste0(ScriptPath, 'm2LLg_simple.R'))
	#undebug(cope_simple)
	#undebug(m2LLg_simple)
	# Jay's covariance estimation function is called cope_simple()
	# use true covariance values as intial values
	cope_out = cope_simple(z ~ X1 + X2, data = d2, x_column = 'x', y_column = 'y', 
		subSampCol = 'grpindx', thetaini = c(100, 1, 1))
	# set spmodel initial values as those estimated by cope_simple()
	spcov_initial = spcov_initial('exponential', 
		range = cope_out['range'], ie = cope_out['ie'], 
		 de = cope_out['de']) 
	# fit regular spmodel
	spmodel_out_nolocal = splm(z ~ X1 + X2, data = d2, xcoord = 'x', ycoord = 'y',
		spcov_initial = spcov_initial, 
#		local = list(index = d2$grpindx),
		control = list(reltol = 1e-12), estmethod = 'reml')
	# fit model using local option in spmodel but only 1 big data group
	spmodel_out = splm(z ~ X1 + X2, data = d2, xcoord = 'x', ycoord = 'y',
		spcov_initial = spcov_initial, 
		local = list(index = d2$grpindx),
		control = list(reltol = 1e-12), estmethod = 'reml')
	# compare covariance parameter estimation using Jay's code vs. spmodel
	summary(spmodel_out_nolocal)
	summary(spmodel_out)
	cope_out

###
###  THEY ALL MATCH VERY CLOSELY
###


  
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#           SIMULATION WITH N = 1000, 25 GROUPS
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
	
	# set seed for reproducibility
	set.seed(1003)
	#set sample size
	n = 1000
	# simulate random spatial x-coordinates
	xcoord = runif(n)
	# simulate random spatial x-coordinates
	ycoord = runif(n)
	# make a data.frame from coordinates
	xyDF = data.frame(x = xcoord, y = ycoord)
	
	# chose a range parameter and simulate pure autocorrelation with variance 1
  range = 1
  # specify covariance values
	spcov_params_val <- spcov_params("spherical", 
		de = 1, ie = .00001, range = range)
	# simulate the data for response variable
	z = sprnorm(spcov_params_val, data = xyDF, 
		xcoord = x, ycoord = y)
	xyz_4y = data.frame(xyDF, z = z)
	# simulate an autocorrelated covariate as well
	z = sprnorm(spcov_params_val, data = xyDF, 
		xcoord = x, ycoord = y)
	xyz_4x = data.frame(xyDF, z = z)

  # design matrix with X1 spatially independent and X1 spatially patterned
	X = as.matrix(cbind(rep(1, times = n), rnorm(n), xyz_4x$z))
  colnames(X) = c('int','X1','X2')

  # make response variable with a specified standard deviation for
  # autocorrelated error and another for independent error
  de = 10
  ie = 1
	Z =  1*X[,1]+ 1*X[,2] + 1*X[,3] + 
    de*xyz_4y$z + ie*rnorm(n, 0, 1)

  # make group sizes random with 25 groups
  ngroups = 25
	gsize = trunc(length(xcoord)/ngroups)
	nxtra = length(xcoord) - gsize*ngroups
	grpindx = kronecker((1:(ngroups-nxtra)),rep(1, times = gsize))
	if(nxtra > 0) grpindx = c(grpindx,
		kronecker(((ngroups - nxtra + 1):ngroups), 
		rep(1, times = gsize + 1)))
	grpindx = grpindx[order(runif(length(grpindx)))]

	# create a data.frame of the data
	d2 = data.frame(x = xcoord, y = ycoord, 
		X1 = X[,2], X2 = X[,3], z = Z, grpindx = grpindx)

	# plot the simulated spatial data
	colorpoints(xcoord, ycoord, Z, nclass = 8)

#-----------------------------------------------------------------------
#                ESTIMATE COVARIANCE PARAMETERS
#-----------------------------------------------------------------------

	source(paste0(ScriptPath, 'cope_simple.R'))
	source(paste0(ScriptPath, 'm2LLg_simple.R'))
	#undebug(cope_simple)
	#undebug(m2LLg_simple)
	# Jay's covariance estimation function is called cope_simple()
	# use true covariance values as intial values
	cope_out = cope_simple(z ~ X1 + X2, data = d2, x_column = 'x', y_column = 'y', 
		subSampCol = 'grpindx', thetaini = c(100, 1, 1))
	# set spmodel initial values as those estimated by cope_simple()
	spcov_initial = spcov_initial('exponential', 
		range = cope_out['range'], ie = cope_out['ie'], 
		 de = cope_out['de']) 
	# compare covariance parameter estimation using Jay's code vs. spmodel
	spmodel_out = splm(z ~ X1 + X2, data = d2, xcoord = 'x', ycoord = 'y',
		spcov_initial = spcov_initial, 
		local = list(index = d2$grpindx),
		control = list(reltol = 1e-12), estmethod = 'reml')
	summary(spmodel_out)
	cope_out

###
###  THEY ARE QUITE A BIT DIFFERENT
###

	# refit using Jay's code using spmodel output as initial values
	spmod_spcov = summary(test)$coefficients$spcov
	cope_out1 = cope_simple(z ~ X1 + X2, data = d2, x_column = 'x', y_column = 'y', 
		subSampCol = 'grpindx', 
		thetaini = c(spmod_spcov['de'], spmod_spcov['range'], spmod_spcov['ie'])
	)
	# any changes in output using Jay's code?
  cbind(cope_out, cope_out1)
  
###
###  NO CHANGE
###

	# refit using spmodel and true initial values
	spcov_initial = spcov_initial('exponential', 
		range = 1, ie = 5, de = 95) 
	spmodel_out1 = splm(z ~ X1 + X2, data = d2, xcoord = 'x', ycoord = 'y',
		spcov_initial = spcov_initial, 
		local = list(index = d2$grpindx),
		control = list(reltol = 1e-12), estmethod = 'reml')
	summary(spmodel_out1)
  cbind(summary(spmodel_out)$coefficients$spcov, 
		summary(spmodel_out1)$coefficients$spcov)
  
###
###  NO CHANGE
###
  
  
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#           REPEAT SIMULATION WITH N = 1000, 25 GROUPS
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
	
	# set seed for reproducibility
	set.seed(1007)
	#set sample size
	n = 1000
	# simulate random spatial x-coordinates
	xcoord = runif(n)
	# simulate random spatial x-coordinates
	ycoord = runif(n)
	# make a data.frame from coordinates
	xyDF = data.frame(x = xcoord, y = ycoord)
	
	# chose a range parameter and simulate pure autocorrelation with variance 1
  range = 1
  # specify covariance values
	spcov_params_val <- spcov_params("spherical", 
		de = 1, ie = .00001, range = range)
	# simulate the data for response variable
	z = sprnorm(spcov_params_val, data = xyDF, 
		xcoord = x, ycoord = y)
	xyz_4y = data.frame(xyDF, z = z)
	# simulate an autocorrelated covariate as well
	z = sprnorm(spcov_params_val, data = xyDF, 
		xcoord = x, ycoord = y)
	xyz_4x = data.frame(xyDF, z = z)

  # design matrix with X1 spatially independent and X1 spatially patterned
	X = as.matrix(cbind(rep(1, times = n), rnorm(n), xyz_4x$z))
  colnames(X) = c('int','X1','X2')

  # make response variable with a specified standard deviation for
  # autocorrelated error and another for independent error
  de = 10
  ie = 1
	Z =  1*X[,1]+ 1*X[,2] + 1*X[,3] + 
    de*xyz_4y$z + ie*rnorm(n, 0, 1)

  # make group sizes random with 25 groups
  ngroups = 25
	gsize = trunc(length(xcoord)/ngroups)
	nxtra = length(xcoord) - gsize*ngroups
	grpindx = kronecker((1:(ngroups-nxtra)),rep(1, times = gsize))
	if(nxtra > 0) grpindx = c(grpindx,
		kronecker(((ngroups - nxtra + 1):ngroups), 
		rep(1, times = gsize + 1)))
	grpindx = grpindx[order(runif(length(grpindx)))]

	# create a data.frame of the data
	d2 = data.frame(x = xcoord, y = ycoord, 
		X1 = X[,2], X2 = X[,3], z = Z, grpindx = grpindx)

	# plot the simulated spatial data
	colorpoints(xcoord, ycoord, Z, nclass = 8)

#-----------------------------------------------------------------------
#                ESTIMATE COVARIANCE PARAMETERS
#-----------------------------------------------------------------------

	source(paste0(ScriptPath, 'cope_simple.R'))
	source(paste0(ScriptPath, 'm2LLg_simple.R'))
	#undebug(cope_simple)
	#undebug(m2LLg_simple)
	# Jay's covariance estimation function is called cope_simple()
	# use true covariance values as intial values
	cope_out = cope_simple(z ~ X1 + X2, data = d2, x_column = 'x', y_column = 'y', 
		subSampCol = 'grpindx', thetaini = c(100, 1, 1))
	# set spmodel initial values as those estimated by cope_simple()
	spcov_initial = spcov_initial('exponential', 
		range = cope_out['range'], ie = cope_out['ie'], 
		 de = cope_out['de']) 
	# compare covariance parameter estimation using Jay's code vs. spmodel
	spmodel_out = splm(z ~ X1 + X2, data = d2, xcoord = 'x', ycoord = 'y',
		spcov_initial = spcov_initial, 
		local = list(index = d2$grpindx),
		control = list(reltol = 1e-12), estmethod = 'reml')
	summary(spmodel_out)
	cope_out

###
###  STILL QUITE A BIT DIFFERENT
###

	# refit using Jay's code using spmodel output as initial values
	spmod_spcov = summary(test)$coefficients$spcov
	cope_out1 = cope_simple(z ~ X1 + X2, data = d2, x_column = 'x', y_column = 'y', 
		subSampCol = 'grpindx', 
		thetaini = c(spmod_spcov['de'], spmod_spcov['range'], spmod_spcov['ie'])
	)
	# any changes in output using Jay's code?
  cbind(cope_out, cope_out1)
  
###
###  NO CHANGE
###

	# refit using spmodel and true initial values
	spcov_initial = spcov_initial('exponential', 
		range = 1, ie = 5, de = 95) 
	spmodel_out1 = splm(z ~ X1 + X2, data = d2, xcoord = 'x', ycoord = 'y',
		spcov_initial = spcov_initial, 
		local = list(index = d2$grpindx),
		control = list(reltol = 1e-12), estmethod = 'reml')
	summary(spmodel_out1)
  cbind(summary(spmodel_out)$coefficients$spcov, 
		summary(spmodel_out1)$coefficients$spcov)
  
###
###  NO CHANGE
###
  
  

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#           SINGLE SIMULATION WITH N = 5000, 100 GROUPS
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
	
	# set seed for reproducibility
	set.seed(1009)
	#set sample size
	n = 5000
	# simulate random spatial x-coordinates
	xcoord = runif(n)
	# simulate random spatial x-coordinates
	ycoord = runif(n)
	# make a data.frame from coordinates
	xyDF = data.frame(x = xcoord, y = ycoord)
	
	# chose a range parameter and simulate pure autocorrelation with variance 1
  range = 1
  # specify covariance values
	spcov_params_val <- spcov_params("spherical", 
		de = 1, ie = .00001, range = range)
	# simulate the data for response variable
	z = sprnorm(spcov_params_val, data = xyDF, 
		xcoord = x, ycoord = y)
	xyz_4y = data.frame(xyDF, z = z)
	# simulate an autocorrelated covariate as well
	z = sprnorm(spcov_params_val, data = xyDF, 
		xcoord = x, ycoord = y)
	xyz_4x = data.frame(xyDF, z = z)

  # design matrix with X1 spatially independent and X1 spatially patterned
	X = as.matrix(cbind(rep(1, times = n), rnorm(n), xyz_4x$z))
  colnames(X) = c('int','X1','X2')

  # make response variable with a specified standard deviation for
  # autocorrelated error and another for independent error
  de = 1
  ie = 1
	Z =  1*X[,1]+ 1*X[,2] + 1*X[,3] + 
    de*xyz_4y$z + ie*rnorm(n, 0, 1)

  # make group sizes random with 25 groups
  ngroups = 100
	gsize = trunc(length(xcoord)/ngroups)
	nxtra = length(xcoord) - gsize*ngroups
	grpindx = kronecker((1:(ngroups-nxtra)),rep(1, times = gsize))
	if(nxtra > 0) grpindx = c(grpindx,
		kronecker(((ngroups - nxtra + 1):ngroups), 
		rep(1, times = gsize + 1)))
	grpindx = grpindx[order(runif(length(grpindx)))]

	# create a data.frame of the data
	d2 = data.frame(x = xcoord, y = ycoord, 
		X1 = X[,2], X2 = X[,3], z = Z, grpindx = grpindx)

	# plots the simulated spatial data
	colorpoints(xcoord, ycoord, xyz_4y$z, nclass = 8)
	colorpoints(xcoord, ycoord, X[,3], nclass = 8)
	colorpoints(xcoord, ycoord, Z, nclass = 8)

#-----------------------------------------------------------------------
#                ESTIMATE COVARIANCE PARAMETERS
#-----------------------------------------------------------------------

	source(paste0(ScriptPath, 'cope_simple.R'))
	source(paste0(ScriptPath, 'm2LLg_simple.R'))
	#undebug(cope_simple)
	#undebug(m2LLg_simple)
	# Jay's covariance estimation function is called cope_simple()
	# use true covariance values as intial values
	cope_out = cope_simple(z ~ X1 + X2, data = d2, 
		x_column = 'x', y_column = 'y', 
		corModel = corModelSpherical,
		subSampCol = 'grpindx', thetaini = c(1, 1, .1))
	# set spmodel initial values as those estimated by cope_simple()
	spcov_initial = spcov_initial('spherical', 
		range = cope_out['range'], ie = cope_out['ie'], 
		 de = cope_out['de']) 
	# compare covariance parameter estimation using Jay's code vs. spmodel
	spmodel_out = splm(z ~ X1 + X2, data = d2, xcoord = 'x', ycoord = 'y',
		spcov_initial = spcov_initial, 
		local = list(index = d2$grpindx),
		control = list(reltol = 1e-12), estmethod = 'reml')
	summary(spmodel_out)
	cope_out

###
###  THEY ARE QUITE A BIT DIFFERENT
###
  
  
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#         SIMULATE IN A LOOP, CHECK FOR BIAS AND MSE
#                BASED ON EXPONENTIAL COVARIANCE
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
	
	# set seed for reproducibility
	set.seed(1011)

	# number of simulations
	nsim = 100

	# empy matrix to store cope_simple results
	store_cope = matrix(0, nrow = 3, ncol = nsim)
	# empy matrix to store spmodel results
	store_spmodel = matrix(0, nrow = 3, ncol = nsim)
for (j in 1:nsim) {
	#set sample size
	n = 1000
	# simulate random spatial x-coordinates
	xcoord = runif(n)
	# simulate random spatial x-coordinates
	ycoord = runif(n)
	# make a data.frame from coordinates
	xyDF = data.frame(x = xcoord, y = ycoord)
	
	# chose a range parameter and simulate pure autocorrelation with variance 1
  range = 1
  # specify covariance values
	spcov_params_val <- spcov_params("exponential", 
		de = 1, ie = .00001, range = range)
	# simulate the data for response variable
	z = sprnorm(spcov_params_val, data = xyDF, 
		xcoord = x, ycoord = y)
	xyz_4y = data.frame(xyDF, z = z)
	# simulate an autocorrelated covariate as well
	spcov_params_val <- spcov_params("exponential", 
		de = 1, ie = .5, range = 0.25*range)
	z = sprnorm(spcov_params_val, data = xyDF, 
		xcoord = x, ycoord = y)
	xyz_4x = data.frame(xyDF, z = z)

  # design matrix with X1 spatially independent and X1 spatially patterned
	X = as.matrix(cbind(rep(1, times = n), rnorm(n), xyz_4x$z))
  colnames(X) = c('int','X1','X2')

  # make response variable with a specified standard deviation for
  # autocorrelated error and another for independent error
  de = 1
  ie = .1
	Z =  1*X[,1]+ 1*X[,2] + 1*X[,3] + 
    de*xyz_4y$z + ie*rnorm(n, 0, 1)

  # make group sizes random with 20 groups
  ngroups = 20
	gsize = trunc(length(xcoord)/ngroups)
	nxtra = length(xcoord) - gsize*ngroups
	grpindx = kronecker((1:(ngroups-nxtra)),rep(1, times = gsize))
	if(nxtra > 0) grpindx = c(grpindx,
		kronecker(((ngroups - nxtra + 1):ngroups), 
		rep(1, times = gsize + 1)))
	grpindx = grpindx[order(runif(length(grpindx)))]

	# create a data.frame of the data
	d2 = data.frame(x = xcoord, y = ycoord, 
		X1 = X[,2], X2 = X[,3], z = Z, grpindx = grpindx)

	# plot the simulated spatial data
	colorpoints(xcoord, ycoord, Z, nclass = 8)

#-----------------------------------------------------------------------
#                ESTIMATE COVARIANCE PARAMETERS
#-----------------------------------------------------------------------

	source(paste0(ScriptPath, 'cope_simple.R'))
	source(paste0(ScriptPath, 'm2LLg_simple.R'))
	#undebug(cope_simple)
	#undebug(m2LLg_simple)
	# Jay's covariance estimation function is called cope_simple()
	# use true covariance values as intial values
	cope_out = cope_simple(z ~ X1 + X2, data = d2, x_column = 'x', y_column = 'y',
		corModel = corModelExponential,
		subSampCol = 'grpindx', thetaini = c(1, 1, 0.1))
	# use true covariance values as intial values
	spcov_ini = spcov_initial('exponential', 
		range = 1, ie = 0.1, de = 1) 
	# compare covariance parameter estimation using Jay's code vs. spmodel
	spmodel_out = splm(z ~ X1 + X2, data = d2, xcoord = 'x', ycoord = 'y',
		spcov_initial = spcov_ini, 
		local = list(index = d2$grpindx),
		control = list(reltol = 1e-12), estmethod = 'reml')
	spmod_spcov = summary(spmodel_out)$coefficients$spcov[1:3]
	store_cope[,j] = cope_out
	store_spmodel[,j] = spmod_spcov
}


# histograms of range parameter estimates
hist(store_cope[3,])
hist(store_spmodel[3,])
median(store_cope[3,])
median(store_spmodel[3,])

# histograms of partial sill parameter estimates
hist(store_cope[1,])
hist(store_spmodel[1,])
median(store_cope[1,])
median(store_spmodel[1,])

####
#### when simulating and estimating exponential covariance
#### Jay's code gets closer to true values (based on median)
#### but Jay's code also has more times with large partial sills and ranges
####


#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#         SIMULATE IN A LOOP, CHECK FOR BIAS AND MSE
#                BASED ON SPHERICAL COVARIANCE
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
	
	# set seed for reproducibility
	set.seed(1011)

	# number of simulations
	nsim = 100

	# empy matrix to store cope_simple results
	store_cope = matrix(0, nrow = 3, ncol = nsim)
	# empy matrix to store spmodel results
	store_spmodel = matrix(0, nrow = 3, ncol = nsim)
for (j in 1:nsim) {
	#set sample size
	n = 1000
	# simulate random spatial x-coordinates
	xcoord = runif(n)
	# simulate random spatial x-coordinates
	ycoord = runif(n)
	# make a data.frame from coordinates
	xyDF = data.frame(x = xcoord, y = ycoord)
	
	# chose a range parameter and simulate pure autocorrelation with variance 1
  range = 1
  # specify covariance values
	spcov_params_val <- spcov_params("spherical", 
		de = 1, ie = .00001, range = range)
	# simulate the data for response variable
	z = sprnorm(spcov_params_val, data = xyDF, 
		xcoord = x, ycoord = y)
	xyz_4y = data.frame(xyDF, z = z)
	# simulate an autocorrelated covariate as well
	spcov_params_val <- spcov_params("spherical", 
		de = .5, ie = .5, range = 0.25*range)
	z = sprnorm(spcov_params_val, data = xyDF, 
		xcoord = x, ycoord = y)
	xyz_4x = data.frame(xyDF, z = z)

  # design matrix with X1 spatially independent and X1 spatially patterned
	X = as.matrix(cbind(rep(1, times = n), rnorm(n), xyz_4x$z))
  colnames(X) = c('int','X1','X2')

  # make response variable with a specified standard deviation for
  # autocorrelated error and another for independent error
  de = 1
  ie = .1
	Z =  1*X[,1]+ 1*X[,2] + 1*X[,3] + 
    de*xyz_4y$z + ie*rnorm(n, 0, 1)

  # make group sizes random with 20 groups
  ngroups = 20
	gsize = trunc(length(xcoord)/ngroups)
	nxtra = length(xcoord) - gsize*ngroups
	grpindx = kronecker((1:(ngroups-nxtra)),rep(1, times = gsize))
	if(nxtra > 0) grpindx = c(grpindx,
		kronecker(((ngroups - nxtra + 1):ngroups), 
		rep(1, times = gsize + 1)))
	grpindx = grpindx[order(runif(length(grpindx)))]

	# create a data.frame of the data
	d2 = data.frame(x = xcoord, y = ycoord, 
		X1 = X[,2], X2 = X[,3], z = Z, grpindx = grpindx)

	# plot the simulated spatial data
	colorpoints(xcoord, ycoord, Z, nclass = 8)

#-----------------------------------------------------------------------
#                ESTIMATE COVARIANCE PARAMETERS
#-----------------------------------------------------------------------

	source(paste0(ScriptPath, 'cope_simple.R'))
	source(paste0(ScriptPath, 'm2LLg_simple.R'))
	#undebug(cope_simple)
	#undebug(m2LLg_simple)
	# Jay's covariance estimation function is called cope_simple()
	# use true covariance values as intial values
	cope_out = cope_simple(z ~ X1 + X2, data = d2, x_column = 'x', y_column = 'y',
		corModel = corModelSpherical,
		subSampCol = 'grpindx', thetaini = c(1, 1, 0.1))
	# use true covariance values as intial values
	spcov_ini = spcov_initial('spherical', 
		range = 1, ie = 0.1, de = 1) 
	# compare covariance parameter estimation using Jay's code vs. spmodel
	spmodel_out = splm(z ~ X1 + X2, data = d2, xcoord = 'x', ycoord = 'y',
		spcov_initial = spcov_ini, 
		local = list(index = d2$grpindx),
		control = list(reltol = 1e-12), estmethod = 'reml')
	spmod_spcov = summary(spmodel_out)$coefficients$spcov[1:3]
	store_cope[,j] = cope_out
	store_spmodel[,j] = spmod_spcov
}


# histograms of range parameter estimates
hist(store_cope[3,])
hist(store_spmodel[3,])
median(store_cope[3,])
median(store_spmodel[3,])

# histograms of partial sill parameter estimates
hist(store_cope[1,])
hist(store_spmodel[1,])
median(store_cope[1,])
median(store_spmodel[1,])

####
#### when simulating and estimating spherical covariance
#### Jay's code gets a little bit closer to true values (based on median)
#### but Jay's code also has more times with large partial sills and ranges
#### and spmodel has nice symmetrical distribution of estimates with a 
#### somewhat downward bias, but not as much as for exponential
####
