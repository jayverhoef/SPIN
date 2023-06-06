library('spmodel')

	set.seed(1001)
  n = 1000

	x = runif(n)
	y = runif(n)
	xyDF = data.frame(x = x, y = y)
	
	spcov_params_val <- spcov_params("spherical", 
		de = 10, ie = .0001, range = .9)
	z = sprnorm(spcov_params_val, data = xyDF, xcoord = x, ycoord = y)
	xyz_4y = data.frame(xyDF, z = z)
	z = sprnorm(spcov_params_val, data = xyDF, xcoord = x, ycoord = y)
	xyz_4x = data.frame(xyDF, z = z)
		
  # design matrix with X1 spatially independent and X2 spatially patterned
	X1 = rnorm(n)
	X2 = xyz_4x$z

  # make response variable with some proportion of nugget
  prop = .95
	Z =  1 + 1*X1 + 1*X2 + 10*prop*xyz_4y$z + 10*(1 - prop)*rnorm(n, 0, 1)


  # make group sizes random with 25 groups
  ngroups = 25
	gsize = trunc(length(x)/ngroups)
	nxtra = length(x) - gsize*ngroups
	grpindx = kronecker((1:(ngroups-nxtra)),rep(1, times = gsize))
	if(nxtra > 0) grpindx = c(grpindx,
		kronecker(((ngroups - nxtra + 1):ngroups), 
		rep(1, times = gsize + 1)))
	grpindx = grpindx[order(runif(length(grpindx)))]

	#make a data.frame of the data
	d2 = data.frame(x = x, y = y, X1 = X1, X2 = X2, z = Z, 
    grpindx.rand = grpindx)

	spcovi = spcov_initial('exponential', range = .30, ie = 30, 
		de = 30, known = c('range','ie','de'))
	test = splm(z ~ X1 + X2, data = d2, xcoord = 'x', ycoord = 'y',
		spcov_initial = spcovi, 
		local = list(index = d2$grpindx.rand),
		control = list(reltol = 1e-8), estmethod = 'reml')

	spcovi = spcov_initial('exponential', range = .30, ie = 30, 
		de = 30, known = c('given'))
	test = splm(z ~ X1 + X2, data = d2, xcoord = 'x', ycoord = 'y',
		spcov_initial = spcovi, 
		local = list(index = d2$grpindx.rand),
		control = list(reltol = 1e-8), estmethod = 'reml')
