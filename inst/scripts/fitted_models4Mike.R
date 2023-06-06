library(gstat)
library(RColorBrewer)
library(viridis)
library(classInt)
library(sp)
library(sf)
library(spmodel)

#-------------------------------------------------------------------------------
#                    Load the Data
#-------------------------------------------------------------------------------
# load the data from disk
basepath = paste0('/mnt/ExtraDrive1/Work/desktop_data/2019_papers',
	'/fastREML/spPart_package/spPart/inst/rawdata/')
load(file = paste0(basepath,'NDVIcali.rda'))
load(file = paste0(basepath,'cali.rda'))


# create a working dataset
d1 = NDVIcali
# get rid of missing values for now
d1 = d1[!is.na(d1$NDVI),]
# rescale some of the variables to make model-fitting more stable
d1$NDVI = d1$NDVI/1000
d1$northing = d1$northing/100000
d1$easting = d1$easting/100000
d1$elev = log(d1$elev/1000 + .3)


# Try anisotropy and random slopes by sectionID
# This is the final model that I want
spcov_ini = spcov_initial('spherical', 
		range = 1, ie = .02, de = 0.1, rotate = 1.5, scale = .5) 
spmodel_out8 = splm(NDVI ~ elev + precip, 
	data = d1, xcoord = 'easting', ycoord = 'northing',
	spcov_initial = spcov_ini, 
	random = ~ sectionID + sectionID:elev,
	local = list(size = 100, var_adjust = 'theoretical'),
	control = list(reltol = 1e-6))
summary(spmodel_out8)
# The variability among the fitted random effects does not look like the
# estimated variance for that effect?
fitted(spmodel_out8, type = 'randcov')
fefit8 = fitted(spmodel_out8)
refit8 = fitted(spmodel_out8, type = 'randcov')
Z1mat = model.matrix( ~ -1 + sectionID, data = d1)
Z2mat = model.matrix( ~ -1 + sectionID:elev, data = d1)
# This one looks OK
plot(fefit8, d1[,'NDVI'], pch = 19)
# but adding the fitted random effects does not look right
plot(fefit8 + Z1mat %*% refit8[[1]] + Z2mat %*% refit8[[2]],
	d1[,'NDVI'], pch = 19)
	
spcov_ini = spcov_initial('spherical', 
		range = 1, ie = .02, de = 0.1, rotate = 1.5, scale = .5) 
spmodel_out9 = splm(NDVI ~ elev + precip, 
	data = d1, xcoord = 'easting', ycoord = 'northing',
	spcov_initial = spcov_ini, 
	random = ~ sectionID + sectionID:elev,
	local = list(size = 100, var_adjust = 'theoretical', method = 'kmeans'),
	control = list(reltol = 1e-6))
summary(spmodel_out9)
# The variability among the fitted random effects does not look like the
# estimated variance for that effect?
fitted(spmodel_out9, type = 'randcov')
fefit9 = fitted(spmodel_out8)
refit8 = fitted(spmodel_out8, type = 'randcov')
Z1mat = model.matrix( ~ -1 + sectionID, data = d1)
Z2mat = model.matrix( ~ -1 + sectionID:elev, data = d1)
# This one looks OK
plot(fefit8, d1[,'NDVI'], pch = 19)
# but adding the fitted random effects does not look right
plot(fefit8 + Z1mat %*% refit8[[1]] + Z2mat %*% refit8[[2]],
	d1[,'NDVI'], pch = 19)


nres9 = residuals(spmodel_out9, type = "normalized")
.04487*Z1mat

d3 = d1[d1$sectionID %in% levels(d1$sectionID)[8:10],]
d3$sectionID = as.factor(as.character(d3$sectionID))

spmodel_out9 = splm(NDVI ~ elev + precip, 
	data = d3, xcoord = 'easting', ycoord = 'northing',
	spcov_type = 'none', 
	random = ~ sectionID + sectionID:elev,
	local = list(size = 100, method = 'kmeans'),
	control = list(reltol = 1e-6))
#summary(spmodel_out9)
fitted(spmodel_out9, type = 'randcov')


n = 800
range = 1.5
set.seed(1001)
xcoord = runif(n)
ycoord = runif(n)
x1 = rnorm(n)
f1 = as.factor(sample(rep(1:5, times = n/5), n))
f2 = as.factor(sample(rep(1:10, times = n/10), n))
spcov = 3*exp(-as.matrix(dist(cbind(xcoord, ycoord)))/range)
X1 = model.matrix( ~ -1 + f1, data.frame(f1))
Z1 = model.matrix( ~ -1 + f2, data.frame(f2))
Z2 = model.matrix( ~ -1 + f2:x1, data.frame(f1 = f2, x1 = x1))
err = t(chol(spcov)) %*% rnorm(n)
re1 = rnorm(10)
re2 = rnorm(10)
y = 3*x1 + X1 %*% 1:5 + Z1%*%re1 + Z2%*%re2 + err

d9 = data.frame(xcoord, ycoord, x1, f1, f2, y)
spmodel_out9 = splm(y ~ x1 + f1, 
	data = d9, xcoord = 'xcoord', ycoord = 'ycoord',
	spcov_type = 'exponential', 
	random = ~ f2 + f2:x1,
#	local = list(size = 100, method = 'kmeans'),
	control = list(reltol = 1e-6))
summary(spmodel_out9)
fitted(spmodel_out9, type = 'randcov')
plot(y, fitted(spmodel_out9))
plot(re1, fitted(spmodel_out9, type = 'randcov')[[1]])
plot(re2, fitted(spmodel_out9, type = 'randcov')[[2]])
plot(y, fitted(spmodel_out9) + 
	Z1 %*%fitted(spmodel_out9, type = 'randcov')[[1]] +
	Z2 %*%fitted(spmodel_out9, type = 'randcov')[[2]] )  
