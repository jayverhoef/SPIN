library(gstat)
library(RColorBrewer)
library(viridis)
library(classInt)
library(sp)
library(sf)

#-------------------------------------------------------------------------------
#                    Load the Data
#-------------------------------------------------------------------------------
# load the data from disk
basepath = paste0('/mnt/ExtraDrive1/Work/desktop_data/2019_papers',
	'/fastREML/spPart_package/spPart/inst/rawdata/')
load(file = paste0(basepath,'NDVIcali.rda'))
load(file = paste0(basepath,'cali.rda'))

X11()
# Look at the response variable (NDVI) and covariates spatially
layout(matrix(1:4, ncol = 2))
# NDVI map
plot(NDVIcali[,c('easting','northing')], pch = 19, 
	col = findColours(classIntervals(NDVIcali$NDVI, n = 20, style = 'fisher'), 
	viridis(20)), cex = .4, main = 'NDVI')
plot(cali, add = TRUE, lwd = 3)
# elev map
plot(NDVIcali[,c('easting','northing')], pch = 19, 
	col = findColours(classIntervals(NDVIcali$elev, n = 20, style = 'fisher'), 
	viridis(20)), cex = .4, main = 'Elevation')
plot(cali, add = TRUE, lwd = 3)
# ecoregions map
pal = brewer.pal(n = 19, name = "Set1")
plot(NDVIcali[,c('easting','northing')], pch = 19, 
	col = findColours(classIntervals(as.integer(NDVIcali$sectionID), 
	n = 19, style = 'fisher'), pal), cex = .4, main = 'Ecoregion')
plot(cali, add = TRUE, lwd = 3)
# precipitation map
plot(NDVIcali[,c('easting','northing')], pch = 19, 
	col = findColours(classIntervals(as.integer(NDVIcali$precip), 
	n = 15, style = 'fisher'), viridis(15)), cex = .4, main = 'Precipitation')
plot(cali, add = TRUE, lwd = 3)
layout(1)

# create a working dataset
d1 = NDVIcali
# get rid of missing values for now
d1 = d1[!is.na(d1$NDVI),]
# rescale some of the variables to make model-fitting more stable
d1$NDVI = d1$NDVI/1000
d1$northing = d1$northing/100000
d1$easting = d1$easting/100000
d1$elev = log(d1$elev/1000 + .3)

#-------------------------------------------------------------------------------
#                    Exploratory Data Analysis
#-------------------------------------------------------------------------------

# look at histograms
layout(matrix(1:4, nrow = 2))
hist(d1$NDVI)
hist(d1$northing)
hist(d1$easting)
hist(d1$elev)
layout(1)

# Response according to covariates
layout(matrix(c(1,1,2,3), nrow = 2, byrow = TRUE))
# scatterplot of response variable to elevation
plot(d1$elev, d1$NDVI, pch = 19, cex = .5)
# boxplot of response variable by ecological section
boxplot(NDVI ~ sectionID, data = d1)
# boxplot of response variable by precipitation categories
boxplot(NDVI ~ precip, data = d1)
layout(1)

# This is the model that we envision, all as fixed effects with no spatial
# autocorrelation
fix_eff_model = lm(NDVI ~ sectionID*elev + precip, data = d1)
summary(fix_eff_model)
d1$resid = resid(fix_eff_model)

# look at the directional semivariogram of the residuals
# it shows most autocorrelation in N/S direction, least in E/W direction
vgm_resid_dir <- variogram(resid~1, loc=~easting+northing, data=d1,
	alpha=c(0,45,90,135), cutoff = 3, width = .3)
       
  old_par = par(mar = c(5,5,3,1))
  plot(vgm_resid_dir$dist, vgm_resid_dir$gamma, type = 'n', xlab = 'Distance (km)',
    ylab = 'Semivariogram', cex.lab = 2, cex.axis = 1.5)
  points(vgm_resid_dir[vgm_resid_dir$dir.hor == 0,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 0,'gamma'], pch = 1, 
    cex = 10*vgm_resid_dir$np/max(vgm_resid_dir$np))
  lines(vgm_resid_dir[vgm_resid_dir$dir.hor == 0,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 0,'gamma'], lwd = 2)
  points(vgm_resid_dir[vgm_resid_dir$dir.hor == 45,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 45,'gamma'], pch = 2, 
    cex = 10*vgm_resid_dir$np/max(vgm_resid_dir$np))
  lines(vgm_resid_dir[vgm_resid_dir$dir.hor == 45,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 45,'gamma'], lwd = 2)
  points(vgm_resid_dir[vgm_resid_dir$dir.hor == 90,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 90,'gamma'], pch = 5, 
    cex = 10*vgm_resid_dir$np/max(vgm_resid_dir$np))
  lines(vgm_resid_dir[vgm_resid_dir$dir.hor == 90,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 90,'gamma'], lwd = 2)
  points(vgm_resid_dir[vgm_resid_dir$dir.hor == 135,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 135,'gamma'], pch = 22, 
    cex = 10*vgm_resid_dir$np/max(vgm_resid_dir$np))
  lines(vgm_resid_dir[vgm_resid_dir$dir.hor == 135,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 135,'gamma'], lwd = 2)
  legend(1, 1.3, legend=(c(0,45,90,135)), pch = c(1,2,5,22), cex = 3)
  par(old_par)

# map of the residuals
d2 = d1
d2$northing = d1$northing*100000
d2$easting = d1$easting*100000
coordinates(d2) = ~ easting + northing
# set the CRS to CONUS Albers
d2 = st_set_crs(st_as_sf(d2), "EPSG:5070")
# convert it back to a SpatialPointsDataFrame
d2 = as(d2, 'Spatial')
# residuals map
plot(d2, pch = 19, 
	col = findColours(classIntervals(d2$resid, n = 20, style = 'fisher'), 
	viridis(20)), cex = .4)
plot(cali, add = TRUE, lwd = 3)

#-------------------------------------------------------------------------------
#                    Fitting Big Data Models
#-------------------------------------------------------------------------------

library(spmodel)

  ngroups = 300
	gsize = trunc(dim(d1)[1]/ngroups)
	nxtra = dim(d1)[1] - gsize*ngroups
	grpindx = kronecker((1:(ngroups-nxtra)),rep(1, times = gsize))
	if(nxtra > 0) grpindx = c(grpindx,
		kronecker(((ngroups - nxtra + 1):ngroups), 
		rep(1, times = gsize + 1)))
	grpindx = grpindx[order(runif(length(grpindx)))]

spcov_ini = spcov_initial('spherical', 
		range = 100, ie = .01, de = 0.1) 
spmodel_out = splm(NDVI ~ sectionID*elev + precip, 
	data = d1, xcoord = 'easting', ycoord = 'northing',
		spcov_initial = spcov_ini, 
		local = TRUE,
		control = list(reltol = 1e-6), estmethod = 'reml')
summary(spmodel_out)

# In this model, I want a model with anisotropy, but it doesn't estimate it?
spcov_ini = spcov_initial('spherical', 
		range = 1, ie = .01, de = 0.1) 
spmodel_out1 = splm(NDVI ~ sectionID*elev + precip, 
	data = d1, xcoord = 'easting', ycoord = 'northing',
	anistropy = TRUE,
	spcov_initial = spcov_ini, 
	local = list(size = 25),
	control = list(reltol = 1e-6), estmethod = 'reml')
summary(spmodel_out1)

# In this model, I give the initial rotate value as 2 (> 0 but < 3.1414), 
# and I get an error, the anisotropy argument is not used.
# Error in optim(par = c(ie_prop_logodds = -2.30258509299405, range_log =
# 4.60517018598809,  : non-finite value supplied by optim
# In addition: Warning message:
# In log(rotate_prop/(1 - rotate_prop)) : NaNs produced
spcov_ini = spcov_initial('spherical', 
		range = 1, ie = .01, de = 0.1, rotate = 2, scale = .5) 
spmodel_out2 = splm(NDVI ~ sectionID*elev + precip, 
	data = d1, xcoord = 'easting', ycoord = 'northing',
	spcov_initial = spcov_ini, 
	local = list(size = 25),
	control = list(reltol = 1e-6), estmethod = 'reml')
summary(spmodel_out2)

# In this model, I give the initial rotate value as 2 (> 0 but < 3.1414)
# I get the same error when I used anisotropy = TRUE 
spcov_ini = spcov_initial('spherical', 
		range = 1, ie = .01, de = 0.1, rotate = 2, scale = .5) 
spmodel_out3 = splm(NDVI ~ sectionID*elev + precip, 
	data = d1, xcoord = 'easting', ycoord = 'northing',
	anisotropy = TRUE,
	spcov_initial = spcov_ini, 
	local = list(size = 25),
	control = list(reltol = 1e-6), estmethod = 'reml')
summary(spmodel_out3)

# In this model, it fits anisotropy as long as I give starting values
spcov_ini = spcov_initial('spherical', 
		range = 1, ie = .01, de = 0.1, rotate = .5, scale = .5) 
spmodel_out4 = splm(NDVI ~ sectionID*elev + precip, 
	data = d1, xcoord = 'easting', ycoord = 'northing',
	anisotropy = TRUE,
	spcov_initial = spcov_ini, 
	local = list(size = 25),
	control = list(reltol = 1e-6), estmethod = 'reml')
summary(spmodel_out4)
fefit4 = fitted(spmodel_out4)
X11()
plot(fefit4, d1$NDVI)
plot(d1$easting, d1$northing)
points(d1[fefit4 > 1.47 & fefit4 < 1.84 & d1$NDVI > .73 & d1$NDVI < 2.88,
	c('easting','northing')], pch = 19, cex = .5, col = 'red')
plot(d1$easting, d1$northing)
points(d1[fefit4 > 0.5 & fefit4 < 4 & d1$NDVI > -1 & d1$NDVI < .5,
	c('easting','northing')], pch = 19, cex = .5, col = 'red')

# In this model, it fits anisotropy even though I don't tell it to do so,
# because I specified starting values
spcov_ini = spcov_initial('spherical', 
		range = 1, ie = .01, de = 0.1, rotate = .5, scale = .5) 
spmodel_out5 = splm(NDVI ~ sectionID*elev + precip, 
	data = d1, xcoord = 'easting', ycoord = 'northing',
	anistropy = FALSE,
	spcov_initial = spcov_ini, 
	local = list(size = 25),
	control = list(reltol = 1e-6), estmethod = 'reml')
summary(spmodel_out5)

# Try random effects: Looks OK
spcov_ini = spcov_initial('spherical', 
		range = 1, ie = .02, de = 0.1) 
spmodel_out6 = splm(NDVI ~ elev + precip, 
	data = d1, xcoord = 'easting', ycoord = 'northing',
	spcov_initial = spcov_ini, 
	random = ~ sectionID,
	local = list(size = 25),
	control = list(reltol = 1e-6), estmethod = 'reml')
summary(spmodel_out6)
fitted(spmodel_out6, type = 'randcov')

# Try random slopes by sectionID: Looks OK
spcov_ini = spcov_initial('spherical', 
		range = 1, ie = .02, de = 0.1) 
spmodel_out7 = splm(NDVI ~ elev + precip, 
	data = d1, xcoord = 'easting', ycoord = 'northing',
	spcov_initial = spcov_ini, 
	random = ~ sectionID + sectionID:elev,
	local = list(size = 25),
	control = list(reltol = 1e-6), estmethod = 'reml')
summary(spmodel_out7)
fitted(spmodel_out7, type = 'randcov')

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
	
plot(d1$easting, d1$northing)
points(d1[d1$NDVI < 1 & fefit8 > 3.5,c('easting','northing')], col = 'red')

# Try anisotropy and random slopes by sectionID
# This is the final model that I want
# Try it with original unscaled data, no starting values
spmodel_out9 = splm(NDVI ~ elev + precip, 
	data = NDVIcali, xcoord = 'easting', ycoord = 'northing',
	anisotropy = TRUE,
	random = ~ sectionID + sectionID:elev,
	local = list(size = 100, var_adjust = 'theoretical'),
	control = list(reltol = 1e-6))
summary(spmodel_out9)
refit = fitted(spmodel_out9, type = 'randcov')
fefit = fitted(spmodel_out9)

Z1mat = model.matrix( ~ -1 + sectionID, data = NDVIcali)[!is.na(NDVIcali$NDVI),]
Z2mat = model.matrix( ~ -1 + sectionID:elev, data = NDVIcali)[!is.na(NDVIcali$NDVI),]
plot(fefit + Z1mat %*% refit[[1]] + Z2mat %*% refit[[2]],
	NDVIcali[!is.na(NDVIcali$NDVI),'NDVI'], pch = 19)
