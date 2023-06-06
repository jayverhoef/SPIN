library(raster)
library(sp)
library(devtools)
library(sf)
library(classInt)
library(viridis)
library(RColorBrewer)
library(nabor)
# only need to install from Github once
#install_github("bmcnellis/ClelandEcoregions")
library(ClelandEcoregions)


#-------------------------------------------------------------------------------
#              California polygon
#-------------------------------------------------------------------------------

filename = '/home/jayverhoef/Desktop/california.shp'
cali = shapefile(filename)

#-------------------------------------------------------------------------------
#              create NDVI data.frame within California
#-------------------------------------------------------------------------------

filename = '/home/jayverhoef/Desktop/NDVIcali.grd'
NDVIcali = raster(filename)
# get the x- y- coordinates from raster and create a data.frame
NDVIcali = as.data.frame(NDVIcali, xy = TRUE)
# turn it into an SpatialPointsDataFrame
coordinates(NDVIcali) <- ~ x + y
# convert to sf to get Albers CONUS projection
NDVIcali = st_set_crs(st_as_sf(NDVIcali), "EPSG:5070")
# convert back to SpatialPointsDataFrame
NDVIcali = as(NDVIcali, 'Spatial')
# get all of the points within the California boundary
NDVIcali = NDVIcali[!is.na(over(NDVIcali, cali)$NAME),]
# check it
plot(NDVIcali, pch = 19, cex = .3)
# turn it into a data.frame
NDVIcali = as.data.frame(NDVIcali)
# change the names
names(NDVIcali) = c('NDVI', 'easting', 'northing')
# check it
str(NDVIcali)
# plot NDVI spatially
plot(NDVIcali[,c('easting','northing')], pch = 19, 
	col = findColours(classIntervals(NDVIcali$NDVI, n = 20, style = 'fisher'), 
	viridis(20)), cex = .4)

#-------------------------------------------------------------------------------
#           create DEM elevation data and add to data.frame
#-------------------------------------------------------------------------------

filename = '/home/jayverhoef/Desktop/testCONUS_DEMraster.grd'
# a raster grid of DEM for CONUS
dem = raster(filename)

# crop dem to California and then aggregate on 20 x 20 grid to reduce size
dem_cali = crop(dem, cali)
dem_cali = aggregate(dem_cali, fact = 20, na.rm = TRUE)

# turn it into a data.frame
dem_caliSPDF = as.data.frame(dem_cali, xy = TRUE)
coordinates(dem_caliSPDF) <- ~ x + y
# convert to sf to get Albers CONUS projection
dem_caliSPDF = st_set_crs(st_as_sf(dem_caliSPDF), "EPSG:5070")
# convert back to SpatialPointsDataFrame
dem_caliSPDF = as(dem_caliSPDF, 'Spatial')
# get all of the points within the California boundary
dem_caliSPDF = dem_caliSPDF[!is.na(over(dem_caliSPDF, cali)$NAME),]
# plot the locations
plot(dem_caliSPDF, pch = 19, cex = .1)
# find the elevation value closest to NDVI points
nn = knn(data = as.data.frame(dem_caliSPDF)[,c('coords.x1','coords.x2')], 
	query = NDVIcali[,c('easting','northing')], k = 1) 
NDVIcali$elev = as.data.frame(
	dem_caliSPDF[as.vector(nn$nn.idx),'testCONUS_DEMraster'])$ testCONUS_DEMraster
# check it
str(NDVIcali)
# plot elev spatially
plot(NDVIcali[,c('easting','northing')], pch = 19, 
	col = findColours(classIntervals(NDVIcali$elev, n = 20, style = 'fisher'), 
	viridis(20)), cex = .4)
# plot locations of missing elevation data
points(NDVIcali[which(is.na(NDVIcali$elev)),c('easting','northing')], pch = 19,
	cex = 1, col = 'red')
# remove location from data set
NDVIcali = NDVIcali[!is.na(NDVIcali$elev),]

#-------------------------------------------------------------------------------
#            create Ecoregion data and add to data.frame
#-------------------------------------------------------------------------------

# get the intermediate resolution ecoregions from the list
ecor = Cleland2007_eco_map[[2]]
# get the single attribute for ecoregion as a factor
DF = as.data.frame(ecor@data[,'SECTION_ID'])
# now the @data has only the ecoregion attribute
ecor@data = DF

# the ecoregion polygons use lat/lon coordinates
ecor = st_set_crs(st_as_sf(ecor), "EPSG:4326")
# change it to a SpatialPolygonsDataFrame
ecor = as(ecor, 'Spatial')
# transform projection to CONUS Albers
ecor = spTransform(ecor, CRS("+init=epsg:5070"))
# plot it
plot(ecor)
# check that it is the same projection as California
plot(cali, add = TRUE, col = 'red')
# get the Ecoregion for each NDVI location
NDVIcaliSPDF = NDVIcali
coordinates(NDVIcaliSPDF) = ~ easting + northing
# set the CRS to CONUS Albers
NDVIcaliSPDF = st_set_crs(st_as_sf(NDVIcaliSPDF), "EPSG:5070")
# convert it back to a SpatialPointsDataFrame
NDVIcaliSPDF = as(NDVIcaliSPDF, 'Spatial')
# now get the ecoregions
ecorFACT <- data.frame(xx=over(NDVIcaliSPDF, ecor))
# add the ecoregions to the data.frame
NDVIcali$sectionID = ecorFACT[,1]
# check it
str(NDVIcali)
# get rid of unused factor levels
NDVIcali$sectionID = as.factor(as.character(NDVIcali$sectionID))
# plot the ecoregions spatially
plot(NDVIcali[,c('easting','northing')], pch = 19, 
	col = findColours(classIntervals(as.integer(NDVIcali$sectionID), 
	n = 19, style = 'fisher'), viridis(19)), cex = .4)
# plot locations of missing ecoregion data
points(NDVIcali[which(is.na(NDVIcali$sectionID)),c('easting','northing')], 
	pch = 19, cex = 1, col = 'red')
# they are all at the edges, so get rid of these points
NDVIcali = NDVIcali[!is.na(NDVIcali$sectionID),]
# check it
str(NDVIcali)

#-------------------------------------------------------------------------------
#              precipitation polygons
#-------------------------------------------------------------------------------

filename = '/home/jayverhoef/Desktop/precip_polys.shp'
precip = shapefile(filename)
precip = st_set_crs(st_as_sf(precip), "EPSG:4326")
# change it to a SpatialPolygonsDataFrame
precip = as(precip, 'Spatial')
# transform projection to CONUS Albers
precip = spTransform(precip, CRS("+init=epsg:5070"))
# plot it
plot(precip)
# check that it is the same projection as California
plot(cali, add = TRUE, col = 'red')
# get the precip category for each NDVI location
NDVIcaliSPDF = NDVIcali
coordinates(NDVIcaliSPDF) = ~ easting + northing
# set the CRS to CONUS Albers
NDVIcaliSPDF = st_set_crs(st_as_sf(NDVIcaliSPDF), "EPSG:5070")
# convert it back to a SpatialPointsDataFrame
NDVIcaliSPDF = as(NDVIcaliSPDF, 'Spatial')
# now get the ecoregions
precipFACT <- data.frame(xx=over(NDVIcaliSPDF, precip))
# add the ecoregions to the data.frame
NDVIcali$precip = as.factor(as.numeric(precipFACT[,2]))
sum(is.na(NDVIcali$precip))
# check it
str(NDVIcali)
# plot precipitation spatially
plot(NDVIcali[,c('easting','northing')], pch = 19, 
	col = findColours(classIntervals(as.integer(NDVIcali$precip), 
	n = 15, style = 'fisher'), viridis(15)), cex = .4)
# plot locations of missing ecoregion data
points(NDVIcali[which(is.na(NDVIcali$precip)),c('easting','northing')], 
	pch = 19, cex = 1, col = 'red')
# they are all at the edges, so get rid of these points
NDVIcali = NDVIcali[!is.na(NDVIcali$precip),]
# check it
str(NDVIcali)

#-------------------------------------------------------------------------------
#              final data.frame
#-------------------------------------------------------------------------------

plot(cali)
plot(precip)
plot(cali, add = TRUE, col = 'red')

# any missing data now?
any(is.na(NDVIcali$elev))
any(is.na(NDVIcali$sectionID))
any(is.na(NDVIcali$precip))
sum(is.na(NDVIcali$NDVI))

# final NDVI map
plot(NDVIcali[,c('easting','northing')], pch = 19, 
	col = findColours(classIntervals(NDVIcali$NDVI, n = 20, style = 'fisher'), 
	viridis(20)), cex = .4)
plot(cali, add = TRUE, lwd = 3)
# final elev map
plot(NDVIcali[,c('easting','northing')], pch = 19, 
	col = findColours(classIntervals(NDVIcali$elev, n = 20, style = 'fisher'), 
	viridis(20)), cex = .4)
plot(cali, add = TRUE, lwd = 3)
# final ecoregions map
pal = brewer.pal(n = 19, name = "Set1")
plot(NDVIcali[,c('easting','northing')], pch = 19, 
	col = findColours(classIntervals(as.integer(NDVIcali$sectionID), 
	n = 19, style = 'fisher'), pal), cex = .4)
plot(cali, add = TRUE, lwd = 3)
# final precipitation map
plot(NDVIcali[,c('easting','northing')], pch = 19, 
	col = findColours(classIntervals(as.integer(NDVIcali$precip), 
	n = 15, style = 'fisher'), viridis(15)), cex = .4)
plot(cali, add = TRUE, lwd = 3)

# write the data to disk
basepath = paste0('/mnt/ExtraDrive1/Work/desktop_data/2019_papers',
	'/fastREML/spPart_package/spPart/inst/rawdata/')
save(NDVIcali, file = paste0(basepath,'NDVIcali.rda'))
save(cali, file = paste0(basepath,'cali.rda'))
