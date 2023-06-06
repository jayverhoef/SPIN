library(SSNbd)
library(doParallel)
library(viridis)
library(rgdal)
datapath = '/home/xverhoef/data/2019_packages/midCol_data_package/inst/'
datapath = '/media/jay/data/desktop_data/2019_packages/midCol_data_package/inst/'

# only need to do this once to create distance matrices
#ssn = importSSN(paste0(packpath,'midcolumbiaLSN.ssn'),
#  o.write = TRUE, predpts = 'pred_601')
#createBigDistMat(ssn, predpts = "pred_601", o.write = TRUE, no.cores = 7)

#ssn = importSSN(paste0(packpath,'midcolumbiaLSN.ssn'),
#  o.write = TRUE, predpts = 'pred_701')
#createBigDistMat(ssn, predpts = "pred_701", o.write = TRUE, no.cores = 7)

#ssn = importSSN(paste0(packpath,'midcolumbiaLSN.ssn'),
#  o.write = TRUE, predpts = 'pred_702')
#createBigDistMat(ssn, predpts = "pred_702", o.write = TRUE, no.cores = 7)

ssn = importSSN(paste0(datapath,'midcolumbiaLSN.ssn'),
  o.write = TRUE, predpts = 'pred_601')
DF = getSSNdata.frame(ssn)
dim(DF)[1]
DFp = getSSNdata.frame(ssn, Name = 'pred_601')
dim(DFp)[1]
ssn = importPredpts(target = ssn, predpts ='pred_701', obj.type = "ssn")
DFp1 = getSSNdata.frame(ssn, 'pred_701')
dim(DFp1)[1]
ssn = importPredpts(target = ssn, predpts ='pred_702', obj.type = "ssn")
DFp2 = getSSNdata.frame(ssn, 'pred_702')
dim(DFp2)[1]
dim(DFp)[1] + dim(DFp1)[1] + dim(DFp2)[1]

Washington <- readOGR(dsn = paste0(datapath,'stateSHP/'), 
	layer = "WA_State_Boundary")
Washington = spTransform(Washington, 
	CRS(proj4string(as.SpatialLinesDataFrame(ssn))))
Oregon <- readOGR(dsn = paste0(datapath,'stateSHP/'), 
	layer = "or_state_boundary")
Oregon = spTransform(Oregon, 
	CRS(proj4string(as.SpatialLinesDataFrame(ssn))))

DFl = as.SpatialLinesDataFrame(ssn)@data
path = paste0('/home/xverhoef/data/2019_papers/fastREML/spPart_package',
	'/spPart/inst/doc/figure/')
linecol = '#377eb8'
obscol = '#984ea3'
predcol = '#fdb462'
png(paste0(path,'StudyArea.png'),
	width = 1920, height = 1920)
	layout(matrix(c(1,1,2,3), nrow = 2, ncol = 2, byrow = TRUE),
		heights = c(2,1))
	par(mar = c(0,0,0,0))
	plot(as.SpatialLinesDataFrame(ssn), col = linecol, lwd = 2)
	points(as.SpatialPoints(ssn, data = 'Obs'), pch = 19, col = obscol)
	lines(c(872706, 872706, 900345, 900345, 872706),
		c(1757134, 1790049, 1790049, 1757134, 1757134), lwd = 8)
	text(820000, 1910000, label = '(a)', cex = 8)

	par(mar = c(0,10,0,0))
	plot(Washington, xlim = c(450000, 1301037), ylim = c(1283766, 2171702),
		lwd = 4)
	plot(Oregon, lwd = 4, add = TRUE)
	plot(as.SpatialLinesDataFrame(ssn), add = TRUE, col = linecol)
	text(1015441, 2000000, label = 'Washington', cex = 5)
	text(976648, 1408001, label = 'Oregon', cex = 5)
	text(500000, 2136499, label = '(b)', cex = 8)

	par(mar = c(4,0,0,15))
	plot(as.SpatialLinesDataFrame(ssn), col = linecol, lwd = 5,
	  xlim = c(872706, 900345), ylim = c(1757134, 1790049))
	points(as.SpatialPoints(ssn, data = 'pred_601'), pch = 19, cex = 2,
		col = predcol)
	points(as.SpatialPoints(ssn, data = 'pred_701'), pch = 19, cex = 2,
		col = predcol)
	points(as.SpatialPoints(ssn, data = 'pred_702'), pch = 19, cex = 2,
		col = predcol)
	points(as.SpatialPoints(ssn, data = 'Obs'), pch = 19, col = obscol,
	cex = 3)
	lines(c(872706, 872706, 900345, 900345, 872706),
		c(1757134, 1790049, 1790049, 1757134, 1757134), lwd = 8)
	text(867522.3, 1789000, label = '(c)', cex = 8)
dev.off()


  registerDoParallel(cores = 7)
  set.seed(101)
	ssnr = mapi(ssn, nIndx = 100)
	start = Sys.time()
	#undebug(cope)
	copeOut = cope(STREAM_AUG ~ I(ELEV/1000) + I(SLOPE*100) + NLCD11PC + PRECIP + I(SNAP_Y/100000) + 
  BFI + I(CUMDRAINAG/10000) + CANOPY + Air_Aug + Flow_Aug, ssnr,
		CorModels = c("Exponential.tailup", "Exponential.taildown", 
		"Exponential.Euclid"), use.nugget = TRUE, 
		partIndxCol = "partIndx", addfunccol = 'afvArea', 
		parallel = TRUE)
	end = Sys.time()
	copeOutTime = difftime(end, start, units = 'mins')
	copeOut
	copeOutTime

	start = Sys.time()
	#undebug(fefe)
	fefeOut = fefe(copeOut)
	end = Sys.time()
	fefeOutTime = difftime(end, start, units = 'mins')
	fefeOut
	fefeOutTime

	# --------------------------------------------------------------------
	# 601 Pred group
	# --------------------------------------------------------------------

#	ssn = importSSN(paste0(packpath,'/midcolumbiaLSN.ssn'),
#		o.write = TRUE, predpts = 'pred_601')
#	createBigDistMat(ssn, predpts = "pred_601", o.write = TRUE, no.cores = 7)

	start = Sys.time()
	pureOut = pure(ecp = copeOut, efe = fefeOut, predsID = 'pred_601', nNN = 25)
	end = Sys.time()
	pureOutTime = difftime(end, start, units = 'mins')
	pureOutTime
	pure601 = pureOut
	path = '/media/jay/Hitachi2GB/00NMML/ActiveRPack/SSNbd_package/'
	save(pure601, file = paste0(path,'SSNbd/data/preds601.rda'))
	
	# --------------------------------------------------------------------
	# 701 Pred group
	# --------------------------------------------------------------------

	ssn = importSSN(paste0(packpath,'/midcolumbiaLSN.ssn'),
		o.write = TRUE, predpts = 'pred_701')
	createBigDistMat(ssn, predpts = "pred_701", o.write = TRUE, no.cores = 7)
	
  registerDoParallel(cores=7)
  set.seed(101)
	ssnr = masi(ssn, nIndx = 50)
	start = Sys.time()
	copeOut = cope(STREAM_AUG ~ I(ELEV/1000) + I(SLOPE*100) + NLCD11PC + PRECIP + I(SNAP_Y/100000) + 
  BFI + I(CUMDRAINAG/10000) + CANOPY + Air_Aug + Flow_Aug, ssnr,
		CorModels = c("Exponential.tailup", "Exponential.taildown", 
		"Exponential.Euclid"), use.nugget = TRUE, 
		subSampIndxCol = "subSampIndx", addfunccol = 'afvArea', 
		parallel = TRUE)
	end = Sys.time()
	copeOutTime = difftime(end, start, units = 'mins')
	copeOut
	copeOutTime

	start = Sys.time()
	fefeOut = fefe(copeOut)
	end = Sys.time()
	fefeOutTime = difftime(end, start, units = 'mins')
	fefeOut
	fefeOutTime

	start = Sys.time()
	pureOut = pure(ecp = copeOut, efe = fefeOut, predsID = 'pred_701', nNN = 25)
	end = Sys.time()
	pureOutTime = difftime(end, start, units = 'mins')
	pure701 = pureOut
	save(pure701, file = paste0(path,'SSNbd/data/preds701.rda'))

	# --------------------------------------------------------------------
	# 702 Pred group
	# --------------------------------------------------------------------

	start = Sys.time()
	pureOut = pure(ecp = copeOut, efe = fefeOut, predsID = 'pred_702', nNN = 25)
	end = Sys.time()
	pureOutTime = difftime(end, start, units = 'mins')
	pure702 = pureOut
	save(pure702, file = paste0(path,'SSNbd/data/preds702.rda'))

	
library(viridis)
class <- classIntervals(c(pure601$pred, pure701$pred, pure702$pred), 
	n = 8, style="fisher")
colcode <- findColours(class, viridis(8))
col601 = colcode[1:length(pure601[,1])]
col701 = colcode[(length(pure601[,1]) + 1):(length(pure601[,1]) + length(pure701[,1]))]
col702 = colcode[(length(pure601[,1]) + length(pure701[,1]) + 1):length(colcode)]

png('/media/jay/Hitachi2GB/00NMML/ActivePresentations/AFS pres 2018/figure/MidColPreds.png',
	width = 1920, height = 1920)
	par(bg=NA) 
	plot(as.SpatialLinesDataFrame(ssn), col = 'grey50', lwd = 3)
	points(as.SpatialPoints(ssn, data = 'pred_601'), pch = 19, cex = .5,
		col = col601)
	points(as.SpatialPoints(ssn, data = 'pred_701'), pch = 19, cex = .5,
		col = col701)
	points(as.SpatialPoints(ssn, data = 'pred_702'), pch = 19, cex = .5,
		col = col702)
dev.off()
hist(c(pure601$pred, pure701$pred, pure702$pred))
