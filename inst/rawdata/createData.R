setwd('/home/jay/Desktop/')
load("RData.USmonthlyMet.bin")
year.id <- 95
loc <- cbind(USpinfo$lon, USpinfo$lat, USpinfo$elev)
xr = c(-111, -99)
yr = c(32, 43)
station.subset<-  (loc[,1]>= xr[1]) & (loc[,1] <= xr[2]) & (loc[,2]>= yr[1]) & (loc[,2]<= yr[2])
ydata<-  USppt[ year.id,8,station.subset]
ydata <- ydata*10 #  cm -> mm conversion
xdata<- loc[station.subset,]
dimnames(xdata)<- list( USpinfo$station.id[station.subset], c( "lon", "lat", "elev"))
xdata<- data.frame( xdata)
good<- !is.na(ydata)
ydata<- ydata[good]
xdata<- xdata[good,]
str(xdata)
plot(xdata[,c('lon','lat')], pch = 19, cex = ydata/100)

d1 = data.frame(prec = ydata, lon = xdata$lon, lat = xdata$lat, elev = xdata$elev)
plot(d1$elev, d1$prec)

year.id <- 97
loc <- cbind(UStinfo$lon, UStinfo$lat, UStinfo$elev)
xr = c(-111, -95)
yr = c(32, 45)
station.subset<-  (loc[,1]>= xr[1]) & (loc[,1] <= xr[2]) & (loc[,2]>= yr[1]) & (loc[,2]<= yr[2])
ydata<-  (UStmax[ year.id,8,station.subset] + UStmin[ year.id,8,station.subset])/2
xdata<- loc[station.subset,]
d2 = data.frame(tempmax = ydata, lon = xdata[,1], lat = xdata[,2], 
  elev = xdata[,3]/1000)
d2 = d2[!is.na(d2$tempmax),]
d2 = d2[d2$elev < 3600/1000,]
str(d2)

path = paste0('/media/jay/data/desktop_data/2019_papers/',
  'fastREML/spPart_package/spPart/R/')
source(paste0(path,'LLtoUTM.R'))
xy = LLtoUTM(cm = mean(d2$lon), lat = d2$lat, lon = d2$lon)$xy
d2$x = xy[,1]/1000
d2$y = xy[,2]/1000

plot(d2[,c('x','y')], 
  pch = 19, cex = 0.5 + 3*(d2$tempmax-min(d2$tempmax))/
  (max(d2$tempmax) - min(d2$tempmax))
)

plot(d2$elev, d2$tempmax)
plot(d2$y, d2$tempmax)

summary(lm(tempmax ~ elev + I(elev^2) + y, data = d2))
d2$resids = residuals(lm(tempmax ~ elev + I(elev^2) + y , data = d2))
plot(d2[,c('x','y')], 
  pch = 19, cex = 0.5 + 3*(d2$resids-min(d2$resids))/
  (max(d2$resids) - min(d2$resids))
)
plot(d2$elev, d2$resids)
plot(d2$y, d2$resids)
d2 = d2[abs(d2$resids) < 6,]

d2$grp = as.integer(
  as.factor(
    paste0(as.character(cut(d2$x, 
    breaks = c(-.001,250,500,750,1000,1250,1600)/1000)),
    as.character(cut(d2$y, 
    breaks = c(-.001,250,500,750,1000,1250,1600)/1000)))
  )
)
ngr_rand = round(dim(d2)[1]/50,0)
set.seed(1001)
d2$grpr = sample(1:ngr_rand, dim(d2)[1], replace = TRUE)
spatial_sd = data.frame(
  grp = aggregate(d2$x, by = list(d2$grp), mean)[,1],
  x = aggregate(d2$x, by = list(d2$grp), mean)[,2], 
  y = aggregate(d2$y, by = list(d2$grp), mean)[,2],
  sd = aggregate(d2$resids, by = list(d2$grp), sd)[,2]
)
plot(spatial_sd[,c('x','y')], cex = 1 +5*(spatial_sd$sd - min(spatial_sd$sd))/
  (max(spatial_sd$sd) - min(spatial_sd$sd)), pch = 19
)

library(gstat)

vgm_resid_dir <- variogram(resids ~ 1, loc=~x+y, data=d2[d2$x < 800/1000,], 
  alpha=c(0,45,90,135), cutoff = 500/1000)
  old_par = par(mar = c(5,5,3,1))
  plot(vgm_resid_dir$dist, vgm_resid_dir$gamma, type = 'n', xlab = 'Distance (km)',
    ylab = 'Semivariogram', cex.lab = 2, cex.axis = 1.5, ylim = c(0,3))
  points(vgm_resid_dir[vgm_resid_dir$dir.hor == 0,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 0,'gamma'], pch = 1, 
    cex = 8*vgm_resid_dir$np/max(vgm_resid_dir$np))
  lines(vgm_resid_dir[vgm_resid_dir$dir.hor == 0,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 0,'gamma'], lwd = 2)
  points(vgm_resid_dir[vgm_resid_dir$dir.hor == 45,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 45,'gamma'], pch = 2, 
    cex = 8*vgm_resid_dir$np/max(vgm_resid_dir$np))
  lines(vgm_resid_dir[vgm_resid_dir$dir.hor == 45,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 45,'gamma'], lwd = 2)
  points(vgm_resid_dir[vgm_resid_dir$dir.hor == 90,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 90,'gamma'], pch = 5, 
    cex = 8*vgm_resid_dir$np/max(vgm_resid_dir$np))
  lines(vgm_resid_dir[vgm_resid_dir$dir.hor == 90,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 90,'gamma'], lwd = 2)
  points(vgm_resid_dir[vgm_resid_dir$dir.hor == 135,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 135,'gamma'], pch = 22, 
    cex = 8*vgm_resid_dir$np/max(vgm_resid_dir$np))
  lines(vgm_resid_dir[vgm_resid_dir$dir.hor == 135,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 135,'gamma'], lwd = 2)
  legend(0.05, 3.3, legend=(c('N-S','NE-SW','E-W','SE-NW')), 
    pch = c(1,2,5,22), cex = 2)
  par(old_par)

vgm_resid_dir <- variogram(resids ~ 1, loc=~x+y, data=d2[d2$x > 800/1000,], 
  alpha=c(0,45,90,135), cutoff = 500/1000)
  old_par = par(mar = c(5,5,3,1))
  plot(vgm_resid_dir$dist, vgm_resid_dir$gamma, type = 'n', xlab = 'Distance (km)',
    ylab = 'Semivariogram', cex.lab = 2, cex.axis = 1.5, ylim = c(0,3.5))
  points(vgm_resid_dir[vgm_resid_dir$dir.hor == 0,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 0,'gamma'], pch = 1, 
    cex = 8*vgm_resid_dir$np/max(vgm_resid_dir$np))
  lines(vgm_resid_dir[vgm_resid_dir$dir.hor == 0,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 0,'gamma'], lwd = 2)
  points(vgm_resid_dir[vgm_resid_dir$dir.hor == 45,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 45,'gamma'], pch = 2, 
    cex = 8*vgm_resid_dir$np/max(vgm_resid_dir$np))
  lines(vgm_resid_dir[vgm_resid_dir$dir.hor == 45,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 45,'gamma'], lwd = 2)
  points(vgm_resid_dir[vgm_resid_dir$dir.hor == 90,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 90,'gamma'], pch = 5, 
    cex = 8*vgm_resid_dir$np/max(vgm_resid_dir$np))
  lines(vgm_resid_dir[vgm_resid_dir$dir.hor == 90,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 90,'gamma'], lwd = 2)
  points(vgm_resid_dir[vgm_resid_dir$dir.hor == 135,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 135,'gamma'], pch = 22, 
    cex = 8*vgm_resid_dir$np/max(vgm_resid_dir$np))
  lines(vgm_resid_dir[vgm_resid_dir$dir.hor == 135,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 135,'gamma'], lwd = 2)
  legend(0.05, 3.5, legend=(c('N-S','NE-SW','E-W','SE-NW')), 
    pch = c(1,2,5,22), cex = 2)
  par(old_par)


  formula = tempmax ~ elev + I(elev^2) + y
	trms <- terms(formula, data = d2)
	respCol <- as.character(as.list(attr(trms,"variables")[-1]))[1]
	covList <- attr(trms,"term.labels")
	#design matrix
	X <- model.matrix(formula, d2)
	z <- as.matrix(d2[, respCol])
	#vectors of spatial coordinates
	xcoords <- d2[,'x']
	ycoords <- d2[,'y']
	n <- length(z)
	p <- sum(svd(X)$d>1e-10)
  grpindx = d2$grpr
  
  xylist = vector(mode = "list", length = max(grpindx)) 
  zlist = vector(mode = "list", length = max(grpindx)) 
  Xlist = vector(mode = "list", length = max(grpindx)) 
  for(i in 1:max(grpindx)) {
    xylist[[i]] = cbind(xcoords[grpindx == i], ycoords[grpindx == i])
    zlist[[i]] = z[grpindx == i]
    Xlist[[i]] = X[grpindx == i,]
	}

distGeoAni(xylist[[1]][,1], xylist[[1]][,2], xylist[[1]][,1], xylist[[1]][,2],
  0, 300, 0.5)


expit = function(x){exp(x)/(1 + exp(x))}
logit = function(x){log(x/(1 - x))}
alph = function(x,mu,scale){exp(scale*(x - mu))/(1 + exp(scale*(x - mu)))}
plot(1:1600, alph(1:1600, 800, -.02))

#'@export
exponential <- function(distance.matrix)
{
	exp(-3*distance.matrix) 
}
spherical <- function(distance.matrix)
{
	CovMat <- (1 - 1.5*distance.matrix + 0.5*distance.matrix^3)
	CovMat[distance.matrix > 1] <- 0
	CovMat
}

theta = c(log(1.5), log(3), log(.4), logit(.1), logit(.7), log(0.5),
  log(3), log(.8), .8, -.02)
theta = c(log(1.5), log(3), log(.4), logit(.4), logit(.7), log(0.5),
  log(3), log(.8), .8, -.02)
  
m2LLnons(theta, zlist, Xlist, xylist)

	parmest <-optim(theta, m2LLnons, zlist = zlist, Xlist = Xlist, 
			xylist = xylist)
m2LLnons(parmest$par, zlist = zlist, Xlist = Xlist, 
			xylist = xylist)
      
exp(parmest$par[1])
exp(parmest$par[2])
exp(parmest$par[3])
expit(parmest$par[4])*90
expit(parmest$par[5])
exp(parmest$par[6])
exp(parmest$par[7])
exp(parmest$par[8])
parmest$par[9]
parmest$par[10]

theta = parmest$par
theta[1] = 1.5
theta[9] = .8
alph((1:1600)/1000, parmest$par[9], parmest$par[10])
plot((1:1600)/1000, alph((1:1600)/1000, parmest$par[9], parmest$par[10]))


m2LLnons <- function(theta, zlist, Xlist, xylist)
{
  # theta[1] is nugW, theta[2] is parsW, theta[3] is rangW, theta[4] is rotaW,
  # theta[5] is minpW, theta[6] is nugE, theta[7] is parsE, theta[8] is rangE,
  # theta[9] is logimu, theta[10] is logiscale
  # if(any(abs(theta) > 10)) return(1e+30)
  if(exp(theta[3]) > 2) return(1e+30)
  if(exp(theta[8]) > 2) return(1e+30)
  p = length(Xlist[[1]][1,])
  ngrps = length(zlist)
  qrlist = vector("list", ngrps)
  Sxx = matrix(0, nrow = p, ncol = p)
  Sxy = matrix(0, nrow = p, ncol = 1)
  logDetV = 0
  logDetXViX = 0
  for(i in 1:ngrps) {
	  dismat <- distGeoAni(xylist[[i]][,1], xylist[[i]][,2], xylist[[i]][,1], 
      xylist[[i]][,2], 90*expit(theta[4]), exp(theta[3]), expit(theta[5]))
	  covMatW <- exp(theta[2])*spherical(dismat) 
    diag(covMatW) = diag(covMatW) + (exp(theta[1]) + 1e-10*exp(theta[1]))
	  dismat <- distGeoAni(xylist[[i]][,1], xylist[[i]][,2], 
      xylist[[i]][,1], xylist[[i]][,2], 0, exp(theta[8]), 1)
	  covMatE <- exp(theta[7])*spherical(dismat) 
    diag(covMatE) = diag(covMatE) + (exp(theta[6]) + 1e-10*exp(theta[6]))
    covMat = outer(alph(xylist[[i]][,1],theta[9], theta[10]),
      alph(xylist[[i]][,1],theta[9], theta[10]))*covMatW + 
      outer(1 - alph(xylist[[i]][,1],theta[9], theta[10]),
      1 - alph(xylist[[i]][,1],theta[9], theta[10]))*covMatE
    qrlist[[i]] = qr(covMat, LAPACK = TRUE)
    ViX = solve(qrlist[[i]], Xlist[[i]])
    XViX = crossprod(Xlist[[i]],ViX)
    Sxx = Sxx + XViX
    Sxy = Sxy + t(crossprod(zlist[[i]],ViX)) 
    logDetV = logDetV + sum(log(abs(diag(qr.R(qrlist[[i]])))))
    logDetXViX = logDetXViX + as.numeric(determinant(XViX, 
      logarithm = TRUE)$modulus)
  }
	
  betaHat = solve(Sxx, Sxy)
  rVir = 0
  for(i in 1:ngrps) 
    rVir = rVir + sum((zlist[[i]] - 
      Xlist[[i]] %*% betaHat) *
      solve(qrlist[[i]], (zlist[[i]] - 
      Xlist[[i]] %*% betaHat)))

	minus2LL1 <- logDetV + rVir + logDetXViX

  p = length(Xlist[[1]][1,])
  ngrps = length(zlist)
  qrlist = vector("list", ngrps)
  Sxx = matrix(0, nrow = p, ncol = p)
  Sxy = matrix(0, nrow = p, ncol = 1)
  logDetV = 0
  logDetXViX = 0
  for(i in 1:ngrps) {
	  dismat <- distGeoAni(xylist[[i]][,1], xylist[[i]][,2], xylist[[i]][,1], 
      xylist[[i]][,2], 180 - 90*expit(theta[4]), exp(theta[3]), expit(theta[5]))
	  covMatW <- exp(theta[2])*spherical(dismat) 
    diag(covMatW) = diag(covMatW) + (exp(theta[1]) + 1e-10*exp(theta[1]))
	  dismat <- distGeoAni(xylist[[i]][,1], xylist[[i]][,2], 
      xylist[[i]][,1], xylist[[i]][,2], 0, exp(theta[8]), 1)
	  covMatE <- exp(theta[7])*spherical(dismat) 
    diag(covMatE) = diag(covMatE) + (exp(theta[6]) + 1e-10*exp(theta[6]))
    covMat = outer(alph(xylist[[i]][,1],theta[9], theta[10]),
      alph(xylist[[i]][,1],theta[9], theta[10]))*covMatW + 
      outer(1 - alph(xylist[[i]][,1],theta[9], theta[10]),
      1 - alph(xylist[[i]][,1],theta[9], theta[10]))*covMatE
    qrlist[[i]] = qr(covMat, LAPACK = TRUE)
    ViX = solve(qrlist[[i]], Xlist[[i]])
    XViX = crossprod(Xlist[[i]],ViX)
    Sxx = Sxx + XViX
    Sxy = Sxy + t(crossprod(zlist[[i]],ViX)) 
    logDetV = logDetV + sum(log(abs(diag(qr.R(qrlist[[i]])))))
    logDetXViX = logDetXViX + as.numeric(determinant(XViX, 
      logarithm = TRUE)$modulus)
  }
	
  betaHat = solve(Sxx, Sxy)
  rVir = 0
  for(i in 1:ngrps) 
    rVir = rVir + sum((zlist[[i]] - 
      Xlist[[i]] %*% betaHat) *
      solve(qrlist[[i]], (zlist[[i]] - 
      Xlist[[i]] %*% betaHat)))

	minus2LL2 <- logDetV + rVir + logDetXViX

  m2LL =  as.numeric(min(minus2LL1,minus2LL2))
  attr(m2LL, 'which') = which(c(minus2LL1,minus2LL2) == m2LL)[1]
  return(m2LL)
}

distGeoAni <- function(xrow, yrow, xcol, ycol, rotate = 0, range = 1, minorp = 1)
{
	# total number of observations for each set of coordinates  
		n.rows <- length(xrow)
		n.cols <- length(xcol)
	# expand all x-coordinates
		sxr <- matrix(xrow, ncol = 1) %*% 
			matrix(rep(1,times = n.cols), nrow = 1)
		sxc <- matrix(rep(1,times = n.rows), ncol = 1) %*% 
			matrix(xcol, nrow = 1)
		syr <- matrix(yrow,ncol = 1) %*% 
			matrix(rep(1,times = n.cols), nrow = 1)
		syc <- matrix(rep(1,times = n.rows), ncol = 1) %*% 
			matrix(ycol, nrow = 1)
	# find difference in coordinates between all pairwise locations
		sxdif <- sxr - sxc
		sydif <- syr - syc
	# rotate coordinates
		newx <- cos(rotate*.0174533)*sxdif - sin(rotate*.0174533)*sydif
		newy <- sin(rotate*.0174533)*sxdif + cos(rotate*.0174533)*sydif
	# scale coordinates by minor and major axes */
		newx <- newx/(range*minorp)
		newy <- newy/range
	# compute distance for the scaled and rotated coordinates */
		sqrt(newx^2 + newy^2)
}

