basepath = '/media/jay/data/2019_papers/fastREML/spPart_package/spPart/data/'
load(paste0(basepath,'effBlockType_est1.rda'))
load(paste0(basepath,'effBlockType_est1se.rda'))
load(paste0(basepath,'effBlockType_est2.rda'))
load(paste0(basepath,'effBlockType_est2se.rda'))
load(paste0(basepath,'effBlockType_RMSPE.rda'))
load(paste0(basepath,'effBlockType_PI90.rda'))

library(spPart)
data(fefealt_est1)
data(fefealt_est2)
data(fefealt_est1_1)
data(fefealt_est2_1)
fefealt_est1 = rbind(fefealt_est1, fefealt_est1_1)
fefealt_est2 = rbind(fefealt_est2, fefealt_est2_1)

fefealt = matrix(NA, nrow = 4, ncol = 6)
fefealt[1,1] = mean(fefealt_est1[,1] - qnorm(.95)*fefealt_est1[,2] < 1 & 
    1 < fefealt_est1[,1] + qnorm(.95)*fefealt_est1[,2])
fefealt[1,2] = mean(fefealt_est1[,1] - qnorm(.95)*fefealt_est1[,3] < 1 & 
    1 < fefealt_est1[,1] + qnorm(.95)*fefealt_est1[,3])
fefealt[1,3] = mean(fefealt_est1[,1] - qnorm(.95)*fefealt_est1[,4] < 1 & 
    1 < fefealt_est1[,1] + qnorm(.95)*fefealt_est1[,4])
fefealt[2,1] = mean(fefealt_est1[,5] - qnorm(.95)*fefealt_est1[,6] < 1 & 
    1 < fefealt_est1[,5] + qnorm(.95)*fefealt_est1[,6])
fefealt[2,2] = mean(fefealt_est1[,5] - qnorm(.95)*fefealt_est1[,7] < 1 & 
    1 < fefealt_est1[,5] + qnorm(.95)*fefealt_est1[,7])
fefealt[2,3] = mean(fefealt_est1[,5] - qnorm(.95)*fefealt_est1[,8] < 1 & 
    1 < fefealt_est1[,5] + qnorm(.95)*fefealt_est1[,8])
fefealt[3,1] = mean(fefealt_est1[,9] - qnorm(.95)*fefealt_est1[,10] < 1 & 
    1 < fefealt_est1[,9] + qnorm(.95)*fefealt_est1[,10])
fefealt[3,2] = mean(fefealt_est1[,9] - qnorm(.95)*fefealt_est1[,11] < 1 & 
    1 < fefealt_est1[,9] + qnorm(.95)*fefealt_est1[,11])
fefealt[3,3] = mean(fefealt_est1[,9] - qnorm(.95)*fefealt_est1[,12] < 1 & 
    1 < fefealt_est1[,9] + qnorm(.95)*fefealt_est1[,12])
fefealt[4,1] = mean(fefealt_est1[,13] - qnorm(.95)*fefealt_est1[,14] < 1 & 
    1 < fefealt_est1[,13] + qnorm(.95)*fefealt_est1[,14])
fefealt[4,2] = mean(fefealt_est1[,13] - qnorm(.95)*fefealt_est1[,15] < 1 & 
    1 < fefealt_est1[,13] + qnorm(.95)*fefealt_est1[,15])
fefealt[4,3] = mean(fefealt_est1[,13] - qnorm(.95)*fefealt_est1[,16] < 1 & 
    1 < fefealt_est1[,13] + qnorm(.95)*fefealt_est1[,16])
    
fefealt[1,4] = mean(fefealt_est2[,1] - qnorm(.95)*fefealt_est2[,2] < 1 & 
    1 < fefealt_est2[,1] + qnorm(.95)*fefealt_est2[,2])
fefealt[1,5] = mean(fefealt_est2[,1] - qnorm(.95)*fefealt_est2[,3] < 1 & 
    1 < fefealt_est2[,1] + qnorm(.95)*fefealt_est2[,3])
fefealt[1,6] = mean(fefealt_est2[,1] - qnorm(.95)*fefealt_est2[,4] < 1 & 
    1 < fefealt_est2[,1] + qnorm(.95)*fefealt_est2[,4])
fefealt[2,4] = mean(fefealt_est2[,5] - qnorm(.95)*fefealt_est2[,6] < 1 & 
    1 < fefealt_est2[,5] + qnorm(.95)*fefealt_est2[,6])
fefealt[2,5] = mean(fefealt_est2[,5] - qnorm(.95)*fefealt_est2[,7] < 1 & 
    1 < fefealt_est2[,5] + qnorm(.95)*fefealt_est2[,7])
fefealt[2,6] = mean(fefealt_est2[,5] - qnorm(.95)*fefealt_est2[,8] < 1 & 
    1 < fefealt_est2[,5] + qnorm(.95)*fefealt_est2[,8])
fefealt[3,4] = mean(fefealt_est2[,9] - qnorm(.95)*fefealt_est2[,10] < 1 & 
    1 < fefealt_est2[,9] + qnorm(.95)*fefealt_est2[,10])
fefealt[3,5] = mean(fefealt_est2[,9] - qnorm(.95)*fefealt_est2[,11] < 1 & 
    1 < fefealt_est2[,9] + qnorm(.95)*fefealt_est2[,11])
fefealt[3,6] = mean(fefealt_est2[,9] - qnorm(.95)*fefealt_est2[,12] < 1 & 
    1 < fefealt_est2[,9] + qnorm(.95)*fefealt_est2[,12])
fefealt[4,4] = mean(fefealt_est2[,13] - qnorm(.95)*fefealt_est2[,14] < 1 & 
    1 < fefealt_est2[,13] + qnorm(.95)*fefealt_est2[,14])
fefealt[4,5] = mean(fefealt_est2[,13] - qnorm(.95)*fefealt_est2[,15] < 1 & 
    1 < fefealt_est2[,13] + qnorm(.95)*fefealt_est2[,15])
fefealt[4,6] = mean(fefealt_est2[,13] - qnorm(.95)*fefealt_est2[,16] < 1 & 
    1 < fefealt_est2[,13] + qnorm(.95)*fefealt_est2[,16])
fefealt = cbind(c(25,50,100,200),fefealt) 
library(xtable)
  print(
    xtable(fefealt, 
      align = c('l',rep('l', times = length(fefealt[1,]))),
      digits = c(0,0,rep(3, times = 6)),
    ),
    size = 'footnotesize',
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
  )
