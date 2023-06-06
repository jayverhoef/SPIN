basepath = '/media/jay/data/2019_papers/fastREML/spPart_package/spPart/data/'
load(paste0(basepath,'effBlockType_est1.rda'))
load(paste0(basepath,'effBlockType_est1se.rda'))
load(paste0(basepath,'effBlockType_est2.rda'))
load(paste0(basepath,'effBlockType_est2se.rda'))
load(paste0(basepath,'effBlockType_RMSPE.rda'))
load(paste0(basepath,'effBlockType_PI90.rda'))

library(spPart)
data(effBlockType_est1)
data(effBlockType_est1se)
data(effBlockType_est2)
data(effBlockType_est2se)
data(effBlockType_RMSPE)
data(effBlockType_PI90)

effBlockType = data.frame(
  COPE = c('RAND','RAND','RAND','COMP','COMP','COMP','MIXD','MIXD','MIXD'),
  FEFE = c('RAND','COMP','MIXD','RAND','COMP','MIXD','RAND','COMP','MIXD'),
  MSE1 = sqrt(apply((effBlockType_est1 - 1)^2,2,mean)),
  MSE2 = sqrt(apply((effBlockType_est2 - 1)^2,2,mean)),
  CI190 = apply(effBlockType_est1 - qnorm(.95)*effBlockType_est1se < 1 & 
    1 < effBlockType_est1 + qnorm(.95)*effBlockType_est1se,2,mean),
  CI290 = apply(effBlockType_est2 - qnorm(.95)*effBlockType_est2se < 1 & 
    1 < effBlockType_est2 + qnorm(.95)*effBlockType_est2se,2,mean),
  RMSPE = apply(effBlockType_RMSPE,2,mean),
  PI90 = apply(effBlockType_PI90,2,mean)
)
