basepath = '/media/jay/data/2019_papers/fastREML/spPart_package/spPart/data/'
load(paste0(basepath,'effBlockType_est1.rda'))
load(paste0(basepath,'effBlockType_est1se.rda'))
load(paste0(basepath,'effBlockType_est2.rda'))
load(paste0(basepath,'effBlockType_est2se.rda'))
load(paste0(basepath,'effBlockType_RMSPE.rda'))
load(paste0(basepath,'effBlockType_PI90.rda'))

library(spPart)
data(partSizes_est1)
data(partSizes_est1se)
data(partSizes_est2)
data(partSizes_est2se)
data(partSizes_RMSPE)
data(partSizes_PI90)
data(partSizes_copetime)
data(partSizes_fefetime)

partSizes = data.frame(
  COPE = kronecker(c(25,50,100,200),rep(1, times = 4)),
  FEFE = kronecker(rep(1, times = 4),c(25,50,100,200)),
  MSE1 = sqrt(apply((partSizes_est1 - 1)^2,2,mean)),
  MSE2 = sqrt(apply((partSizes_est2 - 1)^2,2,mean)),
  CI190 = apply(partSizes_est1 - qnorm(.95)*partSizes_est1se < 1 & 
    1 < partSizes_est1 + qnorm(.95)*partSizes_est1se,2,mean),
  CI290 = apply(partSizes_est2 - qnorm(.95)*partSizes_est2se < 1 & 
    1 < partSizes_est2 + qnorm(.95)*partSizes_est2se,2,mean),
  RMSPE = apply(partSizes_RMSPE,2,mean),
  PI90 = apply(partSizes_PI90,2,mean),
  COPEtime = kronecker(apply(partSizes_copetime,2,mean),rep(1, times = 4)),
  FEFEtime = kronecker(rep(1, times = 4),apply(partSizes_fefetime,2,mean))
)
