basepath = '/media/jay/data/desktop_data/2019_papers/fastREML/spPart_package/spPart/data/'
load(paste0(basepath,'scen1_1000pt_NNGP_est1.rda'))
load(paste0(basepath,'scen1_1000pt_NNGP_est1se.rda'))
load(paste0(basepath,'scen1_1000pt_NNGP_est2.rda'))
load(paste0(basepath,'scen1_1000pt_NNGP_est2se.rda'))
load(paste0(basepath,'scen1_1000pt_NNGP_RMSPE.rda'))
load(paste0(basepath,'scen1_1000pt_NNGP_PI90.rda'))
load(paste0(basepath,'scen1_1000pt_NNGP_time.rda'))

library(spPart)
data(scen1_1000pt_NNGP_est1)
data(scen1_1000pt_NNGP_est1se)
data(scen1_1000pt_NNGP_est2)
data(scen1_1000pt_NNGP_est2se)
data(scen1_1000pt_NNGP_RMSPE)
data(scen1_1000pt_NNGP_PI90)
data(scen1_1000pt_NNGP_time)


FullCompNNGP = data.frame(
  METH = c('Full','spPart','NNGP'),
  MSE1 = sqrt(apply((scen1_1000pt_NNGP_est1 - 1)^2,2,mean)),
  MSE2 = sqrt(apply((scen1_1000pt_NNGP_est2 - 1)^2,2,mean)),
  RMSPE = apply(scen1_1000pt_NNGP_RMSPE,2,mean),
  CI190 = apply(scen1_1000pt_NNGP_est1 - qnorm(.95)*scen1_1000pt_NNGP_est1se < 1 & 
    1 < scen1_1000pt_NNGP_est1 + qnorm(.95)*scen1_1000pt_NNGP_est1se,2,mean),
  CI290 = apply(scen1_1000pt_NNGP_est2 - qnorm(.95)*scen1_1000pt_NNGP_est2se < 1 & 
    1 < scen1_1000pt_NNGP_est2 + qnorm(.95)*scen1_1000pt_NNGP_est2se,2,mean),
  PI90 = apply(scen1_1000pt_NNGP_PI90,2,mean),
  time = apply(scen1_1000pt_NNGP_time,2,mean)
)

  print(
    xtable(FullCompNNGP, 
      align = c('l',rep('l', times = length(FullCompNNGP[1,]))),
      digits = c(0,rep(3, times = length(FullCompNNGP[1,]))),
    ),
    size = 'footnotesize',
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
  )
