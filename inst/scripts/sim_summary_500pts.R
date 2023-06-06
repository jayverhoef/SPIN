library('spPart')
library('xtable')
data(scen1_500pts_est)
data(scen1_500pts_estse)
data(scen1_500pts_RMSPE)
data(scen1_500pts_PI90)
data(scen1_500pts_blk_pred)
data(scen1_500pts_blk_se)


scen1_500pts_est_summ = data.frame(
	RMSE = sqrt(apply(scen1_500pts_est - 1,2,var)),
  CI90 = apply(scen1_500pts_est - qnorm(.95)*scen1_500pts_estse < 1 & 
	1 < scen1_500pts_est + qnorm(.95)*scen1_500pts_estse, 2, mean),
	CovParmEst = c('all', 'rand', NA, 'comp', NA, NA, 'zimm', NA, NA),
	FixEffEst = c('all', 'all', 'rand', 'all', 'rand', 'comp', 'all',
		'rand', 'zimm')
)
DFfefe500sc1 = data.frame(metric = c('RMSE','CI90') ,rbind(scen1_500pts_est_summ$RMSE,
  scen1_500pts_est_summ$CI90))
print(
    xtable(DFfefe500sc1, 
      align = c('l',rep('l', times = length(DFfefe500sc1[1,]))),
      digits = c(0,3,3,3,3,3,3,3,3,3,3),
      caption = 'Fefe results',
      label = 'tab:fefe500sc1'
    ),
    size = 'footnotesize',
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
  )
  
scen1_500pts_pred_summ = data.frame(
	RMSPE = apply(scen1_500pts_RMSPE, 2, mean),
	PI90 = apply(scen1_500pts_PI90, 2, mean),
	CovParmEst = c('all', 'rand', NA, 'comp', NA, 'blnd', NA),
	PredMeth = c('all', 'all', 'nene', 'all', 'nene', 'all', 'nene')
)

scen1_500pts_blk_summ = data.frame(
	RMSPE = sqrt(apply((scen1_500pts_blk_pred[,1:7] - 
	scen1_500pts_blk_pred[,8])^2, 2, mean)),
	PI90 = apply(scen1_500pts_blk_pred[,1:7] - qnorm(.95)*scen1_500pts_blk_se < 
		scen1_500pts_blk_pred[,8] & scen1_500pts_blk_pred[,8] <
		scen1_500pts_blk_pred[,1:7] + qnorm(.95)*scen1_500pts_blk_se, 2, mean),
	CovParmEst = c('all', 'rand', NA, 'comp', NA, 'blnd', NA),
	PredMeth = c('all', 'all', 'nene', 'all', 'nene', 'all', 'nene') 
)

scen1_500pts_est_summ
scen1_500pts_pred_summ 
scen1_500pts_blk_summ

scen2_500pts_est_summ = data.frame(
	RMSE = sqrt(apply(scen2_500pts_est - 1,2,var)),
  CI90 = apply(scen2_500pts_est - qnorm(.95)*scen2_500pts_estse < 1 & 
	1 < scen2_500pts_est + qnorm(.95)*scen2_500pts_estse, 2, mean),
	CovParmEst = c('all', 'rand', NA, 'comp', NA, NA, 'zimm', NA, NA),
	FixEffEst = c('all', 'all', 'rand', 'all', 'rand', 'comp', 'all',
		'rand', 'zimm')
)

scen2_500pts_pred_summ = data.frame(
	RMSPE = apply(scen2_500pts_RMSPE, 2, mean),
	PI90 = apply(scen2_500pts_PI90, 2, mean),
	CovParmEst = c('all', 'rand', NA, 'comp', NA, 'blnd', NA),
	PredMeth = c('all', 'all', 'nene', 'all', 'nene', 'all', 'nene')
)

scen2_500pts_est_summ
scen2_500pts_pred_summ 

scen3_500pts_est_summ = data.frame(
	RMSE = sqrt(apply(scen3_500pts_est - 1,2,var)),
  CI90 = apply(scen3_500pts_est - qnorm(.95)*scen3_500pts_estse < 1 & 
	1 < scen3_500pts_est + qnorm(.95)*scen3_500pts_estse, 2, mean),
	CovParmEst = c('all', 'rand', NA, 'comp', NA, NA, 'zimm', NA, NA),
	FixEffEst = c('all', 'all', 'rand', 'all', 'rand', 'comp', 'all',
		'rand', 'zimm')
)

scen3_500pts_pred_summ = data.frame(
	RMSPE = apply(scen3_500pts_RMSPE, 2, mean),
	PI90 = apply(scen3_500pts_PI90, 2, mean),
	CovParmEst = c('all', 'rand', NA, 'comp', NA, 'blnd', NA),
	PredMeth = c('all', 'all', 'nene', 'all', 'nene', 'all', 'nene')
)

scen3_500pts_est_summ
scen3_500pts_pred_summ 

scen4_500pts_est_summ = data.frame(
	RMSE = sqrt(apply(scen4_500pts_est - 1,2,var)),
  CI90 = apply(scen4_500pts_est - qnorm(.95)*scen4_500pts_estse < 1 & 
	1 < scen4_500pts_est + qnorm(.95)*scen4_500pts_estse, 2, mean),
	CovParmEst = c('all', 'rand', NA, 'comp', NA, NA, 'zimm', NA, NA),
	FixEffEst = c('all', 'all', 'rand', 'all', 'rand', 'comp', 'all',
		'rand', 'zimm')
)

scen4_500pts_pred_summ = data.frame(
	RMSPE = apply(scen4_500pts_RMSPE, 2, mean),
	PI90 = apply(scen4_500pts_PI90, 2, mean),
	CovParmEst = c('all', 'rand', NA, 'comp', NA, 'blnd', NA),
	PredMeth = c('all', 'all', 'nene', 'all', 'nene', 'all', 'nene')
)

scen4_500pts_est_summ
scen4_500pts_pred_summ 

