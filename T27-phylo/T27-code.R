# Chi-squared value calculation for dN/dS analysis in Supplementary Table S4
dt <- data.table(
   segment = 1:8,
   lnL_M7 = c(-382.400741, -267.572990, -296.738780, -1541.406639, -429.560444, -661.580273, -2422.542513, -558.233586),
   lnL_M8 = c(-381.058233, -267.121865, -296.427578, -1521.268204, -427.757214, -657.751507, -2421.319128, -558.149152)
   )
dt[, "lambda" := 2*(-lnL_M7+lnL_M8)] # log-likelihood-ratio statistic
#dt[, "chi2pval" := sprintf("%f",pchisq(lambda, df=2, lower.tail=FALSE))] 
dt[, "chi2pval" := pchisq(lambda, df=2, lower.tail=FALSE)] 
