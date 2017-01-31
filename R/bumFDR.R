### Roxygen based Documentation


#' Fit null and estimate FDR using the BUM model
#'



bumFDR = function(x, statistic=c("normal", "pvalue"), fit=c("smoothing", "native"), a=1e-6, ...)
{
  statistic = match.arg(statistic)
  fit = match.arg(fit)

  ### find null model
 

#  if( fit=="smoothing")
#  {
#     
#     resFDR = log.fdr.fda(x, theo = FALSE, classic = FALSE)
#     eta0 =  resFDR$eta0
#     sd = resFDR$sd
#     pval =resFDR$pval	
#  } 
 # else
  {
     param.out = bum.fitnull(x, statistic=statistic)
     eta0 = param.out$eta0
     sd = param.out$sd

     if(statistic=="normal")
       pval = 2-2*pnorm(abs(x/sd))
     else #statistic="pvalue"
       pval = x

  }

  ### compute FDR values
  Fdr = BUM(1-pval, eta0, a)
  fdr = bum(1-pval, eta0, a)

  return( list(eta0=eta0, sd=sd, fdr=fdr, Fdr=Fdr) )
}

