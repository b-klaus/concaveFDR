
#' Fit null and estimate FDR using the HND model
#' @export




run.hnd = function(x, statistic=c("normal", "pvalue"), fit=c("logfdr", "native"), ...)
{
  statistic = match.arg(statistic)
  fit = match.arg(fit)

  ### find null model

  if( fit=="logfdr")
  {
    resFDR = log.fdr.fda(x, theo = FALSE, classic = FALSE)
     eta0 =  resFDR$eta0
     sd = resFDR$sd
     pval =resFDR$pval

  }
  else
  {
     param.out = hnd.fitnull(x, statistic=statistic)
     eta0 = param.out$eta0
     sd = param.out$sd

     if(statistic=="normal")
       pval = 2-2*pnorm(abs(x/sd))
     else #statistic="pvalue"
       pval = x

  }

cat("DEBUG: eta0=", eta0, "\n")
cat("DEBUG: sd=", sd, "\n")


  ### compute FDR values
  az = qnorm(1-pval/2)
  k = eta02k(eta0)

  Fdr = HND(az, k)
  fdr = hnd(az, k)

  return( list(eta0=eta0, sd=sd, fdr=fdr, Fdr=Fdr) )
}

