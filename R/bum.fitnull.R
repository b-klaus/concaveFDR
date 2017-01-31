

# fit null model using bum model

bum.fitnull = function(x, statistic=c("normal", "pvalue"))
{
  statistic = match.arg(statistic)

  if( statistic == "normal")
  { 
    # null x are normal distributed with mean 0 and sd

    nlogL = function(pp) 
    {
      eta0 = pp[1]
      if (eta0 > 1) eta0=1
      if (eta0 < 0.0001) eta0=0.0001
      sd = pp[2]
      if (sd < 0.0001) sd=0.0001
      y = 2*pnorm(abs(x/sd))-1  # 1-pval

      logfdr = log(bum(y, eta0))
      logf = log(eta0) + dnorm(x, sd=sd, log=TRUE) - logfdr
#cat("DEBUG: ",range(logfdr) , "\n")
#cat("DEBUG: ",sum( logf ) , "\n")
#cat("DEBUG: ",eta0 , " ", sd, "\n")

      return( -sum( logf ) )
    }

    eta0.start=0.5
    sd.start = 1
    param.out = optim(c(eta0.start, sd.start), nlogL)$par
    #param.out = optim(c(eta0.start, sd.start), nlogL, 
     # method="L-BFGS-B", lower=c(0.0001,0.0001), upper=c(1, 1000))$par

    if (param.out[1] > 1) param.out[1]=1

    return(list(eta0=param.out[1], sd=param.out[2]))

  }

  if( statistic == "pvalue")
  {
     # null x are uniform distributed

    nlogL = function(pp) 
    {
      eta0 = pp[1]
      y = 1-x
      logfdr = log(bum(y, eta0))
      logf = log(eta0) + 0 - logfdr

      return( -sum( logf ) )
    }

    param.out = optimize(nlogL, upper=1, lower=0.0001)$minimum

    return(list(eta0=param.out, sd=NULL))
  }
}

