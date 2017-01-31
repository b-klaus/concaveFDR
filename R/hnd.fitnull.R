
#source("hnd/hnd.curves.R")

# fit null model using HND model

hnd.fitnull = function(x, statistic=c("normal", "pvalue"))
{
  statistic = match.arg(statistic)

  if( statistic == "normal")
  { 
    # null x are normal distributed with mean 0 and sd

    nlogL = function(pp) 
    {
      eta0 = pp[1]
      sd = pp[2]
      y = abs(x/sd)
      k = eta02k(eta0)
      logfdr = hnd(y, k, log=TRUE)
      logf = log(eta0) + dnorm(x, sd=sd, log=TRUE) - logfdr

      return( -sum( logf ) )
    }

    eta0.start=0.5
    sd.start = 1
    #param.out = optim(c(eta0.start, sd.start), nlogL)$par
    param.out = optim(c(eta0.start, sd.start), nlogL, 
       method="L-BFGS-B", lower=c(0.0001,0.0001), upper=c(1, 1000))$par

    return(list(eta0=param.out[1], sd=param.out[2]))

  }

  if( statistic == "pvalue")
  {
     # null x are uniform distributed

    nlogL = function(pp) 
    {
      eta0 = pp[1]
      y = qnorm(1-x/2)
      k = eta02k(eta0)
      logfdr = hnd(y, k, log=TRUE)
      logf = log(eta0) + 0 - logfdr

      return( -sum( logf ) )
    }

    param.out = optimize(nlogL, upper=1, lower=0.0001)$minimum

    return(list(eta0=param.out, sd=NULL))
  }
}

