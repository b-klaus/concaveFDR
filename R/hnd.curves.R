

# FDR curves half-normal decay model

###########

# relationship between k and eta0

k2eta0 = function(k) 1/(2*pnorm(k)-1 + sqrt(2/pi)*exp(-k^2/2-log(k)))

eta02k = function(eta0, maxk=1e3)
{
  if(eta0 >= 1) return(Inf)

  findk = function(k) k2eta0(k)-eta0
  kopt = uniroot(findk, c(0,maxk))$root

  return(kopt)
}

###########

# fdr value
hnd = function(y, k, log=FALSE)
{
  if( log==FALSE)  
    ifelse(y < k, 1, exp(-(y-k)^2/2))
  else
   ifelse(y < k, 0, -(y-k)^2/2)
}

###########

# marginal density and distribution

f.hnd = function(y, k) 
{
  e0 = k2eta0(k)
  
  ifelse(y < k, 
      e0*sqrt(2/pi)*exp(-y^2/2),
      e0*sqrt(2/pi)*exp(k^2/2-y*k)
     )
}

F.hnd = function(y, k)
{
  e0 = k2eta0(k)

  e0*ifelse(y < k, 
      2*pnorm(y)-1, 
      2*pnorm(k)-1 + sqrt(2/pi)*exp(k^2/2-log(k))*(exp(-k^2)-exp(-k*y)) 
   )
}

# null density and null distribution

f0.hnd = function(y)
{
  #2*dnorm(y)
  sqrt(2/pi)*exp(-y^2/2)
}

F0.hnd = function(y)
{
  2*pnorm(y)-1
}

####

# Fdr value
HND = function(y, k)
{
  e0 = k2eta0(k)
  Fdr = e0*(1-F0.hnd(y))/(1-F.hnd(y, k))
  return(Fdr)
}




