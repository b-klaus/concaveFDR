### Roxygen based Documentation


#' Simulate z-scores according to a mixture model
#' 
#'
#'\deqn{\eta_0 * N(0,\sigma^2) + (1-\eta_0/2) * Unif(alt.min,alt.max)
#'		+ (1-\eta_0/2)*Unif(-alt.max, alt.min)}
#'
#' also returns functions of the simulated null and alternative densities and
#'CDFs as well as FDR and fdr functions 



#' @export









#Data $z_1, \ldots, z_{200}$ were
#drawn from a  mixture of the normal distribution $N(\mu = 0, \sigma^2=4)$
#with the symmetric uniform alternatives $Unif}(-10, -5)$ and $Unif}(5, 10)$ and a
#null proportion of $\eta_0=0.8$, 
#i.e. $0.8N(0,4) + 0.1Unif}(-5,10) + 0.1Unif}(5,10)$. 





get.random.zscore = function(d= 200, sigma =2, eta0 = 0.8, 
				alt.min = 5, alt.max = 10)
{
  
  round(d0 <- d*eta0)
  d1 <- d-d0
  z <- c( pmin(9.999, abs(rnorm(d0, mean=0, sd=sigma))),
       runif(d1, min=alt.min, max=alt.max))
  z <- sign(rnorm(length(z)))*z

  	f0 = function(x) 2*dnorm(x, sd=sigma) 
	F0 = function(x) 2*pnorm(x, sd=sigma)-1

	fA = function(x) dunif(x, min=alt.min, max=alt.max)
	FA = function(x) punif(x, min=alt.min, max=alt.max)


	f = function(x) eta0*f0(x)+(1-eta0)*fA(x)
	F = function(x) eta0*F0(x)+(1-eta0)*FA(x)


	Fdr.func = function(x) eta0*(1-F0(abs(x)))/(1-F(abs(x)))
	fdr.func = function(x) eta0*f0(abs(x))/ f(abs(x))

  return(list(z = z, f0 = f0, F0 = F0, 
		fA = fA, FA = FA, 
	   	f = f, F = F, 
		Fdr.func = Fdr.func, fdr.func = fdr.func))
}


#test <- get.random.zscore()
