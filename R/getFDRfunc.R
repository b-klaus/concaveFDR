#'Get the FDR functions associated with an FDR analysis on the abs(stat)  scale
#'
#'
#' returns them based on an FDR analysis
#'
#'


#'@export


##############
    # zeta > 0 in the following, gives quantities on absolute value scale

getFDRfunc <- function(result){

	eta0 <- result$param[, "eta0"]
	result$pval[ result$pval == .Machine$double.eps ] = 0
	## quantile(result$pval, 0.9999)## avoid very large values

	if( result$statistic=="pvalue" )
    {
      f0 <- function(zeta) return( result$null.model$f0(zeta, result$param[, "sd"]) ) 
      F0 <- function(zeta) return( result$null.model$F0(zeta, result$param[, "sd"]) )
      get.pval <- function(zeta) return( result$null.model$get.pval(1-zeta, result$param[, "sd"]) )
      x0 = result$param[, "cutoff"]
    } 
    else
    {
      f0 <- function(zeta) return( 2*result$null.model$f0(zeta, result$param[, 5])  ) 
      F0 <- function(zeta) return( 2*result$null.model$F0(zeta, result$param[, 5])-1  )
      get.pval <- function(zeta) return( result$null.model$get.pval(zeta, result$param[, 5]) )
    }



    lfdr.log.p = approxfun(result$pval,  result$lfdr.log, method="constant", rule=1)
	lfdr.log = function(zeta)   lfdr.log.p(get.pval(zeta)) 

    lfdr.gr.p = approxfun(result$pval,  result$lfdr.gr, method="constant", rule=1)	       
	lfdr.gr = function(zeta)   lfdr.gr.p(get.pval(zeta)) 		


	qval.log.p = approxfun(result$pval,  result$qval.log, method="linear", rule=1)
	qval.log = function(zeta)   qval.log.p(get.pval(zeta)) 

    qval.gr.p = approxfun(result$pval,  result$qval.gr, method="linear", rule=1)	       
	qval.gr = function(zeta)   qval.gr.p(get.pval(zeta)) 


	f.log = function(zeta) eta0*(f0(zeta))/lfdr.log(zeta) 	
	F.log = function(zeta) 1-eta0*get.pval(zeta)/qval.log(zeta) 


	f.gr = function(zeta) eta0*(f0(zeta))/lfdr.gr(zeta) 
	F.gr = function(zeta) 1-eta0*get.pval(zeta)/qval.gr(zeta)


	FA.gr = function(zeta)  (F.gr(zeta)-eta0*F0(zeta))/(1-eta0)	
	fA.gr = function(zeta)  (f.gr(zeta)-eta0*f0(zeta))/(1-eta0)

	FA.log = function(zeta)  (F.log(zeta)-eta0*F0(zeta))/(1-eta0)	
	fA.log = function(zeta)  (f.log(zeta)-eta0*f0(zeta))/(1-eta0)


#        tp <- sort(abs(z))
#	plot(tp, FA.gr(tp))

#	plot(tp, FA.log(tp))

#	plot(tp, fA.gr(tp), type = "l")
#	plot(tp, fA.log(tp), type = "l")

#	plot(tp, F.gr(tp))
#	X11()
#	plot(tp, eta0*F0(tp))
#	X11()
#	plot(tp, F.gr(tp) - eta0*F0(tp))

#	plot(tp, F0(tp))

#	plot(tp, f.log(tp))
#	plot(tp, F.log(tp))
#	plot(tp, qval.gr(tp))


	

     return(list(f0 = f0, F0 = F0, 
		fA.gr = fA.gr, fA.log = fA.log,
		FA.gr = FA.gr, FA.log = FA.log,
	   	f.log = f.log, f.gr = f.gr, 
		F.log = F.log, F.gr = F.gr, 
		qval.gr=qval.gr, qval.log = qval.log,
		lfdr.gr=lfdr.gr,  lfdr.log=lfdr.log
		))

}
