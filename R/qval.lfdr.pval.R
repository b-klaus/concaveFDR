## 21.5.2012

#library("logcondens")
#library("R.utils")


### Roxygen based Documentation


#' Estimate the density and distribution function estimation using logcondens
#' on  p-values transformed to z-values
#'constrained such that the known fraction eta0 of null p-values 
#' is taken into account

 	
#' @importFrom logcondens logConDens LocalExtend LocalNormalize LocalF
#' @importFrom R.utils insert


## x are "correct" pvalues (i.e. computed under the right null model) 
qval.lfdr.pval <- function (pv, eta0=1) 
{
      Len = length(pv)

   if (Len != length(unique(pv)) ){
   ### positions of pv-value which are duplicated elsewhere in the 
    ## vector == pvvalues to be removed from the vector
    idxRem = which(duplicated(pv) == TRUE)
   ## index of originally duplicated values which remain in the vector
    idxDup= unique(match(pv[idxRem],pv))
   ### DEBUG: test index sets, recover pvalues belonging to  indexes cut 
   ## recover.idx = idxRem[pv[idxDup[1]] == pv[idxRem]]	
   ### recover indices by each duplicated value separately
   recover.idx= numeric(0)
   for (i in 1:length(idxDup)){
	recover.idx[i] = length(idxRem[pv[idxDup[i]] == pv[idxRem]])
   }
   # remove duplicates
   pv = pv[-idxRem]   
   } 


    ##sort the pvalues	and remove the duplicates	 
    vals = sort(pv)
     
    ##create "artificial" z-Values by mirroring
    zVals = c(-p2z(vals),p2z(vals))
    #DEBUG: duplicated(zVals)  
    #DEBUG: hist(zVals) 	
    #DEBUG: vals2 = sort(z2p(zVals))
    #DEBUG: vals-vals2
    #DEBUG: vals[which(zVals == -Inf)]
    
    #compute log concave density etimator
    res <- logConDens(zVals, smoothed = FALSE, print = FALSE, xs = zVals)


    ##auxillary values for bound corrections

    n <- length(vals)
    x <- zVals[1:length(vals)]
    ## is a zvalue a knot of the log concave estimator?
    ## knots =  "changepoints" of the estimator
    IsKnot <- res$IsKnot[1:length(vals)]
    ## log concave estimator =  exp(phi)
    phi <- res$phi[1:length(vals)]	
    x.knot <- res$knots[1:length(vals)]
    phi.knot <- phi[IsKnot > 0]
    js <- 2:n
    xj <- x[js - 1]
    xj1 <- x[js] 
  
    		 
    target = 1 ##start value of the norm	
    steps = 0 ## factor for modifiying phi
    phiMOD = numeric()
    while((target > 0.05) & (steps < 1)){
	##DEBUG  browser()	
	steps  <-  steps + 0.01
	##Additional values for phi and x (last Data point has to be
	##a Knot => but  the phi values behind that point don't look 
	##trustworthy so they are interpolated!)
	IsKnot <- res$IsKnot[1:length(vals)]	
	K <- (1:n) * IsKnot
    	K <- K[K > 0]
		
		if(IsKnot[length(x)]!=1){	
		K =c(K,length(x)) 
		IsKnot[length(x)]=1
		}

	## adapt phi   
	phi2 = phi[K]*(1 +steps)
	## phi2[length(K)] = 0
	## phi2 = c(phi2, phi2[length(K)])
	
	##DEBUG recompute phi (if phi2 = knots of phi)
	##DEBUG phiTEST = LocalExtend(x, IsKnot, x.knot, phi2) 
	##DEBUG: plot( phi-phiMOD)
 	phiMOD = LocalExtend(x, IsKnot, x.knot, phi2) 
	phiMOD = LocalNormalize(x,phiMOD) 
	### adjust normalising factor to avoid postive phi-values
		if ( phiMOD[length(x)] > 0 ){	
		norm = phiMOD[length(x)]  
		phiMOD = phiMOD - phiMOD[length(x)]		
		}		
	## compute F
	F.raw = LocalF(x, phiMOD)
	
 	##DEBUG plot(vals, F.raw)
	##DEBUG lines(vals, 1-eta0*(1-vals))
   	 ##DEBUG lines(vals, eta0*vals)
	##How many points violate upper bound?	
         s1 <- sum(F.raw > 1-eta0*(1-vals))
	##How many points violate lower bound?
	  s2 <- sum(F.raw <  eta0*vals)	
	## deviations from the boundaries
	dev <- c((F.raw -(  eta0*vals))[which(F.raw -(  eta0*vals) < 0) ],
	(F.raw -(  1-eta0*(1-vals)))[which(F.raw -(  1-eta0*(1-vals)) > 0) ])
	## compute the infintiy norm of the deviations	
	target =  norm(as.matrix(dev), type ="F")
 
	}

	##DEBUG print("end of loop reached")  	

	##trim F (final corrections)
	##upper bound
	tmp = pmin(LocalF(x, phiMOD), (1-eta0*(1-vals)))
	## lower bound	
	F= pmax(tmp ,eta0*vals)
	f = exp(phiMOD)

	##DEBUG plot(vals, F)
	##DEBUG lines(vals, 1-eta0*(1-vals))
   	##DEBUG lines(vals, eta0*vals)

	##compute raw qval	
	qval <-(eta0*vals) / F

	##DEBUG plot(vals, qval)

	##linear interpolate qval[ra]on of the smallest qvalues 
	min1 <- which.min(qval[1:round((1-eta0)*length(qval))])
	
		if(min1 != 1){
	 	qval[1:min1] = approx(x = c(x[1],x[min1]), y = c(0,qval[min1]), 
		xout = x[1:min1] )$y
		}
		
	
	##get index non monotone qvalues
	idxM = which(diff(qval) < 0)
		if(length(idxM) > 0){
		qval[idxM[1] : length(qval)] = qval[idxM[1]]
		}
	## make sure the maximum qval is equal to eta0
	idxEta = which( qval > eta0)
	if(length(idxEta) > 0){
	qval[idxEta[1] : length(qval)] = qval[idxEta[1]]
	}

	

	## compute lfdr and trim values greater than 1
 	lfdr <- pmin( eta0*(dnorm(x)*2/ f),1)
	##DEBUG plot(x, f)
	 ##DEBUG  	hist(x, freq = FALSE);points(x, f)
	##DEBUG plot(x, eta0*(dnorm(x)*2/ f)
	##DEBUG plot(x,f)
	##	
	## makes sure lfdr stays 1 after it has reached 1
	## if it does not reach 1, set index to the length of
	## the fdr vector 
		
		if(sum(lfdr==1) != 0){
		idx1 <- min(length(lfdr),min(which(lfdr==1), na.rm = TRUE))
		}else{
		idx1 <- length(lfdr)		
		}
	
		if(sum(lfdr==1) != 0){
		
		lfdr[idx1:length(lfdr)] = 1
    		}
	### make sure ldfr ist monotone (the inverse 
	# usually happens at this stage if the null model is misspecified!)
	idx.lfdr.e = which(diff(lfdr)<0)	
		if ( length(idx.lfdr.e) && (eta0 < 0.9) ){
	#warning("Non-monotone fdr values
	# computed, null model probably misspecified, try empirical null")
	### make sure ldfr ist monotone by approximation!
	idx.e = max(idx.lfdr.e[1]-1,length(lfdr))
	
		if(idx.e != idx1 ){
		lfdr[idx.e:idx1] = approx(x = c(x[idx.e],x[idx1]), y = c(lfdr[idx.e],1),xout = x[idx.e:idx1] )$y}
		}	
	
	
	##DEBUG plot(vals, lfdr)	

	## bring the results into an order that corresponds to the orginial 
	## ordering
	ra <- rank(pv,ties.method =  "min")	
	

	### DEBUG ra[1:100]
	### DEBUG ra2[1:100]
    	
		
		qval = qval[ra]
		lfdr = lfdr[ra]
		
	##DEBUG plot(pv, lfdr)	
	##DEBUG plot(pv, qval)	
	

	
		if (Len != length(unique(pv)) ){
		### indeces of repeated values to insert depending on 
		### the correseponding index of idxRem
		repVal.idx = rep(idxDup, recover.idx)	
		### order of adding is important to find the right indices!	
		for (i in 1:length(idxRem)){	
		### fill in missing FDR values
		lfdr = insert(x = lfdr, ats=idxRem[i], values = lfdr[repVal.idx[i]])
		qval = insert(x = qval, ats=idxRem[i], values = qval[repVal.idx[i]])
		}}

	##DEBUG pvO = resFDR$pval
	##DEBUG plot(pval.log, lfdr)	
	##DEBUG plot(pval.log, qval)	


   	rval <- list(qval = qval, lfdr = lfdr)
	return(rval)	
	}

	

