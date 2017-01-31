
##  @import splines

###   @importFrom fdrtool censored.fit

#' ### Roxygen based Documentation
#' estimate tail area based FDR and local fdr  via log-concave densities
#' @importFrom pastecs turnpoints extract.turnpoints
#' @importFrom fda create.bspline.basis fdPar smooth.basis eval.fd
#'
#' @export







### find a cutoff by a smoothing technique and return null model parameters
smoothing.cutoff <- function(x, statistic=c("normal", "correlation")){
require(splines)

    if(statistic=="normal") z = x
    if(statistic=="correlation") z = atanh(x)



## DEBUG z = sample(test.stats.cov, sub.sample)
	### increase variance of z to allow stable spline modelling
	sd.mod.fac <- IQR(z[abs(z) <  quantile(z, .75) ]) / 1.349 /4


	z = z / sd.mod.fac
	## DEBUG sd(z)

	try({
 	fit <-censored.fit(z, statistic="normal", cutoff=quantile(abs(z),
	probs=30:99/100)) [ ,c(1,2,3,5)]
	}, silent = TRUE)

	### remove fits with same number of censored observations
	fit = fit[  !duplicated(fit[, "N.cens"]	) ,]


####### eta0  // sd ##########################

			## idxs of very small eta0 values
		which( fit[,3] <= quantile(fit[,3], probs = 0.10))
		## idxs of very large eta0 values
		which( fit[,3] >= quantile(fit[,3], probs = 0.90))

		## union of both
		rm.idx.eta0 = union(which( fit[,3] < quantile(fit[,3], probs = 0.1) ),
		which( fit[,3] >= quantile(fit[,3], probs = 0.9) ))

	if( (quantile(fit[,3], probs = 0.10) < 1)  && ( length(rm.idx.eta0) < 60 ) ){
	fit = fit[- rm.idx.eta0,]
	}

	## plot eta0
	#plot(fit[,1], fit[,3])
	#range(fit[,1])
	eta0.basis = create.bspline.basis( breaks =unique(fit[,1]), norder = 4)
	#argumjents: basis, derivative to penalize, lambda
	eta0fdPar = fdPar(eta0.basis, 2,0.01)
	eta0fd = smooth.basis(fit[,1], fit[,3], eta0fdPar)
	#X11()
	# plot(eta0fd, ylab = expression(eta[0]), xlab = expression(y[c]), cex = 2)
	#X11()
	### alternative spline fun
	#smooth.eta0.func <- splinefun( x = fit[seq(10,70, by = 10),1], y = fit[seq(10,70, by = 10),3] )
	#plot(fit[,1], smooth.eta0.func(fit[,1], deriv = 0), type = "l")

	### first deriv
	 eta0.1 <- function(x){ abs(eval.fd(x, eta0fd$fd,1))}
	 #plot(fit[,1],  eta0.1(fit[,1]), type = "l")
	### second deriv
	 #eta0.2 <- function(x){ abs(eval.fd(x, eta0fd$fd,2))}
	 #plot(fit[,1],  eta0.2(fit[,1]), type = "l")
	### second deriv

	#### test turnpoints

	test.tp = turnpoints( as.vector(  eta0.1(fit[,1]) ) )

	## peaks
	##which(extract.turnpoints(test.tp ) > 0)

	## pits
	tp.pits = extract.turnpoints(test.tp, , level = 0.05 ) < 0
	if(sum(tp.pits) != 0) {
	idx.opt = round(median(which(tp.pits)))
	} else {
	idx.opt = which.min (as.vector(  eta0.1(fit[,1]) ))
	}

	###########





		#eta0.opt = optimize(eta0.1, lower = 2, upper = 4 )$minimum




		#eta0.opt = optimize(eta0.1, lower = min(fit[,1]),
		#upper = max(fit[,1]) )$minimum
		#print("opt-case")
		#eta0 =  censored.fit(z, statistic="normal", cutoff=eta0.opt)[3]
 		#sd = censored.fit(z, statistic="normal", cutoff=eta0.opt)[5]

#		### undo variance inflation
#		z = z*sd.mod.fac

#		eta0 = fit[idx.opt ,3]

# 		sd =  fit[idx.opt ,4] * sd.mod.fac
#
#		## fit[idx.opt ,4] * sd.mod.fac
#		### fit[idx.opt ,4] * sd.mod.fac
#

#
#
#
#
#	### BUG fix if empirical null is grossly wrong!!
#	## (eta0 < 0.65) || (sd < 0.8)
# 	if( FALSE ){
#	pval.log = 2- 2*pnorm(abs(z), sd = 1)
#	sd = 1
#	eta0 = 0.8
#	}else{
#  	 pval.log = 2- 2*pnorm(abs(z), sd = sd)
#	}
#	print(paste("eta0: ", round(eta0, 4), "sd", round(sd,4)))

 	return(fit[idx.opt, "cutoff"] *  sd.mod.fac)
}
