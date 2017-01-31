#'Plot the results of an FDR analysis
#'
#'Plots the results of an FDR analysis using ggplot2 
#' and returns FDR functions
#'
#'

#'@import ggplot2 grid
#'@export





 plotFDR <- function(result){

	
    x <-  result$input.statistics
    funcs <-getFDRfunc(result)
    ##############

    ax = abs(x) 
    if (result$statistic=="pvalue") ax = 1-ax  # reverse p-val plot 
   # xxx = seq(0, max(ax), length.out=500)
    #xxx <- c(xxx, rep(NA, max(0, length(ax) - 500)))
  

	## avoid very large values by thresholding the densities and cdfs
	thresh.plot <- max(hist(ax, breaks = 50, plot=F)$density)

    eta0 <- result$param[, "eta0"]
    dataGG <- data.frame(statistic = ax, 
		f0 = eta0*funcs$f0(ax),
		fA.gr = pmin((1-eta0)*funcs$fA.gr(ax),thresh.plot),
		fA.log = pmin((1-eta0)*funcs$fA.log(ax),thresh.plot),
		f.log = pmin(funcs$f.log(ax), thresh.plot),
		f.gr = pmin(funcs$f.gr(ax), thresh.plot),
		F.log = pmax(funcs$F.log(ax),0),
		F.gr = pmax(funcs$F.gr(ax),0),
		FA.log= pmax((1-eta0)*funcs$FA.log(ax),0),
		FA.gr = pmax((1-eta0)*funcs$FA.gr(ax), 0),
		F0= eta0*funcs$F0(ax),
		qval.gr=funcs$qval.gr(ax), 
		qval.log = funcs$qval.log(ax),
		lfdr.gr=funcs$lfdr.gr(ax),  
		lfdr.log=funcs$lfdr.log(ax)
		#plottingColors = c(plottingColors, 
	#		rep(NA, max(0, length(ax) - length(plottingColors)))) 
	)    
	dataGG$F.log[ax == max(ax) ] = 1 ## remove approximation artifact
	dataGG$FA.log[ax == max(ax) ] = (1-eta0) ## remove approximation artifact

    ll = fdrtool:::pvt.plotlabels(result$statistic, result$param[, 5], eta0)

	### get color values for the  densities and distributions to be plotted
	### to make the ggplot2 object self contained
#	values  <- eval(parse(text = paste0("values = c(",names(plottingColors)[1],  "=", "\"",plottingColors[1],"\"",
#					",", names(plottingColors)[2],  "=", "\"",plottingColors[2],"\"",
#						 ",", names(plottingColors)[3],  "=", "\"",plottingColors[3],"\"",
#						 ",", names(plottingColors)[4],  "=", "\"",plottingColors[4],"\"",
#						 ",", names(plottingColors)[5],  "=", "\"",plottingColors[5],"\"",				
#						 ")")))
# 
# text.tmp = paste0(" c(","\"",names(plottingColors)[1], "\"", "=", "\"",plottingColors[1],"\"",
# 					",", "\"",names(plottingColors)[2], "\"", "=", "\"",plottingColors[2],"\"",
# 						 ",", "\"",names(plottingColors)[3], "\"", "=", "\"",plottingColors[3],"\"",
# 						 ",", "\"",names(plottingColors)[4], "\"", "=", "\"",plottingColors[4],"\"",
# 						 ",", "\"",names(plottingColors)[5], "\"", "=", "\"",plottingColors[5],"\"",				
# 						 ")")
# parse(text = text.tmp)


# text.col <- paste0(names(plottingColors))
################### histogram / density plot #####################################

    den.plot<- ggplot(dataGG, aes(x=statistic))
    den.plot<- den.plot+ geom_histogram(aes(y = ..density..),
		binwidth = round(diff(range(ax)))/50, fill = "white", color = "black")

    den.plot<- den.plot+  xlab(ll$xlab) + ggtitle(paste0("Density estimates\n",ll$main))	
    den.plot <- den.plot+ geom_line(aes(y=f0, colour = "Null"), 
	linetype = 2, size =1.2)
	den.plot <- den.plot + geom_line(aes(y=fA.gr, colour = "Alternative grenander"), 
	 linetype = 1, size =1.2)
    den.plot <- den.plot + geom_line(aes(y=fA.log, colour = "Alternative log concave"), 
	 linetype = 1, size =1.2)
	den.plot <- den.plot  + scale_colour_manual(name = 'Color code', 
         values =  c("Null"="coral3","Alternative grenander"="#7b3294",
			"Alternative log concave"="#008837","Mixture grenander"="#c2a5cf",
		"Mixture log concave"="#a6dba0") )

#pdf("test.pdf")
	 #den.plot
#dev.off()
# save(den.plot, file = "test.RData")
################### CDF plot ############################################



     cdf.plot <- ggplot(dataGG, aes(x=statistic))
	 cdf.plot <- cdf.plot +  geom_line(aes(y=F0, colour = "Null"), 
	linetype = 2, size =1.2)

cdf.plot <- cdf.plot + xlab(ll$xlab) +  ggtitle(paste0("Distribution Functions\n",ll$main))
	cdf.plot <- cdf.plot + ylab("Cummulative distribution function")

	cdf.plot  <- cdf.plot  + geom_line(aes(y=FA.gr, colour = "Alternative grenander"), 
	 linetype = 1, size =1.2)
    cdf.plot  <- cdf.plot  + geom_line(aes(y=FA.log, colour = "Alternative log concave"), 
	 linetype = 1, size =1.2)
	 
		   cdf.plot  <- cdf.plot  + geom_line(aes(y=F.gr, colour = "Mixture grenander"), 
	 linetype = 1, size =1.2)

    cdf.plot  <- cdf.plot  + geom_line(aes(y=F.log, colour = "Mixture log concave"), 
	 linetype = 1, size =1.2)


	cdf.plot  <- cdf.plot + scale_colour_manual(name = 'Color code', 
         values =  c("Null"="coral3",
                     "Alternative grenander"="#7b3294",
                    "Alternative log concave"="#008837",
                    "Mixture grenander"="#c2a5cf",
                     "Mixture log concave"="#a6dba0"))
#pdf("test.pdf")
#cdf.plot 
#dev.off()


################### FDR / lfdr plot ############################################


	  FDR.plot <- ggplot(dataGG, aes(x=statistic))
	 FDR.plot <- FDR.plot +  geom_line(aes(y=qval.gr, colour = "Fdr grenander" ), 
	linetype = 1, size =1.2) 
  	 FDR.plot <- FDR.plot + xlab(ll$xlab) + ylab("(Local) False Discovery Rate")
	 FDR.plot <- FDR.plot +  ggtitle(paste0("Local and tail area based false discovery rates\n",ll$main)) 
	FDR.plot  <- FDR.plot  + geom_line(aes(y=qval.log, colour = "Fdr log concave"), 
	 linetype = 1, size =1.2)
    FDR.plot  <- FDR.plot  + geom_line(aes(y=lfdr.gr, colour = "local fdr grenander"), 
	 linetype = 1, size =1.2)
	FDR.plot  <- FDR.plot  + geom_line(aes(y=lfdr.log, colour = "local fdr log concave"), 
	 linetype = 1, size =1.2)
	FDR.plot  <- FDR.plot + scale_colour_manual(name = 'Color code', 
	             values =  c("Fdr grenander"="#7b3294",
	                         "Fdr log concave"="#008837",
	                         "local fdr grenander"="#c2a5cf",
	                         "local fdr log concave"="#a6dba0"))

# 
# pdf("test.pdf")
# FDR.plot
# dev.off() 


################### Combine plots in a single one ##################################


	
    ret <- list(densities = den.plot, 
				CDFs = cdf.plot, 
				FDRs = FDR.plot, 
				plottingColors = plottingColors)

  	return(invisible(ret))

}


#####


  
