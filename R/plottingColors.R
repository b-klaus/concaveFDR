#'Named vectors of colors used to plot the FDRs, densities etc.
#'
#'Colors for plotting
#' 
#'
#'@export



	plottingColors <- c("coral3", "#7b3294","#008837", "#c2a5cf","#a6dba0",
				"#7b3294","#008837", "#c2a5cf","#a6dba0" ) 
    names(plottingColors) = c("Null", "Alternative grenander", 
		"Alternative log concave","Mixture grenander", 
		 "Mixture log concave", "Fdr grenander", "Fdr log concave", 
		"local fdr grenander", "local fdr log concave")
