###    Estimate (Local) False Discovery Rates For Diverse Test Statistics
###    Using log--concave density estimation, main function. It is based
###    on code by Korbinian Strimmer.
###
###
### This file is part of the `concaveFDR' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 3, or at your option, any later version,
### incorporated herein by reference.
###
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
###
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA

### Roxygen based Documentation

#' Estimate (Local) False Discovery Rates For Diverse Test Statistics
#' Using Log-Concave Densities
#'
#' \code{concaveFDR} takes a vector of z-scores (or of correlations, p-values, or
#' t-statistics), and estimates for each case both the tail area-based Fdr as
#' well as the density-based fdr (=q-value resp. local false discovery rate).
#' The parameters of the null distribution are estimated adaptively from the
#' data (except for the case of p-values where this is not necessary).
#'
#' The algorithm implemented in this function proceeds as follows:
#'
#' \enumerate{ \item A suitable cutoff point is determined.  If
#' \code{cutoff.method} is "fndr" then first an approximate null model is
#' fitted and subsequently a cutoff point is sought with false nondiscovery
#' rate as small as possible (see \code{\link{fndr.cutoff}}).  If
#' \code{cutoff.method} is "pct0" then a specified quantile (default value:
#' 0.75) of the data is used as the cutoff point.  If \code{cutoff.method}
#' equals "locfdr" then the heuristic of the "locfdr" package (version 1.1-6)
#' is employed to find the cutoff (z-scores and correlations only).
#' If \code{cutoff.method} is "smoothing" the a smoothed b-spline is
#' fitted to the eta0 estimates at various cutoff points and flat regions
#' are found by inspecting local minima of the smoothed b-spline curve.
#'
#' \item The
#' parameters of the null model are estimated from the data using
#' \code{\link{censored.fit}}. This results in estimates for scale parameters
#' und and proportion of null values (\code{eta0}). Not that choosing theo = TRUE
#' will results in using the theoretical null model without a scale parameter.
#'
#'  \item Subsequently the
#' corresponding p-values are computed, and a modified \code{\link{grenander}} or
#' log-concave density estimation
#' algorithm is employed to obtain the overall density and distribution
#' function (note that this respects the estimated \code{eta0}). The choice
#' of the density estimation algorithm can be made by using the \code{alternative}
#' argument.
#'
#'   \item
#' Finally, q-values and local fdr values are computed for each case.  }
#'
#' The assumed null models all have (except for p-values) one free scale
#' parameter.  Note that the z-scores and the correlations are assumed to have
#' zero mean.
#'
#' @param x vector of the observed test statistics.
#' @param statistic one of "normal" (default), "correlation", "pvalue".  This
#' species the null model.
#' @param plot plot a figure with estimated densities, distribution functions,
#' and (local) false discovery rates.
#' @param verbose print out status messages.
#' @param theo Should the theoretical null model be used? In this case no null
#' model scale parameters (e.g. sigma for z-scores) are estimated. Note that the df and
#' the kappa parameters for t-statistics and the correlation statistics have to
#' be set manually as they depeden on the original sample sizes.
#' @param scale_param scale parameter used when calculating the statistics, only needs to
#' be specified if \code{theo} is set to TRUE and the statistics are t-scores
#' or correlations. See \code{\link{dcor0}} for details on how to set this for correlations.
#' For t-statistics, this depedends on the type of t-test performed.
#' @param alternative The estimation algorithm used for estimation of the alternative
#' densitiy. "Grenander" corresponds to the method implemented in fdrtool, while
#' "log-concave" is new method introduced in the package. At the moment the FDR values
#' are computed using both approaches.
#' @param cutoff.method one of "fndr" (default), "pct0", "locfdr", or smoothing.
#' @param pct0 fraction of data used for fitting null model - only if
#' \code{cutoff.method}="pct0"
#' @param color.figure determines whether a color figure or a black and white
#' figure is produced (defaults to "TRUE", i.e. to color figure).
#' @return A list with the following components:
#' \item{ pval }{a vector with p-values for each case.}
#' \item{ qval}{a vector with q-values (Fdr) for each case..}
#' \item{ lfdr }{ a vector with local fdr values for each case.}
#' \item{ statistic }{the specified type of null model.}
#' \item{ param }{ a vector containing the estimated parameters (the null
#' proportion \code{eta0} and the free parameter of the null model).}
#' @author Bernd Klaus, Korbinian Strimmer (\url{http://strimmerlab.org}).
#' @seealso \code{\link{pval.estimate.eta0}}, \code{\link{censored.fit}}.
#' @references
#'
#' Strimmer, K. (2008a).  A unified approach to false discovery rate
#' estimation. BMC Bioinformatics 9: 303.  Available from
#' \url{http://www.biomedcentral.com/1471-2105/9/303/}.
#'
#' Strimmer, K. (2008b). fdrtool: a versatile R package for estimating local
#' and tail area- based false discovery rates.  Bioinformatics 24: 1461-1462.
#' Available from
#' \url{http://bioinformatics.oxfordjournals.org/cgi/content/abstract/24/12/1461}.
#'
#'
#' @examples
#'
#' # load "fdrtool" library and p-values
#' library("fdrtool")
#' data(pvalues)
#'
#'
#' # estimate fdr and Fdr from p-values
#'
#' data(pvalues)
#' fdr <- concaveFDR(pvalues, statistic="pvalue")
#' fdr$qval.log # estimated Fdr values using log-concave density estimation
#' fdr$qval.gr # estimated local fdr using Grenander density estimation
#'
#' # the same but with black and white figure
#' fdr <- concaveFDR(pvalues, statistic="pvalue", color.figure=FALSE)
#'
#'
#' # estimate fdr and Fdr from z-scores
#'
#' sd.true = 2.232
#' n = 500
#' z = rnorm(n, sd=sd.true)
#' z = c(z, runif(30, 5, 10)) # add some contamination
#' fdr =  concaveFDR(z)
#'
#' # you may change some parameters of the underlying functions
#' fdr =  concaveFDR(z, cutoff.method="pct0", pct0=0.9)
#'
#'
#'
#'
#' @export
#' @import fdrtool


# @importFrom fdrtool  get.nullmodel


concaveFDR = function(x,
  statistic=c("normal", "correlation", "pvalue"),
  plot=TRUE, color.figure=TRUE, verbose=TRUE, theo = FALSE, scale_param = NULL,
  alternative = c("grenander", "log-concave"),
  cutoff.method=c("fndr", "pct0", "locfdr", "smoothing"),
  pct0=0.75)
{
  statistic = match.arg(statistic)
  cutoff.method = match.arg(cutoff.method)

  if ( is.vector(x) == FALSE )
  	stop("input test statistics must be given as a vector!")

  if ( length(x) < 200 ) warning("There may be too few input test statistics for reliable FDR calculations!")
  if (statistic=="pvalue")
  {
    if (max(x) > 1 | min(x) < 0)
      stop("input p-values must all be in the range 0 to 1!")
  }

#### step 1 ####

  if(verbose) cat("Step 1... determine cutoff point\n")

  # determine cutoff point for censoring

  if (cutoff.method=="pct0")
  {
    # use specified quantile

    if(statistic=="pvalue") x0 = quantile(x, probs=1-pct0)
    else x0 = quantile(abs(x), probs=pct0)
  }
  else if ( cutoff.method=="locfdr" & (statistic=="normal" | statistic=="correlation") )
  {
    # use same procedure as in locfdr R package (due to Brit Katzen-Turnbull)

    if(statistic=="normal") z = x
    if(statistic=="correlation") z = atanh(x)

    iqr = as.double(diff(quantile(z, probs=c(.25, .75))))
    sdhat = iqr/(2*qnorm(.75))
    N = length(z)
    # b = 3.55-.44*log(N, 10)                               # locfdr 1.1-3
    b = ifelse(N > 500000, 1,  4.3 * exp(-0.26*log(N,10)) ) # locfdr 1.1-6
    z0 = b*sdhat

    if(statistic=="normal") x0 = z0
    if(statistic=="correlation") x0 = tanh(z0)
  }

  else
  {
    if(cutoff.method=="locfdr")
    warning("cutoff.method=\"locfdr\" only available for normal and correlation statistic.")

    # control false nondiscovery rate
    x0 = fndr.cutoff(x, statistic)

  }


   if ( cutoff.method=="smoothing"   & (statistic=="normal" | statistic=="correlation")  ){
   x0 = smoothing.cutoff(x, statistic)
   }

	if  ( cutoff.method=="smoothing"   & !(statistic=="normal" | statistic=="correlation") ){
    warning("cutoff.method=\"smoothing\" only available for normal and correlation statistic.")
  	}



#### step 2 ####

if(!theo){

  if(verbose) cat("Step 2... estimate parameters of null distribution and eta0\n")

  try({
  cf.out <- censored.fit(x=x, cutoff=x0, statistic=statistic)},
  silent = TRUE)

  if (statistic=="pvalue")
    scale.param = NULL
  else
    scale.param <- cf.out[1,5] # variance parameter

  eta0 = cf.out[1,3]

}

#### step 3 ####

  if(verbose) cat("Step 3... compute p-values and estimate empirical PDF/CDF\n")

  nm <- fdrtool:::get.nullmodel(statistic)

  if(!theo){

  pval <- nm$get.pval(x, scale.param)

  } else {


   if( (statistic %in% c("correlation", "studentt")) & is.null(scale_param) ) {

     stop(paste0("statistic is ", statistic, ", which needs a user specified scale_param for the theoretical null model"))

     }
  scale_param <- switch(statistic,
           "normal" = 1,
           "correlation" = scale_param,
           "pvalue" = NULL,
           "studentt" = scale_param)

  pval <- nm$get.pval(x, scale_param)
  x0 <- quantile(pval, probs=1-pct0)
  try({
    cf.out <- censored.fit(x=pval, cutoff=x0, statistic="pvalue")},
    silent = TRUE)

  eta0 = cf.out[1,3]

  }


  # determine cumulative empirical distribution function (pvalues)
  ee <- fdrtool:::ecdf.pval(pval, eta0=eta0)

  g.pval <- grenander(ee)

  #cat("DEBUG: Grenander eta0=", g.pval$f.knots[length(g.pval$f.knots)], "\n")
  #cat("DEBUG: estimated eta0=", eta0 , "\n\n")

  # mixture density and CDF
  f.pval = approxfun( g.pval$x.knots,  g.pval$f.knots, method="constant", rule=2)
  f0.pval = function(x) return( ifelse(x > 1 | x < 0, 0, rep(1, length(x))) )

  F.pval = approxfun( g.pval$x.knots,  g.pval$F.knots, method="linear",
           yleft=0, yright=g.pval$F.knots[length(g.pval$F.knots)])
  F0.pval = function(x) return( ifelse(x > 1, 1, ifelse(x < 0, 0, x )) )

  #fdr.pval = function(p) pmin( eta0   / f.pval(p), 1) # eta0*f0/ f
  fdr.pval = function(p)
  {
    p[ p == .Machine$double.eps ] = 0
    pmin( eta0   / f.pval(p), 1) # eta0*f0/ f
  }

  Fdr.pval = function(p) pmin( eta0*p / F.pval(p), 1) # eta0*F0/ F


#### step 4 ####

  if(verbose) cat("Step 4... compute q-values and local fdr\n")


  qval.gr <- Fdr.pval(pval)
  lfdr.gr <- fdr.pval(pval)

  resLog = qval.lfdr.pval(pval, eta0=eta0)

  lfdr.log = resLog$lfdr
  qval.log = resLog$qval

  rm(resLog)


   # if()

#### return results ####


  result = list(input.statistics = x, pval=pval,
		qval.gr=qval.gr, qval.log = qval.log,
		lfdr.gr=lfdr.gr,  lfdr.log=lfdr.log,
             statistic=statistic, param=cf.out,
		null.model=nm[c("f0",	"F0","get.pval")],
		alternative = alternative)

  if(theo){

    result$param <- cbind(result$param[, 1:4, drop = FALSE], scale_param)
    }

  if (plot)
  {
    if(verbose) cat("Step 5... prepare for plotting\n")

    ##############
    # zeta > 0 in the following

	  plots <- plotFDR(result)

	  #dev.new(width = 15, height = 8)
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(3, 1)))

    vplayout <- function(x, y){
      viewport(layout.pos.row = x, layout.pos.col = y)
    }
    print(plots$densities, vp = vplayout(1, 1))
    print(plots$CDFs, vp = vplayout(2, 1))
    print(plots$FDRs, vp = vplayout(3, 1))

	result$plots = plots

  }

  if(verbose) cat("\n")

  return(invisible(result))
}


#####

## create labels for plots
pvt.plotlabels <- function(statistic, scale.param, eta0)
{
   if (statistic=="pvalue")
   {
     main = paste("Type of Statistic: p-Value (eta0 = ", round(eta0, 4), ")", sep="")
     xlab ="1-pval"
   }

   if (statistic=="studentt")
   {
     df = scale.param
     main = paste("Type of Statistic: t-Score (df = ", round(df,3),
                       ", eta0 = ", round(eta0, 4), ")", sep="")
     xlab = "abs(t)"
   }

   if (statistic=="normal")
   {
     sd = scale.param
     main = paste("Type of Statistic: z-Score  (sd = ", round(sd,3),
                       ", eta0 = ", round(eta0, 4), ")", sep="")
     xlab = "abs(z)"
   }

   if (statistic=="correlation")
   {
     kappa =scale.param
     main = paste("Type of Statistic: Correlation (kappa = ", round(kappa,1),
                       ", eta0 = ", round(eta0, 4), ")", sep="")
     xlab = "abs(r)"
   }

   return(list(main=main, xlab=xlab))
}

