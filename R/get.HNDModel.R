#'Draws statistics from the half normal dacy (HND) model 
#'
#'
#' 
#'
#'



#Fdr.HNDscore = function(x) HND(abs(x),k) 
#fdr.HNDscore = function(x) hnd(abs(x),k)



get.random.HNDscore = function(m=1000, eta0=0.8)
{

  k = eta02k(0.8) 	
  HNDs = runif(round(m/2),0,1)
  score = c(-F.hndInv(HNDs,k ), F.hndInv(HNDs,k ))

  return(score)
}


F.hndInv = function(pval, k)
{
 

  findroot <- function(pv){
  	tmpfun = function(y)  F.hnd(y,k)-pv
  	uniroot(tmpfun, interval = c(0, 100))$root
  }
  
vapply(pval, FUN = findroot, FUN.VALUE = numeric(1)) 


}


##Test
#m = 1000
#HNDs = runif(m/2,0,1)
#test = c(-F.hndInv(HNDs,k ), F.hndInv(HNDs,k ))
#ylim =c(0,0.5) 
#hist(test, freq = FALSE)
#lines(abs(test), 0.5*f.hnd(abs(test),k), col = "darkgreen", type = "p")

#hist(F.hndInv(HNDs,k ), freq = FALSE, ylim =c(0,0.8) )
#lines(F.hndInv(HNDs,k ), f.hnd(F.hndInv(HNDs,k ),k), col = "darkgreen", type = "p")

#test = get.random.HNDscore(1000)
