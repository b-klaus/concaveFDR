   ## conversion of p to  two sided z-values
   ### cut off at 20 to avoid computational problems
   p2z = function(pv) pmax(pmin(qnorm(1-pv/2),20),-20) 
 
   # #conversion of two sided z to p-values
   z2p = function(z) 2- 2*pnorm( pmax(pmin(abs(z),20),-20) )
