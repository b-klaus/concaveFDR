

# FDR curves BUM model


# y = 1-pval

# fdr value
bum = function(y, lambda, a=1e-3)
{
  out = lambda/(lambda + a*(1-lambda)*(1-y)^(a-1))

  return( ifelse (out < 0.00000001, 0.00000001, out) )
}

# Fdr value
BUM = function(y, lambda, a=1e-3)
{
   lambda/(lambda+(1-lambda)*(1-y)^(a-1))
}


