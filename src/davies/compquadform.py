import numpy as np
from generalized_chi_squared import davies_method
from typing import List

def davies(
        x: List[float],
        coeff: List[float],
        nc: List[float],
        df: List[int],
        sigma: float = 0,
        lim: int = 10000,
        acc: float = 1e-4) -> float:
    _x = np.asarray(x, dtype=np.float64)
    _coeff = np.asarray(coeff, dtype=np.float64)
    _nc = np.asarray(nc, dtype=np.float64)
    _df = np.asarray(df, dtype=np.int64)

    r = len(_df)
    if len(_x) <= 0:
        raise Exception('x cannot be empty')
    if np.any(_nc < 0):
        raise Exception('Non centrality parameters must be non-negative')
    if np.any(_df < 0):
        raise Exception('df should be non-negative')
    if r != len(_nc):
        raise Exception('df and nc should have the same length')
    if r != len(_coeff):
        raise Exception('df and coeff should have the same length')

    return davies_method(_coeff, _nc, _df, _x, sigma, lim, acc)

# davies <- function(q, lambda, h = rep(1, length(lambda)), delta = rep(0, length(lambda)), sigma = 0, lim = 10000, acc = 0.0001) {

#   r <- length(lambda)
#   if (any(delta < 0)) stop("All non centrality parameters in 'delta' should be positive!") 
#   if (length(h) != r) stop("lambda and h should have the same length!")
#   if (length(delta) != r) stop("lambda and delta should have the same length!")
  
#   out <- .C("qfc", lambdas = as.double(lambda), noncentral = as.double(delta), df = as.integer(h), r = as.integer(r), sigma = as.double(sigma), q = as.double(q), lim = as.integer(lim), acc = as.double(acc), trace = as.double(rep(0, 7)), ifault = as.integer(0), res = as.double(0), PACKAGE = "CompQuadForm")

#   out$res <- 1 - out$res

#   if (out$res > 1) warning("Consider playing with 'lim' or 'acc'.")
  
#   return(list(trace = out$trace, ifault = out$ifault, Qq = out$res))
  
# }