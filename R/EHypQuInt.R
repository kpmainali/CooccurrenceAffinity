#' Largest interval of alpha values
#'
#' This function outputs the largest interval of alpha values within which the Extended Hypergeometric distribution function
#' at x is >= q  and the complementary distribution function 1 - F(x-) is >= 1-q.
#'
#' @details  This function outputs the largest interval of alpha values within which the Extended Hypergeometric distribution
#' function at x is >= q  and the complementary distribution function 1 - F(x-) is >= 1-q.
#' This is the set of alpha values "compatible" with x being a q'th quantile for the Extended Hypergeometric.
#'
#' @param marg a 3-entry integer vector containing (mA,mB,N)
#' @param x nteger co-occurrence counts that should properly fall within the closed interval  [max(0,mA+mB-N), min(mA,mB)]
#' @param q a quantile falling strictly between 0 and 1
#' @param scal an integer parameter (default 10) that should fall somewhere between 2 and 10
#'
#' @return This function returns the largest interval of alpha values
#'
#' @author Eric Slud
#'
#' @references to be added


EHypQuInt <-
function(marg, x, q, scal=log(2*marg[3]^2)) {
  #  marg = c(mA, mB, N),
  # x an observed co-occurrence count, q the desired quantile
  require(BiasedUrn)
  mA=marg[1]; mB=marg[2]; N=marg[3]
  xmin = max(mA+mB-N,0);  xmax = min(mA,mB)
  maxabs = log(2*N^2)
  if(length(intersect(c(mA,mB), c(0,N))))
    return("Degenerate co-occurrence distribution!")
  null.mean = mA*mB/N
  newint = if(x >= null.mean) c(-1,scal) else c(-scal,1)

  # --------------
  # The function EHypCent() is a utility forming part of the EHypQuInt function, and would not be called as a stand-alone.
  # Its output is a vector of centered Extended-Hypergeometric distribution function values with the quantile q subtracted.
  # The function is used in root-finding to locate alpha-value intervals for which the corresponding
  # Extended Hypergeometric distribution assigns distribution function close to q at t.

  EHypCent = function(alp,t) {
    # t an observed co-occurrence count, alp a sequence of alphas
    L = length(alp)
    out = numeric(L)
    for(i in 1:L) out[i]= pFNCHypergeo(t,mA,N-mA,mB,exp(alp[i]))-q
    out
  }
  #-------------

  if(x < xmin) c(-Inf,-Inf) else {
    if(x == xmin) c(-maxabs,
      uniroot(EHypCent, newint, t=x, extendInt="y")$root) else {
        if(x > xmax) c(Inf,Inf) else {
          if(x == xmax) c(uniroot(EHypCent,
            newint, t=x-1, extendInt="yes")$root,maxabs) else
              c(uniroot(EHypCent, newint, extendInt="y", t=x-1)$root,
                uniroot(EHypCent, newint, extendInt="y", t=x)$root)
        }
      }
  }
}
