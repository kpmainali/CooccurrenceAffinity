#' Interval of alpha values for which X is a specified q'th quantile
#'
#' This function outputs the largest interval of log-odds parameter values alpha for which
#' the Extended Hypergeometric distribution function at x is >= q and the complementary distribution function 1 - F(x-) is >= 1-q.
#'
#' @details  This function outputs the endpoints  a1, a2 defined by
#'
#' F(x, a1) = q      and    F(x-1, a2) = q
#'
#' where  F(z, a) = F(z, mA,mB,N, exp(a)) is the extended Hypergeometric distribution function.
#'
#' The interval of alpha values with these endpoints a1, a2 is viewed as the set of alpha values "compatible" with x being a q'th quantile for the Extended Hypergeometric.

#'
#' @param x integer co-occurrence count that should properly fall within the closed interval from  max(0,mA+mB-N) to  min(mA,mB)
#' @param marg a 3-entry integer vector (mA,mB,N) consisting of the first row and column totals and the table total for a 2x2 contingency table
#' @param q a quantile falling strictly between 0 and 1
#' @param scal an integer parameter (default 2*N^2, capped at 10 within the function) that should be 2 or greater
#'
#' @return This function returns the vector (a1, a2) defined above, the endpoints of the set of alpha values for which x is a q'th quantile of the Extended Hypergeometric distribution.
#'
#' @author Eric Slud
#'
#' @references
#'
#' @export


EHypQuInt <-
  function(x, marg, q, scal=log(2*marg[3]^2)) {
    #  marg = c(mA, mB, N),
    # x an observed co-occurrence count, q the desired quantile
    require(BiasedUrn)
    # cap scal at 10 to avoid error in pFNCHypergeo
    scal = min(scal,10)
    mA=marg[1]; mB=marg[2]; N=marg[3]
    xmin = max(mA+mB-N,0);  xmax = min(mA,mB)
    maxabs = log(2*N^2)
    if(length(intersect(c(mA,mB), c(0,N))))
      return("Degenerate co-occurrence distribution!")
    null.mean = mA*mB/N
    newint = if(x >= null.mean) c(-1,scal) else c(-scal,1)

    #--------------
    # The function EHypCent() is a utility forming part of the EHypQuInt function, and would not be called as a stand-alone.
    # Its output is a vector of centered Extended-Hypergeometric distribution function values with the quantile q subtracted.
    # The function is used in root-finding to locate alpha-value intervals for which the corresponding
    # Extended Hypergeometric distribution assigns distribution function close to q at t.

    EHypCent = function(alp,t) {
      # t an observed co-occurrence count, alp a sequence of alphas
      L = length(alp)
      out = numeric(L)
      for(i in 1:L) out[i]= pFNCHypergeo(t,mA,N-mA,mB,exp(alp[i]))-q
      out }
    #-------------

    if(x < xmin) c(-Inf,-Inf) else {
      if(x == xmin) c(-maxabs,
                      uniroot(EHypCent, newint, t=x, extendInt="y")$root) else {
                        if(x > xmax) c(Inf,Inf) else {
                          if(x == xmax) c(uniroot(EHypCent,
                                                  newint, t=x-1, extendInt="yes")$root,maxabs) else
                                                    c(uniroot(EHypCent, newint, extendInt="y", t=x-1)$root,
                                                      uniroot(EHypCent, newint, extendInt="y", t=x)$root) }}}}
