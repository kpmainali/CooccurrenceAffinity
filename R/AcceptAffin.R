#' Calculate the "Acceptability Function"
#'
#' This function calculates the "Acceptability Function" of Blaker (2000, Thm.1, p.785) for the Extended Hypergeometric distribution.
#'
#' @details  This function calculates the "Acceptability Function" of Blaker (2000, Thm.1, p.785) for the Extended Hypergeometric distribution,
#' a function from which the "Acceptability Interval: is calculated by another
#' Affinity package function  AcceptAffCI = function(x, marg, lev, CPint).
#' Here "CPint" is the exact conservative ("Clopper-Pearson-type") interval calculated as $Int2 in the function AlphInts below.
#'
#' @param marg a 3-entry integer vector containing (mA,mB,N)
#' @param x nteger co-occurrence counts that should properly fall within the closed interval  [max(0,mA+mB-N), min(mA,mB)]
#' @param alp a vector of (one or more) real "alpha" values
#' @param lev a confidence level, generally somewhere from 0.8 to 0.95  (default 0.95)
#'
#' @return This function returns the "Acceptability Function" that is later used by another function AcceptAffCI() to compute "Acceptability Interval"
#'
#' @author Eric Slud
#'
#' @references to be added
#'
#' @example
#' to be added

AcceptAffin <-
  function(x, marg, alph, lev) {
    mA=marg[1]; mB=marg[2]; N=marg[3]
    K = length(alph)
    out=numeric(K)
    for(i in 1:K) {
      ealp = exp(alph[i])
      p1 = 1 - pFNCHypergeo(x - 1, mA,N-mA,mB,ealp)
      p2 = pFNCHypergeo(x, mA,N-mA,mB,ealp)
      a1 = p1 + pFNCHypergeo(qFNCHypergeo(p1,mA,N-mA,mB,ealp) - 1,
                          mA,N-mA,mB,ealp)
      a2 = p2+1 - pFNCHypergeo(qFNCHypergeo(1-p2, mA,N-mA,mB,ealp),
                          mA,N-mA,mB,ealp)
      out[i] = pmin(a1,a2) }
    out
  }
