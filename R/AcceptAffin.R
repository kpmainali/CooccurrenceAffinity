#' Calculates the "Acceptability Function" used in defining Blaker's (2000) Acceptability Interval and computing the latter in the function AcceptAffCI().
#'
#' This function calculates the "Acceptability Function" of Blaker (2000, Thm.1, p.785) for the log-odds parameter alpha in the Extended Hypergeometric distribution.
#'
#' @details  This function calculates the "Acceptability Function" of Blaker (2000, Thm.1, p.785) for the log-odds parameter alpha
#' in the Extended Hypergeometric distribution, a function from which the "Acceptability Interval" is calculated by another CooccurrenceAffinity package function AcceptAffCI().
#'
#' @param x integer co-occurrence count that should properly fall within the closed interval \[max(0,mA+mB-N), min(mA,mB)\]
#' @param marg a 3-entry integer vector (mA,mB,N) consisting of the first row and column totals and the table total for a 2x2 contingency table
#' @param alph a vector of (one or more) real-valued "alpha" values, where alpha  is the log-odds parameter in the Extended Hypergeometric distribution
#'
#' @return This function returns the "Acceptability Function" that is later used by another function AcceptAffCI() to compute "Acceptability Interval".
#'
#' @author Eric Slud
#'
#' @references
#' Blaker, H. (2000), “Confidence curves and improved exact confidence intervals for discrete distributions", Canadian Journal of Statistics 28, 783-798.
#'
#' @example
#' to be added
#'
#' @export


AcceptAffin <-
  function(x, marg, alph) {
    mA=marg[1]; mB=marg[2]; N=marg[3]
    if(length(intersect(c(mA,mB), c(0,N))))
      return("Degenerate co-occurrence distribution!")
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
    out   }
