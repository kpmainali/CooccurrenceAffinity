#' Acceptability Interval
#'
#' This function calculates the "Acceptability Interval" of Blaker for the log-odds parameter alpha in the Extended Hypergeometric distribution.
#'
#' @details  This function calculates the "Acceptability Interval" based on "Acceptability Function" computed by AcceptAffin().
#' This interval, developed by Blaker (2000), was proved in that paper's Theorem 1 in a more general class of estimation problems
#' to have three essential properties: it falls within the CI.CP confidence interval; it maintains the property of being conservative,
#' i.e., of having coverage probability under the Extended Hypergeometric (mA,mB,N, alpha) distribution at least as large as the nominal level;
#' and it is larger when the confidence level is larger.
#'
#' @param x integer co-occurrence count that should properly fall within the closed interval \[max(0,mA+mB-N), min(mA,mB)\]
#' @param marg a 3-entry integer vector (mA,mB,N) consisting of the first row and column totals and the table total for a 2x2 contingency table
#' @param lev a confidence level, generally somewhere from 0.8 to 0.95  (default 0.95)
#' @param CPint the exact conservative ("Clopper-Pearson-type") interval CI.CP calculated in the function AlphInts()
#'
#' @return This function returns the "Acceptability Interval" of Blaker (2000). The code is adapted from Blaker's Splus code for the case of an unknown binomial proportion.
#'
#' @author Eric Slud
#'
#' @references
#' Blaker, H. (2000), “Confidence curves and improved exact confidence intervals for discrete distributions", Canadian Journal of Statistics 28, 783-798.
#'
#' @example
#' inst/examples/AcceptAffCI_example.R
#'
#' @export

AcceptAffCI <-
  function(x, marg, lev, CPint) {
    # falls within CP-type interval; enpdts given by
    # roots of equation AcceptAffin(x,marg, ., lev) = 1-lev
    alph=1-lev
    mA=marg[1]; mB=marg[2]; N=marg[3]
    if(length(intersect(c(mA,mB), c(0,N))))
      return("Degenerate co-occurrence distribution!")
    maxabs = min(log(2*N^2), 10)
    lower= -maxabs; upper= maxabs
    # AcceptAffin not smooth, so "invert" via bisections
    minx = max(mA+mB-N,0)
    maxx = min(mA,mB)
    ini.est = if(x==minx) -maxabs else {
      if(x==maxx) maxabs else
        log(x*(N-mA-mB+x)/((mA-x)*(mB-x)))  }
    if(x > minx) lower = Bisect(function(alph)
      AcceptAffin(x,marg,alph)-1+lev, c(CPint[1],ini.est))[1]
    if(x < maxx) upper = Bisect(function(alph)
      -AcceptAffin(x,marg,alph)+1-lev, c(ini.est, CPint[2]))[1]
    as.numeric(c(lower,upper))  }
