#' Coverage Probabilities for Confidence Intervals about alpha, for fixed true alpha
#'
#' This function calculates the coverage probability at the true value alpha of the four types of Cpnfidence Intervals (CI.CP, CI.Blaker, CI.midQ, CI.midP) computed in AlphInts().
#'
#' @details   See AlphInts() documentation for details of computation of the four confidence intervals CI.CP, CI.Blaker, CI.midQ, CI.midP. The confidence intervals are calculated for each x in the allowed range from max(mA+mB-N,0) to min(mA,mB), and the probability that X=x times the indicator of alph falling in each of them is summed.
#'
#' @param marg a 3-entry integer vector (mA,mB,N) consisting of the first row and column totals and the table total for a 2x2 contingency table
#' @param alph True log-odds-ratio value alpha at which coverage probabilities (under Extended Hypergeometric with parameters mA,mB,N, exp(alp)) are to be calculated
#' @param scal an integer parameter (default 2*N^2, capped at 10 within the function) that should be 2 or greater
#' @param lev a confidence level, generally somewhere from 0.8 to 0.95 (default 0.95)
#'
#' @return A vector covPrb containing the coverage propbabilities for the four Confidence Intervals
#'
#' @author Eric Slud
#'
#' @references
#' to be added
#'
#' @example
#' inst/examples/Covrg_example.R
#'
#' @export


Covrg <-
  function(marg, alph, scal=log(2*marg[3]^2), lev=0.95) {
    # require(BiasedUrn)
    mA=marg[1]; mB=marg[2]; N=marg[3]
    if(length(intersect(c(mA,mB), c(0,N))))
      return("Degenerate co-occurrence distribution!")
    xrng = c(max(0, mA+mB-N), min(mA,mB))
    xn = xrng[2]-xrng[1]+1
    arrCI = array(0, c(xn,4,2))
    minx = xrng[1]-1
    for(k in 2:(xn-1)) {
      tmp = AlphInts(k+minx,marg,scal=scal, lev=lev)
      arrCI[k,,] = rbind(tmp$CI.CP, tmp$CI.Blaker,
                         tmp$CI.midQ, tmp$CI.midQ) }
    ## special case for X = xrng[1] or xrng[2]
    arrCI[1,,] = rep(1,4) %o% MinX.Int(marg, scal=scal, lev=(1+lev)/2)
    arrCI[xn,,] = rep(1,4) %o% MaxX.Int(marg, scal=scal, lev=(1+lev)/2)
    ## arrCI[xn,2,] = 0.5*(log(2*marg[3]^2) + MaxX.Int(marg, scal=scal, lev=lev))
    ## step 2 array of coverage indicators times prob's
    covCI = array(0, c(xn,4))
    for(j in 1:4) for(b in 1:xn) covCI[b,j] =
      (alph >= arrCI[b,j,1] & alph <= arrCI[b,j,2])*
      BiasedUrn::dFNCHypergeo(minx+b,mA,N-mA,mB,exp(alph))
    ## step 3 coverage prob's
    c(covPrb = t(covCI) %*% rep(1,xn))   }
