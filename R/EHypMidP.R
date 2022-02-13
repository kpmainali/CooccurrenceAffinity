#' Quantile of the Extended Hypergeometric distribution approximated by the midP distribution function
#'
#' This function does the analogous calculation to that of EHypQuInt, but with the Extended Hypergeometric distribution
#' function F(x) = F(x,mA,mB,N, exp(alpha)) replaced by (F(x) + F(x-1))/2.
#'
#' @details  This function does the analogous calculation to that of CI.CP, but with the Extended Hypergeometric distribution
#' function F(z, alpha) = F(z,mA,mB,N, exp(alpha)) replaced by (F(z,alpha) + F(z-1,alpha))/2.
#'
#' @param x integer co-occurrence count that should properly fall within the closed interval from  max(0,mA+mB,N) to  min(mA,mB)
#' @param marg a 3-entry integer vector (mA,mB,N) consisting of the first row and column totals and the table total for a 2x2 contingency table
#' @param lev a confidence level, generally somewhere from 0.8 to 0.95  (default 0.95)
#'
#' @return This function returns the interval of alpha values with endpoints
#' (F(x,alpha)+F(x-1,alpha))/2 = (1+lev)/2   and  (F(x,alpha)+F(x+1,alpha))/2 = (1-lev)/2.
#'
#' The idea of calculating a Confidence Interval this way is analogous to the midP CI used for unknown binomial proportions (Agresti 2013, p.605).
#'
#' @author Eric Slud
#'
#' @references
#' Agresti, A. (2013) Categorical Data Analysis, 3rd edition, Wiley.
#'
#' @example
#' to be added
#'
#' @export

EHypMidP <-
  function(x, marg, lev) {
    mA=marg[1]; mB=marg[2]; N=marg[3]
    if(length(intersect(c(mA,mB), c(0,N))))
      return("Degenerate co-occurrence distribution!")
    ## cap absmax at 10 to avoid error in pFNCHypergeo
    absmax = min(log(2*N^2),10)
    ## NB. This function can fail to find interval when the
    #   marg  numbers are very large and the x value too extreme.
    #   In that case, the midQ interval is used in place of midP.
    midP.EHyp = function(alp)  pFNCHypergeo(x,mA,N-mA,mB,exp(alp))-0.5*
      dFNCHypergeo(x,mA,N-mA,mB,exp(alp))
    lower = if(x==max(mA+mB-N,0)) absmax else  {
      tmp = try(uniroot(function(alp) midP.EHyp(alp) - (1+lev)/2,
                        c(-absmax,absmax)), silent=T)
      if(class(tmp)!="try-error") tmp$root else NA}
    upper = if(x==min(mA,mB)) absmax else   {
      tmp = try(uniroot(function(alp) midP.EHyp(alp) - (1-lev)/2,
                        c(-absmax,absmax)), silent=T)
      if(class(tmp)!="try-error") tmp$root else NA}
    if(is.na(lower+upper)) NA else c(lower,upper) }
