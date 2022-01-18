#' Largest interval of alpha values-2
#'
#' This function does the analogous calculation to that of EHypQuInt,
#' but with the Extended Hypergeometric distribution function  F(x) = F(x,mA,mB,N, exp(alpha)) replaced by  (F(x) + F(x-1))/2.
#'
#' @details  This function does the analogous calculation to that of EHypQuInt,
#' but with the Extended Hypergeometric distribution function  F(x) = F(x,mA,mB,N, exp(alpha)) replaced by  (F(x) + F(x-1))/2.
#'
#' @param marg a 3-entry integer vector containing (mA,mB,N)
#' @param x nteger co-occurrence counts that should properly fall within the closed interval  [max(0,mA+mB-N), min(mA,mB)]
#' @param lev a confidence level, generally somewhere from 0.8 to 0.95  (default 0.95)
#'
#' @return This function returns the largest interval of alpha values
#'
#' @author Eric Slud
#'
#' @references to be added
#'
#' @example
#' to be added
#'
#' @export

EHypMidP <-
  function(marg, x, lev) {
    mA=marg[1]; mB=marg[2]; N=marg[3]
    absmax = log(2*N^2)
    midP.EHyp = function(alp)  pFNCHypergeo(x,mA,N-mA,mB,exp(alp))-0.5*
      dFNCHypergeo(x,mA,N-mA,mB,exp(alp))
    lower = if(x==max(mA+mB-N,0)) absmax else
      uniroot(function(alp) midP.EHyp(alp) - (1+lev)/2,
        c(-absmax,absmax))$root
    upper = if(x==min(mA,mB)) absmax else
      uniroot(function(alp) midP.EHyp(alp) - (1-lev)/2,
        c(-absmax,absmax))$root
    c(lower,upper)
  }
