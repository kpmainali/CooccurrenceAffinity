#' Repeat as stand-alone calculations of part of ML.Alpha()
#'
#' This function calculates ... to be added
#'
#' @details   to be added
#'
#' @param x integer co-occurrence count that should properly fall within the closed interval  [max(0,mA+mB-N), min(mA,mB)]
#' @param marg a 3-entry integer vector (mA,mB,N) consisting of the first row and column totals and the table total for a 2x2 contingency table
#'
#' @return to be added
#'
#' @author Eric Slud
#'
#' @references
#' to be added
#'
#' @example
#' to be added
#'
#' @export


logLikExtHyp <-
  function(x,marg,alpha)  {
    mA=marg[1]; mB=marg[2]; N=marg[3]
    out = numeric(length(alpha))
    require(BiasedUrn)
    for(i in 1:length(alpha))   out[i] =
      log(dFNCHypergeo(x,mA,N-mA,mB,exp(alpha[i])))
    out }
