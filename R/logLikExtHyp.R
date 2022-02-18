#' log of Extended Hypergeometric Likelihiood at (X, mA,mB,N, alpha)
#'
#' This function calculates the logarithm of the Extended Hypergeometric likelihood at specified x and alpha, with marginal totals mA, mB, N fixed.
#'
#' @details   This is simply the logarithm of the Extended Hypergeometric (Harkness 1965) or Fisher noncentral Hypergeometric, as calculated by the R package BiasedUrn. The formula is  log(pFNCHypergeo(x,mA,N-mA,mB,exp(alpha))
#'
#' @param x integer co-occurrence count that should properly fall within the closed interval \[max(0,mA+mB-N), min(mA,mB)\]
#' @param marg a 3-entry integer vector (mA,mB,N) consisting of the first row and column totals and the table total for a 2x2 contingency table
#' @param alpha a real number, the log odds ratio or affinity parameter for the 2x2 contingency table
#'
#' @return scalar loglikelihood value
#'
#' @author Eric Slud
#'
#' @references
#' Fog, A. (2015), BiasedUrn: Biased Urn Model Distributions. R package version 1.07.
#'
#' Harkness, W. (1965), “Properties of the extended hypergeometric distribution“, Annals of Mathematical Statistics, 36, 938-945.
#'
#' @example
#' examples/logLikExtHyp_example.R
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
