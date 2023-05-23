#' midP.EHyp computation
#'
#' Helper function
#'
#' @details   This is a helper function.
#'
#' param x integer co-occurrence count that should properly fall within the closed interval \[max(0,mA+mB-N), min(mA,mB)\]
#' @param alp "alpha" parameter, the log-odds parameter in the Extended Hypergeometric distribution
#'
#' @return helper function for midP CI computation with EHypMidP
#'
#' @author Eric Slud
#'
#'
#' @export

midP.EHyp <-
  function(alp)  {
    x <- mA <- N <- mB <- NULL
    BiasedUrn::pFNCHypergeo(x,mA,N-mA,mB,exp(alp))-0.5*BiasedUrn::dFNCHypergeo(x,mA,N-mA,mB,exp(alp)) }
