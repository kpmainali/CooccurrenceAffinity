#' midP.EHyp computation
#'
#' Helper function
#'
#' @details   This is a helper function.
#'
#' @param x integer co-occurrence count that should properly fall within the closed interval \[max(0,mA+mB-N), min(mA,mB)\]
#'
#' @return
#'
#' @author Eric Slud
#'
#' @references
#'
#' @example
#' to be added
#'
#' @export

midP.EHyp <-
  function(alp)  pFNCHypergeo(x,mA,N-mA,mB,exp(alp))-0.5*dFNCHypergeo(x,mA,N-mA,mB,exp(alp))

