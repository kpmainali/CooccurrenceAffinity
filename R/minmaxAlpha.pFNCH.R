#' integer-endpoint of range for which BiasedUrn::pFNCHHypergeo() works without error
#'
#' This function calculates an integer-endpoint of range for which BiasedUrn::pFNCHHypergeo() works without error.
#'
#' @details   Without this function, BiasedUrn::pFNCHHypergeo() returns inconsistency message for extreme examples like: AlphInts(20,c(204,269,2016), lev=0.9, scal=10).
#' This problem is solved within our package by restricting the range of allowed alpha to the computed (alphmin, alphmax) range.
#'
#' @param x integer co-occurrence count that should properly fall within the closed interval \[max(0,mA+mB-N), min(mA,mB)\]
#' @param marg a 3-entry integer vector (mA,mB,N) consisting of the first row and column totals and the table total for a 2x2 contingency table
#'
#' @return minimum and maximum of Alpha
#'
#' @author Eric Slud
#'
#' @references
#' Fog, A. (2015), BiasedUrn: Biased Urn Model Distributions. R package version 1.07.
#'
#' Harkness, W. (1965), “Properties of the extended hypergeometric distribution“, Annals of Mathematical Statistics, 36, 938-945.
#'
#' @example
#' examples/minmaxAlpha.pFNCH_example.R
#'
#' @export


minmaxAlpha.pFNCH <- function(x, marg) {
  require(BiasedUrn)
  ## use pFNCHyper in BiasedUrn to establish range of alphas
  ## within (-10,10) over which that function works properly.
  mA=marg[1]; mB=marg[2]; N=marg[3]
  if(length(intersect(c(mA,mB), c(0,N))))
    return("Degenerate co-occurrence distribution!")
  pFNCH = function(alp) pFNCHypergeo(x,mA,N-mA,mB,exp(alp))
  tmp = try(pFNCH(10), silent=T)
  alphmax = if(class(tmp)!="try-error") 10 else NULL
  tmp = try(pFNCH(-10), silent=T)
  alphmin = if(class(tmp)!="try-error") -10 else NULL
  if(is.null(alphmax)) {
    tmp=NULL
    a0 = 0
    while(a0 < 9 & class(tmp)!="try-error") {
      a0 = a0+1
      tmp = try(pFNCH(a0), silent=T) }
    alphmax = a0-1 }
  if(is.null(alphmin)) {
    tmp=NULL
    a0 = 0
    while(a0 < 10 & class(tmp)!="try-error") {
      a0 = a0+1
      tmp = try(pFNCH(-a0), silent=T) }
    alphmin = -a0+1 }
  c(alphmin=alphmin, alphmax=alphmax) }
