#' Bisect
#'
#' This function bisects...
#'
#' @details   This function bisects...
#'
#' @param x nteger co-occurrence counts that should properly fall within the closed interval  [max(0,mA+mB-N), min(mA,mB)]
#' @param marg a 3-entry integer vector containing (mA,mB,N)
#' @param scal an integer parameter (default 10) that should fall somewhere between 2 and 10
#' @param lev a confidence level, generally somewhere from 0.8 to 0.95  (default 0.95)
#'
#' @return This function returns
#'
#' @author Eric Slud
#'
#' @references to be added
#'
#' @example
#' to be added
#'
#' @export

Bisect <-
  function(ffnc, intrv, tol=1e-8) {
    #  root of increasing possibly non-smooth function by bisections
    if(ffnc(intrv[1])>= tol | ffnc(intrv[2])<= tol)
      return("Same sign at endpts!")
    while(abs(intrv[2]-intrv[1]) > tol) {
      newp = (intrv[2]+intrv[1])/2
      newv = ffnc(newp)
      if(newv<0) intrv[1]=newp else intrv[2]=newp }
    newv=ffnc((intrv[2]+intrv[1])/2)
    c(root=(intrv[1]+intrv[2])/2, fval=newv)
  }
