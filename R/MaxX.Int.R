#' MaxX.Int computation
#'
#' Helper function
#'
#' @details   This is a helper function.
#'
#' @param marg a 3-entry integer vector (mA,mB,N) consisting of the first row and column totals and the table total for a 2x2 contingency table
#' @param scal an integer parameter (default 2*N^2, capped at 10 within the function) that should be 2 or greater
#' @param lev a confidence level, generally somewhere from 0.8 to 0.95  (default 0.95)
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

MaxX.Int <-
  function(marg, scal=log(2*marg[3]^2), lev=0.95) {
    ## special one-sided interval used when X= max possible value
    # cap scal at 10 to avoid error in pFNCHypergeo
    if(length(intersect(marg[1:2], c(0,marg[3]))))
      return("Degenerate co-occurrence distribution!")
    scal = min(scal,10)
    maxx = min(marg[1:2])
    Upper = scal
    Lower = uniroot( function(xa)
      dFNCHypergeo(maxx,marg[1],marg[3]-marg[1],marg[2],
                   exp(xa)) - (1-lev), c(-1,scal), extendInt="yes")$root
    c(Lower,Upper)  }

