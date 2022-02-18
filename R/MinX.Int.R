#' MinX.Int
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

MinX.Int <-
  function(marg, scal=log(2*marg[3]^2), lev=0.95) {
    ## special one-sided interval used when X= min possible value
    # cap scal at 10 to avod error in pFNCHypergeo
    scal=min(scal,10)
    if(length(intersect(marg[1:2], c(0,marg[3]))))
      return("Degenerate co-occurrence distribution!")
    minx = max(marg[1]+marg[2]-marg[3],0)
    Lower = -scal
    Upper = uniroot( function(xa)
      dFNCHypergeo(minx,marg[1],marg[3]-marg[1],marg[2],
                   exp(xa)) - (1-lev), c(-scal,1), extendInt="yes")$root
    c(Lower,Upper)  }
