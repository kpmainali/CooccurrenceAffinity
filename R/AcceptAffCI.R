#' Acceptability Interval
#'
#' This function calculates the "Acceptability Interval"
#'
#' @details  This function calculates the "Acceptability Interval" based on "Acceptability Function" computed by AcceptAffin().
#'
#' @param marg a 3-entry integer vector containing (mA,mB,N)
#' @param x nteger co-occurrence counts that should properly fall within the closed interval  [max(0,mA+mB-N), min(mA,mB)]
#' @param lev a confidence level, generally somewhere from 0.8 to 0.95  (default 0.95)
#' @param CPint the exact conservative ("Clopper-Pearson-type") interval calculated as $Int2 in the function AlphInts()
#'
#' @return This function returns the "Acceptability Interval"
#'
#' @author Eric Slud
#'
#' @references to be added
#'
#' @example
#' to be added

AcceptAffCI <-
  function(x, marg, lev, CPint){
    # falls within CP-type interval; enpdts given by
    # roots of equation AcceptAffin(x,marg, ., lev) = 1-lev
    alph=1-lev
    mA=marg[1]; mB=marg[2]; N=marg[3]
    maxabs = log(2*N^2)
    lower= -maxabs; upper= maxabs
    # AcceptAffin not smooth, so "invert" via bisections
    minx = max(mA+mB-N,0)
    maxx = min(mA,mB)
    ini.est = if(x==minx) -maxabs else {
      if(x==maxx) maxabs else
        log(x*(N-mA-mB+x)/((mA-x)*(mB-x)))  }
    if(x > minx) lower = Bisect(function(alph)
      AcceptAffin(x,marg,alph,lev)-1+lev, c(CPint[1],ini.est))[1]
    if(x < maxx) upper = Bisect(function(alph)
      -AcceptAffin(x,marg,alph,lev)+1-lev, c(ini.est, CPint[2]))[1]
    c(lower,upper)
  }
