#' Maximum likelihood estimate and intervals of alpha
#'
#' This function calculates the maximum likelihood estimate of alpha, which is the "affinity" index of co-occurrence
#' championed in our paper Mainali et al as an index of cooccurrence-based similarity, along with the intervals computed in AlphInts.
#'
#' @details   This function calculates the maximum likelihood estimate of alpha, which is the "affinity" index of co-occurrence
#' championed in our paper Mainali et al as an index of cooccurrence-based similarity, along with the intervals computed in AlphInts,
#' now called CI2 through CI5. The boolean "bound" parameter is an option to prevent the intervals containing alpha-estimates
#' to extend to plus or minus infinity, based on a Bayesian argument. The bound substituted for the Infinite endpoints is
#' provably larger than the largest value the MLE can take whenever x avoids the enpoints max(mA+mB-N,0) and min(mA,mB) of its logical range.
#' The recommended confidence interval for alpha is CI2.
#'
#' @param x nteger co-occurrence counts that should properly fall within the closed interval  [max(0,mA+mB-N), min(mA,mB)]
#' @param marg a 3-entry integer vector containing (mA,mB,N)
#' @param scal an integer parameter (default 10) that should fall somewhere between 2 and 10
#' @param lev a confidence level, generally somewhere from 0.8 to 0.95  (default 0.95)
#'
#' @return This function returns maximum likelihood estimate and various confidence intervals of alpha
#'
#' @author Eric Slud
#'
#' @references to be added
#'
#' @example
#' to be added

ML.Alpha <-
  function( x, marg, bound=T, scal=log(2*marg[3]^2), lev=0.95) {
    require(BiasedUrn)
    mA = marg[1]; mB = marg[2]; N = marg[3]
    if(x<0 | x< mA+mB-N | x > min(mA,mB)) return("Impossible x!")
    upbd = if(bound) log(2*N^2) else Inf
    ModInts = AlphInts(x, marg, lev=lev, scal=scal)
    for(k in 1:3) ModInts[[k]] = pmax(-upbd, pmin(upbd, ModInts[[k]]))
    tmp = optimize( function(t) log(dFNCHypergeo(x,mA,N-mA,mB,exp(t))),
                    ModInts[[1]]+c(-1,1), maximum=T)
    ### Outputs are: likelihood maximizing alpha, maximized logLik, and
    ###    flag=T if MLE falls within interval AlphInt Int1,
    ###    and two ConfInts (Int2, Int3 from AlphInts)
    list(est = max(-upbd,min(upbd,tmp$max)), LLK = -tmp$obj,
         Flag = (ModInts[[1]][1] <= tmp$max & ModInts[[1]][2] >= tmp$max),
         AlphInt = ModInts[[1]],
         lev = lev, CI2 = ModInts[[2]], CI3 = ModInts[[3]],
         CI4 = ModInts[[4]], CI5 = ModInts[[5]])
  }
