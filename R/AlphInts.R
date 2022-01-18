#' Five intervals, three using EHypQuInt(), one using EHypMidP(), and one using AcceptAffCI()
#'
#' This function calculates five intervals, one of them should be interpreted as an interval of alpha that is inevitable in discrete data.
#' The other intervals are various types of Confidence Intervals.
#'
#' @details  This function calculates five intervals, three of them using EHypQuInt, one using EHypMidP,
#' and one using AcceptAffCI. First ("Int1") is the interval of alpha values compatible with x as median for the
#' Extended Hypergeometric distribution (with fixed margins and alpha);
#' second ("Int2") an "exact" conservative test-based confidence intrval (in a sense analogous to the Clopper-Pearson (1934)
#' confidence interval for unknown binomial proportion) for alpha based on data (x,mA,mB,N);
#' third the Acceptability Confidence Interval ("Int3") of Blaker (2000)
#' which is a better confidence interval than the CP-type interval "Int2" in the sense of being contained within "Int2"
#' but still provably conservative.
#' The fourth confidence interval ("Int4") is the one given in formula (2) above of the Introduction to this documentation,
#' with endpoints obtained as the midpoints of quantile intervals respectively to the (1+lev)/2 and (1-lev)/2
#' quantiles of the Extended Hypergeometric distribution; and the fifth ("Int5")
#' which behaves very similarly to "Int4" is defined by the midP approach analogous to the midP confidence interval
#' for binomial proportions (Agresti 2013, p.605), and is calculated from EHypMidP.
#'
#' The first of these intervals provides quantification of the underlying discreteness of the
#' Extended Hypergeometric and the impact of that discreteness on the estimation of alpha:
#' Int1 is an interval that will contain the MLE alpha-hat, and the mid-point of that interval is another reasonable estimator of alpha from the data.
#' The recommended (slightly conservative) confidence interval is Int3,
#' while the very similar intervals Int4 and Int5 have coverage typically closer than
#' Int2 or Int3 to the nominal level f coverage, at the cost of occasionally under-covering by as much as 0.04 or 0.05 for confidence levels 0.90 or 0.95.
#'
#' @param x integer co-occurrence counts that should properly fall within the closed interval  [max(0,mA+mB-N), min(mA,mB)]
#' @param marg a 3-entry integer vector containing (mA,mB,N)
#' @param scal an integer parameter (default 10) that should fall somewhere between 2 and 10
#' @param lev a confidence level, generally somewhere from 0.8 to 0.95  (default 0.95)
#'
#' @return This function returns various types of intervals including confidence intervals
#'
#' @author Eric Slud
#'
#' @references to be added
#'
#' @example
#' to be added
#'
#' @export

AlphInts <-
  function(x, marg, scal=log(2*marg[3]^2), lev=0.95) {
    # function to calculate intervals for alpha = log(odds)  values
    # Interval 1 = EHypQuInt interval for q=0.5
    # Interval 2 = Exact (conservative, CP-style) CI giving outer interval
    #     from EHypQuInt evaluations at q=(1+lev)/2 and (1-lev)/2
    # Interval 3 = Conservative Blaker (2000) style (better than Intrv.2)
    # Interval 4 = (Avg of q=(1+lev)/2 val's, Avg of q=(1-lev)/2  val's)
    # Interval 5 = midP style interval
    require(BiasedUrn)
    alph = 1-lev
    Int2 = c(EHypQuInt(marg, x, (1+lev)/2)[1],
             EHypQuInt(marg, x, (1-lev)/2)[2])
    list(Int1 = EHypQuInt(marg, x, 0.5),
         Int2 = Int2,
         Int3 = AcceptAffCI(x,marg,lev, Int2),
         Int4 = c(mean(EHypQuInt(marg, x, (1+lev)/2)),
                  mean(EHypQuInt(marg, x, (1-lev)/2))),
         Int5 = EHypMidP(marg,x,lev))
  }
