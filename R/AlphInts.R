#' Median interval, four confidence intervals, null expectation of cooccurrence count, and p-value
#'
#' This function calculates
#' (i) MedianIntrvl, the interval of alpha values for which the co-occurrence count is a median,
#' (ii) four Confidence Intervals, two using EHypQuInt(), one using EHypMidP(), and one using AcceptAffCI(),
#' (iii) the Expected Co-occurrence count under the Null distribution, and
#' (iv) the p-value for the observed co-occurrence count.
#'
#' @details  This function calculates five intervals, three of them using EHypQuInt, one using EHypMidP, and one using AcceptAffCI.
#' First ("MedianIntrvl") is the interval of alpha values compatible with x as median for the Extended Hypergeometric distribution (Harkness 1965)
#' with fixed margins and alpha; second ("CI.CP") an "exact" conservative test-based 2-sided confidence interval (analogous to the Clopper-Pearson (1934)
#' confidence interval for unknown binomial proportion) for alpha based on data (x,mA,mB,N); third the Acceptability Confidence Interval ("CI.Blaker")
#' of Blaker (2000, Theorem 1) which is a better confidence interval than the CP-type interval "CI.CP" in the sense of being contained within "CI.CP"
#' but still provably conservative, i.e., with coverage probability always at least as large as the nominal level.
#' The fourth confidence interval ("CI.midQ") is the one given in formula (2) above of the Introduction to this documentation,
#' with endpoints obtained as the midpoints of quantile intervals respectively to the (1+lev)/2 and (1-lev)/2 quantiles of the Extended Hypergeometric distribution;
#' and the fifth ("CI.midP") which behaves very similarly to "CI.midQ" is defined by the midP approach analogous
#' to the midP confidence interval for binomial proportions (Agresti 2013, p.605), and is calculated from EHypMidP.
#'
#' The first of these intervals quantifies the underlying discreteness of the Extended Hypergeometric and its impact on the estimation of alpha.
#' MedianIntrvl is an interval that will contain the MLE alpha-hat, and the mid-point of that interval is another reasonable estimator of alpha from the data.
#' The recommended (slightly conservative) confidence interval is CI.Blaker, while the very similar intervals CI.midQ and CI.midP have
#' coverage generally closer than CI.CP or CI.Blaker to the nominal level of coverage, at the cost of occasionally under-covering
#' by as much as 0.04 or 0.05 for confidence levels 0.90 or 0.95. The comparison among intervals, and different possible goals that CIs of
#' conservative or close-to-nominal coverage can serve, are similar to those compared by  Brown et al. (2001) for interval estimation of an unknown binomial proportion.
#'
#' Two other output list components are computed. First is Null.Exp, the expected co-occurrence count under the null (hypergeometric, corresponding to alpha=0)
#' distribution, and second is the two-sided p-value for the equal-tailed test of the null hypothesis alpha=0. This p-value is calculated when pval="Blaker"
#' according to Blaker's (2000) "Acceptability" function; if the input parameter pval is anything else, the p-value is calculated using the same idea as the midP confidence interval.

#'
#' @param x integer co-occurrence count that should properly fall within the closed interval \[max(0,mA+mB-N), min(mA,mB)\]
#' @param marg a 3-entry integer vector (mA,mB,N) consisting of the first row and column totals and the table total for a 2x2 contingency table
#' @param scal an integer parameter (default 2*N^2, capped at 10 within the function) that should be 2 or greater
#' @param lev a confidence level, generally somewhere from 0.8 to 0.95  (default 0.95)
#' @param pvalType a character string telling what kind of p-value to calculate. ‘Blaker’ or “midP’.
#' If ‘pvalType=Blaker” (the default value), the p-value is calculated according to "Acceptability" function of Blaker (2000).
#' If ‘pvalType=midP’, the p-value is calculated using the same idea as the midP confidence interval.
#'
#' @return A list of seven components: the median interval MedianIntrvl; the four two-sided Confidence Intervals described above,
#' two (CI.CP and CI.Blaker) conservative and two (CI.midQ and CI.midP) with coverage probabilities generally closer to the nominal level;
#' the null expectation Null.Exp of the co-occurrence count associated with alpha=0; and pval, the two-sided p-value for the hypothesis test of alpha=0,
#' calculated by the method selectied, which is the Blaker acceptability-function method if pvalType="Blaker" and otherwise
#' the "midP" p-value associated with the midP confidence-interval type.
#'
#' Of the four Confidence intervals produced, CI.Blaker is the recommended conservative interval and CI.midP the interval to use if coverage close to the nominal is desired.
#'
#' @author Eric Slud
#'
#' @references
#' Agresti, A. (2013) Categorical Data Analysis, 3rd edition, Wiley.
#'
#' Blaker, H. (2000), “Confidence curves and improved exact confidence intervals for discrete distributions", Canadian Journal of Statistics 28, 783-798.
#'
#' Brown, L., T. Cai, and A. DasGupta (2001), “Interval Estimation for a Binomial Proportion,” Statistical Science, 16, 101–117.
#'
#' Clopper, C., and E. Pearson (1934), “The Use of Confidence or Fiducial Limits Illustrated in the Case of the Binomial,” Biometrika, 26, 404–413.
#'
#' Fog, A. (2015), BiasedUrn: Biased Urn Model Distributions. R package version 1.07.
#'
#' Harkness, W. (1965), “Properties of the extended hypergeometric distribution“, Annals of Mathematical Statistics, 36, 938-945.
#'
#' @example
#' examples/AlphInts_example.R
#'
#' @export

AlphInts <-
  function(x, marg, scal=log(2*marg[3]^2), lev=0.95, pvalType="Blaker") {
    # function to calculate intervals for alpha = log(odds)  values
    # p-val either from Baker or midP CI
    # Interval 1 = MedianIntrvl =   EHypQuInt interval for q=0.5
    # Interval 2 = CI.CP =  Exact (conservative, CP-style) CI giving
    #     from EHypQuInt evaluations at q=(1+lev)/2 and (1-lev)/2
    # Interval 3 = CI.Blaker = Conservative Blaker (2000) style (better than Intrv.2)
    # Interval 4 = CI.midQ = (Avg of q=(1+lev)/2 val's, Avg of q=(1-lev)/2  val's)
    # Interval 5 = CI.midP =  midP style interval
    # require(BiasedUrn)
    # to avoid errors in pFNCHypergeo, cap  scal  at 10
    mA=marg[1]; mB=marg[2]; N=marg[3]
    if(x<0 | x< mA+mB-N | x > min(mA,mB)) return("Impossible x!")
    if(length(intersect(c(mA,mB), c(0,N)))){
      warning(paste0("If the mA or mB value is equal to 0 or N,","\n",
                "then the corresponding co-occurrence distribution is degenerate at min(mA,mB).", "\n",
                "This means that the co-occurrence count X will always be min(mA,mB) regardless of alpha.", "\n",
                "In this case alpha is undefined, and no computations are done."),
            call. = TRUE, immediate. = TRUE)
      # warning("If the mA or mB value is equal to 0 or N, then the corresponding co-occurrence distribution is degenerate at min(mA,mB). This means that the co-occurrence count X will always be min(mA,mB) regardless of alpha. In this case alpha is undefined, and no computations are done.", call. = TRUE, immediate.=TRUE)
      return("Degenerate co-occurrence distribution!")}
    if(x==max(mA+mB-N,0))
      warning(paste("MLE = -Infty is capped, along with lower confidence limits","\n"), call. = TRUE, immediate.=TRUE)
    if(x==min(mA,mB))
      warning(paste("MLE = Infty is capped, along with upper confidence limits","\n"), call. = TRUE, immediate.=TRUE)

    scal = min(scal,10)
    alph = 1-lev
    Int2 = c(EHypQuInt(x, marg, (1+lev)/2)[1],
             EHypQuInt(x, marg, (1-lev)/2)[2])
    Int4 = c(mean(EHypQuInt(x, marg, (1+lev)/2)),
             mean(EHypQuInt(x, marg, (1-lev)/2)))
    Int5 = EHypMidP(x, marg,lev)
    list(MedianIntrvl = EHypQuInt(x,marg, 0.5),
         CI.CP = Int2, CI.Blaker = AcceptAffCI(x,marg,lev, Int2),
         CI.midQ = Int4, CI.midP = if(is.na(sum(Int5))) Int4 else Int5,
         Null.Exp = mA*mB/N, pval = if(pvalType=="Blaker")
           AcceptAffin(x,marg,0) else {
             pr = phyper(x,mA,N-mA,mB)+phyper(x-1,mA,N-mA,mB)
             min(pr,2-pr) }) }
