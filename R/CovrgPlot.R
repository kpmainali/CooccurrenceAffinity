#' Coverage probabilities of the confidence intervals
#'
#' This function calculates and plots the coverage probabilities of the confidence intervals CI1 and CI2, which is useful in showing how conservative
#' and more nomimal recommended confidence inervals map against each other.
#'
#' @details   This function calculates and plots the coverage probabilities of the confidence intervals CI1 and CI2,
#' for fixed marginals "marg" = (mA,mB,N) and confidence level "lev" as a function of alpha.
#' The plots serve as a useful exhibit to show the range of always conservative (and sometimes quite conservative)
#' coverage probabilities for CI1, and the much closer-to-nominal coverage probabilities for the
#' recommended Confidence Interval CI2. The outputted vector "alphvec" of alpha values is the set of
#' all discontinuities in the slope of the coverag-probability function together with the corresponding exact (not simulated)
#' coverage probabilities of each of the two confidence intervals at the alphvec points in the two-column array  covPrb.
#'
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
#'
#' @export

CovrgPlot <-
  function(marg, scal=log(2*marg[3]^2), lev=0.95,
           alim = 4, minp = 0.7, intpts = 1, plotCI=1:4)  {
    ## plotCI= set of integer indices of 4 CI coverage curves to plot
    #  plot is restricted to alpha in interval in (-alim, alim)
    #    and shows only coverage probability values >= minp
    #  parameter intpts controls the number of points used
    #    to plot curves (more, for larger intpts)
    ## step 0 the vector of alpha values to use as evaluation points
    require(BiasedUrn)
    mA=marg[1]; mB=marg[2]; N=marg[3]
    if(length(intersect(c(mA,mB), c(0,N))))
      return("Degenerate co-occurrence distribution!")
    xrng = c(max(0, mA+mB-N), min(mA,mB))
    xn = xrng[2]-xrng[1]+1
    avec0L = avec0U = numeric(xn)
    avec0L[1] = -log(2*marg[3]^2) ; avec0U[xn] = log(2*marg[3]^2)
    avec0U[1] = MinX.Int(marg, scal=scal, lev=lev)[2]
    avec0L[xn] = MinX.Int(marg, scal=scal, lev=lev)[1]
    for(x in 2:(xn-1)) {
      avec0L[x] = EHypQuInt(x+xrng[1]-1, marg, (1+lev)/2, scal=scal)[1]
      avec0U[x] = EHypQuInt(x+xrng[1]-1, marg, (1-lev)/2, scal=scal)[1] }
    avec = sort(union(round(avec0L,4), round(avec0U,4)))
    ## add at least intpts extra equally spaced points between all the points
    avec = bvec = avec[abs(avec)<alim+0.5]
    for(j in 2:length(avec)) bvec = c(bvec, seq(avec[j-1],avec[j],
                                                length=intpts+2))
    avec = sort(unique(bvec))
    an = length(avec)
    ## step 1 the array of CIs
    arrCI = array(0, c(xn,4,2))
    minx = xrng[1]-1
    for(k in 2:(xn-1)) {
      tmp = AlphInts(k+minx,marg,scal=scal, lev=lev)
      arrCI[k,,] = rbind(tmp$CI.CP, tmp$CI.Blaker,
                         tmp$CI.midQ, tmp$CI.midP) }
    ## special case for X = xrng[1] or xrng[2]
    arrCI[1,,] = rep(1,4) %o% MinX.Int(marg, scal=scal, lev=(1+lev)/2)
    arrCI[xn,,] = rep(1,4) %o% MaxX.Int(marg, scal=scal, lev=(1+lev)/2)
    ## step 2 array of coverage indicators times prob's
    covCI = array(0, c(an, xn,4))
    for(j in 1:4) for(k in 1:an) for(b in 1:xn) covCI[k,b,j] =
      (avec[k] >= arrCI[b,j,1] & avec[k] <= arrCI[b,j,2])*
      dFNCHypergeo(minx+b,mA,N-mA,mB,exp(avec[k]))
    ## step 3 coverage prob's
    covPrb = array(0, c(an,4))
    for(j in 1:4) covPrb[,j] = c(covCI[,,j] %*% rep(1,xn))
    ## graph preparation
    colors = c("black","blue","green","orange")
    plot(avec, c(minp,rep(1,an-1)), xlab="alpha", ylab="Cov Prob", type="n",
         ylim=c(minp,1), xlim=c(-alim,alim), main=paste0("Coverage Prob's",
                                                         " for 4 Types of EHyp ",100*lev,"% CIs, ","(",mA,",",mB,",",N,")"))
    for(j in plotCI) lines(avec,covPrb[,j], type="b",lty=j, col=colors[j])
    legend(-alim,minp+0.20*(1-minp), legend=c("Consrv.Exact",
                                              "Consrv.Accept","MidQ.Exact", "MidP.Exact"),
           lty=1:2, lwd=1.6, col=colors)
    abline(h=lev, col="brown", lwd=2, lty=5)
    list(covPrb=covPrb, arrCI=arrCI, avec=avec)}
