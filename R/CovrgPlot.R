#' Coverage probabilities of the confidence intervals, calculated and plotted
#'
#' This function calculates and plots the coverage probabilities of the four confidence intervals CI.CP, CI.Blaker, CI.midQ and CI.midP
#' produced by functions AlphInts() and ML.Alpha(). The coverage plots are useful in showing how the conservative confidence intervals (CI.CP and CI.Blaker)
#' tend to have above-nominal coverage probabilities versus the intervals CI.midQ and CI.midP that have coverage much closer to the nominal confidence level lev.
#'
#' @details   This function calculates and plots the coverage probabilities of the confidence intervals CI.CP, CI.Blaker, CI.midQ and CI.midQ produced by function AlphInts(),
#' for fixed marginals "marg" = (mA,mB,N) and confidence level "lev" as a function of alpha. The plots serve as a useful exhibit to show the range of always conservative
#' (and sometimes quite conservative) coverage probabilities for CI1, and the much closer-to-nominal coverage probabilities for the recommended Confidence Interval CI2.
#' The outputted vector "avec" of alpha values is the set of all discontinuities in the slope of the coverage-probability function together with
#' "intpts: equally spaced values between them. The plotted coverage probabilities are the corresponding exact (not simulated) coverage probabilities
#' of each of the selected (up to 4) confidence intervals at the alphvec points in the four-column array "covPrb".
#'
#' @param marg a 3-entry integer vector containing (mA,mB,N)
#' @param scal an integer parameter (default 2*N^2, capped at 10 within the function) that should be 2 or greater
#' @param lev a confidence level, generally somewhere from 0.8 to 0.95 (default 0.95)
#' @param alim the absolute value of alpha at which the plotted figure is cutoff; so the x-axis limits are (-alim, alim)
#' @param intpts an integer controlling how many points are used to interpolate between values alpha for which
#' the Extended Hypergeometric distribution function is equal exactly to (1+lev)/2 or (1-lev)/2 at some integer co-occurrence value x
#' @param plotCI a vector of up to 4 numbers out of the set {1,2,3,4}, corresponding to which
#' Confidence Intervals (1=CP, 2=Blaker, 3=midQ, 4=midP) should have their coverage probability values plotted simultaneously (with different colors) in the figure
#'
#' @return This function plots a figure of overlaid coverage probabilities with differently colored lines, explained in the figure legend.
#' In addition, it returns "avec", a vector of alpha values at which coverage probabilities are calculated,
#' an xn x 4 x 2 array "arrCI" where xn is the length of "avec" and for each index i in 1:xn the array arrCI[i,,] is the 4x2 matrix of 4 Confidence Interval
#' lower and upper alpha-value endpoints; and "covPrb", the xn x 4 matrix of coverage probabilities for the four CI's at each of the alpha values avec[i] for i =1,...,xn.
#'
#' @author Eric Slud and Kumar Mainali
#'
#' @references to be added
#'
#' @example
#' examples/CovrgPlot_example.R
#'
#' @export

CovrgPlot <-
  function(marg, scal=log(2*marg[3]^2), lev=0.95,
           alim = 4, intpts = 1, plotCI=1:4)  {
    ## plotCI= set of integer indices of 4 CI coverage curves to plot
    #  plot is restricted to alpha in interval in (-alim, alim)
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




    covPrbDF <- data.frame(cbind(avec, covPrb))
    colnames(covPrbDF) <- c("avec","Clopper-Pearson CI","Blaker CI","MidQ CI", "MidP CI")

    # subset the data for the selected CIs
    head(covPrbDF)
    covPrbDF <- covPrbDF[, c(1,(plotCI+1))]
    head(covPrbDF)

    covPrbmelt <- reshape::melt(covPrbDF, id=c("avec"))

    # prepare line plot
    # ----------------------
    lp <- ggplot(covPrbmelt, aes(x=avec, y=value, group=variable)) +
      geom_line(aes(color=variable)) +
      geom_point(aes(color=variable)) +
      geom_hline(yintercept = lev, color="black", linetype="dashed") +
      facet_grid(variable ~ ., switch = "y") +
      theme(legend.position = "none") +
      xlab("Alpha MLE") + ylab(paste0("True coverage probability by the ", lev*100, "% confidence interval"))



    # prepare histogram plot
    # ---------------------------

    tmp <- covPrbmelt
    tmp$abovethr[tmp$value >= lev] <- "yes"
    tmp$belowthr[tmp$value < lev] <- "yes"

    count <- plyr::ddply(tmp, .(variable), summarize,
                         above = length(abovethr[!is.na(abovethr)]),
                         below = length(belowthr[!is.na(belowthr)]))

    count$all <- count$above+count$below
    count$aboveperc <- count$above/count$all*100
    count$belowperc <- count$below/count$all*100

    hp <- ggplot(covPrbmelt, aes(x=value))+
      geom_histogram(aes(fill=variable), color="white")+
      facet_grid(variable ~ .) +
      geom_vline(xintercept = lev, color="black", linetype="dashed") +
      theme(legend.position = "none") +
      xlab("Alpha MLE") + ylab("Count")

    # sometimes, e.g., when requesting just Blaker, all histogram mass can fall above threshold (lev) and histogram print of 0% for less than lev
    # is cut off. for that situation, extend the xlim by a bit to show the printing of 0%
    if(min(covPrbmelt$value) >= lev) {
      hp <- hp + xlim((lev-0.005), max(ggplot_build(hp)$data[[1]]$xmax))
    }



    # find the maximum of the count values plotted on y axis of histogram, to find location to print percentage
    maxcount <- max(ggplot_build(hp)$data[[1]]$count)
    hp <- hp +
      geom_text(data = count, aes(x = lev,  y = maxcount*0.95, label = paste0(round(aboveperc, 1), "%")), vjust = 0.5, hjust = -0.25, size=3) +
      geom_text(data = count, aes(x = lev,  y = maxcount*0.95, label = paste0(round(belowperc, 1), "%")), vjust = 0.5, hjust = 1.25, size=3)

    # add threshold in the first lineplot
    thresholdlegend <- data.frame(cbind(count[,1, drop=FALSE], legendtext=NA))
    thresholdlegend$legendtext[1] <- paste0("true coverage probability of ", round(lev*100, 2), "% threshold")

    lp <- lp +
      geom_text(data = thresholdlegend, aes(x = 0,  y = lev, label = legendtext), hjust = 0.5, vjust = 1.5, size=3)


    list(covPrbDF, covPrbmelt, lev)

    # plot grid of both line plot and histogram
    finalplot <- plot_grid(lp, hp, align = "h", ncol = 2, rel_widths = c(2/3, 1/3))
    print(finalplot)

  }

