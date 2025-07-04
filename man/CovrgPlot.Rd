% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CovrgPlot.R
\name{CovrgPlot}
\alias{CovrgPlot}
\title{Coverage probabilities of the confidence intervals, calculated and plotted}
\usage{
CovrgPlot(
  marg,
  scal = log(2 * marg[3]^2),
  lev = 0.95,
  alim = 4,
  intpts = 1,
  plotCI = 1:4
)
}
\arguments{
\item{marg}{a 3-entry integer vector containing (mA,mB,N)}

\item{scal}{an integer parameter (default 2*N^2, capped at 10 within the function) that should be 2 or greater}

\item{lev}{a confidence level, generally somewhere from 0.8 to 0.95 (default 0.95)}

\item{alim}{the absolute value of alpha at which the plotted figure is cutoff; so the x-axis limits are (-alim, alim)}

\item{intpts}{an integer controlling how many points are used to interpolate between values alpha for which
the Extended Hypergeometric distribution function is equal exactly to (1+lev)/2 or (1-lev)/2 at some integer co-occurrence value x}

\item{plotCI}{a vector of up to 4 numbers out of the set \{1,2,3,4\}, corresponding to which
Confidence Intervals (1=CP, 2=Blaker, 3=midQ, 4=midP) should have their coverage probability values plotted simultaneously (with different colors) in the figure}
}
\value{
This function plots a figure of overlaid coverage probabilities with differently colored lines, explained in the figure legend.
In addition, it returns "avec", a vector of alpha values at which coverage probabilities are calculated,
an xn x 4 x 2 array "arrCI" where xn is the length of "avec" and for each index i in 1:xn the array arrCI[i,,] is the 4x2 matrix of 4 Confidence Interval
lower and upper alpha-value endpoints; and "covPrb", the xn x 4 matrix of coverage probabilities for the four CI's at each of the alpha values avec[i] for i =1,...,xn.

\describe{
\item{variable}{Column name of the dataframe produced by reshape::melt}
\item{value}{Column name of the dataframe produced by reshape::melt}
\item{abovethr}{Column name of the dataframe produced by reshape::melt}
\item{belowthr}{Column name of the dataframe produced by reshape::melt}
\item{aboveperc}{Column name of the dataframe produced by plyr::ddply}
\item{belowperc}{Column name of the dataframe produced by plyr::ddply}
\item{legendtext}{Column name of the dataframe}
}
}
\description{
This function calculates and plots the coverage probabilities of the four confidence intervals CI.CP, CI.Blaker, CI.midQ and CI.midP
produced by functions AlphInts() and ML.Alpha(). The coverage plots are useful in showing how the conservative confidence intervals (CI.CP and CI.Blaker)
tend to have above-nominal coverage probabilities versus the intervals CI.midQ and CI.midP that have coverage much closer to the nominal confidence level lev.
}
\details{
This function calculates and plots the coverage probabilities of the confidence intervals CI.CP, CI.Blaker, CI.midQ and CI.midQ produced by function AlphInts(),
for fixed marginals "marg" = (mA,mB,N) and confidence level "lev" as a function of alpha. The plots serve as a useful exhibit to show the range of always conservative
(and sometimes quite conservative) coverage probabilities for CI1, and the much closer-to-nominal coverage probabilities for the recommended Confidence Interval CI2.
The outputted vector "avec" of alpha values is the set of all discontinuities in the slope of the coverage-probability function together with
"intpts: equally spaced values between them. The plotted coverage probabilities are the corresponding exact (not simulated) coverage probabilities
of each of the selected (up to 4) confidence intervals at the alphvec points in the four-column array "covPrb".
}
\examples{
# the script has been enclosed under \donttest{}
# to bypass the CRAN's 5 second limit on example files
# --------------------------------------------------------------

\donttest{

  # True coverage probability by the 95\% Confidence Interval
  CovrgPlot(marg = c(50,70,150), lev = 0.95)

  # If confidence level is not specified, it shows results for 95\% CI
  CovrgPlot(marg = c(50,70,150))

  # True coverage probability by the 90\% Confidence Interval
  CovrgPlot(marg = c(50,70,150), lev = 0.90)

  # Select some of the confidence intervals to plot
  CovrgPlot(marg = c(50,70,150), lev = 0.95, plotCI=1:3)
  CovrgPlot(marg = c(50,70,150), lev = 0.95, plotCI=c(1,3,4))
  CovrgPlot(marg = c(50,70,150), lev = 0.95, plotCI=1)

} #end of \donttest{}
}
\references{
to be added
}
\author{
Eric Slud and Kumar Mainali
}
\keyword{internal}
