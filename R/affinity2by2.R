#' Maximum likelihood estimate and intervals of alpha, null expectation, p-value and traditional indices from a 2x2 table
#'
#' This function uses ML.Alpha() and supplements to the outcome with traditional indices of Jaccard, Sorenson, and Simpson.
#' ML.Alpha() calculates the maximum likelihood estimate and other quantities computed in AlphInts(),
#' for the log-odds parameter alpha in the Extended Hypergeometric distribution with fixed margins (mA,mB) and
#' table-total N, which is the "log-affinity" index of co-occurrence championed in a paper by Mainali et al. (2022) as an index of co-occurrence-based similarity.
#'
#' @details   See the details of ML.Alpha(). In addition to the output of ML.Alpha, this function also computes Jaccard, Sorenson and Simpson indices.
#'
#' @param x integer co-occurrence count that should properly fall within the closed interval \[max(0,mA+mB-N), min(mA,mB)\]
#' @param marg a 3-entry integer vector (mA,mB,N) consisting of the first row and column totals and the table total for a 2x2 contingency table
#' @param bound a boolean parameter which when TRUE replaces the MLE of "+/-Infinity", applicable when x is respectively at the upper extreme min(mA,mB)
#' or the lower extreme max(mA+mB-N,0) of its possible range, by a finite value with absolute value upper-bounding the value of
#' MLEs attainable for values of x not equal to its extremes
#' @param scal an integer parameter (default 2*N^2, capped at 10 within the function) that should be 2 or greater
#' @param lev a confidence level, generally somewhere from 0.8 to 0.95  (default 0.95)
#' @param pvalType a character string telling what kind of p-value to calculate. ‘Blaker’ or “midP’.
#' If ‘pvalType=Blaker” (the default value), the p-value is calculated according to "Acceptability" function of Blaker (2000).
#' If ‘pvalType=midP’, the p-value is calculated using the same idea as the midP confidence interval.
#'
#' @return This function returns maximum likelihood estimate of alpha, the interval-endpoints of alpha values for which x is a median,
#' four confidence intervals for alpha, described in detail under documentation for AlphInts(), and traditional indices of Jaccard, Sorenson and Simpson.
#' In addition there are two output list-components for the null-distribution expected co-occurrence count and the p-value
#' for the test of the null hypothesis alpha=0, calculated as in AlphInts.
#'
#' @author Kumar Mainali and Eric Slud
#'
#' @references
#' Fog, A. (2015), BiasedUrn: Biased Urn Model Distributions. R package version 1.07.
#'
#' Harkness, W. (1965), “Properties of the extended hypergeometric distribution“, Annals of Mathematical Statistics, 36, 938-945.
#'
#' Mainali, K., Slud, E., Singer, M. and Fagan, W. (2022), "A better index for analysis of co-occurrence and similarity", Science Advances, to appear.

#'
#' @example
#' inst/examples/affinity2by2_example.R
#'
#' @export


affinity2by2 <-
  function(x, marg, bound=TRUE, scal=log(2*marg[3]^2), lev=0.95, pvalType="Blaker") {

    mlout <- ML.Alpha(x=x, marg=marg, bound=bound, scal=scal, lev=lev, pvalType=pvalType)

    newlist <- list(
         jaccard = as.numeric(x/(marg[1]+marg[2]-x)),
         sorensen = as.numeric(2*x/(marg[1]+marg[2])),
         simpson = as.numeric(x/min(marg[1], marg[2]))
    )

    c(mlout, newlist)
  }


