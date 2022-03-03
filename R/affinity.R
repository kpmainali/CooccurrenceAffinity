#' Computes alpha, probability, expected co-occurence, median interval, various confidence intervals, other indices of affinity, etc.
#'
#' This is the principal function of "CooccurrenceAffinity" package that analyzes occurrence or abundance data (e.g., species by site)
#' using other functions of this package and returns several quantities in one main dataframe and (optionally) up to 11 square matrices.
#' This function processes data using dataprep() function and then feeds the data to analytical pipeline which includes ML.Alpha() and AlphInts().
#' The outputs of the function in the $all dataframe include the following: \cr
#' * alpha_mle: maximum likelihood estimate of the log-odds parameter alpha in the Extended Hypergeometric distribution with fixed margins (mA,mB) and
#' table-total N, which is the "log-affinity" index of co-occurrence championed in a paper by
#' Mainali et al. (2022) as an index of co-occurrence-based similarity; computed in ML.Alpha() \cr
#' * exp_cooccur: expected co-occurrence count under the null (hypergeometric, corresponding to alpha=0) distribution; computed as ML.Alpha()$Null.Exp \cr
#' * p_value: the commonly reported P-value of the observed co-occurrences; computed by AlphInts()$pval \cr
#' * alpha_medianInt: the interval of alpha values compatible with x as median for the Extended Hypergeometric distribution (Harkness 1965) \cr
#' with fixed margins and alpha; computed in AlphInts() as $MedianIntrvl \cr
#' * conf_level: confidence level for estimating the various types of confidence intervals \cr
#' * ci_: fout types of confidence intervals (see details below) \cr
#' * jaccard: Jaccard index \cr
#' * sorensen: Sørensen-Dice index \cr
#' * simpson: Simpson index \cr
#'
#'
#'
#'
#'
#' @details   This function calculates "alpha_mle", which is the maximum likelihood estimate of the log-odds paramater alpha within the Extended Hypergeometric distribition (Harkness 1965)
#' based on the count x and fixed table margins (mA,mB) and total N, which is the "affinity" index of co-occurrence championed in the paper of Mainali et al. (2022)
#' as an index of cooccurrence-based similarity.
#'
#' This function calculates five intervals, three of them using EHypQuInt, one using EHypMidP, and one using AcceptAffCI.
#' First ("alpha_medianInt") is the interval of alpha values compatible with x as median for the Extended Hypergeometric distribution (Harkness 1965)
#' with fixed margins and alpha. Computed as AlphInts()$MedianIntrvl.
#' This interval quantifies the underlying discreteness of the Extended Hypergeometric and its impact on the estimation of alpha.
#' MedianIntrvl is an interval that will contain the MLE alpha-hat, and the mid-point of that interval is another reasonable estimator of alpha from the data.
#'
#' There are four confidence intervals computed in AlphInts(), called ci_cp, ci_blaker, ci_midQ, ci_midP, matching name of the outputs in standalone outputs of AlphInts(),
#' except the differences in capital/small letters. The boolean "bound" parameter is an option to prevent the intervals containing alpha-estimates
#' to extend to plus or minus infinity, based on a Bayesian argument.
#' The bound substituted for the Infinite endpoints is provably larger than the largest value the MLE can take whenever x avoids the enpoints max(mA+mB-N,0) and min(mA,mB)
#' of its logical range. The recommended confidence interval for alpha is CI.Blaker if a reliably conservative (over-large) coverage probability is desired, and CI.midP otherwise.
#'
#'
#' "ci_cp", computed as AlphInts()$CI.CP is an "exact" conservative test-based 2-sided confidence interval (analogous to the Clopper-Pearson (1934)
#' confidence interval for unknown binomial proportion) for alpha based on data (x,mA,mB,N)
#'
#' "ci_blaker", computed as AlphInts()$CI.Blaker is the Acceptability Confidence Interval of Blaker (2000, Theorem 1)
#' which is a better confidence interval than the CP-type interval "CI.CP" in the sense of being contained within "CI.CP"
#' but still probably conservative, i.e., with coverage probability always at least as large as the nominal level.
#'
#' "ci_midQ", computed as AlphInts()$CI.midQ has the endpoints obtained as the midpoints of quantile intervals
#' respectively to the (1+lev)/2 and (1-lev)/2 quantiles of the Extended Hypergeometric distribution.
#'
#' "ci_midP", computed as AlphInts()$CI.midQ, behaves very similarly to "CI.midQ" and is defined by the midP approach analogous
#' to the midP confidence interval for binomial proportions (Agresti 2013, p.605), and is calculated from EHypMidP().
#'
#' The recommended (slightly conservative) confidence interval is CI.Blaker, while the very similar intervals CI.midQ and CI.midP have
#' coverage generally closer than CI.CP or CI.Blaker to the nominal level of coverage, at the cost of occasionally under-covering
#' by as much as 0.04 or 0.05 for confidence levels 0.90 or 0.95. The comparison among intervals, and different possible goals that CIs of
#' conservative or close-to-nominal coverage can serve, are similar to those compared by  Brown et al. (2001) for interval estimation of an unknown binomial proportion.
#'
#' "p_value" is the two-sided p-value for the equal-tailed test of the null hypothesis alpha=0. This p-value is calculated when pval="Blaker"
#' according to Blaker's (2000) "Acceptability" function; if the input parameter pvalType of AlphInts() is anything else,
#' the p-value is calculated using the same idea as the midP confidence interval.
#'
#' ADDITIONAL ARGUMENTS can be supplied from ML.Alpha() and AlphInts().
#'
#'
#'
#' @param data occurrence matrix (binary or abundance) in matrix or dataframe format
#' @param row.or.col specify if the pairs of rows or columns are analyzed for affinity. 'row' or 'column'.
#' @param which.row.or.col a vector of name or the number of row/column if a subset of the data is intended to be analyzed; optional argument with default of all rows/columns.
#' @param datatype specify if the datatype is 'abundance' or 'binary'; optional argument with default 'binary'.
#' @param threshold cutoff for converting an abundance data to binary; needed if datatype is 'abundance'
#' @param class0.rule 'less.or.equal' or 'less'. 'less.or.equal' converts a threshold or lower values to zero and all the others to 1. 'less' converts a threshold and higher values to 1.
#' @param sigPval acceptable rate of false positives or probability of rejecting the null hypothesis when it is true, commonly known as alpha in hypothesis testing
#' @param sigdigit the number of decimals for rouding alpha mle, its intervals, expected cooccurence under the null, jaccard, sorensen and simpson indices
#' @param squarematrix a vector of quantities so that a square matrix for each of them is generated on the top of the main long matrix of all outputs
#'
#'
#' @return This function returns one main long dataframe ($all) with various outputs in columns (a list given under "details") for each of the pairs of the entities in row.
#' This function also outputs optionally upto 11 square matrices of NxN entities.
#'
#' @author Kumar Mainali
#'
#' @references
#'
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
#' Mainali, K., Slud, E., Singer, M. and Fagan, W. (2022), "A better index for analysis of co-occurrence and similarity", Science Advances, to appear.

#'
#' @example
#' examples/affinity_example.R
#'
#' @export


affinity <-

  function(data, row.or.col, which.row.or.col=NULL, datatype=NULL, threshold=NULL, class0.rule=NULL, sigPval=NULL, sigdigit=NULL, squarematrix=NULL, ...) {

    # use another function dataprep() to process data and throw errors if needed
    occur.mat <- dataprep(data = data, row.or.col = row.or.col, which.row.or.col=which.row.or.col, datatype=datatype, threshold=threshold, class0.rule=class0.rule)


    # if a user asks for a squarematrix that does not exists in the list, throw an error and stop
    if(!is.null(squarematrix)) {
      indexvec <- c("alpha_mle", "alpha_mle_sig", "p_value", "cooccur.null", "cooccur.obs", "jaccard", "jaccard_sig", "sorensen", "sorensen_sig", "simpson", "simpson_sig", "all")
      if(sum(is.na(match(squarematrix, indexvec))) > 0) {
        squarematrix[!squarematrix %in% indexvec]
        stop(paste(c("squarematrix cannot take the following value(s) you supplied: ", squarematrix[!squarematrix %in% indexvec]), collapse=" "))
      }
    }


    # optional square multiple output dataframes - create empty matrices
    # ------------------------------------------------------------------
    # make a square matrix for inserting estimates of every pair of columns

    if(!is.null(squarematrix)) {

      # if 'all' added with or without some other selections, all of the square matrix will be computed
      if("all" %in% squarematrix) {
        squarematrix <- c("alpha_mle", "alpha_mle_sig", "p_value", "cooccur.null", "cooccur.obs", "jaccard", "jaccard_sig", "sorensen", "sorensen_sig", "simpson", "simpson_sig")
      }

      dummy.df <- data.frame(matrix(ncol = ncol(occur.mat), nrow = ncol(occur.mat)))
      colnames(dummy.df) <- rownames(dummy.df) <- colnames(occur.mat)

      for(sqmat in squarematrix) {
        new.df <- dummy.df
        assign(paste0(sqmat, ".df"), new.df)
      }

    }

    # main single dataframe with output of all analysis - create empty matrices
    # -------------------------------------------------------------------------
    # this output dataframe will have rows = number of cells in top right half of square matrix of columns that are being compared
    sum(1:(ncol(occur.mat)-1))   # the number of cells in the upper or lower half of triangle in NxN square matrix, excluding the diagonal
    output.long <- data.frame(matrix(ncol = 19, nrow = sum(1:(ncol(occur.mat)-1))))
    colnames(output.long) <- c("entity_1", "entity_2", "entity_1_count_mA", "entity_2_count_mB", "obs_cooccur_X", "sites_total_N",
                               "p_value", "exp_cooccur", "alpha_mle", "alpha_medianInt",
                               "conf_level", "ci_blaker", "ci_cp", "ci_midQ", "ci_midP",
                               "jaccard", "sorensen", "simpson", "errornote")



    # perform pairwise analysis and fill in the empty matrices created above
    # ----------------------------------------------------------------------
    for(col_main in 1:(ncol(occur.mat)-1)) {
      for(col_pair in (col_main+1):ncol(occur.mat)) {

        sub <- occur.mat[, c(col_main, col_pair)]
        colnames(sub) <- c("entity_main", "entity_pair")
        # remove all rows with missing data at least in one of the two columns to compare
        sub <- na.omit(sub)

        sub$cooccur <- sub$entity_main + sub$entity_pair
        sub

        # entities count
        mA <- sub$entity_main
        mB <- sub$entity_pair

        # co-occurrence calculation
        t <- as.data.frame(table(sub$cooccur))
        # sometimes there is no co-occurrence and no value of 2
        if(match("2", t$Var1, 0) != 0) {
          X <- t$Freq[t$Var1==2]
        } else {
          X = 0
        }
        X

        serialnum <- (col_main-1)*ncol(occur.mat) + col_pair - sum(1:col_main)

        # how many digits to round to a number in the output
        if(is.null(sigdigit)) sigdigit <- 3


        # ----------- SAVE ALL OUTPUTS TO A SINGLE LONG DATAFRAME --------------

        output.long$entity_1[serialnum] <- colnames(occur.mat)[col_main]
        output.long$entity_2[serialnum] <- colnames(occur.mat)[col_pair]
        output.long$entity_1_count_mA[serialnum] <- sum(mA)
        output.long$entity_2_count_mB[serialnum] <- sum(mB)
        output.long$obs_cooccur_X[serialnum] <- X
        output.long$sites_total_N[serialnum] <- nrow(sub)

        myjaccard <- round( X/(sum(mA)+sum(mB)-X), sigdigit )
        mysorensen <- round( 2*X/(sum(mB)+sum(mB)), sigdigit )
        mysimpson <- round( X/min(sum(mA), sum(mB)), sigdigit )

        output.long$jaccard[serialnum] <- myjaccard
        output.long$sorensen[serialnum] <- mysorensen
        output.long$simpson[serialnum] <- mysimpson


        # ---------- SAVE OUTPUTS TO REQUESTED SQUARE MATRICES -------------

        # Jaccard
        if(exists("jaccard.df")) jaccard.df[col_pair, col_main] <- myjaccard

        # Sorensen
        if(exists("sorensen.df")) sorensen.df[col_pair, col_main] <- mysorensen

        # Simpson
        if(exists("simpson.df")) simpson.df[col_pair, col_main] <- mysimpson


        # ---------- Additional Inputs only if Distribution is NOT Degenerate -------------

        # condition to bypass a failed alpha
        if(ML.Alpha(x = X, marg = c(sum(mA), sum(mB), nrow(sub)))[1] == "Degenerate co-occurrence distribution!") {
          output.long$errornote[serialnum] <- "Degenerate co-occurrence distribution!"
          next
        } else {

          mle <- ML.Alpha(x = X, marg = c(sum(mA), sum(mB), nrow(sub)), ...)


          # ----------- CONTINUE LONG DATAFRAME OUTPUT --------------

          # round P-value to 4 digits if it is >0.0001. If not, save it in scientific notation with 4 digits
          if(mle$pval > 0.0001) {
            output.long$p_value[serialnum] <- round(mle$pval, 4)
          }
          if(mle$pval <= 0.0001) {
            output.long$p_value[serialnum] <- formatC(mle$pval, format = "e", digits = 4)
          }

          output.long$exp_cooccur[serialnum] <- round(mle$Null.Exp, sigdigit)
          output.long$alpha_mle[serialnum] <- round(mle$est, sigdigit)

          output.long$alpha_medianInt[serialnum] <- paste0("[", round(mle$MedianIntrv[1], sigdigit), ", ", round(mle$MedianIntrv[2], sigdigit), "]")

          output.long$conf_level[serialnum] <- mle$lev
          output.long$ci_blaker[serialnum] <- paste0("[", round(mle$CI.Blaker[1], sigdigit), ", ", round(mle$CI.Blaker[2], sigdigit), "]")
          output.long$ci_cp[serialnum] <- paste0("[", round(mle$CI.CP[1], sigdigit), ", ", round(mle$CI.CP[2], sigdigit), "]")
          output.long$ci_midQ[serialnum] <- paste0("[", round(mle$CI.midQ[1], sigdigit), ", ", round(mle$CI.midQ[2], sigdigit), "]")
          output.long$ci_midP[serialnum] <- paste0("[", round(mle$CI.midP[1], sigdigit), ", ", round(mle$CI.midP[2], sigdigit), "]")


          # ---------- CONTINUE SQUARE OUTPUT -------------

          if(exists("alpha_mle.df")) alpha_mle.df[col_pair, col_main] <- round(mle$est, sigdigit)
          if(exists("cooccur.null.df")) cooccur.null.df[col_pair, col_main] <- round(mle$Null.Exp, sigdigit)
          if(exists("cooccur.obs.df")) cooccur.obs.df[col_pair, col_main] <- X

          if(exists("p_value.df")) {
            if(mle$pval > 0.0001) {
              p_value.df[col_pair, col_main] <- round(mle$pval, 4)
            }
            if(mle$pval <= 0.0001) {
              p_value.df[col_pair, col_main] <- formatC(mle$pval, format = "e", digits = 4)
            }
          }


          # significant Alpha, Jaccard, Sonrenson and Simpson - if not significant, the cell keeps carrying NA
          if(is.null(sigPval)) sigPval <- 0.05

          if(mle$pval <= sigPval) {
            if(exists("alpha_mle_sig.df")) alpha_mle_sig.df[col_pair, col_main] <- round(mle$est, sigdigit)
            if(exists("jaccard_sig.df")) jaccard_sig.df[col_pair, col_main] <- myjaccard
            if(exists("sorensen_sig.df")) sorensen_sig.df[col_pair, col_main] <- mysorensen
            if(exists("simpson_sig.df")) simpson_sig.df[col_pair, col_main] <- mysimpson
          }

        }

      }
    }


    # create a list of all objects that a user wants, including the main dataframe that we provide without asking

    # main single output of everything
    finalout <- list(
      all = output.long,
      occur_mat = occur.mat
    )

    # optional square matrices
    if(!is.null(squarematrix)) {
      for(sqmat in squarematrix) {
        finalout[[sqmat]] <- get(paste0(sqmat, ".df"))
      }
    }

    message("~~~~~~~~~~ printing head of all elements of the output list ~~~~~~~~~~")
    print(sapply(finalout, head))

    # to do - the complete data is not printed if the output is assigned to a variable. If not assigned, it still prints the complete output.
    # for now, turn off the message and explore later.
    # message("*************** the returned output has NOT been printed completely; if you have not already done so, please assign the return output of affinity() to a variable and inspect the variable ***************")
    invisible(return(finalout))

  }
