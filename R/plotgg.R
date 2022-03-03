#' Heatmap plot of affinity() output
#'
#' This function works on the output of affinity() and gives a heatmap plot for the numeric columns of $all dataframe except
#' the interval columns (median interval and confidence intervals) and confidence level (which is a constant for all pairs in one run of the code)
#'
#' @details   This function is really a ggplot behind the scene where I have taken care of the default value of many arguments for generating a useful plot.
#' It generates a plot for the lower triangle of an NxN square matrix where both row and columns carry the same set of entities,
#' such that all pairwise analyses are shown in the plot (upper triangle is the mirror image of the lower triange
#' and diagonals are the relation to the self which are excluded).
#'
#' The plots can be requested using the column names of $all of the main output of affinity(). The function can include additional arguments
#' either inside plot.gg() or by appending to it with a "+" that is characteristic of ggplot().
#'
#' "legendlimit" centers the legend color to white by default at null expectation in the case of alpha_mle,
#' and negative and positive values stretch between pastel blue and pastel red colors, respectively
#' such that the color spectrum is applied NOT to the range of data but to the same extent of values
#' on both sides of zero, which is max(abs(valrange)) and -(max(abs(valrange))). For example, if alpha_mle
#' ranges between -1.25 to 2.0, then the color spectrum always ranges between -2.0 and 2.0 but the legend can be printed
#' to span between -1.25 and 2.0 with "dataframe" and -2.0 and 2.0 with "balanced".
#'
#' For "entity_1_count_mA", "entity_2_count_mB", and "sites_total_N", there is no natural midpoint.
#' So, "balanced" and "datarange" both use the natural behavior of ggplot in creating the color spectrum that spans between the extremes of the data.
#'
#' For "obs_cooccur_X", and "exp_cooccur" also, there is no natural midpoint.
#' To make the two plots of observed and expected cooccurrence counts comparable visually, one color scale has been applied in these two plots
#' such that the spectrum ranges between the extremes of "obs_cooccur_X", and "exp_cooccur" collectively.
#'
#'
#'
#' @param data the output of affinity()
#' @param variable a column name in $all dataframe in affinity() output; should be a quantitative column;
#' can be one of the following: "entity_1_count_mA", "entity_2_count_mB", "obs_cooccur_X", "sites_total_N",
#' "p_value", "exp_cooccur", "alpha_mle", "jaccard", "sorensen", "simpson"
#' @param legendlimit "datarange" or "balanced"; if "datarange", the legend spans to the range of data,
#' if "balanced, the legend spans in equal magnitude from the center (= white by default) in both directions;
#' note that, irrespective of the value-span in legend, the color spectrum of the plot and legend always goes from the center (=white by default)
#' to two directions in equal magnitude. See details for more information.
#' @param col a set of three colors c("#87beff", "white", "#fd6a6c") by default to represent low, mid and high values in the plot;
#' these colors are applied with ggplot::scale_fill_gradient2()
#' @param show.value a boolean to show ("TRUE") or hide ("FALSE") values in the plot; "TRUE" by default if <=20 entities to compare, otherwise "FALSE" by default
#' @param value.digit the number of digits in values if they are printed; default 2
#' @param text.size the size of values if they are printed; default 2.5
#' @param plot.margin same as ggplot's plot.margin which includes top, right, bottom and left margins as "margin(1,1,5,2, "cm")"
#'
#' @return This function returns a heatmap plot generated with ggplot() behind the scene.
#'
#' @author Kumar Mainali
#'
#' @references
#'
#' @example
#' examples/plot.gg_example.R
#'
#' @export


plotgg <-
  function(data, variable, legendlimit, col=NULL, show.value=NULL, value.digit=NULL, text.size=NULL, plot.margin=NULL, ...) {

    if(!variable %in% colnames(data$all)) {
      stop("the variable does not exist in the data")
    }

    if(variable %in% c("alpha_medianInt", "conf_level", "ci_blaker", "ci_cp", "ci_midQ", "ci_midP")) {
      stop("honestly, we do not like to plot intervals and confidence level... at least for now")
    }

    require(ggplot2)

    gp <- ggplot(data$all, aes(x = entity_1, y = entity_2, fill = get(variable))) +
      geom_tile(color = "gray") + coord_fixed() + labs(fill = variable) +
      ylim(rev(colnames(data$occur_mat)[-1])) + xlim(colnames(data$occur_mat)[-length(colnames(data$occur_mat))]) +
      theme(panel.background = element_blank(), axis.title = element_blank(),
            axis.text.x = element_text(angle = 35, vjust = 0.85, hjust=1), axis.text.y = element_text(vjust = 0.5, hjust = 0.1),
            axis.ticks.length=unit(.25, "cm"))

    # plot with specified margin if supplied
    if(!is.null(plot.margin)) {
      gp <- gp + theme(plot.margin = plot.margin)
    }

    if(is.null(col)) col <- c("#87beff", "white", "#fd6a6c")


    # -------- midpoint and color range for "p_value", "jaccard", "sorensen", "simpson", "alpha_mle" -----------
    # find the legend range for balanced scaling of color
    if(variable %in% c("p_value", "jaccard", "sorensen", "simpson")) {
      upperlimit <- 1
      lowerlimit <- 0
      midpoint <- 0.5
      message(paste0("...for the balanced stretch of color for ", variable, ", we've used the conventional range of 0-1 for the limits and 0.5 for the midpoint..."))
    }
    if(variable %in% c("alpha_mle")) {
      valrange <- range(data$all[[variable]][!is.na(data$all[[variable]])])
      upperlimit <- max(abs(valrange))
      lowerlimit <- -(upperlimit)
      midpoint <- 0
    }

    # plot with legend midpoint and limits calculated above
    if(variable %in% c("p_value", "jaccard", "sorensen", "simpson", "alpha_mle")) {
      if(legendlimit == "datarange") {
        gp <- gp + scale_fill_gradient2(midpoint = midpoint, low = col[1], mid = col[2], high = col[3], space ="Lab", na.value = "grey50")
      } else if(legendlimit == "balanced") {
        gp <- gp + scale_fill_gradient2(midpoint = midpoint, low = col[1], mid = col[2], high = col[3], space ="Lab", na.value = "grey50", limits = c(lowerlimit, upperlimit))
      } else {
        stop("legendlimit should be either datarange or balanced")
      }
    }

    # -------- midpoint and color range for OTHER variables -----------
    # midpoint and limits for other variables such as cooccurrence counts are not defined
    if(variable %in% c("entity_1_count_mA", "entity_2_count_mB", "obs_cooccur_X", "sites_total_N", "exp_cooccur")) {

      # create one color scale for observed and expected cooccurrence so that the images can be compared
      if(variable %in% c("obs_cooccur_X", "exp_cooccur")) {
        merge <- c(data$all$obs_cooccur_X[!is.na(data$all$obs_cooccur_X)], data$all$exp_cooccur[!is.na(data$all$exp_cooccur)])
        midpoint <- mean(merge)
        upperlimit <- max(merge)
        lowerlimit <- min(merge)
      } else {
        midpoint <- mean(data$all[[variable]][!is.na(data$all[[variable]])])
        upperlimit <- max(data$all[[variable]][!is.na(data$all[[variable]])])
        lowerlimit <- min(data$all[[variable]][!is.na(data$all[[variable]])])
      }


      if(legendlimit == "datarange") {
        gp <- gp + scale_fill_gradient2(midpoint = midpoint, low = col[1], mid = col[2], high = col[3], space ="Lab", na.value = "grey50")
      } else if(legendlimit == "balanced") {
        gp <- gp + scale_fill_gradient2(midpoint = midpoint, low = col[1], mid = col[2], high = col[3], space ="Lab", na.value = "grey50", limits = c(lowerlimit, upperlimit))
        message("there is no natural balanced range of values on the two sides of midpoint for the selected variable")
        if(variable %in% c("obs_cooccur_X", "exp_cooccur")) {
          message("however, one color scale has been applied in the plots of obseved and expected cooccurrences so that the colors across the plots can be compared")
        }
        if(variable %in% c("entity_1_count_mA", "entity_2_count_mB", "sites_total_N")) {
          message("limits on legend color == datarange")
        }
      } else {
        stop("legendlimit should be either datarange or balanced")
      }
    }


    # how many digits to print
    if(is.null(value.digit)) value.digit <- 2

    # text size for value printing
    if(is.null(text.size)) text.size <- 2.5

    # print the value by default if total elements to compare are <=20 OR if show.value=T
    if(is.null(show.value) & ncol(data$occur_mat) <= 20 | isTRUE(show.value)) {
      gp <- gp + geom_text(aes(label = round(get(variable), value.digit)), size = text.size)
      message("you can hide the printed values with show.value=F")
      message("use the argument value.digit to change number of digits and text.size to adjust the text size")
    }

    gp

  }
