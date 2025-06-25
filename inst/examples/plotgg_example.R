if(requireNamespace("cooccur", quietly = TRUE)) {

  data(finches, package = "cooccur")
  head(finches)
  require(ggplot2)


  # the remainder of the script has been enclosed under \donttest{}
  # to bypass the CRAN's 5 second limit on example files
  # --------------------------------------------------------------

  \donttest{

    # plotting various variables
    # ---------------------------------------------
    # compute alpha and other quantities for island-pair affinity (beta diversity)
    # the square matrices are not used for plotting
    myout <- affinity(data = finches, row.or.col = "col")
    # myout

    plotgg(data = myout, variable = "alpha_mle", legendlimit = "datarange")
    # in the example above, null expectation of the alpha_mle (=0) has white color,
    # and negative and positive values stretch between "#87beff" and "#fd6a6c", respectively
    # so that the color spectrum is applied NOT to the range of data but to the same extent of values
    # on both sides of zero, which is max(abs(valrange)) and -(max(abs(valrange))).
    # however, the legend can be printed to show the extent of data with "datarange"
    # or the entire spectrum where the color is applied with "balanced".
    plotgg(data = myout, variable = "alpha_mle", legendlimit = "balanced")
    # notice that the two plots above are identical but the legend has
    # different range with the same color scale.


    plotgg(data = myout, variable = "sorensen", legendlimit = "balanced")
    plotgg(data = myout, variable = "jaccard", legendlimit = "balanced")

    # in the case of observed and expected cooccurrences, one color scale is applied for both plots
    # so that the shades of color across plots can be visually compared
    plotgg(data = myout, variable = "exp_cooccur", legendlimit = "datarange")
    plotgg(data = myout, variable = "exp_cooccur", legendlimit = "balanced")
    plotgg(data = myout, variable = "obs_cooccur_X", legendlimit = "balanced")

    plotgg(data = myout, variable = "entity_1_count_mA", legendlimit = "datarange")
    plotgg(data = myout, variable = "entity_2_count_mB", legendlimit = "datarange")
    plotgg(data = myout, variable = "total_N", legendlimit = "datarange")
    # for "entity_1_count_mA", "entity_2_count_mB", "sites_total_N",
    # if legendlimit is set to "balanced", it will be changed to "datarange"
    plotgg(data = myout, variable = "entity_2_count_mB", legendlimit = "balanced")


    # change color of the plot
    # ---------------------------------------------
    plotgg(data = myout, variable = "alpha_mle", legendlimit = "balanced")
    plotgg(data = myout, variable = "alpha_mle", legendlimit = "balanced",
           col = c('#99cc33', 'white', '#ff9933'))

    plotgg(data = myout, variable = "obs_cooccur_X", legendlimit = "balanced")
    plotgg(data = myout, variable = "obs_cooccur_X", legendlimit = "balanced",
           col = c('#99cc33', 'white', '#ff9933'))


    # change the characteristics of text printed in the plot
    # ------------------------------------------------------
    plotgg(data = myout, variable = "alpha_mle", legendlimit = "balanced")

    # change the number of digits; the default is 2
    plotgg(data = myout, variable = "alpha_mle", legendlimit = "balanced", value.digit = 3)

    # make the fonts bigger; the default is 2.5
    plotgg(data = myout, variable = "alpha_mle", legendlimit = "balanced", text.size = 3.5)



    # hide values from the plot
    plotgg(data = myout, variable = "alpha_mle", legendlimit = "balanced", show.value = FALSE)


    # increase or decrease margin
    # ---------------------------------------------

    myout <- affinity(data = finches, row.or.col = "row")
    # myout

    plotgg(data = myout, variable = "alpha_mle", legendlimit = "balanced")
    plotgg(data = myout, variable = "alpha_mle", legendlimit = "balanced",
           plot.margin = ggplot2::margin(1,1,5,2, "cm"))


    # change angle of x-axis tick label; the default is 35 degrees
    # ------------------------------------------------------------
    plotgg(data = myout, variable = "alpha_mle", legendlimit = "balanced")
    plotgg(data = myout, variable = "alpha_mle", legendlimit = "balanced") +
      ggplot2::theme(axis.text.x = element_text(angle = 45))

    # to change to 90 degrees, adjust vjust
    plotgg(data = myout, variable = "alpha_mle", legendlimit = "balanced") +
      ggplot2::theme(axis.text.x = element_text(angle = 90))
    plotgg(data = myout, variable = "alpha_mle", legendlimit = "balanced") +
      ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5))


    # additional elements in the plot
    # ----------------------------------
    # because it is ggplot output, you can use the arguments of ggplot() to make changes

    # add plot title and change legend title
    plotgg(data = myout, variable = "alpha_mle", legendlimit = "balanced") +
      ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
      ggplot2::ggtitle("Affinity of island pairs measured with Alpha MLE") +
      ggplot2::labs(fill = 'Alpha MLE')


    # an example of much bigger dataset
    # -------------------------------------
    data("beetles", package = "cooccur")
    dim(beetles)
    myout <- affinity(data = beetles, row.or.col = "row")

    plotgg(data = myout, variable = "alpha_mle", legendlimit = "balanced") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
    plotgg(data = myout, variable = "alpha_mle", legendlimit = "balanced",
           show.value = TRUE, text.size = 1.5, value.digit = 1) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

  } #end of \donttest{}


} else {

  message("The cooccur package is not installed.
          You can install it with install.packages('cooccur').")

}
