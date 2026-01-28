data(finches)
head(finches)

library(ggplot2)


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
  # so that the color spectrum is applied NOT to the range of data
  # but to the same extent of values
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


  # plot only statistically significant tiles (based on p_value)
  # -----------------------------------------------------------
  # sig.only = TRUE masks non-significant tiles (p_value > 0.05) to NA
  plotgg(data = myout, variable = "alpha_mle", legendlimit = "balanced", sig.only = TRUE)

  # you can also supply a stricter p-value cutoff (e.g., 0.01)
  plotgg(data = myout, variable = "alpha_mle", legendlimit = "balanced", sig.only = 0.01)


  # change color of the plot and text
  # ---------------------------------------------
  plotgg(data = myout, variable = "alpha_mle", legendlimit = "balanced")
  plotgg(data = myout, variable = "alpha_mle", legendlimit = "balanced",
         col = c('#99cc33', 'black', '#ff9933'), text.col = "white")
  plotgg(data = myout, variable = "alpha_mle", legendlimit = "balanced",
         col = c('#99cc33', '#ff9933'), text.col = "white")

  plotgg(data = myout, variable = "obs_cooccur_X", legendlimit = "balanced")
  plotgg(data = myout, variable = "obs_cooccur_X", legendlimit = "balanced",
         col = c('black', 'red'), text.col = "white")


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
  # bad ->
  plotgg(data = myout, variable = "alpha_mle", legendlimit = "balanced") +
    ggplot2::theme(axis.text.x = element_text(angle = 90))
  # good ->
  plotgg(data = myout, variable = "alpha_mle", legendlimit = "balanced") +
    ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5))


  # additional elements in the plot
  # ----------------------------------
  # because it is ggplot output, you can use the arguments of ggplot() to make changes

  # add plot title and change legend title
  plotgg(data = myout, variable = "alpha_mle", legendlimit = "balanced") +
    ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    ggplot2::ggtitle("Affinity of island pairs measured with Alpha MLE") +
    ggplot2::labs(fill = 'My Personal\nTitle')



  # show/hide entities that are entirely empty (all-NA tiles)
  # --------------------------------------------------------
  # Here we create an artificial "empty" entity by setting one column to NA.
  # This guarantees that all pairwise comparisons involving that entity have no usable data,
  # so the corresponding tiles become NA for variables such as alpha_mle.

  finches2 <- as.matrix(finches)
  storage.mode(finches2) <- "numeric"
  finches2[, 3] <- NA_real_  # make the first entity entirely missing (choose any column)
  myout2 <- affinity(data = finches2, row.or.col = "col")

  # Default behavior: drop.empty = TRUE (empty entity removed from the axes)
  plotgg(data = myout2, variable = "alpha_mle", legendlimit = "balanced")

  # Keep empty entities (legacy/full grid): shows the empty row/column
  plotgg(data = myout2, variable = "alpha_mle", legendlimit = "balanced", drop.empty = FALSE)

  # keep empty entities even after masking (shows rows/columns with all-NA tiles)
  plotgg(data = myout2, variable = "alpha_mle", legendlimit = "balanced",
         sig.only = TRUE, drop.empty = FALSE)


  # automatic suppression of numeric values on tiles
  # ------------------------------------------------
  # By default, numeric values are printed on tiles only when the number of
  # plotted entities is reasonably small (<= 20). This avoids severe visual
  # clutter when the heatmap becomes large.

  finches_big <- finches

  # duplicate columns to artificially inflate the number of entities
  finches_big <- cbind(finches_big, finches_big[, 1:5])
  colnames(finches_big)[(ncol(finches) + 1):ncol(finches_big)] <-
    paste0(colnames(finches)[1:5], "_dup")

  myout_big <- affinity(data = finches_big, row.or.col = "col")

  # Numeric values are NOT printed because the number of entities exceeds 20
  plotgg(data = myout_big, variable = "alpha_mle", legendlimit = "balanced")

  # To force printing numeric values despite the large number of entities:
  plotgg(data = myout_big, variable = "alpha_mle", legendlimit = "balanced", show.value = TRUE)

} #end of \donttest{}


