#' Heatmap plot of affinity() output
#'
#' This function works on the output of \code{\link{affinity}} and uses
#' \code{ggplot2::ggplot()} to generate a heatmap for numeric columns of the
#' \code{$all} dataframe, excluding interval columns (median interval and
#' confidence intervals) and the confidence level (which is constant across
#' pairs in a single run).
#'
#' @details
#' This function is a wrapper around \code{ggplot2} with carefully chosen
#' defaults to generate an interpretable heatmap of pairwise associations.
#' The plot shows the lower triangle of an \eqn{N \times N} matrix (diagonal
#' excluded), where both rows and columns represent the same set of entities.
#' The upper triangle is omitted because it is a mirror image of the lower
#' triangle.
#'
#' By default (\code{drop.empty = TRUE}), entities whose values are entirely
#' \code{NA} for the selected \code{variable} are removed from both axes. This
#' avoids plotting empty rows and columns when an entity has no usable values
#' (e.g., due to degenerate distributions or missing data). Set
#' \code{drop.empty = FALSE} to retain all entities and reproduce the full grid,
#' including empty rows or columns.
#'
#' If \code{sig.only} is enabled, values of the selected \code{variable} are
#' masked to \code{NA} wherever \code{p_value} exceeds the specified cutoff, so
#' only statistically significant tiles are shown. Use \code{sig.only = TRUE}
#' to apply the default cutoff (0.05), or supply a numeric cutoff (e.g.,
#' \code{sig.only = 0.01}). This option cannot be used when
#' \code{variable = "p_value"}.
#'
#' The plot can be requested using column names from the \code{$all} dataframe
#' returned by \code{affinity}. Additional \code{ggplot2} layers or theme
#' modifications can be added by appending them with \code{+}, as in standard
#' \code{ggplot2} usage.
#'
#' The \code{legendlimit} argument controls how the color scale is defined.
#' For \code{alpha_mle}, the default midpoint is 0 (null expectation), and the
#' color scale can be either data-driven (\code{"datarange"}) or symmetrically
#' balanced around zero (\code{"balanced"}), using the maximum absolute value
#' observed. For indices bounded in \eqn{[0,1]} (\code{p_value}, \code{jaccard},
#' \code{sorensen}, \code{simpson}), the balanced scale uses limits \eqn{[0,1]} with
#' midpoint 0.5. For count-based variables, no natural midpoint exists; the
#' color scale spans the observed range. For \code{obs_cooccur_X} and
#' \code{exp_cooccur}, a shared color scale is applied so the two plots are
#' visually comparable.
#'
#' @param data Output list returned by \code{\link{affinity}}.
#' @param variable Name of a numeric column in \code{data$all} to plot.
#' @param legendlimit Either \code{"datarange"} or \code{"balanced"}.
#' @param col Vector of three colors (low, mid, high) passed to
#'   \code{ggplot2::scale_fill_gradient2()}.
#' @param show.value Logical; if \code{TRUE}, values are printed on tiles. If
#'   \code{NULL}, values are printed automatically when the number of plotted
#'   entities is \eqn{\le 20}.
#' @param value.digit Number of digits used when printing values; default is 2.
#' @param text.size Size of printed values; default is 2.5.
#' @param plot.margin Plot margin passed to \code{ggplot2::theme()}.
#' @param drop.empty Logical; if \code{TRUE} (default), entities whose values are
#'   all \code{NA} for the selected \code{variable} are removed from the plot.
#'   Set to \code{FALSE} to keep all entities.
#'   #' @param sig.only Logical or numeric. If \code{FALSE} (default), plot all values.
#'   If \code{TRUE}, show only significant tiles using cutoff 0.05 by masking
#'   non-significant values to \code{NA}. If numeric, the value is used as the
#'   p-value cutoff (e.g., \code{sig.only = 0.01}). Requires a \code{p_value}
#'   column in \code{data$all}. Not allowed when \code{variable = "p_value"}.
#' @param sig.only Logical or numeric. If \code{FALSE} (default), all values are
#'   plotted. If \code{TRUE}, only statistically significant tiles are shown
#'   using a default cutoff of 0.05 by masking non-significant values to
#'   \code{NA}. If numeric, the value is used as the p-value cutoff (e.g.,
#'   \code{sig.only = 0.01}). Requires a \code{p_value} column in
#'   \code{data$all} and is not allowed when \code{variable = "p_value"}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A heatmap plot generated with \code{ggplot2}.
#'
#' @author Kumar Mainali
#'
#' @example
#' inst/examples/plotgg_example.R
#'
#' @export





plotgg <-
  function(data, variable, legendlimit, col=NULL, show.value=NULL, value.digit=NULL, text.size=NULL, plot.margin=NULL, drop.empty=TRUE, sig.only=FALSE, ...) {

    if(!variable %in% colnames(data$all)) {
      stop("the variable does not exist in the data")
    }

    if(variable %in% c("alpha_medianInt", "conf_level", "ci_blaker", "ci_cp", "ci_midQ", "ci_midP")) {
      stop("honestly, we do not like to plot intervals and confidence level... at least for now")
    }

    # ---------------- drop entities with all-NA values for selected variable ----------------
    plot_data <- data$all

    # keep an unmasked copy for legend scaling (used when sig.only is active)
    legend_data <- plot_data


    # ---------------- optionally mask values by significance ----------------
    if (!identical(sig.only, FALSE)) {

      if (!"p_value" %in% colnames(plot_data)) {
        stop("sig.only requires a p_value column in data$all")
      }

      if (variable == "p_value") {
        stop("sig.only cannot be used when variable = 'p_value'")
      }

      # interpret sig.only: TRUE -> 0.05, numeric -> user cutoff
      if (isTRUE(sig.only)) {
        p_cut <- 0.05
      } else if (is.numeric(sig.only) && length(sig.only) == 1 && is.finite(sig.only)) {
        p_cut <- sig.only
      } else {
        stop("sig.only must be FALSE, TRUE, or a single numeric p-value cutoff (e.g., 0.05)")
      }

      # mask the selected variable where p_value is not significant
      plot_data[[variable]][plot_data$p_value > p_cut] <- NA
    }



    # Full ordered entity list as used by current xlim/ylim logic
    entities_full <- colnames(data$occur_mat)

    if (isTRUE(drop.empty)) {

      # rows where the plotted variable exists (non-NA)
      non_na <- !is.na(plot_data[[variable]])

      # entities that participate in at least one non-NA tile
      keep <- unique(c(plot_data$entity_1[non_na], plot_data$entity_2[non_na]))
      keep <- keep[!is.na(keep)]

      # preserve original order from occur_mat (critical for consistent axes/triangle)
      keep <- entities_full[entities_full %in% keep]

      # need at least 2 entities to form pairs
      if (length(keep) < 2) {
        stop(paste0("drop.empty=TRUE removed all entities for variable='", variable,
                    "' (no non-NA values to plot)."))
      }

      # filter plot data to kept entities only
      plot_data <- plot_data[plot_data$entity_1 %in% keep & plot_data$entity_2 %in% keep, , drop = FALSE]
    } else {
      keep <- entities_full
    }

    # Axis vectors for the lower-triangle layout you already use
    x_entities <- keep[-length(keep)]
    y_entities <- rev(keep[-1])



    entity_1 <- entity_2 <- labs <- ylim <- unit <- NULL
    # gp <- ggplot2::ggplot(plot_data, ggplot2::aes(x = entity_1, y = entity_2, fill = get(variable))) +
    gp <- ggplot2::ggplot(plot_data, ggplot2::aes(x = entity_1, y = entity_2, fill = .data[[variable]])) +
      ggplot2::geom_tile(color = "gray") + ggplot2::coord_fixed() + ggplot2::labs(fill = variable) +
      ggplot2::ylim(y_entities) + ggplot2::xlim(x_entities) +
      ggplot2::theme(panel.background = element_blank(), axis.title = element_blank(),
                     axis.text.x = element_text(angle = 35, vjust = 0.85, hjust=1),
                     axis.text.y = element_text(vjust = 0.5, hjust = 0.1),
                     axis.ticks.length=ggplot2::unit(.25, "cm"))


    # plot with specified margin if supplied
    if(!is.null(plot.margin)) {
      gp <- gp + ggplot2::theme(plot.margin = plot.margin)
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
    # if(variable %in% c("alpha_mle")) {
    #   valrange <- range(plot_data[[variable]][!is.na(plot_data[[variable]])])
    #   upperlimit <- max(abs(valrange))
    #   lowerlimit <- -(upperlimit)
    #   midpoint <- 0
    # }
    if(variable %in% c("alpha_mle")) {

      # legend scaling: use original data if sig.only is active, otherwise use plotted data
      scale_source <- if (!identical(sig.only, FALSE)) legend_data else plot_data

      valrange <- range(scale_source[[variable]][!is.na(scale_source[[variable]])])
      upperlimit <- max(abs(valrange))
      lowerlimit <- -(upperlimit)
      midpoint <- 0
    }


    # plot with legend midpoint and limits calculated above
    # if(variable %in% c("p_value", "jaccard", "sorensen", "simpson", "alpha_mle")) {
    #   if(legendlimit == "datarange") {
    #     gp <- gp + ggplot2::scale_fill_gradient2(midpoint = midpoint, low = col[1],
    #                                              mid = col[2], high = col[3], space ="Lab",
    #                                              na.value = "grey50")
    #   } else if(legendlimit == "balanced") {
    #     gp <- gp + ggplot2::scale_fill_gradient2(midpoint = midpoint, low = col[1],
    #                                              mid = col[2], high = col[3], space ="Lab",
    #                                              na.value = "grey50", limits = c(lowerlimit, upperlimit))
    #   } else {
    #     stop("legendlimit should be either datarange or balanced")
    #   }
    # }

    # plot with legend midpoint and limits calculated above
    if (variable %in% c("p_value", "jaccard", "sorensen", "simpson", "alpha_mle")) {

      if (legendlimit == "datarange") {
        # If sig.only is active, keep the datarange scale from the ORIGINAL (unmasked) data.
        # This matters mainly for alpha_mle; otherwise the scale can collapse after masking.
        if (!identical(sig.only, FALSE) && variable == "alpha_mle") {
          lim <- range(legend_data[[variable]], na.rm = TRUE)
          gp <- gp + ggplot2::scale_fill_gradient2(
            midpoint = midpoint, low = col[1], mid = col[2], high = col[3],
            space = "Lab", na.value = "grey50", limits = lim
          )
        } else {
          gp <- gp + ggplot2::scale_fill_gradient2(
            midpoint = midpoint, low = col[1], mid = col[2], high = col[3],
            space = "Lab", na.value = "grey50"
          )
        }
      } else if (legendlimit == "balanced") {
        gp <- gp + ggplot2::scale_fill_gradient2(
          midpoint = midpoint, low = col[1], mid = col[2], high = col[3],
          space = "Lab", na.value = "grey50", limits = c(lowerlimit, upperlimit)
        )
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
        gp <- gp + ggplot2::scale_fill_gradient2(midpoint = midpoint, low = col[1], mid = col[2], high = col[3], space ="Lab", na.value = "grey50")
      } else if(legendlimit == "balanced") {
        gp <- gp + ggplot2::scale_fill_gradient2(midpoint = midpoint, low = col[1], mid = col[2], high = col[3], space ="Lab", na.value = "grey50", limits = c(lowerlimit, upperlimit))
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
    n_entities_plot <- length(keep)

    if ((is.null(show.value) && n_entities_plot <= 20) || isTRUE(show.value)) {
      # gp <- gp + ggplot2::geom_text(aes(label = round(get(variable), value.digit)), size = text.size)
      gp <- gp + ggplot2::geom_text(ggplot2::aes(label = round(.data[[variable]], value.digit)), size = text.size)
      message(paste0(
        "------------------------------ plotgg() notes ------------------------------\n",
        "- Values are printed on tiles. To hide them: show.value = FALSE\n",
        "- To change rounding of printed values: value.digit = <integer>\n",
        "- To change label size: text.size = <numeric>\n",
        "- To keep entities whose tiles are all NA: drop.empty = FALSE\n",
        "---------------------------------------------------------------------------"
      ))
    }

    gp

  }
