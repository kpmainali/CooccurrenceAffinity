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
#' \code{sig.only = 0.01}). Requires a \code{p_value} column in \code{data$all}.
#' When \code{variable = "p_value"}, p-values above the cutoff are masked to
#' \code{NA}.
#'
#' Legend titles are mapped to human-readable labels (some shown on two lines),
#' rather than using raw column names from \code{data$all}.
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
#' observed.  For indices bounded in \eqn{[0,1]} (\code{p_value}, \code{jaccard},
#' \code{sorensen}, \code{simpson}), the balanced scale uses fixed limits
#' \eqn{[0,1]}. For \code{p_value}, the color mapping is reversed so smaller
#' p-values appear more intense. For count-based variables, no natural midpoint exists; the
#' color scale spans the observed range. For \code{obs_cooccur_X} and
#' \code{exp_cooccur}, a shared color scale is applied so the two plots are
#' visually comparable.
#'
#' When \code{show.value = TRUE}, numeric values are printed on each tile using
#' \code{ggplot2::geom_text()}. If \code{show.value = NULL} (default), values are
#' printed automatically when the number of plotted entities is \eqn{\le 20}.
#' Rounding and text appearance are controlled by \code{value.digit},
#' \code{text.size}, and \code{text.col}.
#'
#' @param data Output list returned by \code{\link{affinity}}.
#' @param variable Name of a numeric column in \code{data$all} to plot.
#' @param legendlimit Either \code{"datarange"} or \code{"balanced"}.
#' @param show.value Logical; if \code{TRUE}, values are printed on tiles. If
#'   \code{NULL}, values are printed automatically when the number of plotted
#'   entities is \eqn{\le 20}.
#' @param value.digit Number of digits used when printing values; default is 2.
#' @param text.size Size of printed values; default is 2.5.
#' @param text.col Color of printed values on tiles (used when values are shown).
#' @param plot.margin Plot margin passed to \code{ggplot2::theme(plot.margin = ...)}.
#'   Typically a \code{ggplot2::margin(t, r, b, l, unit)} object.
#' @param drop.empty Logical; if \code{TRUE} (default), entities whose values are
#'   all \code{NA} for the selected \code{variable} are removed from the plot.
#'   Set to \code{FALSE} to keep all entities.
#' @param sig.only Logical or numeric. If \code{FALSE} (default), all values are
#'   plotted. If \code{TRUE}, tiles are masked to \code{NA} wherever
#'   \code{p_value > 0.05}. If numeric, the value is used as the p-value cutoff
#'   (e.g., \code{sig.only = 0.01}). Requires a \code{p_value} column in
#'   \code{data$all}. When \code{variable = "p_value"}, p-values above the cutoff
#'   are masked to \code{NA}.
#' @param col Color specification for the fill scale. For \code{alpha_mle}
#'   (diverging), supply \code{c(low, high)} or \code{c(low, mid, high)}.
#'   For all other variables (sequential), supply \code{c(low, high)}.
#'   If \code{NULL}, defaults are used (including an auto-generated low color for
#'   sequential scales).
#' @param ... Additional arguments (currently unused).
#'
#' @return A heatmap plot generated with \code{ggplot2}.
#'
#' @author Kumar Mainali
#'
#' @seealso \code{\link{affinity}}
#'
#' @example
#' inst/examples/plotgg_example.R
#'
#' @export



plotgg <- function (data,
                     variable,
                     legendlimit,
                     col = NULL,
                     show.value = NULL,
                     value.digit = NULL,
                     text.size = NULL,
                     text.col = NULL,
                     plot.margin = NULL,
                     drop.empty = TRUE,
                     sig.only = FALSE,
                     ...) {

  # ----------------------------- input validation -----------------------------
  if (is.null(data) || is.null(data$all) || !is.data.frame(data$all)) {
    stop("data must be the output list returned by affinity() and contain a data.frame in data$all")
  }

  # occur_mat is used ONLY to determine entity ordering for axes.
  # It MUST be present; we coerce data.frame -> matrix but do not silently guess order.
  if (is.null(data$occur_mat)) {
    stop("data must contain occur_mat to determine entity ordering")
  }
  if (!is.matrix(data$occur_mat)) {
    if (is.data.frame(data$occur_mat)) {
      data$occur_mat <- as.matrix(data$occur_mat)
    } else {
      stop("occur_mat must be a matrix or coercible data.frame")
    }
  }

  if (!variable %in% colnames(data$all)) {
    stop("the variable does not exist in data$all")
  }

  # variables we explicitly refuse to plot
  banned_vars <- c("alpha_medianInt", "conf_level", "ci_blaker", "ci_cp", "ci_midQ", "ci_midP")
  if (variable %in% banned_vars) {
    stop("honestly, we do not like to plot intervals and confidence level... at least for now")
  }

  # legendlimit validation
  if (!legendlimit %in% c("datarange", "balanced")) {
    stop("legendlimit should be either 'datarange' or 'balanced'")
  }

  # ----------------------------- legend title mapping -------------------------
  # Two-line titles use "\n". For alpha-hat use plotmath via expression().
  fill_title_tbl <- data.frame(
    variable   = c("entity_1_count_mA", "entity_2_count_mB", "obs_cooccur_X", "exp_cooccur",
                   "total_N", "p_value", "alpha_mle", "jaccard", "sorensen", "simpson"),
    fill_title = I(list(
      "Entity A\nCount",
      "Entity B\nCount",
      "Observed\nCo-Occurrences",
      "Expected\nCo-Occurrences",
      "Total N",
      "P-Value",
      "Affinity\n(Alpha MLE)",
      "Jaccard",
      "S\u00f8rensen",
      "Simpson"
    )),
    stringsAsFactors = FALSE
  )

  idx <- match(variable, fill_title_tbl$variable)
  if (is.na(idx)) {
    stop(paste0(
      "unsupported variable for plotgg(): '", variable, "'.\n",
      "Supported variables: ", paste(fill_title_tbl$variable, collapse = ", ")
    ))
  }
  fill_title <- fill_title_tbl$fill_title[[idx]]

  # defaults
  # Color policy:
  # - Defaults (when col is NULL):
  #   * alpha_mle (diverging): blue -> white -> red  (negative -> 0 -> positive)
  #   * all sequential variables: single-hue sequential lightened(high) -> high
  # - User input:
  #   * Sequential variables: col must be c(low, high)
  #   * alpha_mle: col may be c(low, high) (mid defaults to white) or c(low, mid, high)
  user_supplied_col <- !is.null(col)
  default_seq_high <- "#fd6a6c"
  default_div_cols <- c("#87beff", "white", "#fd6a6c")
  if (is.null(col)) {
    col <- if (variable == "alpha_mle") default_div_cols else default_seq_high
  }
  if (is.null(value.digit)) value.digit <- 2
  if (is.null(text.size)) text.size <- 2.5
  if (is.null(text.col))  text.col  <- "black"

  # ----------------------------- variable groups -----------------------------
  vars_01     <- c("p_value", "jaccard", "sorensen", "simpson")
  vars_signed <- c("alpha_mle")
  vars_count  <- c("entity_1_count_mA", "entity_2_count_mB",
                   "obs_cooccur_X", "exp_cooccur", "total_N")

  supported_vars <- c(vars_01, vars_signed, vars_count)
  if (!variable %in% supported_vars) {
    stop(paste0(
      "unsupported variable for plotgg(): '", variable, "'.\n",
      "Supported variables: ", paste(supported_vars, collapse = ", ")
    ))
  }

  # ----------------------------- color normalization --------------------------
  # Helper: create a lighter shade from a given hex color.
  make_lighter <- function(hex, t = 0.20) {
    rgb <- grDevices::col2rgb(hex) / 255
    out <- 1 - (1 - rgb) * t
    grDevices::rgb(out[1], out[2], out[3])
  }

  if (variable %in% vars_signed) {
    # alpha_mle (diverging):
    # - Default (col was NULL): blue -> white -> red
    # - User may supply:
    #     col=c(low,high)     -> uses low -> white -> high
    #     col=c(low,mid,high) -> uses low -> mid  -> high
    if (!user_supplied_col && length(col) >= 3) {
      # keep default as-is (blue/white/red)
      col <- col[1:3]
    } else if (length(col) == 2) {
      message("* Note (alpha_mle): using white as the default midpoint color. ",
              "To override, supply three colors: col = c(low, mid, high).\n")
      col <- c(col[1], "white", col[2])
    } else if (length(col) >= 3) {
      col <- col[1:3]
    } else {
      stop("alpha_mle requires col=c(low,high) or col=c(low,mid,high)")
    }
  } else {
    # sequential variables:
    # - Default (col was NULL): single-hue sequential lightened(high) -> high
    # - User must supply col=c(low,high)
    if (!user_supplied_col && length(col) == 1) {
      high_col <- col[1]
      col <- c(make_lighter(high_col, t = 0.20), high_col)
    } else {
      if (length(col) != 2) {
        stop("Sequential variables require col=c(low, high).")
      }
    }
  }

  # ----------------------------- data preparation ----------------------------
  # plot_data: what we actually draw (may be masked and/or subset by drop.empty)
  plot_data <- data$all

  # legend_data: ALWAYS the original, unmasked data used for legend scaling
  legend_data <- data$all

  # Full ordered entity list (critical for consistent axes / triangle layout)
  entities_full <- colnames(data$occur_mat)

  if (is.null(entities_full) || length(entities_full) < 2) {
    stop("data$occur_mat must have at least 2 column names representing entities")
  }

  # ---------------------- optional masking by significance --------------------
  # This only masks tiles; legend scaling remains anchored to legend_data.
  if (!identical(sig.only, FALSE)) {

    if (!"p_value" %in% colnames(plot_data)) {
      stop("sig.only requires a p_value column in data$all")
    }

    # interpret sig.only: TRUE -> 0.05, numeric -> user cutoff
    if (isTRUE(sig.only)) {
      p_cut <- 0.05
    } else if (is.numeric(sig.only) && length(sig.only) == 1 && is.finite(sig.only)) {
      p_cut <- sig.only
    } else {
      stop("sig.only must be FALSE, TRUE, or a single numeric p-value cutoff (e.g., 0.05)")
    }

    # Apply masking
    if (variable == "p_value") {
      # Mask p-values above cutoff (and NA) when plotting p_value itself
      mask_idx <- is.na(plot_data$p_value) | plot_data$p_value > p_cut
      plot_data$p_value[mask_idx] <- NA
    } else {
      # Mask the selected variable where p_value is not significant OR missing
      mask_idx <- is.na(plot_data$p_value) | plot_data$p_value > p_cut
      plot_data[[variable]][mask_idx] <- NA
    }
  }

  # ----------------------------- drop.empty logic ----------------------------
  # If drop.empty=TRUE: remove entities that have all-NA values for the selected variable
  # *after* any sig.only masking (expected behavior).
  if (isTRUE(drop.empty)) {

    non_na <- !is.na(plot_data[[variable]])

    keep <- unique(c(plot_data$entity_1[non_na], plot_data$entity_2[non_na]))
    keep <- keep[!is.na(keep)]

    # preserve original order from occur_mat
    keep <- entities_full[entities_full %in% keep]

    if (length(keep) < 2) {
      stop(paste0(
        "drop.empty=TRUE removed all entities for variable='", variable,
        "' (no non-NA values to plot)."
      ))
    }

    plot_data <- plot_data[
      plot_data$entity_1 %in% keep & plot_data$entity_2 %in% keep,
      , drop = FALSE
    ]

  } else {
    keep <- entities_full
  }

  # Axis vectors for your lower-triangle layout
  x_entities <- keep[-length(keep)]
  y_entities <- rev(keep[-1])

  # ----------------------------- plot skeleton -------------------------------
  entity_1 <- entity_2 <- NULL

  gp <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = entity_1, y = entity_2, fill = .data[[variable]])
  ) +
    ggplot2::geom_tile(color = "gray") +
    ggplot2::coord_fixed() +
    ggplot2::labs(fill = fill_title) +
    ggplot2::ylim(y_entities) +
    ggplot2::xlim(x_entities) +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 35, vjust = 0.85, hjust = 1),
      axis.text.y = ggplot2::element_text(vjust = 0.5, hjust = 0.1),
      axis.ticks.length = ggplot2::unit(.25, "cm"),
      legend.title = ggplot2::element_text(
        face = "bold",
        size = ggplot2::rel(0.9)
      ),
      legend.text = ggplot2::element_text(
        size = ggplot2::rel(0.85)
      )
    )

  if (!is.null(plot.margin)) {
    gp <- gp + ggplot2::theme(plot.margin = plot.margin)
  }

  # ---------------------- legend scaling: compute limits ----------------------
  # All legend computations use legend_data (original, unmasked).

  get_scale_vector <- function(var) {
    if (var %in% c("obs_cooccur_X", "exp_cooccur")) {
      c(legend_data$obs_cooccur_X, legend_data$exp_cooccur)
    } else {
      legend_data[[var]]
    }
  }

  scale_vec <- get_scale_vector(variable)
  scale_vec <- scale_vec[is.finite(scale_vec)]

  if (length(scale_vec) == 0) {
    stop(paste0("cannot compute legend limits: no finite values found for variable='", variable, "'"))
  }

  lim_data <- range(scale_vec, na.rm = TRUE)

  legend_mode_requested <- legendlimit
  legend_mode_effective <- legendlimit

  if (variable %in% vars_01) {

    if (legendlimit == "balanced") {
      lim_final <- c(0, 1)
    } else {
      lim_final <- lim_data
    }

  } else if (variable %in% vars_signed) {

    if (legendlimit == "balanced") {
      A <- max(abs(lim_data))
      lim_final <- c(-A, A)
    } else {
      lim_final <- lim_data
    }

  } else {

    lim_final <- lim_data

    if (legendlimit == "balanced") {
      legend_mode_effective <- "datarange"
    }
  }

  # ---------------------- legend notes: standardized output -------------------
  fmt_range <- function(x) {
    if (any(!is.finite(x)) || length(x) != 2) return("<unavailable>")
    paste0("[", formatC(x[1], format = "fg", digits = 4), ", ",
           formatC(x[2], format = "fg", digits = 4), "]")
  }

  legend_lines <- character(0)

  if (variable %in% vars_signed) {
    legend_lines <- c(legend_lines,
                      "- Palette: diverging (centered at 0 for alpha_mle).")
  } else {
    legend_lines <- c(legend_lines,
                      "- Palette: sequential (light -> high).")
  }

  if (legend_mode_requested == legend_mode_effective) {
    legend_lines <- c(legend_lines,
                      paste0("- legendlimit='", legend_mode_requested, "' applied as requested."))
  } else {
    legend_lines <- c(legend_lines,
                      paste0("- No natural balanced scale exists for this variable; legendlimit='",
                             legend_mode_effective, "' was used instead."))
  }

  if (variable %in% vars_01 && legend_mode_effective == "balanced") {
    legend_lines <- c(legend_lines,
                      "- Limits used: fixed [0, 1].")
  } else if (variable %in% vars_signed && legend_mode_effective == "balanced") {
    legend_lines <- c(legend_lines,
                      paste0("- Limits used: symmetric around 0 ", fmt_range(lim_final), "."))
  } else {
    legend_lines <- c(legend_lines,
                      paste0("- Limits used: observed data range ", fmt_range(lim_final), "."))
  }

  if (variable %in% c("obs_cooccur_X", "exp_cooccur")) {
    legend_lines <- c(legend_lines,
                      "- A shared legend range is used for obs_cooccur_X and exp_cooccur to enable visual comparison.")
  }

  if (variable == "p_value") {
    legend_lines <- c(legend_lines,
                      "- Note: p_value colors are reversed so smaller p-values appear more intense.")
  }


  legend_lines <- c(
    legend_lines,
    "- Legend limits are computed from the original (unmasked) values; sig.only masks tiles but does not change the color scale."
  )

  message(paste0(
    "------------------------------ legend notes --------------------------------\n",
    paste0(legend_lines, collapse = "\n"),
    "\n----------------------------------------------------------------------------\n"
  ))

  # ---------------------- apply fill scale (palette-aware) --------------------
  if (variable %in% vars_signed) {

    gp <- gp + ggplot2::scale_fill_gradient2(
      # name = fill_title,
      midpoint = 0,
      low = col[1], mid = col[2], high = col[3],
      space = "Lab",
      na.value = "grey50",
      limits = lim_final
    )

  } else {

    seq_low  <- col[1]
    seq_high <- col[2]

    # Reverse color mapping ONLY for p_value (legend numeric direction stays low->high)
    if (variable == "p_value") {
      gp <- gp + ggplot2::scale_fill_gradient(
        # name = fill_title,
        low = seq_high,
        high = seq_low,
        na.value = "grey50",
        limits = lim_final
      )
    } else {
      gp <- gp + ggplot2::scale_fill_gradient(
        # name = fill_title,
        low = seq_low,
        high = seq_high,
        na.value = "grey50",
        limits = lim_final
      )
    }
  }

  # --------------------------- value label handling ---------------------------
  # n_entities_plot <- length(keep)
  #
  # if ((is.null(show.value) && n_entities_plot <= 20) || isTRUE(show.value)) {
  #
  #   gp <- gp + ggplot2::geom_text(
  #     ggplot2::aes(label = round(.data[[variable]], value.digit)),
  #     size = text.size,
  #     color = text.col
  #   )
  #
  #   message(paste0(
  #     "------------------------------ plotgg() notes ------------------------------\n",
  #     "- Values are printed on tiles. To hide them: show.value = FALSE\n",
  #     "- To change rounding of printed values: value.digit = <integer>\n",
  #     "- To change label size: text.size = <numeric>\n",
  #     "- To change text color: text.col = <color>\n",
  #     "- To keep entities whose tiles are all NA: drop.empty = FALSE\n",
  #     "---------------------------------------------------------------------------"
  #   ))
  # }

  n_entities_plot <- length(keep)

  values_drawn <- (is.null(show.value) && n_entities_plot <= 20) || isTRUE(show.value)

  if (values_drawn) {
    gp <- gp + ggplot2::geom_text(
      ggplot2::aes(label = round(.data[[variable]], value.digit)),
      size = text.size,
      color = text.col
    )
  }

  # --- ALWAYS print plotgg notes ---
  message(paste0(
    "------------------------------ plotgg() notes ------------------------------\n",
    if (values_drawn) {
      "- Values are printed on tiles.\n"
    } else {
      "- Values are NOT printed on tiles.\n"
    },
    "- To show values explicitly: show.value = TRUE\n",
    "- To hide values explicitly: show.value = FALSE\n",
    "- To change rounding of printed values: value.digit = <integer>\n",
    "- To change label size: text.size = <numeric>\n",
    "- To change text color: text.col = <color>\n",
    "- To keep entities whose tiles are all NA: drop.empty = FALSE\n",
    "---------------------------------------------------------------------------\n"
  ))

  # --- EXTRA message ONLY for automatic suppression ---
  if (is.null(show.value) && n_entities_plot > 20) {
    message(paste0(
      "------------------- note on automatic label suppression -------------------\n",
      "- Tile values are not printed automatically because the number of entities exceeds 20.\n",
      "- To force printing values: show.value = TRUE\n",
      "---------------------------------------------------------------------------"
    ))
  }

  gp

}
