#' Occurrence matrix (e.g., species by site) data preparation for affinity() function
#'
#' This function checks the format of the data for its appropriateness, converts abundance to binary and subsets the data for the selected columns or rows.
#' Note that the affinity can be computed between columns or between rows. In the latter case, the dataset is transposed to bring rows into the columns.
#'
#' @details  This function does the following:
#' 1) checks if the supplied data is in matrix or dataframe formats which are the acceptable formats
#' 2) if rows are selected for affinity analysis, it transposes the dataframe
#' 3) subsets the data if specific columns or rows are selected for analysis; the selection can be made with number or name of the rows/columns
#' 4) checks if the selected cols/rows are in numeric or integer format or not
#' 5) checks if the selected cols/rows have data in binary 1/0 format or not; if datatype is specified as abundance, it converts it to binary format following the supplied rule
#
#'
#' @param data occurrence matrix (binary or abundance) in matrix or dataframe format
#' @param row.or.col specify if the pairs of rows or columns are analyzed for affinity. 'row' or 'column'.
#' @param which.row.or.col a vector of name or the number of row/column if a subset of the data is intended to be analyzed; optional argument with default of all rows/columns.
#' @param datatype specify if the datatype is 'abundance' or 'binary'; optional argument with default 'binary'.
#' @param threshold cutoff for converting an abundance data to binary; needed if datatype is 'abundance'
#' @param class0.rule 'less.or.equal' or 'less'. 'less.or.equal' converts a threshold or lower values to zero and all the others to 1. 'less' converts a threshold and higher values to 1.
#'
#' @return A dataframe in binary 1/0 format ready to be analyzed by affinity(). Abundance data is converted to binary.
#' A subset of the input data is returned if certain rows or columns selected.
#' If rows are being analyzed for affinity between pairs, they are brought to columns by transposing the data.
#'
#' @author Kumar Mainali
#'
#' @references
#'
#' @example
#' examples/dataprep_example.R
#'
#' @export



dataprep <- function(data, row.or.col, which.row.or.col=NULL, datatype=NULL, threshold=NULL, class0.rule=NULL) {


  #  ============================ SOME HELPER FUNCTIONS FIRST ============================

  dataformatcheck <- function(data) {

    # This function is a utility forming part of the dataprep() function, and would not be called as a stand-alone.
    # it checks the format of all columns of a dataframe and stops the code if it is not either numeric or integer.

    cols.improper <- c()
    for(cols in colnames(data)) {
      class(data[[cols]])
      if(!is.integer(data[[cols]]) & !is.numeric(data[[cols]])) {
        cols.improper <- c(cols.improper, cols)
      }
    }
    if(length(cols.improper) > 0) {
      msg.a <- paste("the following columns/rows are NOT in numeric or integer format: ")
      msg.b <- paste(cols.improper, collapse = ", ")
      msg.c <- paste(" --- ***the data should be in numeric or integer formats only***")
      stop(paste0(msg.a, msg.b, msg.c))
    }

  }


  abun2binary <- function(data, threshold, class0.rule) {

    # This function is a utility forming part of the dataprep() function, and would not be called as a stand-alone.
    # it converts each column of a dataframe from numeric or integer format to binary 1/0
    # for this function to work, the object should first pass through dataformatcheck() function, ensuring it is either numeric or integer.

    if(is.null(threshold) & is.null(class0.rule)) {
      stop("arguments threshold and class0.rule  missing for converting abundance data to binary format")
    }
    if(is.null(threshold)) {
      stop("argument threshold missing for converting abundance data to binary format")
    }
    if(is.null(class0.rule)) {
      stop("argument class0.rule missing for converting abundance data to binary format")
    }

    # class0.rule = "less.or.equal" for "<= threshold" and "less" for "< threshold".
    if(class0.rule == "less.or.equal") {
      data <- data.frame(ifelse(data <= threshold, yes = 0, no = 1))
    }
    if(class0.rule == "less") {
      data <- data.frame(ifelse(data < threshold, yes = 0, no = 1))
    }
    return(data)

  }


  databinarycheck <- function(data) {

    # This function is a utility forming part of the dataprep() function, and would not be called as a stand-alone.
    # it checks if each column of a dataframe is in binary 1/0 format or not; it throws error message if not in binary format.
    # for this function to work, the object should first pass through dataformatcheck() function, ensuring it is either numeric or integer.

    cols.improper <- c()
    for(cols in colnames(data)) {
      class(data[[cols]])
      # get unique values in the data and check if they include only 1 and 0
      uniqvals <- unique(data[[cols]])
      uniqvals %in% c(0,1,NA)
      all(uniqvals %in% c(0,1,NA))

      if(all(uniqvals %in% c(0,1,NA)) != TRUE) {
        cols.improper <- c(cols.improper, cols)
      }
    }

    if(length(cols.improper) > 0) {
      msg.a <- paste("the following columns/rows are NOT in binary format: ")
      msg.b <- paste(cols.improper, collapse = ", ")
      msg.c <- paste(" --- ***the data should include 1 and 0 only; if you have abundance data, add the arguments datatype, threshold and class0.rule***")
      stop(paste0(msg.a, msg.b, msg.c))
    } else {
      print("------------ as expected, the data ready for analysis has only 1 and 0... 1 = present, 0 = absent used for the interpretation")
    }

  }

  #  ============================ END OF THE HELPER FUNCTIONS ============================


  # check if the input data is either dataframe or a matrix
  if(!is.matrix(data) & !is.data.frame(data)) {
    stop("the input data must be a dataframe or a matrix")
  }

  # check if "data" is a dataframe; if it is a matrix, convert to dataframe
  if(is.matrix(data)) {
    data <- as.data.frame(data)
  }

  # if rows selected, transpose the dataframe
  if(row.or.col == "row") {
    data <- as.data.frame(t(data))
  }


  # subset data if certain rows or columns are selected for analysis
  # ----------------------------------------------------------------
  if(!is.null(which.row.or.col)) {

    # check if at least two columns/rows have been selected for analysis
    if(length(which.row.or.col) < 2) {
      stop("you need to select at least two entities for co-occurrence analysis")
    }

    # check if all selected rows or columns are in the data
    # selected row or columns can be either character, numeric (e.g., c(1,3)) or integer (e.g., 1:3)
    if(is.character(which.row.or.col)) {
      # check which provided column name as character is missing in data, if any
      missingcolposition <- which(!which.row.or.col %in% colnames(data))
      missingcols <- paste(which.row.or.col[missingcolposition], collapse = ", ")
      if(length(missingcolposition) > 0) {
        stop(paste("the following", row.or.col, "are missing in data:", missingcols))
      }
    } else {
      missingcolposition <- which(!which.row.or.col %in% 1:ncol(data))
      missingcols <- paste(which.row.or.col[missingcolposition], collapse = ", ")
      if(length(missingcolposition) > 0) {
        stop(paste("the following", row.or.col, "numbers are missing in data:", missingcols))
      }
    }

    # subset the data for selected rows or columns
    data <- subset(data, select = which.row.or.col)

  }


  # a missing data is not an absent data. for example, in an analysis of affinity between two species based on their occurrences/co-occurrences in multiple sites,
  # a site that is unknown for a species occurrence should be omitted from total number of sites studied, affecting N.
  # as missing data can be species-specific, such records are removed on species-pairwise basis.
  # that means, if a site has missing data for at least one of the pair of species being analyzed for co-occurrence, the site is eliminated from analysis.
  if(sum(is.na(data) >0)) {
    print("***** Missing data detected and will be eliminated on a pairwise basis; please read the details section of function help for an explanation *****")
  }

  # check if the input data is only numeric or integer only
  dataformatcheck(data)

  # convert abundance data to binary data if datatype = abundance
  # datatype will be binary by default
  if(!is.null(datatype)) {
    if(!datatype %in% c("abundance", "binary")) {
      stop("datatype can be either abundance or binary")
    }
    if(datatype == "abundance") {
      data <- abun2binary(data, threshold = threshold, class0.rule = class0.rule)
    } else {
      data <- data
    }
  }

  # check if the input data is binary 1/0 only
  databinarycheck(data)

  return((data))

}
