# the following function checks the format of all columns of a dataframe and stops the code if it is not either numeric or integer.
dataformatcheck <- function(data) {
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


# the following function converts each column of a dataframe from numeric or integer format to binary 1/0
# for this function to work, the object should first pass through dataformatcheck() function, ensuring it is either numeric or integer.
abun2binary <- function(data, threshold, class0.rule) {
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


# the following function checks if each column of a dataframe is in binary 1/0 format or not; it throws error message if not in binary format.
# for this function to work, the object should first pass through dataformatcheck() function, ensuring it is either numeric or integer.
databinarycheck <- function(data) {
  cols.improper <- c()
  for(cols in colnames(data)) {
    class(data[[cols]])
    # get unique values in the data and check if they include only 1 and 0
    uniqvals <- unique(data[[cols]])
    uniqvals %in% c(0,1)
    all(uniqvals %in% c(0,1))

    if(all(uniqvals %in% c(0,1)) != TRUE) {
      cols.improper <- c(cols.improper, cols)
    }
  }

  if(length(cols.improper) > 0) {
    msg.a <- paste("the following columns/rows are NOT in binary format: ")
    msg.b <- paste(cols.improper, collapse = ", ")
    msg.c <- paste(" --- ***the data should include 1 and 0 only; if you have abundance data, add the arguments datatype, threshold and class0.rule***")
    stop(paste0(msg.a, msg.b, msg.c))
  } else {
    print("------------ only 1 and 0 found in the data... 1 = present, 0 = absent used for the interpretation")
  }
}


# the following function does the following:
# 1) checks if the supplied data is in matrix or dataframe format
# 2) if the data is in matrix format, it converts it to dataframe
# 3) if rows are selected for affinity analysis, it transposes the dataframe
# 4) takes the columns or rows for analysis with the number or name and subsets the matrix/dataframe only for those selected cols/rows
# 5) checks if the selected cols/rows are in numeric or integer format or not
# 6) checks if the selected cols/rows have data in binary 1/0 format or not

dataprep <- function(data, row.or.col, which.row.or.col, datatype=NULL, threshold=NULL, class0.rule=NULL) {

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

  datasub <- subset(data, select = which.row.or.col)

  # check if the input data is only numeric or integer only
  dataformatcheck(datasub)

  # convert abundance data to binary data if datatype = abundance
  # datatype will be binary by default
  if(!is.null(datatype)) {
    if(datatype == "abundance") {
      datasub <- abun2binary(datasub, threshold = threshold, class0.rule = class0.rule)
    }
  }


  # check if the input data is binary 1/0 only
  databinarycheck(datasub)

  return((datasub))

}


# after dataprep() function, run the affinity() function. affinity() will have dataprep() at the beginning.
# affinity should have an argument something like "indices" with the eighty options. pvalue and cumulative prob and a few others listed in working space should always be included.
# if a user selects jaccard and alpha, then give square dataframes of these as well as p value, expected cooccurrence, observed cooccurrence and then also a combined dataframe
# of all selected indices. allow them to save all square matrices in a list. is there a way to save a list of dataframes as a single file??

