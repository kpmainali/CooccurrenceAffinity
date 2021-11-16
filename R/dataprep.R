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
    # print(paste(length(cols.improper), "column(s) are not in proper format; the columns should be in either numeric or integer format"))
    # print(paste("columns in improper format are the following:"))
    # print(cols.improper)
    # stopifnot("at least one column is not in proper format, see above for details"= length(cols.improper) == 0)

    msg.a <- paste("the following columns/rows are NOT in numeric or integer format: ")
    msg.b <- paste(cols.improper, collapse = ", ")
    msg.c <- paste(" --- ***the data should be in numeric or integer formats only***")
    stop(paste0(msg.a, msg.b, msg.c))

  }
}

# the following function checks if each column of a dataframe is in binary 1/0 format or not; it throws error message if not in binary format.
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
    msg.c <- paste(" --- ***the data should include 1 and 0 only***")
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

dataprep <- function(data, row.or.col, which.row.or.col) {
  # data = df
  # head(data)

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

  # print(head(datasub))

  sapply(datasub, class)
  dataformatcheck(datasub)
  databinarycheck(datasub)

  return(head(datasub))

}

# after dataprep() function, run the affinity() function. affinity() will have dataprep() at the beginning.
# affinity should have an argument something like "indices" with the eighty options. pvalue and cumulative prob and a few others listed in working space should always be included.
# if a user selects jaccard and alpha, then give square dataframes of these as well as p value, expected cooccurrence, observed cooccurrence and then also a combined dataframe
# of all selected indices. allow them to save all square matrices in a list. is there a way to save a list of dataframes as a single file??


matrix.data <- matrix(1:40, nrow = 10, ncol = 4)

matrix.data <- matrix(sample(c(0,1), size = 40, replace = T), nrow = 10, ncol = 4)

row.names(matrix.data) <- paste0("id_", 1:nrow(matrix.data))
colnames(matrix.data) <- paste0("variable_", 1:ncol(matrix.data))
class(matrix.data)
matrix.data

matrix.data[,'variable_1']
class(matrix.data)

df <- as.data.frame(matrix.data)
class(df)
df["variable_1"]


# bad example
dataprep(data = matrix.data, row.or.col = "row", which.row.or.col = c("id_1", "id_4", "id_11", "id_39"))
dataprep(data = matrix.data, row.or.col = "row", which.row.or.col = c(4,7,17))

dataprep(data = matrix.data, row.or.col = "col", which.row.or.col = 2:12)
dataprep(data = matrix.data, row.or.col = "col", which.row.or.col = c("variable_1", "variable_9", "variable_6"))

# good example
dataprep(data = matrix.data, row.or.col = "row", which.row.or.col = c("id_1", "id_4"))
dataprep(data = matrix.data, row.or.col = "col", which.row.or.col = c("variable_1", "variable_4"))





# an example of temporal beta diversity
df <- read.csv("/Users/kpmainali/RPackages/CooccurrenceAffinity/data/argentario.csv")

head(df)
dim(df)

# ----------
dataformatcheck(df[["life.form"]])
dataformatcheck(df$period1)

# --------

sub$cooccur <- sub$period1 + sub$period2

head(sub)
print(sum(is.na(sub)))
A <- sub$period1
A <- A[!is.na(A)]
B <- sub$period2
B <- B[!is.na(B)]
t <- as.data.frame(table(sub$cooccur)); print(t)
X <- t$Freq[t$Var1==2]; X

# affinity mle -- make sure R package "BiasedUrn" has been installed to your computer
alpha_mle_list <- NewAlph(x=X, mA=sum(A), mB=sum(B), N=nrow(sub)); alpha_mle_list
# alpha MLE
alpha_mle_list$AlphMLE

