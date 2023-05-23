matrix.data <- matrix(1:40, nrow = 10, ncol = 4)

row.names(matrix.data) <- paste0("id_", 1:nrow(matrix.data))
colnames(matrix.data) <- paste0("variable_", 1:ncol(matrix.data))

# add some missing data and zero abundance
matrix.data[1,1] <- matrix.data[2,3] <- matrix.data[1,4] <- matrix.data[1,2] <- NA
matrix.data[10,4] <- 0
matrix.data
# abundance data with some missing and some zero occurrences

# some good examples
dataprep(data = matrix.data, row.or.col = "col", datatype = "abundance",
         threshold = 9, class0.rule = "less")
dataprep(data = matrix.data, row.or.col = "row", which.row.or.col = c("id_2", "id_4"),
         datatype = "abundance", threshold = 10, class0.rule = "less")
dataprep(data = matrix.data, row.or.col = "col", which.row.or.col = c("variable_1", "variable_4"),
         datatype = "abundance", threshold = 8, class0.rule = "less")
dataprep(data = matrix.data, row.or.col = "col",
         which.row.or.col = c("variable_1", "variable_3", "variable_4"),
         datatype = "abundance", threshold = 8, class0.rule = "less.or.equal")
dataprep(data = matrix.data, row.or.col = "row", datatype = "abundance",
         threshold = 10, class0.rule = "less")
dataprep(data = matrix.data, row.or.col = "col", datatype = "abundance",
         threshold = 10, class0.rule = "less")

# bad examples of specifying the rows or cols that are not in the data
\dontrun{
  dataprep(data = matrix.data, row.or.col = "row",
           which.row.or.col = c("id_1", "id_4", "id_11", "id_39"), datatype = "abundance",
           threshold = 10, class0.rule = "less")
  dataprep(data = matrix.data, row.or.col = "row", which.row.or.col = c(4,7,17),
           datatype = "abundance", threshold = 10, class0.rule = "less")
  dataprep(data = matrix.data, row.or.col = "col", which.row.or.col = 2:12, datatype = "abundance",
           threshold = 10, class0.rule = "less")
  dataprep(data = matrix.data, row.or.col = "col",
           which.row.or.col = c("variable_1", "variable_9", "variable_6"), datatype = "abundance",
           threshold = 10, class0.rule = "less")
}


# what if you pick just one column or row
\dontrun{
  dataprep(data = matrix.data, row.or.col = "row", which.row.or.col = c("id_4"),
           datatype = "abundance", threshold = 10, class0.rule = "less")
}

# the function fails when a required argument is missing
\dontrun{
  dataprep(data = matrix.data, row.or.col = "col", which.row.or.col = c("variable_1", "variable_4"),
           datatype = "abundance", threshold = 10)
  dataprep(data = matrix.data, row.or.col = "col", which.row.or.col = c("variable_1", "variable_4"),
           datatype = "abundance", class0.rule = "less.or.equal")
  dataprep(data = matrix.data, row.or.col = "col", which.row.or.col = c("variable_1", "variable_4"),
           datatype = "abundance")
}

# what if you have abundance data but do not specify the datatype
\dontrun{
  dataprep(data = matrix.data, row.or.col = "col", which.row.or.col = c("variable_1", "variable_4"))
}

# however, if it is a binary data, it's okay to not specify the datatype
# although specifying is a good practice
matrix.bindata <- dataprep(data = matrix.data, row.or.col = "col", datatype = "abundance",
                           threshold = 9, class0.rule = "less")
matrix.bindata
dataprep(data = matrix.bindata, row.or.col = "col")
dataprep(data = matrix.bindata, row.or.col = "row")
