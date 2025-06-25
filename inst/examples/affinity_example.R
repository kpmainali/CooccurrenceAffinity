# when you have a binary presence absence occurrence data
# -------------------------------------------------------

if(requireNamespace("cooccur", quietly = TRUE)) {

  data(finches, package = "cooccur")
  head(finches)
  # this dataset carries the occurrence records of 13 species in row in 17 islands in columns
  dim(finches)


  # the remainder of the script has been enclosed under \donttest{}
  # to bypass the CRAN's 5 second limit on example files
  # --------------------------------------------------------------

  \donttest{

    # compute alpha and other quantities for island-pair affinity (beta diversity)
    myout <- affinity(data = finches, row.or.col = "col")
    myout

    # you can simply flip the analysis to rows to compute species-pair affinity
    myout <- affinity(data = finches, row.or.col = "row")
    myout


    # the rows of the outputs above include every single pair of the entities,
    # producing many columns for various quantities.
    # # can output an NxN square matrix for selected columns.
    # an example is here
    myout <- affinity(data = finches, row.or.col = "col", squarematrix = c("alpha_mle", "jaccard"))
    # it is a list of three elements: one main dataframe and two square matrices
    length(myout)
    myout
    head(myout)

    # you can also compute all the square matrices with an "all"
    myout <- affinity(data = finches, row.or.col = "col", squarematrix = c("all"))
    # this one has 12 elements
    length(myout)
    myout

    # when you want to compute for only certain pairs
    myout <- affinity(data = finches, row.or.col = "col", which.row.or.col = 4:6,
                      squarematrix = c("alpha_mle"))
    myout

    myout <- affinity(data = finches, row.or.col = "col",
                      which.row.or.col = c("Isabella", "Espanola"), squarematrix = c("alpha_mle"))
    print(myout)

  } #end of \donttest{}



  # if you select only one column, the computation stops
  \dontrun{
    myout <- affinity(data = finches, row.or.col = "col",
                      which.row.or.col = c("Darwin"), squarematrix = c("alpha_mle"))
  }



  \donttest{

    # you can also add additional arguments bringing them from ML.Alpha() or AlphInts()
    myout1 <- affinity(data = finches, row.or.col = "col",
                       which.row.or.col = c("Isabella", "Espanola"), lev=0.95, pvalType="Blaker")
    myout1
    myout2 <- affinity(data = finches, row.or.col = "col",
                       which.row.or.col = c("Isabella", "Espanola"), lev=0.90, pvalType="Blaker")
    myout2
    identical(myout1, myout2)
    # myout1 and myout2 were generated with identical arguments except a difference in "lev",
    # which gave different confidence intervals

    myout3 <- affinity(data = finches, row.or.col = "col",
                       which.row.or.col = 4:6, lev=0.95, pvalType="Blaker")
    myout3
    myout4 <- affinity(data = finches, row.or.col = "col",
                       which.row.or.col = 4:6, lev=0.95, pvalType="midP")
    myout4
    myout3$all$p_value
    myout4$all$p_value
    # the p values are (or, can be) different

    # when you have abundance data requiring conversion to binary
    # -----------------------------------------------------------
    # abundance data is converted to binary based on a threshold supplied.
    # it might be a good idea to explore dataprep() function and its examples
    # first before workign on affinity() for abundance data.
    matrix.data <- matrix(runif(400, 0, 10), nrow = 100, ncol = 4)
    row.names(matrix.data) <- paste0("id_", 1:nrow(matrix.data))
    colnames(matrix.data) <- paste0("variable_", 1:ncol(matrix.data))

    # add some missing data and zero abundance
    matrix.data[1,1] <- matrix.data[2,3] <- matrix.data[1,4] <- matrix.data[1,2] <- NA
    matrix.data[10,4] <- 0
    head(matrix.data)
    # now this is an abundance data with some missing and some zero occurrences

    # inspecting how the abundance is converted to binary first
    dataprep(data = matrix.data, row.or.col = "col", datatype = "abundance",
             threshold = 5, class0.rule = "less")
    myout10 <- affinity(data = matrix.data, row.or.col = "col",
                        datatype = "abundance", threshold = 5, class0.rule = "less")
    myout10

    # you can also feed the output of dataprep() to affinity()
    myinput <- dataprep(data = matrix.data, row.or.col = "col",
                        datatype = "abundance", threshold = 5, class0.rule = "less")
    myout11 <- affinity(data = myinput, row.or.col = "col", datatype = "binary")
    myout11
    # myout 10 and myout11 are identical
    identical(myout10, myout11)

  } # end of \donttest{}


  } else {

  message("The cooccur package is not installed.
          You can install it with install.packages('cooccur').")
}
