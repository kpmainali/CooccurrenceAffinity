# when you have a binary presence absence occurrence data
require(cooccur)
data(finches)

head(finches)
# this dataset carries the occurrence records of 13 species in row in 17 islands in columns
dim(finches)

# compute alpha and other quantities for island-pair affinity (beta diversity)
myout <- affinity(data = finches, row.or.col = "col")
myout

# you can simply flip the analysis to rows to compute species-pair affinity
myout <- affinity(data = finches, row.or.col = "row")
myout

# the rows of the outputs above include every single pair of the entities, producing many columns for various quantities.
# for several of these columns, the function allows outputing the square NxN matrix of the entities
# an example is here

myout <- affinity(data = finches, row.or.col = "col", squarematrix = c("alpha_mle", "jaccard"))
# it is a list of three elements: one main dataframe and two square matrices
length(myout)
myout

# you can also compute all the square matrices with an "all"
myout <- affinity(data = finches, row.or.col = "col", squarematrix = c("all"))
# this one has 12 elements
length(myout)
myout

help(AlphInts)
dataprep(data = finches, row.or.col = "row")
