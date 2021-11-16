library(cooccur)
data(finches)
finches
dim(finches)

data(beetles)
beetles

data(finches)
N_matrix <- matrix(data = rbinom(n = nrow(finches)*ncol(finches),1,prob = 0.75),
                       nrow = nrow(finches),
                       ncol = ncol(finches)
                       ,byrow = T)
create.N.matrix(N_matrix)


data(finches)
colnames(finches)
row.names(finches)
class(finches)

cooccur.finches <- cooccur(mat=finches,
   type="spp_site",
   thresh=TRUE,
   spp_names=TRUE)
obs.v.exp(cooccur.finches)


pair(cooccur.finches,"Geospiza fortis",all=TRUE)
# so this one computes obs cooccur, prob and exp cooccur of one species to all others. So, it is not NxN solution but 1xN solution. Such outputs are created for every single species.

pair.attributes(cooccur.finches)


# Mainali notes:
# ======================

# 1) i have heatmap of NxN square output where I show affinity of each entity to all other entities. diagonally, I turn off affinity to the self.
# I should give options for plotting entire heatmap like I did. I should also allow with an argument upper triangle or lower traingle for the plotting.
# For the empty triangle, I shouuld plan to show either the number of affinity, or observed and exp cooccurrence or probability or something.

# 2) I should give option to save multiple outputs: alpha MLE, alpha median, alpha interval, jaccard, sorenson, simpson, p-value, observed cooccur, expected cooccur.
# these outputs can each be saved as NxN square matrix. Or, a simple dataframe of all outputs like in pair() example of cooccur.
# Perhaps save in both format so that when this output is take to plot(), the function should give arguments like variable = "alpha" or "alphaMLE". Then it should pick
# the square matrix of that variable, and not the single big dataframe.

# 3) The input data should be data frame or matrix. Check if it has numeric or integers only. And check if it has just 1 and 0. If other numbers, what to do?
# if I find >1, give warning, does this data carry abundance rather than binary presence/absence recoded as 1/0? Because observed numbers greater than 1.
# If it has abundance data that you want to convert to binary presence/absence, please use the function .....
# If the input is not numeric or integer (e.g., characters or factors), fail the code with error message
# argument "compute.for" shoudl have options of row and columns.

# 4) a function to convert abundance to binary presence absence -- check in R universe, there must already be such a function.

# 5) the input data should have a few options. First and default is species by site like balls in boxes. Second, "2x2 contigency table". If this second selected, then row vs column option
# for analysis is ignored. Should we allow raster as input? Perhaps add a raster2matrix() function. For this one, have example binary raster of 3-5 species. Same number of rasters.
# allow adding study area as a polygon. study area should be optional. if not provided, it considers all pixels in the raster. when a raster is imported, first check the bitdepth
# to minimize usage of ram. then replace the original name to import raster to important another raster. then stack and check resolution and extent match. theoretically, if extent do not match, okay.
# but resolutions must match. extracting values from each raster and aligning them might get challenging. check later.
# if this test passes, then crop rasters for the boundary of study area so that a smaller raster is created.
# if the study area does not entirely fall inside a raster, fail and give error message. if this step passes,
# stack all rasters and erase the rasters by the study area so that outside of it is NA.
# search for R functions elsewhere for this to take care of. Next, utilize parallel computing for just reading the values entirely in a raster.
# do not use extract().
#
#
#
# 6) CooccurrenceIndices() - make a function for computing all 80 cooccurrence indices based on the list of functions in Keil et al
# and write a Medium article comparing with alpha
#


# Example datasets
# -------------------
# species affinity: my hypothetical one - save output for 100 sites, 0.6 Jaccard example with various prevalences. So, two species data perhaps?
# site affinity: island spatial beta diversity example from Scientific Report
# temporal affinity: same
# raster

#


# functions
# ----------

# raster2binarymatrix() with argument "input.type" with options: probability, abundance, binary and argument "probability.threshold".
# argument "input.data" has options: rasterstack or individual raster. then options: raster1, raster2, raster3 point to various previously loaded rasters.

# abundance2binary


# install.packages("remotes")
# remotes::install_github("RS-eco/rasterSp")
library(rasterSp)
data(package = "rasterSp")
