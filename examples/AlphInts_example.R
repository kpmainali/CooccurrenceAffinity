unlist(AlphInts(30,c(50,80,120), lev=0.9))

AlphInts(30,c(50,80,120), lev=0.9)$CI.CP
AlphInts(30,c(50,80,120), lev=0.9)$MedianIntrvl

EHypMidP(30,c(50,80,120), 0.9)
AlphInts(30,c(50,80,120), lev=0.9)$CI.midP
# NB the third argument of AlphInts is "scal" if not named,
#   so must use "lev=0.9" to define the confidence level.
