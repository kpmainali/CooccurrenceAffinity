unlist(AlphInts(30,c(50,80,120), lev=0.9))

AlphInts(30,c(50,80,120), lev=0.9)$CI.CP
AlphInts(30,c(50,80,120), lev=0.9)$MedianIntrvl

EHypMidP(30,c(50,80,120), 0.9)
AlphInts(30,c(50,80,120), lev=0.9)$CI.midP
# NB the third argument of AlphInts is "scal" if not named,
# so must use "lev=0.9" to define the confidence level.

EHypQuInt(30,c(50,80,120), 0.5)
AlphInts(30,c(50,80,120), lev=0.9)$MedianIntrvl

# Alpha capped warning examples
AlphInts(60,c(80,80,100), lev=0.9)
ML.Alpha(60,c(80,80,100), lev=0.9)

AlphInts(80,c(80,80,100), lev=0.9)
ML.Alpha(80,c(80,80,100), lev=0.9)

# impossible x warning examples
AlphInts(81,c(80,80,100), lev=0.9)
ML.Alpha(81,c(80,80,100), lev=0.9)

# Degenerate distribution warning example
AlphInts(80,c(80,100,100), lev=0.9)
ML.Alpha(80,c(80,100,100), lev=0.9)

