unlist(ML.Alpha(30,c(50,80,120), lev=0.9))
AlphInts(30,c(50,80,120), lev=0.9)

AlphInts(61,c(80,80,100), lev=0.9)
ML.Alpha(61,c(80,80,100), lev=0.9)

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
