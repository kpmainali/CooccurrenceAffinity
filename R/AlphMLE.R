AlphMLE  =  function(x,mA,mB,N, bd=50) {
    aint = pmax(-bd,pmin(bd,AlphInt(x,mA,mB,N)))
    tmp = optimize( function(t) log(dFNCHypergeo(x,mA,N-mA,mB,exp(t))),
                    aint+c(-1,1), maximum=T)
    ### Output is:  alpha value maximizing  likelihood
    ### Function prints "*" if this MLE falls outside AlphInt interval 
    if(aint[1] > tmp$max | aint[2] < tmp$max) cat("*")
    tmp$max  }