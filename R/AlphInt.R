AlphInt = function(x,mA,mB,N, scal=10) {
  ### function to calculate interval of alpha = log(odds)  values
  ### for which  x>=1 is a median of pFNCHypergeo(x,mA,N-mA,mB,exp(alpha))
  ### NB. scal=10 seems to be largest ever needed.
  if(length(intersect(c(mA,mB), c(0,N))))
    return("Degenerate co-occurrence distribution!")
  
  require(BiasedUrn)
  EHypCent = function(alp,t) {
    L = length(alp)
    out = numeric(L)
    for(i in 1:L) out[i] =
      pFNCHypergeo(t,mA,N-mA,mB,exp(alp[i]))-0.5
    out }
  null.mean = mA*mB/N
  newint = if(x >= null.mean) c(-1,scal) else c(-scal,1)
  if(x<max(mA+mB-N,0)) c(-Inf,-Inf) else {
    if(x==max(mA+mB-N,0)) c(-Inf,
                            uniroot(EHypCent, newint, t=x, extendInt="yes")$root) else {
                              if(x>min(mA,mB)) c(Inf,Inf) else {
                                if(x==min(mA,mB)) c(uniroot(EHypCent,
                                                            newint, t=x-1, extendInt="yes")$root,Inf) else
                                                              c(uniroot(EHypCent, newint, t=x-1, extendInt="yes")$root,
                                                                uniroot(EHypCent, newint, t=x, extendInt="yes")$root)
                              }}}}