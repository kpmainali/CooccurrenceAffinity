NewAlph = function(x,mA,mB,N, cap=TRUE) {
   ### covers exceptional cases but always provides ML alpha 
   ### capped at maximum absolue value log(4*N^2) is cap=T,
   ### which is relevant only to avoid value infinity 
   ### when x is at extreme of null-hypothesis, as of  11/3/2018
    maxabs = if(cap) log(4*N^2) else Inf
    if(x< max(0,mA+mB-N) | x > min(mA,mB)) {
         cat("Impossible co-occurrence number! \n")
         return(list(AlphMLE=NULL)) }
    if(length(intersect(c(mA,mB), c(0,N)))) {
         cat("Degenerate co-occurrence distribution!")
         return(list(AlphMLE=NULL)) }
    pcdf.null = phyper(x,max(mA,mB),N-max(mA,mB),min(mA,mB))
    aux.null = dhyper(x,max(mA,mB),N-max(mA,mB),min(mA,mB)) 
    pmid.null = pcdf.null - 0.5*aux.null
    pval.null = c(Lower=pcdf.null, Upper=1-pcdf.null+aux.null)
    midP.null = c(Lower=pmid.null, Upper=1-pmid.null)
    aint = AlphInt(x,mA,mB,N, scal=1)


    if(x==max(0, mA+mB-N)) {
         cat("Affinity negatively infinite! \n")
         return(list(AlphInt=aint, AlphMLE=-maxabs, 
             pval.null=pval.null, midP.null=midP.null)) }
    if(x==min(mA,mB)) {
         cat("Affinity infinite! \n")
         return(list(AlphInt=aint, AlphMLE=maxabs, 
             pval.null=pval.null, midP.null=midP.null)) }
    require(BiasedUrn)
    amle = AlphMLE(x,mA,mB,N)
    list(AlphInt=aint, AlphMLE=amle, 
         pval.null=pval.null, midP.null=midP.null)  }