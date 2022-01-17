# CooccurrenceAffinity

This package computes affinity between two entities based on their co-occurrence (using binary presence/absence data). 

The package refers to and requires an existing package called BiasedUrn, the primary functions in which calculate distributional characteristics of the Fisher Noncentral Hypergeometric distribution (pFNCHypergeo) otherwise known as Extended Hypergeometric (Harkness 1965), which is the way we refer to it in these notes.

The applications served in the present "Affinity" package are primarily Ecology and other biological science data analyses in which associations of species occurrences are of interest, although the same kinds of hypothesis and estimates of statistical association among pairs of entities arise in other settings in biological and social science. (See Mainali, Slud, Singer and Fagan 2021 and Supplements for full discussion. In that paper, the connection between the scientific problems of quantifying association are related explicitly to the classic balls-in-boxes formulation underlying the hypergeometric and extended-hypergeometric distribution families.)

The statistical content of our package "Affinity" is the likelihood-based (frequentist) point and interval estimation of the log-odds-ratio parameter in the Extended Hypergeometric distribution, for fixed values of the prevalence (mA, mB, respectively the total numbers of boxes containing type-A balls and type-B balls) and 2x2 table-total (N) parameters and a value (X) of the co-occurrence count, i.e. of the number of boxes containing both a type-A and a type-B ball. We call this parameter "alpha" or, interchangeably, the natural logarithm of "affinity". Its exponential, the odds or "affinity", is understood intuitively as the ratio of the odds of a site (a "box") being occupied by a type-A ball when it is already occupied by a type-B balls over the odds of type-A occupancy when no type-B ball is in that box.

Our primary functions ML.Alpha and AlphInts calculate the maximum likelihood estimate (MLE) of alpha as well as intervals (a1(x,q),a2(x,q)) of alpha values for which F(x,mA,mB,N,exp(alpha)) >= q and 1-F(x-1,mA,mB,N,exp(alpha)) >= 1-q, for specified choices of the quantile q, where F(x,mA,mB,N,exp(alpha)) denotes the Extended Hypergeometric distribution function for the co-occurance count X. The mid-point of this interval, for q=1/2, is a second reasonable statistical estimate of alpha. Furthermore, "test-based" confidence intervals for alpha are also immediately obtained from this function. For example, a two-sided 90% confidence interval would be reported either as:

    ( a1(x,0.95), a2(x,0.05) )                                         (1)

which is probably a conservative confidence interval for alpha, analogous to the Clopper-Pearson (1934) interval for binomial proportions, or
    ( (a1(x,0.95)+a2(x,0.95)/2, (a1(x,0.05)+a2(x,0.05)/2 )             (2)
which (as we will see below) has coverage much closer to its nominal level of 90%.  All these estimates and confidence intervals are viewed as functions of the co-occurrence count x for a 2x2 table with fixed marginal counts mA, mB and table-total N. 

Two other confidence intervals for alpha are calculated in the package functions  AlphaInts  and ML.Alpha. One is another conservative confidence interval based on theoretical results of Blaker (2000) (that is, an interval whose coverage probability is provably at least as large as the nominal  confidence level) using the so-called "Acceptability Function" in that paper's Theorem 1. This confidence interval also provably lies within the first interval (1) above, so nothing is long in using it in preference to (1) except that it is somewhat less direct to explain. The last confidence interval we calculate is similar in performance to (2) defined above, but is close in spirit to the "mid-P" confidence interval defined in standard references like Agresti (2013) for the unknown probability of success in Binomial triala. This interval expressed for the Extended Hypergeometric distribution is
   (b(x,mA,mB,N,0.95), b(x,mA,mB,N,0.05))                                     (3)

where b = b(x,mA,mB,N, gamma)   solves
    (F(x,mA,mB,N,exp(b))+F(x-1,mA,mB,N,exp(b)))/2 = gamma

The test-based confidence intervals for alpha described in the previous paragraphs have more reliable moderate-sample coverage than Confidence Intervals based on a normal-distribution approximation to the MLE of alpha. This will be established in a separate small simulation study. The situation is closely related to that of confidence intervals for an unknown binomial-dstribution success probability p (Brown, Cai and DasGupta 2001). The test-based interval (a1(x,0.05), a2(x,0.95)) is analogous to the Clopper-Pearson (1934) confidence interval for the binomial p. The famous Wald interval for binomial p would correspond here to a symmetric confidence interval round the MLE based on the approximate normal distribution of the MLE of alpha. 


References

Agresti, A. (2013) Categorical Data Analysis, 3rd edition, Wiley.

Blaker, H. (2000), “Confidence curves and im[proved exact confudence intervals for discrete distributions, Canadian Journal of Statistics 28, 783-798.

Brown, L., T. Cai, and A. DasGupta (2001), “Interval Estimation for a Binomial Proportion,” Statistical Science, 16, 101–117.

Clopper, C., and E. Pearson (1934), “The Use of Confidence or Fiducial Limits Illustrated in the Case of the Binomial,” Biometrika, 26, 404–413.

Mainali, K., Slud, E., Singer, M. and Fagan, B. (2021), “A better index for analysis of co-occurrence and similarity”, Science Advances [in press].

Wilson, E. (1927), “Probable Inference, the Law of Succession, and Statistical Inference,” Journal of the American Statistical Association, 22, 209–212.
