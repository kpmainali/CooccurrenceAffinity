# True coverage probability by the 95% Confidence Interval
CovrgPlot(marg = c(50,70,150), lev = 0.95)

# If confidence level is not specified, it shows results for 95% CI
CovrgPlot(marg = c(50,70,150))

# True coverage probability by the 90% Confidence Interval
CovrgPlot(marg = c(50,70,150), lev = 0.90)

# Select some of the confidence intervals to plot
CovrgPlot(marg = c(50,70,150), lev = 0.95, plotCI=1:3)
CovrgPlot(marg = c(50,70,150), lev = 0.95, plotCI=c(1,3,4))
CovrgPlot(marg = c(50,70,150), lev = 0.95, plotCI=1)
