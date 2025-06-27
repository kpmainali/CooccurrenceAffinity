1.0.2 (2025-06-26)
------------------
- Dropped **cooccur** from Suggests; the `finches` dataset now ships internally and can be loaded via `data(finches)` without requiring **cooccur**.

1.0.1 (2025-06-24)
------------------
- Guarded all cooccur examples so they only run when cooccur is installed
- Escaped braces in the CovrgPlot parameter documentation
- Wrapped long-running examples in \donttest{} blocks


1.0 (2023-05-02)
----------------
- Initial release with:
  - `affinity()` for computing co-occurrence affinity metrics and their intervals  
  - `dataprep()` to turn abundance data into presence/absence  
  - `CovrgPlot()` to visualize coverage probability of confidence intervals  
  - `plotgg()` ggplot2 wrapper for affinity outputs  
  - several core internal routines (e.g., `ML.Alpha()`, `AlphInts()`)
  - Example datasets pulled from the **cooccur** package


### Development timeline at GitHub, before CRAN publication

#### May 3, 2023

Package published on CRAN; available at https://cran.r-project.org/web/packages/CooccurrenceAffinity/index.html.

#### Feb 6, 2023

We are pleased to inform you that the interaction issue between our package and BiasedUrn v2.0.8 has been resolved in the latest version, BiasedUrn v2.0.9. We kindly request that you remove any previous versions of BiasedUrn and install v2.0.9 to ensure proper operation of the CooccurrenceAffinity package. We extend our heartfelt thanks to Agner Fox for promptly updating BiasedUrn and addressing these important issues.

#### Jan 27, 2023

After inspecting this issue <https://github.com/kpmainali/CooccurrenceAffinity/issues/6>, we have discovered that the recent revision of our dependency package BiasedUrn is causing R to crash occasionally while running CooccurrenceAffinity. We are actively working to resolve this issue from within our package. In the meantime, we strongly advise against updating BiasedUrn to version 2.0.8. If you have already upgraded, we recommend removing this version and installing the prior version 1.07 as a temporary solution. We will provide updates as soon as the issue is resolved.

#### Nov 4, 2022

Completed writing the package manuscript serving as the vignette. Preprint available: https://www.biorxiv.org/content/10.1101/2022.11.01.514801v1

#### Sept 23, 2022

`CovrgPlot()` revised to generate multipanel plots.

#### Jul 12, 2022

Added function `minmaxAlpha.pFNCH()` to address range inconsistency in BiasedUrn::pFNCHypergeo() for extreme examples. Without this function, BiasedUrn::pFNCHHypergeo() returns inconsistency message for extreme examples like: AlphInts(20,c(204,269,2016), lev=0.9, scal=10). This problem is solved within our package by restricting the range of allowed alpha to the computed (alphmin, alphmax) range.

#### Mar 18, 2022

Fixed a bug in `ML.Alpha()` affecting computation of the upper bound of the alpha MLE. Users before this date should reinstall.

#### Feb–Mar 2022

Substantial revision of existing functions’ descriptions, examples, and arguments. Added new functions for binary matrix analysis and plotting.

#### Dec 2021

This package was released with a basic set of functions in Dec 2021. 
