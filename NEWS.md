2.0.0 (2026-01-28)
------------------

### plotgg()

- Changed default behavior: entities whose values are entirely `NA` for the selected variable are now removed from both axes by default (`drop.empty = TRUE`).
  This avoids plotting empty rows or columns caused by degenerate distributions or missing data.
  Set `drop.empty = FALSE` to retain all entities and reproduce the legacy full-grid layout.
  If fewer than two entities remain after dropping empty rows/columns, `plotgg()` now stops with an informative error rather than producing an empty or misleading plot.

- New argument `sig.only`:
  When enabled, values of the selected variable are masked to `NA` wherever `p_value` exceeds a specified cutoff, so only statistically significant tiles are displayed.

  - `sig.only = TRUE` applies a default cutoff of `0.05`
  - `sig.only = <numeric>` uses a user-specified p-value cutoff
    This option requires a `p_value` column in `data$all`.
    When `variable = "p_value"`, p-values above the cutoff are masked to `NA`, allowing direct visualization of significant p-values.

- Stable legend scaling under significance filtering:
  When `sig.only` is used, legend limits and color scaling (`"datarange"` or `"balanced"`) are computed from the original, unfiltered data rather than the masked tiles.
  This prevents legend collapse and ensures visual comparability between full and significance-filtered plots.

- Stricter input validation and supported-variable enforcement:
  `plotgg()` now requires `data$occur_mat` to be present and uses it as the authoritative source for entity ordering on both axes.
  Only explicitly supported variables can be plotted; interval and confidence-related columns are rejected, and unsupported or legacy variable names are no longer accepted.

- Revised color-scale logic:
  Color scales are now variable-aware rather than uniform:

  - `alpha_mle` uses a diverging color scale centered at zero (null expectation)
  - All other variables use sequential color scales
  - Variables bounded in `[0,1]` (`p_value`, `jaccard`, `sorensen`, `simpson`) use fixed `[0,1]` limits under `"balanced"` scaling
  - For `p_value`, color mapping is reversed so smaller p-values appear more intense
  - `obs_cooccur_X` and `exp_cooccur` share a common color scale to enable direct visual comparison
    When `col` is not supplied, appropriate defaults are generated automatically.

- Improved legend titles and diagnostics:
  Legend titles are now mapped to human-readable labels (including multi-line titles and mathematical notation where appropriate) rather than raw column names.
  Standardized legend notes are printed describing palette type, effective scaling mode, applied limits, and whether masking affects the plot.

- Improved value labeling:
  A new argument `text.col` controls the color of numeric labels printed on tiles.
  Values continue to be printed automatically when the number of plotted entities is ≤ 20, or when `show.value = TRUE`.



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
