## CooccurrenceAffinity 1.0.1 (2025-06-24)

- Guarded all cooccur examples so they only run when cooccur is installed
- Escaped braces in the CovrgPlot parameter documentation
- Wrapped long-running examples in \donttest{} blocks


## CooccurrenceAffinity 1.0 (2023-05-02)

- Initial release with:
  - `affinity()` for computing co-occurrence affinity metrics and their intervals  
  - `dataprep()` to turn abundance data into presence/absence  
  - `CovrgPlot()` to visualize coverage probability of confidence intervals  
  - `plotgg()` ggplot2 wrapper for affinity outputs  
  - several core internal routines (e.g., `ML.Alpha()`, `AlphInts()`)
  - Example datasets pulled from the **cooccur** package
