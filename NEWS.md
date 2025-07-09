# markeR 0.9.4 (09 Jul, 2025)

* Minor bug fix: corrected p-value rounding in `PlotScores`

# markeR 0.9.3 (03 Jul, 2025)

* Updated documentation and internal code to meet Bioconductor submission guidelines.
* Fixed minor bugs across multiple functions.
* Added unit tests using `testthat` for all exported functions.
* Reduced size of demo data to improve package loading time and final size. 

# markeR 0.9.2 (25 Jun, 2025)

* Fixed broken links in README and vignettes
* Added GitHub Actions workflows:
  - `R-CMD-check` 
  - Matrix-based check for minimal supported `R` versions
* Unified `VariableAssociation()` function by modularly integrating `GSEA_VariableAssociation()` and `Score_VariableAssociation()`
* Added scripts to fully reproduce all analyses from the original `markeR` manuscript (`inst/Paper`)

# markeR 0.9.1 (20 Jun, 2025)

* Added package logo  
* Updated and simplified README file with concise installation instructions and main usage workflow
* Creation of dedicated tutorials:
  - **Benchmarking mode** 
  - **Discovery mode** 
  - **Gene set similarity** 
* Improved function documentation  
* Minor bug fixes and internal cleanup 
* Published full codebase for reproducing analyses shown in markeR's paper

# markeR 0.9.0 (21 Apr, 2025)

* Initial release of the package.
* Implementation of score-based and enrichment-based methods to evaluate gene signatures as phenotype markers.
* Visualization of individual genes' expression, scores, and enrichment results
* Add pkgdown documentation site: https://diseasetranscriptomicslab.github.io/markeR/
