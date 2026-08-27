# markeR 1.3.1 (8 Jul 2026)

## Minor Changes
- Updated the citation to the markeR paper published in NAR Genomics and Bioinformatics.
- Updated the contact email address.
- Fixed `runGSEA` function.

# markeR 1.2 (29 Apr, 2026) 

- Bioconductor Release 3.23

# markeR 1.1.2 (12 Mar, 2026) 

## Minor Changes
- Moved Python bridge scripts from `inst/python/` to a top-level `python/` 
  directory, as these are supplementary scripts not part of the R package itself.
- Added `requirements.txt` to the `python/` directory listing all needed 
  Python dependencies (`rpy2`, `pandas`, `numpy`, and optionally 
  `ipython` and `jupyter`) for easier environment setup. 
- Removed redundant code snippets from the Python bridge scripts.

# markeR 1.1.1 (11 Mar, 2026) 
 
  - Added `p.adjust.method` parameter across all functions performing or
  depending on multiple testing correction, allowing users to specify
  any correction method supported by `stats::p.adjust()`, beyond the default 
  Benjamini-Hochberg FDR.
- Added Python bridge scripts in `inst/python/` for users who wish to call
  markeR from a Python environment via `rpy2`. Includes a tutorial workflow
  script and a generic command-line wrapper capable of invoking any exported
  markeR function. See `inst/python/README.md` for installation and usage.

# markeR 1.0.0 (31 Oct, 2025)

- Official Bioconductor Release.

# markeR 0.99.5 (17 Sep, 2025)

- Minor fix in `.onAttach()` to avoid errors when checking `ggplot2` version and ensure the startup warning works correctly.

# markeR 0.99.4 (17 Sep, 2025)

## General
- Addressed feedback from the Bioconductor review process with updates to documentation and vignette style.  

## Documentation and vignette
- Updated vignette style to **Bioconductor’s BiocStyle** with automatic table of contents.  
- Improved vignette content with small corrections.
- Revised dataset documentation by adding explicit `usage: data(object)` entries.  

## Functions
- Updated `geneset_similarity()` color handling: replaced the single `color_values` parameter with three new parameters — `color`, `neutral_color`, and `cold_color`, for more interpretable visualization.  

# markeR 0.99.3 (21 Aug, 2025)

## Package size and structure
- Reduced package size below the 5 MB limit by converting long vignettes into `pkgdown` articles and keeping only a shorter vignette in the package.
- Moved `inst/Paper` to a dedicated `paper` branch for better repository organization.
- Removed unnecessary `LICENSE` file (already declared in `DESCRIPTION`).

## Documentation
- Added a concise main vignette (`markeR`) with installation, introduction, and a basic workflow.
- Converted three longer vignettes into `pkgdown` articles (linked at the end of the main vignette).
- Added runnable examples for `VariableAssociation`. 

## NAMESPACE and dependencies
- Replaced broad imports with `importFrom()` for most packages (except `ggplot2`, retained as full import).
- Removed unused `patchwork` import.
- Added missing imports from `stats` and `grDevices` to resolve `R CMD check` notes.

## Code quality
- Replaced all `sapply()` calls with `vapply()`.
- Replaced `1:...` usage with `seq_len()` or `seq_along()`.
- Standardized assignment to `<-` instead of `=`.
- Fixed some redundant `stop()`/`warning()` conditions to provide clearer input validation.
- Addressed “no visible binding” notes by using `.data$` or `utils::globalVariables()`.

# markeR 0.99.2 (23 Jul, 2025)

* Minor fixes in documentation

# markeR 0.99.1 (23 Jul, 2025)

* Fix documentation (invalid characters, deep nesting, missing value in data)
* Remove citation, given that a DOI is not yet available
* Removed unwanted files from the repository 

# markeR 0.99.0 (18 Jul, 2025)

* First submission to Bioconductor

# markeR 0.9.5 (18 Jul, 2025)

* Added `VisualiseIndividualGenes()` wrapper to unify individual gene visualisation functions (`ExpressionHeatmap`, `ROCandAUCplot`, etc.) under a single, user-friendly interface.
* Ensured all data arguments are data frames for consistency across functions.
* Minor bug fix: corrected p-value rounding in `PlotScores`

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
