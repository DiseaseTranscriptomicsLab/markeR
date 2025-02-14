#' markeR: A Toolkit for Evaluating Gene Signatures as Phenotype Markers
#'
#' The **markeR** package provides methods to analyze gene signatures and assess their
#' relevance in marking specific phenotypes. It includes tools for score-based evaluation,
#' enrichment analysis, and machine-learning classification.
#'
#' @section Methods:
#' The package implements three main approaches for gene signature analysis:
#'
#' \strong{1. Score-based Methods:} These approaches compute a numerical score representing
#' the expression of a gene signature per sample.
#' - **Log2 Median-Centered Scoring:** Measures signature expression by centering expression
#'   values relative to their median.
#' - **Single-Sample Gene Set Enrichment Analysis (ssGSEA):** Calculates an enrichment score
#'   for a given gene set in each sample.
#'
#' \strong{2. Enrichment-based Methods:} These methods evaluate whether a gene signature
#' is significantly enriched in a ranked gene list.
#' - **Gene Set Enrichment Analysis (GSEA):** A rank-based method for testing enrichment
#'   significance in predefined gene sets.
#'
#' \strong{3. Classification-based Methods:} These approaches classify samples into
#' phenotypic groups based on gene signature expression.
#' - **Random Forest Classification:** Uses decision trees to predict sample labels based
#'   on gene signature scores.
#'
#'
#' @section Main Functions and Future Modules:
#' The current release of **markeR** includes two primary functions for score-based analysis:
#' \itemize{
#'   \item \code{\link{CalculateScores}}: Calculates gene signature scores for each sample using either the ssGSEA or log2 median-centered method.
#'   \item \code{\link{PlotScores}}: Visualizes the calculated scores across conditions using violin plots.
#' }
#' Future updates will expand the package with similar pairs of functions for:
#' \itemize{
#'   \item Enrichment-based analysis (to calculate and visualize GSEA statistics).
#'   \item Classification-based analysis (to train classifiers, for example, using Random Forests, and to evaluate performance, e.g., via ROC curves).
#'   \item Gene-level visualization modules that display the expression patterns of individual genes within the signatures, for example, using heatmaps.
#' }
#'
#' @docType package
#' @name markeR
#' @aliases markeR-package
"_PACKAGE"
