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
#' @section Features:
#' - Flexible tools for computing gene signature scores.
#' - Visualization functions to compare signatures across conditions.
#' - Machine-learning models for phenotype classification.
#' - Tools for validating and refining gene signatures.
#'
#'
#' @keywords internal
#' @docType package
#' @name markeR
#' @aliases markeR-package
"_PACKAGE"
