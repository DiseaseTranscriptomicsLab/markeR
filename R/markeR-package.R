#' markeR: A Toolkit for Evaluating Gene Signatures as Phenotype Markers
#'
#' The **markeR** package provides a suite of methods for analyzing gene signatures and evaluating
#' their relevance in marking specific phenotypes. It includes tools for score-based evaluation,
#' enrichment analysis, and machine-learning classification. The package offers a variety of visualization
#' functions to assist in interpreting gene expression patterns and their association with phenotypes.
#'
#' @section Methods:
#' The package implements three main approaches for gene signature analysis:
#'
#' \strong{1. Score-based Methods:} These approaches compute a numerical score representing the
#' expression of a gene signature per sample.
#' - **Log2 Median-Centered Scoring:** Measures signature expression by centering expression
#'   values relative to their median.
#' - **Single-Sample Gene Set Enrichment Analysis (ssGSEA):** Calculates an enrichment score
#'   for a given gene set in each sample.
#' - **Ranking Method:** Computes scores based on gene rankings (e.g., using statistical methods).
#'
#' \strong{2. Enrichment-based Methods:} These methods evaluate whether a gene signature is
#' significantly enriched in a ranked gene list.
#' - **Gene Set Enrichment Analysis (GSEA):** A rank-based method for testing enrichment
#'   significance in predefined gene sets.
#'
#' \strong{3. Classification-based Methods:} These approaches classify samples into phenotypic
#' groups based on gene signature expression.
#' - **Random Forest Classification:** Uses decision trees to predict sample labels based on gene
#'   signature scores.
#'
#' @section Main Functions:
#' The current release of **markeR** includes a variety of functions, organized into four main categories:
#'
#' \subsection{1. Score-based Methods:}
#' \itemize{
#'   \item \code{\link{CalculateScores}}: Calculates gene signature scores for each sample using either the ssGSEA,
#'         log2 median-centered, or ranking method.
#'   \item \code{\link{PlotScores}}: Visualizes the calculated scores across conditions using violin plots, density plots,
#'         or heatmaps, depending on the chosen parameters.
#'   \item \code{\link{FDR_Simulation}}: Computes false discovery rates using random gene sets for each signature.
#' }
#'
#' \subsection{2. Enrichment-based Methods:}
#' \itemize{
#'   \item \code{\link{calculateDE}}: Performs differential expression analysis to identify genes associated with phenotypes.
#'   \item \code{\link{plotVolcano}}: Generates volcano plots to visualize differentially expressed genes.
#'   \item \code{\link{runGSEA}}: Performs Gene Set Enrichment Analysis (GSEA) using \code{fgsea} for each contrast in a list of differential expression results.
#'   \item \code{\link{plotGSEAenrichment}}: Generates enrichment plots for gene sets using the \code{fgsea::plotEnrichment()} function.
#'   \item \code{\link{plotNESlollipop}}: Generates a lollipop plot to visualize Gene Set Enrichment Analysis (GSEA) results.
#'   \item \code{\link{plotCombinedGSEA}}: Creates a scatter plot visualizing multiple GSEA (Gene Set Enrichment Analysis) results across different contrasts.
#' }
#'
#' \subsection{3. Visualization Functions:}
#' \itemize{
#'   \item \code{\link{IndividualGenes_Violins}}: Creates violin plots for individual gene expression across sample groups.
#'   \item \code{\link{CorrelationHeatmap}}: Computes and visualizes a correlation heatmap for a set of genes.
#'   \item \code{\link{ExpressionHeatmap}}: Generates a heatmap for the expression of a set of genes.
#'   \item \code{\link{ROCandAUCplot}}: Computes ROC curves and AUC values for each gene, visualizing the results in plots.
#'   \item \code{\link{CohenDH_IndividualGenes}}: Computes and visualizes Cohen’s d for each gene to assess effect sizes.
#'   \item \code{\link{plotPCA}}: Performs PCA and visualizes the results, highlighting important genes that explain the variance.
#' }
#'
#' \subsection{4. Future Updates:}
#' Future updates will expand the package with additional functions for:
#' \itemize{
#'   \item Classification-based analysis (e.g., training Random Forest classifiers and evaluating performance via ROC curves).
#'   \item Additional gene-level visualization modules to display expression patterns of individual genes within signatures.
#' }
#'
#' @docType package
#' @name markeR
#' @aliases markeR-package
"_PACKAGE"
