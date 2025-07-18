#' VisualiseIndividualGenes: Wrapper for Visualising Individual Genes in Gene Sets
#'
#' This wrapper function helps explore individual gene behavior within a gene set used in \code{markeR}.
#' It dispatches to specific visualisation functions based on \code{plot_type}, supporting various plot types:
#' heatmaps, violin plots, correlation analysis, PCA, ROC/AUC, and effect size heatmaps.
#'
#' @param type Character. Specifies the type of plot to generate.
#'   Must be one of:
#'   \describe{
#'     \item{"violin"}{Violin plots of individual gene expression by group.}
#'     \item{"correlation"}{Correlation heatmap of selected genes.}
#'     \item{"expression"}{Expression heatmap of selected genes.}
#'     \item{"roc"}{ROC plots for classification performance of individual genes.}
#'     \item{"auc"}{AUC plots for classification performance of individual genes.}
#'     \item{"rocauc"}{ROC and AUC plots for classification performance of individual genes.}
#'     \item{"cohend"}{Effect size (Cohen's D) heatmap of individual genes.}
#'     \item{"pca"}{PCA plot of selected genes.}
#'   }
#' @param data Required. Expression data matrix or data frame, with samples as rows and genes as columns.
#' @param genes Required. Character vector of gene names to include in the visualisation.
#' @param metadata Optional. Data frame with sample metadata, required for some plot types (e.g., violin, roc, cohend).
#' @param ... Additional arguments passed to the specific plotting function.
#'
#' @section Additional required arguments (passed via \code{...}) per \code{plot_type}:
#' \describe{
#'   \item{violin}{Requires \code{GroupingVariable} (column name in \code{metadata} for grouping).}
#'   \item{roc, auc, rocauc}{Requires \code{condition_var} (metadata column with condition labels) and \code{class} (positive class label).}
#'   \item{cohend}{Requires \code{condition_var} and \code{class} (same as roc).}
#'   \item{correlation, expression, pca}{No additional mandatory arguments required.}
#' }
#'
#' @return The output of the specific plotting function called, usually a ggplot or ComplexHeatmap object,
#'   often wrapped in a list with additional data.
#'
#' @seealso \code{\link{IndividualGenes_Violins}}, \code{\link{CorrelationHeatmap}},
#'   \code{\link{ExpressionHeatmap}}, \code{\link{ROCandAUCplot}},
#'   \code{\link{CohenD_IndividualGenes}}, \code{\link{plotPCA}}
#'
#' @examples
#' # Example data
#' set.seed(123)
#' expr_data <- matrix(rexp(1000, rate = 1), nrow = 50, ncol = 20)
#' rownames(expr_data) <- paste0("Gene", 1:50)
#' colnames(expr_data) <- paste0("Sample", 1:20)
#'
#' sample_info <- data.frame(
#'   SampleID = colnames(expr_data),
#'   Condition = rep(c("A", "B"), each = 10),
#'   Diagnosis = rep(c("Disease", "Control"), each = 10),
#'   stringsAsFactors = FALSE
#' )
#' rownames(sample_info) <- sample_info$SampleID
#'
#' selected_genes <- row.names(expr_data)[1:5]
#'
#' # Violin plot
#' VisualiseIndividualGenes(
#'   type = "violin",
#'   data = expr_data,
#'   metadata = sample_info,
#'   genes = selected_genes,
#'   GroupingVariable = "Condition",
#'   nrow=1
#' )

# Correlation heatmap
#' VisualiseIndividualGenes(
#'   type = "correlation",
#'   data = expr_data,
#'   genes = selected_genes
#' )
#'
#' # Expression heatmap
#' VisualiseIndividualGenes(
#'   type = "expression",
#'   data = expr_data,
#'   genes = selected_genes
#' )
#'
#' # PCA plot
#' VisualiseIndividualGenes(
#'   type = "pca",
#'   data = expr_data,
#'   genes = selected_genes,
#'   metadata = sample_info,
#'   ColorVariable="Condition"
#'
#' )
#'
#' # ROC plot
#' VisualiseIndividualGenes(
#'   type = "roc",
#'   data = expr_data,
#'   metadata = sample_info,
#'   genes = selected_genes,
#'   condition_var = "Diagnosis",
#'   class = "Disease"
#'
#' )
#'
#' # AUC plot
#' VisualiseIndividualGenes(
#'   type = "auc",
#'   data = expr_data,
#'   metadata = sample_info,
#'   genes = selected_genes,
#'   condition_var = "Diagnosis",
#'   class = "Disease"
#' )
#'
#' # ROC&AUC plot
#' VisualiseIndividualGenes(
#'   type = "rocauc",
#'   data = expr_data,
#'   metadata = sample_info,
#'   genes = selected_genes,
#'   condition_var = "Diagnosis",
#'   class = "Disease"
#' )
#'
#' # Cohen's D plot
#' VisualiseIndividualGenes(
#'   type = "cohend",
#'   data = expr_data,
#'   metadata = sample_info,
#'   genes = selected_genes,
#'   condition_var = "Diagnosis",
#'   class = "Disease"
#' )
#'
#'
#'
#' @export
VisualiseIndividualGenes <- function(type, data, genes, metadata = NULL, ...) {
  args <- list(...)
  data <- as.data.frame(data) # Ensure data is a data frame
  # Define required additional arguments for each plot type (besides data, genes, metadata)
  required_args <- list(
    violin      = c("GroupingVariable"),
    correlation = character(0),
    expression  = character(0),
    roc         = c("condition_var", "class"),
    auc         = c("condition_var", "class"),
    rocauc      = c("condition_var", "class"),
    cohend      = c("condition_var", "class"),
    pca         = character(0)
  )

  type <- tolower(type)
  if (!type %in% names(required_args)) {
    stop("Invalid type. Choose from: ",
         paste(names(required_args), collapse = ", "))
  }

  # Check mandatory arguments
  if (missing(data) || is.null(data)) stop("Argument 'data' is required.")
  if (missing(genes) || is.null(genes)) stop("Argument 'genes' is required.")

  # For plot types that require metadata, check presence
  if (type %in% c("violin", "roc", "cohend") && (missing(metadata) || is.null(metadata))) {
    stop("Argument 'metadata' is required for type '", type, "'.")
  }

  # Check for additional required args passed via ...
  missing_args <- setdiff(required_args[[type]], names(args))
  if (length(missing_args) > 0) {
    stop("Missing required argument(s) for type '", type, "': ",
         paste(missing_args, collapse = ", "))
  }

  # Prepare argument list to pass to specific plotting function
  call_args <- c(list(data = data, genes = genes))
  if (!is.null(metadata)) call_args$metadata <- metadata
  call_args <- c(call_args, args)

  # Dispatch to correct internal function
  switch(type,
         violin      = do.call(IndividualGenes_Violins, call_args),
         correlation = do.call(CorrelationHeatmap, call_args),
         expression  = do.call(ExpressionHeatmap, call_args),
         roc         = do.call(ROCandAUCplot, c(call_args,plot_type="roc")),
         auc         = do.call(ROCandAUCplot, c(call_args,plot_type="auc")),
         rocauc      = do.call(ROCandAUCplot, c(call_args,plot_type="all")),
         cohend      = do.call(CohenD_IndividualGenes, call_args),
         pca         = do.call(plotPCA, call_args))
}
