#' Visualize Individual Genes' Expression Using Various Methods
#'
#' This function provides multiple visualization methods for analyzing individual genes
#' from an expression dataset. Users can generate violin plots, expression heatmaps,
#' correlation heatmaps, ROC curves, and effect size plots.
#'
#' @param data A data frame where rows correspond to genes and columns correspond to samples.
#'   The row names should be gene identifiers.
#' @param metadata A data frame containing sample metadata. The first column should contain
#'   sample identifiers matching the column names of `data`.
#' @param genes A character vector of gene names to be plotted. If more than 20 genes are provided, a warning is issued.
#' @param method A character string specifying which visualization method to use.
#'   Options are `"violins"`, `"expr_heatmap"`, `"corr_heatmap"`, `"roc"`, `"effectsize"`, or `"all"`.
#'   Default is `"violins"`. If `"all"` is selected, all methods will be applied.
#' @param plot A logical value indicating whether to display the generated plots (`TRUE`) or not (`FALSE`).
#'   Default is `TRUE`.
#' @param GroupingVariable A character string specifying the column name in `metadata` to be used
#'   for grouping samples on the x-axis in violin plots.
#' @param ColorVariable A character string specifying the column name in `metadata` used
#'   to color the points in violin plots.
#' @param ColorValues An optional named vector mapping unique values of `ColorVariable` to specific colors.
#' @param pointSize A numeric value specifying the size of jittered points in violin plots. Default is `2`.
#' @param title An optional character string specifying the title of the plot.
#' @param widthTitle An optional integer specifying the maximum width of the title before inserting line breaks. Default is `16`.
#' @param y_limits A numeric vector of length 2 specifying the limits of the y-axis in violin plots.
#'   If `NULL`, the y-axis is adjusted automatically. Default is `NULL`.
#' @param legend_nrow A numeric value specifying the number of rows in the legend for violin plots.
#'   If `NULL`, the default number of rows is determined automatically.
#' @param ncol An optional numeric value specifying the number of columns in the grid layout for violin plots.
#'   If `NULL`, it is automatically determined.
#' @param nrow An optional numeric value specifying the number of rows in the grid layout for violin plots.
#'   If `NULL`, it is automatically determined.
#' @param divide A character string specifying a column in `metadata` used to divide the violin plots into facets.
#'   Default is `NULL`.
#' @param invert_divide A logical value indicating whether to invert the facet layout of violin plots.
#'   Default is `FALSE`.
#' @param xlab An optional character string specifying the label for the x-axis in violin plots.
#'   Default is the name of `GroupingVariable`.
#' @param colorlab An optional character string specifying the legend title for the color variable.
#'   Default is `NULL`.
#'
#' @return A named list where each selected method (e.g., `"violins"`, `"expr_heatmap"`, etc.) is a key.
#'   Each element in this list is itself a named list with the following structure:
#'   \describe{
#'     \item{`data`}{A data frame or matrix containing the transformed data used for plotting, formatted according to the selected method.}
#'     \item{`plot`}{The generated plot object (e.g., a ggplot2 object) corresponding to the visualization.}
#'     \item{`aux`}{Additional parameters, specific of each method.}
#'   }
#'   If `plot = TRUE`, the plots are also displayed in the output.
#'
#' @details
#' This function calls multiple sub-functions based on the selected visualization method:
#'
#' - `IndividualGenes_Violins()`: Generates violin plots for gene expression (log2 transformed).
#' - `IndividualGenes_ExprHeatmap()`: Creates a heatmap of expression values.
#' - `IndividualGenes_CorrHeatmap()`: Computes and visualizes a correlation heatmap.
#' - `IndividualGenes_ROC()`: Plots receiver operating characteristic (ROC) curves.
#' - `IndividualGenes_EffectSize()`: Computes and visualizes effect sizes.
#'
#' The function returns a list of plots. If `method = "all"`, all available visualizations are included.
#'
#' @export
VisualiseIndividualGenes <- function(data,
                                     metadata,
                                     genes,
                                     method = c("violins","expr_heatmap","corr_heatmap","roc","effectsize","all"),
                                     plot = TRUE,
                                     GroupingVariable = NULL,
                                     ColorVariable = NULL,
                                     ColorValues = NULL,
                                     pointSize = 2,
                                     title = NULL,
                                     widthTitle = 16,
                                     y_limits = NULL,
                                     legend_nrow = NULL,
                                     ncol = NULL,
                                     nrow = NULL,
                                     divide = NULL,
                                     invert_divide = FALSE,
                                     xlab = NULL,
                                     colorlab = NULL) {

  method <- match.arg(method)  # Validate method input

  if (!is.data.frame(data)) stop("Error: data must be a data-frame")
  if (!is.null(metadata) && !is.data.frame(metadata)) stop("Error: metadata must be a data-frame")
  if (!is.list(gene_sets)) stop("Error: gene_sets must be a list")

  if (length(genes) >20) warning("Number of genes is too high. Consider a smaller set of genes.")


  ResultsList <- list()

  if (method == "violins" | method == "all") {
    ResultsList[["violins"]] <- IndividualGenes_Violins(data = data, metadata = metadata,
                                                        genes = genes, GroupingVariable = GroupingVariable,
                                                        ColorVariable = ColorVariable, ColorValues = ColorValues,
                                                        pointSize = pointSize, title = title,
                                                        widthTitle = widthTitle, y_limits = y_limits,
                                                        legend_nrow = legend_nrow, ncol = ncol,
                                                        nrow = nrow, divide = divide, invert_divide = invert_divide,
                                                        xlab = xlab, colorlab = colorlab, plot = plot)
  }

  if (method == "expr_heatmap" | method == "all") {
    ResultsList[["expr_heatmap"]] <- IndividualGenes_ExprHeatmap(data = data, metadata = metadata,
                                                                 genes = genes, plot = plot)
  }

  if (method == "corr_heatmap" | method == "all") {
    ResultsList[["corr_heatmap"]] <- IndividualGenes_CorrHeatmap(data = data, metadata = metadata,
                                                                 genes = genes, plot = plot)
  }

  if (method == "roc" | method == "all") {
    ResultsList[["roc"]] <- IndividualGenes_ROC(data = data, metadata = metadata,
                                                genes = genes, plot = plot)
  }

  if (method == "effectsize" | method == "all") {
    ResultsList[["effectsize"]] <- IndividualGenes_EffectSize(data = data, metadata = metadata,
                                                              genes = genes, plot = plot)
  }

  return(ResultsList)
}
