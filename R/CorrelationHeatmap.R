#' CorrelationHeatmap: Generate correlation heatmaps with optional grouping
#'
#' This function generates correlation heatmaps using the `ComplexHeatmap`
#' package. It allows users
#' to compute correlation matrices for a set of genes and visualize them in a
#' heatmap. If a grouping
#' variable is provided (`separate.by`), multiple heatmaps are created, each
#' corresponding to a different
#' level of the grouping variable.
#'
#' @param data A numeric counts data frame where rows correspond to genes and
#' columns to samples.
#' @param metadata A data frame containing metadata. Required if `separate.by`
#' is specified.
#' @param genes A character vector of gene names to be included in the
#' correlation analysis.
#' @param separate.by A character string specifying a column in `metadata` to
#' separate heatmaps by (e.g., "Condition").
#' @param method Correlation method: `"pearson"` (default), `"spearman"`, or
#' `"kendall"`.
#' @param colorlist A named list specifying the colors for the heatmap (`low`,
#' `mid`, `high`), corresponding to the limits of the colorscale.
#' @param limits_colorscale A numeric vector of length 3 defining the limits
#' for the color scale (default: min, 0, max).
#' @param widthTitle Numeric value controlling the width of the plot title.
#' Default is `16`.
#' @param title A string specifying the main title of the heatmap(s).
#' @param cluster_rows Logical; whether to cluster rows (default = `TRUE`).
#' @param cluster_columns Logical; whether to cluster columns (default =
#' `TRUE`).
#' @param detailedresults Logical; if `TRUE`, additional analysis results are
#' stored in the output list (default = `FALSE`).
#' @param legend_position Character; position of the legend (`"right"` - default -
#' or `"top"`).
#' @param titlesize Numeric; font size of the heatmap title (default = `20`).
#' @param show_row_names A character string specifying whether row names (genes) should be displayed.
#' @param show_column_names A character string specifying whether column names (samples) should be displayed.
#'
#' @return A list containing:
#'   \describe{
#'     \item{`data`}{Correlation matrices for each condition (or a single matrix
#'      if `separate.by = NULL`).}
#'     \item{`plot`}{The generated heatmap object(s).}
#'     \item{`aux`}{A list containing additional analysis results if
#'     `detailedresults = TRUE`.
#'       \describe{
#'         \item{If `separate.by` is specified:}{
#'           A list where each element corresponds to a different condition.
#'           Each sublist contains:
#'           \itemize{
#'             \item `method`: The correlation method used.
#'             \item `corrmatrix`: The computed correlation matrix for that
#'             condition.
#'             \item `metadata`: The subset of metadata corresponding to the
#'             condition.
#'             \item `heatmap`: The `ComplexHeatmap` object before being drawn.
#'           }
#'         }
#'         \item{If `separate.by = NULL` (single heatmap case):}{
#'           A list containing:
#'           \itemize{
#'             \item `method`: The correlation method.
#'             \item `corrmatrix`: The computed correlation matrix.
#'           }
#'         }
#'       }
#'     }
#'   }
#'
#' @examples
#' \dontrun{
#' data_matrix <- matrix(rnorm(100), nrow = 10, ncol = 10)
#' rownames(data_matrix) <- paste0("Gene", 1:10)
#' colnames(data_matrix) <- paste0("Sample", 1:10)
#'
#' # Basic usage
#' result <- CorrelationHeatmap2(data_matrix, genes = rownames(data_matrix))
#'
#' # Using metadata to separate by condition
#' metadata <- data.frame(Sample = colnames(data_matrix),
#'                        Condition = rep(c("A", "B"), each = 5))
#' result <- CorrelationHeatmap2(data_matrix, metadata, genes =
#' rownames(data_matrix), separate.by = "Condition")
#' }
#'
#' @importFrom grid gpar
#' @importFrom grid grid.text
#' @importFrom grid unit
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom ComplexHeatmap draw
#' @importFrom circlize colorRamp2
#' @importFrom stats cor
#'
#' @export
CorrelationHeatmap <- function(data, metadata = NULL, genes, separate.by = NULL,
                                method = c("pearson","spearman","kendall"), colorlist = list(low = "blue", mid = "white", high = "red"),
                                limits_colorscale = NULL, widthTitle = 16, title = NULL,
                                cluster_rows = TRUE, cluster_columns = TRUE,
                                detailedresults = FALSE, legend_position = c("right", "top"), titlesize=20,
                               show_row_names=TRUE,
                               show_column_names=TRUE) {


  # Choose legend position: "side" (vertical) or "top" (horizontal)
  legend_position <- match.arg(legend_position)
  method <- match.arg(method)
  if (legend_position == "right") {
    leg_side <- "right"
    leg_direction <- "vertical"
    if (method == "spearman"){
      title_leg <- "Spearman’s \nCorrelation \nCoefficient (ρ)"
    } else if (method == "pearson"){
      title_leg <- "Pearson's  \nCorrelation \nCoefficient (r)"
    } else if (method == "kendall"){
      title_leg <- "Kendall's  \nCorrelation \nCoefficient (τ)"
    }

  } else {
    leg_side <- "top"
    leg_direction <- "horizontal"
    title_leg <- paste0(method, "'s coefficient")
  }

  resultsList <- list()
  resultsList[["data"]] <- list()
  resultsList[["plot"]] <- list()
  resultsList[["aux"]] <- list()

  # Subset data to selected genes
  #data <- data[rownames(data) %in% genes, , drop = FALSE]
  data <- na.omit(as.data.frame(data[genes,])) # to keep input order

  if (!is.null(separate.by) && is.null(metadata)) {
    stop("separate.by is not NULL but metadata is missing. Please specify metadata.")
  }

  # Helper function: create a heatmap using ComplexHeatmap
  create_heatmap <- function(corrmat, annot_title = NULL, direct="horizontal", titleleg="") {
    col_fun <- circlize::colorRamp2(
      if (is.null(limits_colorscale)) c(min(corrmat), 0, max(corrmat)) else limits_colorscale,
      unname(colorlist)
    )
    ht <- ComplexHeatmap::Heatmap(
      corrmat,
      name = titleleg,  # Shared legend name
      col = col_fun,
      cluster_rows = cluster_rows,
      cluster_columns = cluster_columns,
      show_row_names = show_row_names,
      show_column_names = show_column_names,
      heatmap_legend_param = list(
        title = titleleg,
        direction = direct,
        title_gp = grid::gpar(fontsize = 12),   # Legend title size
        labels_gp = grid::gpar(fontsize = 10),  # Legend label size
        legend_height = grid::unit(2, "cm"),          # Adjust legend height
        legend_width = grid::unit(4, "cm")            # Adjust legend width
      ),
      column_title = annot_title
    )
    return(ht)
  }


  if (!is.null(separate.by)) {
    df_data_merge <- list()
    heatmap_list <- list()
    conditions <- unique(metadata[[separate.by]])

    for (cond in conditions) {
      metadata_subset <- subset(metadata, get(separate.by) == cond)
      # Assume the first column of metadata contains sample IDs
      samples <- metadata_subset[,1]
      data_subset <- data[, samples, drop = FALSE]
      # Check if there are at least 2 samples to compute correlation
      if (ncol(data_subset) < 2) {
        warning(paste("Not enough samples for condition", cond, "to compute correlation. Skipping."))
        next
      }
      data_subset <- log2(data_subset)
      corrmat <- stats::cor(t(data_subset), method = method)
      df_data_merge[[cond]] <- corrmat

      ht <- create_heatmap(corrmat, annot_title = cond, direct=leg_direction, titleleg=title_leg)
      heatmap_list[[cond]] <- ht

      if (detailedresults) {
        resultsList[["aux"]][[cond]] <- list(
          method = method,
          corrmatrix = corrmat,
          metadata = metadata_subset,
          heatmap = ht
        )
      }
    }

    if (length(heatmap_list) == 0) {
      stop("No valid conditions with enough samples found.")
    }

    # Combine heatmaps using the '+' operator to merge legends
    combined_ht <- heatmap_list[[conditions[1]]]
    if (length(heatmap_list) > 1) {
      for (i in 2:length(conditions)) {
        if (!is.null(heatmap_list[[conditions[i]]])) {
          combined_ht <- combined_ht + heatmap_list[[conditions[i]]]
        }
      }
    }

    combined_ht <- ComplexHeatmap::draw(combined_ht,
                                        merge_legend = TRUE,
                                        heatmap_legend_side = leg_side,
                                        annotation_legend_side = leg_side,
                                        column_title = if (!is.null(title)) wrap_title(title, width = widthTitle) else NULL,
                                        column_title_gp = grid::gpar(fontsize = titlesize))



    resultsList[["data"]] <- df_data_merge
    resultsList[["plot"]] <- combined_ht

  } else {
    # Single heatmap (no separate.by)
    data <- log2(data)
    corrmat <- stats::cor(t(data), method = method)
    ht <- create_heatmap(corrmat, direct=leg_direction, titleleg=title_leg)
    ht_drawn <- ComplexHeatmap::draw(ht,
                                     heatmap_legend_side = leg_side,
                                     annotation_legend_side = leg_side,
                                     column_title = if (!is.null(title)) wrap_title(title, width = widthTitle) else NULL,
                                     column_title_gp = grid::gpar(fontsize = titlesize))

    resultsList[["data"]] <- corrmat
    resultsList[["plot"]] <- ht_drawn

    if (detailedresults) {
      resultsList[["aux"]] <- list(method = method, corrmatrix = corrmat)
    }

  }

  invisible(resultsList)
}

