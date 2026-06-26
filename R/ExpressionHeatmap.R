#' ExpressionHeatmap: Generate an expression heatmap with customizable sample
#' annotations and separate legend positions
#'
#' This function creates a heatmap of Z-score scaled gene expression using the
#' `ComplexHeatmap` package. Genes are displayed as rows and samples as columns.
#' A color annotation bar is added on top based on specified metadata columns.
#' The user can control the position of the heatmap color scale (scale_position)
#' and the annotation legend (legend_position) independently.
#'
#' @param data A numeric expression matrix where rows correspond to genes and
#'   columns to samples.
#' @param metadata A data frame containing metadata for the samples. It must
#'   contain a column named `"Sample"` with sample IDs matching the column names
#'   of `data`.
#' @param genes A character vector of gene names to include in the heatmap.
#' @param annotate.by A character vector of metadata column names to be used for
#'   sample annotations (e.g., \code{c("Condition", "Batch")}). If provided, a
#'   color bar is added on top.
#' @param annotation_colors Optional. A named list where each element
#'   corresponds to an annotation variable and provides a named vector mapping
#'   each unique level to a color. If not provided, default Brewer palettes are
#'   used.
#' @param colorlist A named list specifying the colors for the heatmap (for
#'   scaled expression) with elements `low`, `mid`, and `high`. Default is
#'   \code{list(low = "blue", mid = "white", high = "red")}.
#' @param cluster_rows Logical; whether to cluster rows (default = \code{TRUE}).
#' @param cluster_columns Logical; whether to cluster columns (default =
#'   \code{TRUE}). If \code{FALSE}, the columns are reordered based on the
#'   values in `annotate.by`.
#' @param title A string specifying the main title of the heatmap.
#' @param titlesize Numeric; font size of the heatmap title (default = 20).
#' @param scale_position A character string specifying the position of the
#'   heatmap color scale. Options are \code{"right"} (default), \code{"top"}, or
#'   \code{"bottom"}. The scale legend will adopt a vertical orientation if on
#'   the right and horizontal if on top or bottom.
#' @param legend_position A character string specifying the position of the
#'   annotation legend. Options are \code{"top"} (default), \code{"right"}, or
#'   \code{"bottom"}.
#' @param show_row_names A character string specifying whether row names (genes)
#'   should be displayed.
#' @param show_column_names A character string specifying whether column names
#'   (samples) should be displayed.
#'
#'
#' @return Invisibly returns a list with:
#' \describe{
#'   \item{data}{Scaled expression matrix (Z-scores).}
#'   \item{plot}{Generated ComplexHeatmap object.}
#' }
#'
#' @examples
#' # Simulate gene expression data (genes as rows, samples as columns)
#' set.seed(1)
#' expr <- matrix(rnorm(25), nrow = 5, ncol = 5)
#' rownames(expr) <- paste0("Gene", 1:5)
#' colnames(expr) <- paste0("Sample", 1:5)
#'
#' # Simulate metadata for samples
#' metadata <- data.frame(
#'   Sample = colnames(expr),
#'   Condition = rep(c("A", "B"), length.out = 5),
#'   Batch = rep(c("X", "Y"), length.out = 5),
#'   stringsAsFactors = FALSE
#' )
#'
#' # Define annotation colors for the metadata variables
#' annotation_colors <- list(
#'   Condition = c(A = "orange", B = "purple"),
#'   Batch = c(X = "green", Y = "blue")
#' )
#'
#' # Generate the expression heatmap
#' ExpressionHeatmap(
#'   data = expr,
#'   metadata = metadata,
#'   genes = rownames(expr),
#'   annotate.by = c("Condition", "Batch"),
#'   annotation_colors = annotation_colors,
#'   cluster_columns = FALSE,
#'   title = "Demo Expression Heatmap",
#'   scale_position = "right",
#'   legend_position = "top",
#'   titlesize = 14
#' )
#'
#' @importFrom grid gpar
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom ComplexHeatmap draw
#' @importFrom ComplexHeatmap HeatmapAnnotation
#' @importFrom circlize colorRamp2
#' @importFrom RColorBrewer brewer.pal
#'
#' @export
ExpressionHeatmap <- function(data, metadata = NULL, genes, annotate.by = NULL,
                              annotation_colors = NULL,
                              colorlist = list(low = "blue", mid = "white", high = "red"),
                              cluster_rows = TRUE, cluster_columns = TRUE,
                              title = NULL, titlesize = 20,
                              fontsize = 10,
                              scale_position = c("right", "top", "bottom"),
                              legend_position = c("top", "right", "bottom"),
                              show_row_names=TRUE,
                              show_column_names=FALSE) {
  data <- as.data.frame(data) # Ensure data is a data frame
  # Ensure the scale_position and legend_position arguments are matched correctly
  scale_position <- match.arg(scale_position)
  legend_position <- match.arg(legend_position)

  # Define scale legend parameters
  if (scale_position == "right") {
    scale_side <- "right"
    scale_direct <- "vertical"
  } else {
    scale_side <- scale_position  # "top" or "bottom"
    scale_direct <- "horizontal"
  }

  # Define annotation legend parameters
  if (legend_position == "top") {
    annot_side <- "top"
    annot_direct <- "horizontal"# horizontal currently not working
  } else if (legend_position == "bottom") {
    annot_side <- "bottom"
    annot_direct <- "horizontal"# horizontal currently not working
  } else {
    annot_side <- "right"
    annot_direct <- "vertical"
  }

  # Subset data to specified genes
  data <- data[rownames(data) %in% genes, , drop = FALSE]

  # Ensure metadata is provided if annotate.by is specified
  if (!is.null(annotate.by) && is.null(metadata)) stop("annotate.by is specified but metadata is NULL.")

  # Ensure metadata matches sample order if provided
  if (!is.null(metadata)) {
    colnames(metadata)[1] <- "Sample"
    rownames(metadata) <- metadata$Sample
    metadata <- metadata[colnames(data), , drop = FALSE]
  }


  # If column clustering is disabled and annotations are provided, reorder
  # columns based on annotation values
  if (!cluster_columns && !is.null(metadata) && !is.null(annotate.by)) {
    annotate.by <- annotate.by[annotate.by %in% colnames(metadata)]
    if (length(annotate.by) > 0) {
      order_idx <- do.call(order, metadata[, annotate.by, drop = FALSE])
      data <- data[, order_idx, drop = FALSE]
      metadata <- metadata[order_idx, , drop = FALSE]
    }
  }

  # Scale expression data (Z-score per gene)
  data_scaled <- t(scale(t(data)))

  # Create sample annotations if annotate.by is provided
  sample_annotation <- NULL
  if (!is.null(metadata) && !is.null(annotate.by)) {
    ann_data <- metadata[, annotate.by, drop = FALSE]
    ann_colors <- list()
    for (var in annotate.by) {
      if (!is.null(annotation_colors) && var %in% names(annotation_colors)) {
        ann_colors[[var]] <- annotation_colors[[var]]
      } else {
        unique_vals <- unique(ann_data[[var]])
        palette_choices <- c("Set1", "Set2", "Set3", "Paired", "Dark2",
                             "Set1", "Set2", "Set3", "Paired", "Dark2")

        n_colors <- length(unique_vals)
        n_colors <- ifelse(n_colors < 3, 3, min(n_colors, 9))  # enforce min = 3

        pal <- RColorBrewer::brewer.pal(
          n_colors,
          palette_choices[(which(annotate.by == var) - 1) %% length(palette_choices) + 1]
        )

        # Optionally subset colors if n_colors > length(unique_vals)
        pal <- pal[seq_along(unique_vals)]


        names(pal) <- unique_vals
        ann_colors[[var]] <- pal
      }
    }
    sample_annotation <- ComplexHeatmap::HeatmapAnnotation(df = ann_data, col = ann_colors,
                                           annotation_legend_param = list(
                                             direction = annot_direct,
                                             title_gp = grid::gpar(fontsize = 12),
                                             labels_gp = grid::gpar(fontsize = 10)
                                           ))
  }

  # Define color scale for the heatmap (for scaled expression)
  col_fun <- circlize::colorRamp2(c(-2, 0, 2), unname(colorlist))

  # Build legend parameters for the heatmap color scale using scale_position
  scale_legend_param <- list(
    direction = scale_direct,
    title_gp = grid::gpar(fontsize = 12),
    labels_gp = grid::gpar(fontsize = 10)
  )

  # Create the heatmap using ComplexHeatmap
  ht <- ComplexHeatmap::Heatmap(data_scaled,
                                name = "Expression \nZ-score",
                                col = col_fun,
                                cluster_rows = cluster_rows,
                                cluster_columns = cluster_columns,
                                show_row_names = show_row_names,
                                show_column_names = show_column_names,
                                row_names_gp    = grid::gpar(fontsize = fontsize),
                                column_names_gp = grid::gpar(fontsize = fontsize),
                                top_annotation = sample_annotation,
                                column_title = title,
                                column_title_gp = grid::gpar(fontsize = titlesize),
                                heatmap_legend_param = scale_legend_param)

  # Return the Heatmap object (not yet drawn) together with the legend-side
  # parameters so the caller (renderPlot) can call draw() fresh inside a
  # proper graphics device, avoiding stale viewport dimensions.
  invisible(list(data = data_scaled, plot = ht,
                 scale_side = scale_side, annot_side = annot_side))
}
