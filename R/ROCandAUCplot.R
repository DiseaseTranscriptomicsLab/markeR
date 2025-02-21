#' ROC and AUC Plot Function
#'
#' This function computes ROC curves and AUC values for each gene based on gene expression data and sample metadata.
#' It can generate ROC plots, an AUC heatmap, or both arranged side‐by‐side.
#'
#' @param data A data frame or matrix containing gene expression data, with genes as rows and samples as columns.
#' @param metadata A data frame containing sample metadata. The first column should contain sample identifiers that match the column names of \code{data}.
#' @param genes A character vector specifying which genes to plot. If \code{NULL} (default), all genes in \code{data} are used.
#'   A warning is issued if more than 30 genes are selected.
#' @param condition_var A character string specifying the column name in \code{metadata} representing the condition of interest.
#'   (Mandatory; no default.)
#' @param class A character string specifying the positive class label for the condition.
#'   (Mandatory; no default.)
#' @param group_var An optional character string specifying the column name in \code{metadata} used for grouping samples (e.g., cell types).
#'   If not provided (\code{NULL}), all samples are treated as a single group.
#' @param plot_type A character string indicating which plot(s) to generate.
#'   Accepted values are \code{"roc"} (only ROC curves), \code{"heatmap"} (only the AUC heatmap), or \code{"all"} (both arranged side-by-side).
#'   Default is \code{"roc"}.
#' @param title An optional character string specifying a custom title for the plot(s). If not provided, default titles are generated.
#' @param widthTitle A numeric value used when wrapping the title. Default is \code{16}.
#' @param roc_params A list of additional parameters for customizing the ROC plot. Possible elements include:
#'   \describe{
#'     \item{\code{nrow}}{An integer specifying the number of rows in the ROC plot grid. If \code{NULL} (default), it is calculated automatically.}
#'     \item{\code{ncol}}{An integer specifying the number of columns in the ROC plot grid. If \code{NULL} (default), it is calculated automatically.}
#'     \item{\code{colors}}{A named vector of colors for the different groups. If \code{NULL} (default), a default color palette is generated.}
#'   }
#' @param heatmap_params A list of additional parameters for customizing the AUC heatmap. Possible elements include:
#'   \describe{
#'     \item{\code{cluster_rows}}{Logical; if \code{TRUE} (default), rows are clustered.}
#'     \item{\code{cluster_columns}}{Logical; if \code{TRUE} (default), columns are clustered.}
#'     \item{\code{col}}{A vector of length 2 of colors to be used for the minimum and maximum values of the color scale. Defaults to \code{c("#FFFFFF", "#21975C")}.}
#'     \item{\code{limits}}{A numeric vector of length 2 specifying the minimum and maximum values for the color scale.
#'       If not provided, defaults to \code{c(0.5, 1)}.}
#'     \item{\code{name}}{A character string for the legend title of the color scale. Default is \code{"AUC"}.}
#'     \item{\code{row_names_gp}}{Optional graphical parameters for row names (passed to \pkg{ComplexHeatmap}).}
#'     \item{\code{column_names_gp}}{Optional graphical parameters for column names (passed to \pkg{ComplexHeatmap}).}
#'   }
#' @param commomplot_params A list of parameters for customizing the layout of the combined plot when \code{plot_type = "all"}.
#'   Possible elements include:
#'   \describe{
#'     \item{\code{widths}}{A numeric vector specifying the relative widths of the ROC and heatmap panels.}
#'     \item{\code{heights}}{A numeric vector specifying the relative heights of the panels.}
#'   }
#'
#' @return Invisibly returns a list containing:
#'   \describe{
#'     \item{\code{roc_plot}}{The \pkg{ggplot2} object of the ROC curves (if generated).}
#'     \item{\code{heatmap}}{The \pkg{ComplexHeatmap} object (if generated).}
#'     \item{\code{combined}}{The combined grid arrangement (if \code{plot_type = "all"}).}
#'     \item{\code{auc_values}}{A data frame with the calculated AUC values.}
#'   }
#'
#' @details
#' The function processes gene expression data and metadata to compute ROC curves and AUC values for each gene.
#' Depending on the value of \code{plot_type}, it produces ROC plots (using \pkg{ggplot2}), an AUC heatmap (using \pkg{ComplexHeatmap}),
#' or both arranged side-by-side (using \pkg{gridExtra}). When \code{group_var} is not provided, all samples are treated as a single group.
#'
#' @examples
#' \dontrun{
#'   # Example: Generate both ROC plots and heatmap with custom parameters
#'   result <- ROCandAUCplot(data = myData,
#'                           metadata = myMetadata,
#'                           genes = c("Gene1", "Gene2", "Gene3"),
#'                           condition_var = "Genotype",
#'                           class = "Mutant",
#'                           group_var = "CellType",
#'                           plot_type = "all",
#'                           title = "My Custom Title",
#'                           roc_params = list(ncol = 3),
#'                           heatmap_params = list(limits = c(0.4, 1)),
#'                           commomplot_params = list(widths = c(1, 1)))
#' }
#'
#' @export
ROCandAUCplot <- function(data, metadata,
                          genes = NULL,
                          condition_var,
                          class,
                          group_var = NULL,
                          plot_type = "roc",  # "roc", "heatmap", or "all"
                          title = NULL,       # Custom plot title (optional)
                          widthTitle = 16,
                          roc_params = list(),     # ROC-specific parameters
                          heatmap_params = list(),
                          commomplot_params = list()) {  # Additional parameters for combined plots

  # If genes are not specified, use all available genes
  if (is.null(genes)) {
    genes <- rownames(data)
  } else {
    genes <- genes[genes %in% row.names(data)]
  }

  if (!all(colnames(data) %in% metadata[, 1]))
    stop("Not all samples in the data are described in the metadata.")

  # Subset and order metadata based on data
  row.names(metadata) <- metadata[, 1]
  colnames(metadata)[1] <- "sample"
  metadata <- metadata[colnames(data), ]

  data <- data[genes,]

  # Convert count data to log2 scale (avoid log(0) issues by adding 1)
  count_data_log2 <- as.data.frame(t(log2(data + 1)))
  count_data_log2 <- cbind(sample = row.names(count_data_log2), count_data_log2)
  row.names(count_data_log2) <- NULL

  # Merge metadata with count data
  data_roc <- merge(metadata, count_data_log2, by = "sample")

  if (length(genes) > 30)
    warning("Too many genes selected. Consider reducing the number.")

  # If group_var is provided, ensure it exists; otherwise, treat all samples as one group.
  if (!is.null(group_var)) {
    if (!group_var %in% colnames(metadata)) {
      stop(paste("Error: The specified group_var", group_var, "is not found in the metadata."))
    }
    groups <- unique(data_roc[[group_var]])
  } else {
    data_roc$Group <- "All"
    groups <- "All"
    group_var <- "Group"
  }

  roc_df <- data.frame()
  auc_values <- data.frame(Gene = character(), Group = character(), AUC = numeric())

  for (gene in genes) {
    for (group in groups) {

      # Subset data for the current group
      subset_data <- data_roc[data_roc[[group_var]] == group, ]

      # Convert condition to binary labels (1 for class, 0 for others)
      subset_data$Binary_Label <- as.numeric(subset_data[[condition_var]] == class)

      # Compute ROC curve
      roc_obj <- pROC::roc(subset_data$Binary_Label, subset_data[[gene]], direction = "<", quiet = TRUE)

      # Adjust AUC if needed (ensure AUC is always ≥ 0.5)
      auc_value <- pROC::auc(roc_obj)
      if (auc_value < 0.5) {
        auc_value <- 1 - auc_value
        roc_obj <- pROC::roc(subset_data$Binary_Label, subset_data[[gene]], direction = ">", quiet = TRUE)
      }

      # Store ROC curve points
      roc_points <- data.frame(
        FPR = rev(1 - roc_obj$specificities),
        TPR = rev(roc_obj$sensitivities),
        Gene = gene,
        Group = group
      )
      roc_df <- rbind(roc_df, roc_points)

      # Store AUC values with facet labels
      auc_values <- rbind(auc_values, data.frame(Gene = gene, Group = group, AUC = auc_value))
    }
  }

  # Check for NAs in auc_values
  if (any(is.na(auc_values$AUC))) {
    stop("AUC values contain NAs. Please check the data.")
  }

  ### Function to generate the ROC plot
  makeROCplot <- function() {
    # Set default ROC parameters if missing
    roc_defaults <- list(
      nrow = NULL,
      ncol = NULL,
      colors = NULL  # Custom colors
    )
    roc_params_local <- utils::modifyList(roc_defaults, roc_params)

    # Determine grid layout
    n <- length(unique(roc_df$Gene))
    if (is.null(roc_params_local$ncol) && is.null(roc_params_local$nrow)) {
      roc_params_local$ncol <- ceiling(sqrt(n))
      roc_params_local$nrow <- ceiling(n / roc_params_local$ncol)
    } else if (is.null(roc_params_local$ncol)) {
      roc_params_local$ncol <- ceiling(n / roc_params_local$nrow)
    } else if (is.null(roc_params_local$nrow)) {
      roc_params_local$nrow <- ceiling(n / roc_params_local$ncol)
    }

    # Set colors
    if (is.null(roc_params_local$colors)) {
      roc_params_local$colors <- scales::hue_pal()(length(groups))
      base::names(roc_params_local$colors) <- groups
    }

    # Position for AUC text
    auc_values$FPR <- 0.2
    auc_values$TPR <- 0.9

    # Determine plot title for ROC. If group_var is provided, include it; otherwise, indicate "All Samples".
    final_title <- if (!is.null(title)) {
      wrap_title(title, width = widthTitle)
    } else {
      if (group_var != "Group") paste("ROC Curves for Each Gene Across", group_var) else "ROC Curves for Each Gene Across All Samples"
    }

    roc_plot_local <- ggplot2::ggplot(roc_df, ggplot2::aes(x = FPR, y = TPR, color = Group, group = Group)) +
      ggplot2::geom_line(size = 1) +
      ggplot2::facet_wrap(~ Gene, scales = "free", ncol = roc_params_local$ncol, nrow = roc_params_local$nrow) +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = final_title,
                    x = "False Positive Rate (1 - Specificity)",
                    y = "True Positive Rate (Sensitivity)") +
      ggplot2::theme(legend.position = "bottom",
                     plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")) +
      ggplot2::geom_abline(linetype = "dashed", color = "gray") +
      ggplot2::scale_color_manual(values = roc_params_local$colors) +
      ggplot2::labs(color = group_var)
    return(roc_plot_local)
  }

  ### Function to generate the heatmap
  makeHeatmap <- function() {
    # Default heatmap parameters
    heatmap_defaults <- list(
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      col = circlize::colorRamp2(c(0.5, 1), c("#FFFFFF", "#21975C")),  # Default AUC scale: white to green
      name = "AUC"
    )

    # Validate user parameters
    allowed_heatmap_params <- c("cluster_rows", "cluster_columns", "col", "name", "row_names_gp", "column_names_gp", "limits")
    invalid_params <- base::setdiff(base::names(heatmap_params), allowed_heatmap_params)
    if (base::length(invalid_params) > 0) {
      stop(paste("Invalid heatmap parameter(s):", base::paste(invalid_params, collapse = ", ")))
    }

    # Merge user-provided heatmap params with defaults
    heatmap_params_local <- utils::modifyList(heatmap_defaults, heatmap_params)

    # Use user-provided limits or default ones
    if (is.null(heatmap_params_local$limits)) {
      limits <- c(0.5, 1)
    } else {
      limits <- heatmap_params_local$limits
    }

    # Apply the colorRamp2 using the (possibly user-defined) color vector from heatmap_params_local$col.
    heatmap_params_local$col <- circlize::colorRamp2(limits, heatmap_params_local$col)
    heatmap_params_local$limits <- NULL  # Remove 'limits'

    # Prepare AUC matrix for heatmap
    auc_matrix <- stats::reshape(auc_values, idvar = "Gene", timevar = "Group", direction = "wide")
    colnames(auc_matrix) <- base::gsub("AUC\\.", "", colnames(auc_matrix))  # Remove "AUC." prefix
    rownames(auc_matrix) <- auc_matrix$Gene
    auc_matrix <- auc_matrix[, -1, drop=F]  # Remove gene column

    # Convert to a numeric matrix
    auc_matrix <- as.matrix(auc_matrix)
    auc_matrix <- matrix(as.numeric(auc_matrix),
                         nrow = base::nrow(auc_matrix),
                         dimnames = list(rownames(auc_matrix), colnames(auc_matrix)))

    # Determine plot title for heatmap
    final_title <- if (!is.null(title)) {
      wrap_title(title, width = widthTitle)
    } else {
      if (group_var != "Group") paste("AUC Heatmap for Each Gene Across", group_var) else "AUC Heatmap for Each Gene Across All Samples"
    }

    heatmap_obj_local <- ComplexHeatmap::Heatmap(
      auc_matrix,
      show_row_names=TRUE,
      col = heatmap_params_local$col,
      name = heatmap_params_local$name,
      cluster_rows = heatmap_params_local$cluster_rows,
      cluster_columns = heatmap_params_local$cluster_columns,
      column_title = final_title,
      column_title_gp = grid::gpar(fontsize = 13, fontface = "bold"),
      cell_fun = function(j, i, x, y, width, height, fill) {
        grid::grid.text(sprintf("%.2f", auc_matrix[i, j]), x, y,
                        gp = grid::gpar(fontsize = 10, col = "black"))
      }
    )
    return(heatmap_obj_local)
  }

  ### Generate the requested plot(s)
  if (plot_type == "roc") {
    roc_plot <- makeROCplot()
    print(roc_plot)
    return(invisible(list(roc_plot = roc_plot, auc_values = auc_values)))
  } else if (plot_type == "heatmap") {
    heatmap_obj <- makeHeatmap()
    ComplexHeatmap::draw(heatmap_obj)
    return(invisible(list(heatmap = heatmap_obj, auc_values = auc_values)))
  } else if (plot_type == "all") {
    # Generate both plots
    roc_plot <- makeROCplot()
    heatmap_obj <- makeHeatmap()

    # Convert the ROC plot to a grob
    roc_grob <- ggplot2::ggplotGrob(roc_plot)
    # Capture the heatmap as a grob
    heatmap_grob <- invisible(grid::grid.grabExpr(ComplexHeatmap::draw(heatmap_obj)))

    # Arrange them side-by-side
    if (!requireNamespace("gridExtra", quietly = TRUE)) {
      stop("Package 'gridExtra' is required for combining plots. Please install it.")
    }
    combined <- gridExtra::grid.arrange(roc_grob, heatmap_grob, ncol = 2,
                                        widths = commomplot_params$widths,
                                        heights = commomplot_params$heights)

    return(invisible(list(combined = combined, roc_plot = roc_plot, heatmap = heatmap_obj, auc_values = auc_values)))
  } else {
    stop("Invalid plot_type. Use 'roc', 'heatmap', or 'all'.")
  }
}
