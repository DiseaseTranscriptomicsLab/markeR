#' ROC and AUC Plot Function
#'
#' This function computes ROC curves and AUC values for each gene based on gene expression data and sample metadata.
#' It can generate ROC plots, an AUC heatmap / barplot, or both arranged side‐by‐side.
#'
#' @param data A data frame or matrix containing gene expression data, with genes as rows and samples as columns.
#' @param metadata A data frame containing sample metadata. The first column should contain sample identifiers that match the column names of \code{data}.
#' @param genes A character vector specifying which genes to plot. If \code{NULL} (default), all genes in \code{data} are used.
#'   A warning is issued if more than 30 genes are selected.
#' @param condition_var A character string specifying the column name in \code{metadata} representing the condition of interest.
#'   (Mandatory; no default.)
#' @param class A character string or vector specifying the positive class label for the condition.
#'   (Mandatory; no default.)
#' @param group_var An optional character string specifying the column name in \code{metadata} used for grouping samples (e.g., cell types).
#'   If not provided (\code{NULL}), all samples are treated as a single group. Should be a categorical variable.
#' @param plot_type A character string indicating which plot(s) to generate.
#'   Accepted values are \code{"roc"} (only ROC curves), \code{"auc"} (only the AUC heatmap/barplot), or \code{"all"} (both arranged side-by-side).
#'   Default is \code{"roc"}.
#' @param title An optional character string specifying the main title of the plot.
#' @param titlesize A numeric value specifying the size of the title. Default is \code{14}.
#' @param roc_params A list of additional parameters for customizing the ROC plot. Possible elements include:
#'   \describe{
#'     \item{\code{nrow}}{An integer specifying the number of rows in the ROC plot grid. If \code{NULL} (default), it is calculated automatically.}
#'     \item{\code{ncol}}{An integer specifying the number of columns in the ROC plot grid. If \code{NULL} (default), it is calculated automatically.}
#'     \item{\code{colors}}{A named vector of colors for the different groups. If \code{NULL} (default), a default color palette is generated.}
#'   }
#' @param auc_params A list of additional parameters for customizing the AUC heatmap or AUC barplot. Possible elements include:
#'   \describe{
#'     \item{\code{cluster_rows}}{Logical; if \code{TRUE} (default), rows are clustered.}
#'     \item{\code{cluster_columns}}{Logical; if \code{TRUE} (default), columns are clustered.}
#'     \item{\code{colors}}{If `group_var` is used, should be a vector of length 2 of colors to be used for the minimum and maximum values of the color scale. Defaults to \code{c("#FFFFFF", "#21975C")}.
#'     If `group_var` is `NULL`, then should be a single color to fill the barplot. If `NULL`, defaults to \code{"#3B415B"}. If a vector is provided, only the first color will be used.}
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
#' Depending on the value of \code{plot_type}, it produces ROC plots (using \pkg{ggplot2}), an AUC heatmap (using \pkg{ComplexHeatmap})
#' or AUC barplot (if \code{group_var} is `NULL`), or both arranged side-by-side (using \pkg{gridExtra}).
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
#'                           auc_params= list(limits = c(0.4, 1)),
#'                           commomplot_params = list(widths = c(1, 1)))
#' }
#'
#' @importFrom pROC roc
#' @importFrom pROC auc
#' @importFrom scales hue_pal
#' @import ggplot2
#' @importFrom circlize colorRamp2
#' @importFrom utils modifyList
#' @importFrom stats reshape
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom ComplexHeatmap draw
#' @importFrom grid gpar
#' @importFrom grid grid.text
#' @importFrom grid grid.grabExpr
#' @importFrom gridExtra grid.arrange
#'
#' @export
ROCandAUCplot <- function(data, metadata,
                          genes = NULL,
                          condition_var,
                          class,
                          group_var = NULL,
                          plot_type = "roc",  # "roc", "heatmap", or "all"
                          title = NULL,       # Custom plot title (optional)
                          titlesize = 14,
                          roc_params = list(),     # ROC-specific parameters
                          auc_params= list(),
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

  # Convert count data to log2 scale
  count_data_log2 <- as.data.frame(t(log2(data )))
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
      subset_data$Binary_Label <- as.numeric(subset_data[[condition_var]] %in% class)

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

  auc_values$AUC <- as.numeric(auc_values$AUC)

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
    # auc_values$FPR <- 0.2
    # auc_values$TPR <- 0.9

    # Determine plot title for ROC. If group_var is provided, include it; otherwise, indicate "All Samples".
    final_title <-  "ROC curves"


    legend_position <- if (length(unique(roc_df[[group_var]])) > 1 && group_var != "All") "bottom" else "none"

    roc_plot_local <- ggplot2::ggplot(roc_df, ggplot2::aes(x = FPR, y = TPR, color = Group, group = Group)) +
      ggplot2::geom_line(size = 1) +
      ggplot2::facet_wrap(~ Gene, scales = "free", ncol = roc_params_local$ncol, nrow = roc_params_local$nrow) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = final_title,
        x = "False Positive Rate (1 - Specificity)",
        y = "True Positive Rate (Sensitivity)",
        color = group_var
      ) +
      ggplot2::theme(
        legend.position = legend_position,
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
      ) +
      ggplot2::geom_abline(linetype = "dashed", color = "gray") +
      ggplot2::scale_color_manual(values = roc_params_local$colors)

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
    allowed_heatmap_params <- c("cluster_rows", "cluster_columns", "colors", "name", "row_names_gp", "column_names_gp", "limits")
    invalid_params <- base::setdiff(base::names(auc_params), allowed_heatmap_params)
    if (base::length(invalid_params) > 0) {
      stop(paste("Invalid heatmap parameter(s):", base::paste(invalid_params, collapse = ", ")))
    }

    # Merge user-provided heatmap params with defaults
    heatmap_params_local <- utils::modifyList(heatmap_defaults, auc_params)

    # Use user-provided limits or default ones
    if (is.null(heatmap_params_local$limits)) {
      limits <- c(0.5, 1)
    } else {
      limits <- heatmap_params_local$limits
    }

    # Apply the colorRamp2 using the (possibly user-defined) color vector from heatmap_params_local$col.
    heatmap_params_local$col <- circlize::colorRamp2(limits, heatmap_params_local$colors)
    heatmap_params_local$limits <- NULL  # Remove 'limits'
    # keep only columns Gene, Group and AUC
    #auc_values <- auc_values[, c("Gene", "Group", "AUC")]
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
    final_title <-  "AUC values"


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


  makeBarplot <- function() {
    # Sort AUCs in descending order
    auc_sorted <- auc_values[order(-auc_values$AUC), ]

    # Determine title
    final_title <-  "AUC values"


    fillcolor <- ifelse(is.null(auc_params$colors), "#3B415B", auc_params$colors[1])

    barplot <- ggplot2::ggplot(auc_sorted, ggplot2::aes(y = reorder(Gene, AUC), x = AUC)) +
      ggplot2::geom_bar(stat = "identity", fill = fillcolor) +
      #ggplot2::coord_flip()  +
      coord_cartesian(xlim = c(0.5, 1))+
      ggplot2::labs(y = "Gene", x = "AUC", title = final_title) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 14, face = "bold",hjust = 0.5),
        axis.text = ggplot2::element_text(size = 10),
      )

    return(barplot)
  }

  plotAUC <- function() {
    # Check the number of unique conditions
    n_groups <- length(unique(auc_values[["Group"]]))

    if (n_groups == 1) {
      # For a single condition, use barplot
      return(makeBarplot())
    } else {
      # Multiple groups — use heatmap
      return(makeHeatmap())  # assumes makeHeatmap() is in the same environment and uses auc_values
    }
  }




  ### Generate the requested plot(s)
  # create default title specifying the groups for the comparison, with the "class" as positive class
  # condition_var = "Condition",
  # class = "Senescent",

  if (is.null(title)) title <- paste("ROC and AUC for variable", condition_var, "(", paste(class, collapse = ", "), " vs others)")

  if (plot_type == "roc") {
    roc_plot <- makeROCplot()

    roc_plot <- ggpubr::annotate_figure(roc_plot, top = grid::textGrob(title, gp = grid::gpar(cex = 1.3, fontsize = titlesize)))

    print(roc_plot)
    return(invisible(list(roc_plot = roc_plot, auc_values = auc_values)))
  } else if (plot_type == "auc") {

    auc_obj <- plotAUC()

    if (length(unique(auc_values[["Group"]])) > 1){

      # add title to the heatmap
      ComplexHeatmap::draw(auc_obj, column_title = title, column_title_gp = grid::gpar(fontsize = titlesize, fontface = "bold"))

    } else {
      auc_obj <- ggpubr::annotate_figure(auc_obj, top = grid::textGrob(title, gp = grid::gpar(cex = 1.3, fontsize = titlesize)))
      print(auc_obj)
    }

    return(invisible(list(auc_plot = auc_obj, auc_values = auc_values)))

  } else if (plot_type == "all") {
    # Generate both plots
    roc_plot <- makeROCplot()
    auc_obj <- plotAUC()

    # Convert the ROC plot to a grob
    roc_grob <- ggplot2::ggplotGrob(roc_plot)

    if (length(unique(auc_values[["Group"]])) > 1){
      a <- 1
      auc_grob <- invisible(grid::grid.grabExpr(ComplexHeatmap::draw(auc_obj)))
    } else {
      auc_grob <- ggplot2::ggplotGrob(auc_obj)
    }

    # Capture the heatmap as a grob


    # Arrange them side-by-side
    combined <- gridExtra::grid.arrange(roc_grob, auc_grob, ncol = 2,
                                        widths = commomplot_params$widths,
                                        heights = commomplot_params$heights,
                                        top = grid::textGrob(title, gp = grid::gpar(cex = 1.3, fontsize = titlesize)))

    return(invisible(list(combined = combined, roc_plot = roc_plot, auc_plot = auc_obj, auc_values = auc_values)))
  } else {
    stop("Invalid plot_type. Use 'roc', 'auc', or 'all'.")
  }
}
