#' Cohen's d Heatmap Function
#'
#' This function computes Cohen's d for each gene based on gene expression data and sample metadata.
#' For each gene, it compares the expression values between samples where \code{condition_var} equals
#' \code{class} (the positive class) versus the remaining samples. The resulting effect sizes are then
#' visualized as a heatmap.
#'
#' @param data A data frame or matrix containing gene expression data, with genes as rows and samples as columns.
#' @param metadata A data frame containing sample metadata. The first column should contain sample identifiers that match the column names of \code{data}.
#' @param genes A character vector specifying which genes to include. If \code{NULL} (default), all genes in \code{data} are used.
#'   A warning is issued if more than 30 genes are selected.
#' @param condition_var A character string specifying the column name in \code{metadata} representing the condition of interest.
#'   (Mandatory; no default.)
#' @param class A character string specifying the positive class label for the condition.
#'   (Mandatory; no default.)
#' @param group_var An optional character string specifying the column name in \code{metadata} used for grouping samples.
#'   If not provided (\code{NULL}), all samples are treated as a single group.
#' @param title An optional character string specifying a custom title for the heatmap.
#'   If not provided, a default title is generated.
#' @param widthTitle A numeric value used when wrapping the title. Default is \code{16}.
#' @param heatmap_params A list of additional parameters for customizing the heatmap. Possible elements include:
#'   \describe{
#'     \item{\code{cluster_rows}}{Logical; if \code{TRUE} (default), rows are clustered.}
#'     \item{\code{cluster_columns}}{Logical; if \code{TRUE} (default), columns are clustered.}
#'     \item{\code{col}}{A vector of length 2 of colors to be used for the minimum and maximum values of the color scale.
#'       Defaults to \code{c("#FFFFFF", "#21975C")}, but note that the default mapping for Cohen's d is set to a divergent scale.}
#'     \item{\code{limits}}{A numeric vector of length 2 specifying the minimum and maximum values for the color scale.
#'       If not provided, defaults to \code{c(-2, 2)}.}
#'     \item{\code{name}}{A character string for the legend title of the color scale. Default is \code{"Cohen's d"}.}
#'     \item{\code{row_names_gp}}{Optional graphical parameters for row names (passed to \pkg{ComplexHeatmap}).}
#'     \item{\code{column_names_gp}}{Optional graphical parameters for column names (passed to \pkg{ComplexHeatmap}).}
#'   }
#'
#' @return Invisibly returns a list containing:
#'   \describe{
#'     \item{\code{plot}}{The \pkg{ComplexHeatmap} object of the Cohen's d heatmap.}
#'     \item{\code{data}}{A data frame with the calculated Cohen's d values for each gene and group.}
#'   }
#'
#' @details
#' This function computes Cohen's d for each gene by comparing the expression values between samples with
#' \code{condition_var == class} and those that do not. The effect sizes are then visualized as a heatmap using
#' \pkg{ComplexHeatmap}. When \code{group_var} is not provided, all samples are treated as a single group.
#'
#' @examples
#' \dontrun{
#'   result <- CohenDHeatmap(data = myData,
#'                           metadata = myMetadata,
#'                           genes = c("Gene1", "Gene2", "Gene3"),
#'                           condition_var = "Genotype",
#'                           class = "Mutant",
#'                           group_var = "CellType",
#'                           title = "Cohen's d Heatmap",
#'                           heatmap_params = list(limits = c(-2, 2)))
#' }
#'
#' @importFrom grid gpar
#' @importFrom grid grid.text
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom ComplexHeatmap draw
#' @importFrom circlize colorRamp2
#'
#' @export
CohenDHeatmap <- function(data, metadata,
                          genes = NULL,
                          condition_var,
                          class,
                          group_var = NULL,
                          title = NULL,
                          widthTitle = 16,
                          heatmap_params = list()) {

  # If genes are not specified, use all available genes
  if (is.null(genes)) {
    genes <- rownames(data)
  } else {
    genes <- genes[genes %in% rownames(data)]
  }

  if (!all(colnames(data) %in% metadata[, 1]))
    stop("Not all samples in the data are described in the metadata.")

  # Subset and order metadata based on data
  rownames(metadata) <- metadata[, 1]
  colnames(metadata)[1] <- "sample"
  metadata <- metadata[colnames(data), ]

  # Keep only the selected genes
  data <- data[genes, ]

  # Convert count data to log2 scale
  count_data_log2 <- as.data.frame(t(log2(data)))
  count_data_log2 <- cbind(sample = rownames(count_data_log2), count_data_log2)
  rownames(count_data_log2) <- NULL

  # Merge metadata with count data
  data_merge <- merge(metadata, count_data_log2, by = "sample")

  if (length(genes) > 30)
    warning("Too many genes selected. Consider reducing the number.")

  # If group_var is provided, ensure it exists; otherwise, treat all samples as one group.
  if (!is.null(group_var)) {
    if (!group_var %in% colnames(metadata)) {
      stop(paste("Error: The specified group_var", group_var, "is not found in the metadata."))
    }
    groups <- unique(data_merge[[group_var]])
  } else {
    data_merge$Group <- "All"
    groups <- "All"
    group_var <- "Group"
  }

  # Helper function to compute Cohen's d
  compute_cohens_d <- function(x, y) {
    n1 <- length(x)
    n2 <- length(y)
    if(n1 < 2 || n2 < 2) return(NA)
    m1 <- mean(x)
    m2 <- mean(y)
    s1 <- sd(x)
    s2 <- sd(y)
    pooled_sd <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
    if (pooled_sd == 0) return(NA)
    d <- (m1 - m2) / pooled_sd
    d <- abs(d) # dor this application we don't need the direction of the classification, just if it works
    return(d)
  }

  effect_values <- data.frame(Gene = character(), Group = character(), CohensD = numeric())

  # Compute Cohen's d for each gene and group
  for (gene in genes) {
    for (group in groups) {
      subset_data <- data_merge[data_merge[[group_var]] == group, ]
      pos <- subset_data[[gene]][subset_data[[condition_var]] == class]
      neg <- subset_data[[gene]][subset_data[[condition_var]] != class]
      d_val <- compute_cohens_d(pos, neg)
      effect_values <- rbind(effect_values, data.frame(Gene = gene, Group = group, CohensD = d_val))
    }
  }

  if (any(is.na(effect_values$CohensD))) {
    warning("Some Cohen's d values are NA. Check if there are enough samples in each group.")
  }

  # Prepare matrix for heatmap
  effect_matrix <- reshape(effect_values, idvar = "Gene", timevar = "Group", direction = "wide")
  colnames(effect_matrix) <- gsub("CohensD\\.", "", colnames(effect_matrix))
  rownames(effect_matrix) <- effect_matrix$Gene
  effect_matrix <- effect_matrix[, -1, drop = FALSE]
  effect_matrix <- as.matrix(effect_matrix)
  effect_matrix <- matrix(as.numeric(effect_matrix),
                          nrow = nrow(effect_matrix),
                          dimnames = list(rownames(effect_matrix), colnames(effect_matrix)))

  # Determine heatmap title
  final_title <- if (!is.null(title)) {
    wrap_title(title, width = widthTitle)
  } else {
    if (group_var != "Group") paste("Cohen's d Heatmap for Each Gene Across", group_var) else "Cohen's d Heatmap for Each Gene Across All Samples"
  }

  # Set default heatmap parameters if missing
  heatmap_defaults <- list(
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    col = c("#F9F4AE" ,"#B44141"),
    name = "Cohen's d"
  )
  heatmap_params_local <- modifyList(heatmap_defaults, heatmap_params)

  # Use user-provided limits or default ones
  if (is.null(heatmap_params_local$limits)) {
    limits <- c(min(effect_matrix), max(effect_matrix))
  } else {
    limits <- heatmap_params_local$limits
  }

  # Apply the colorRamp2 using the (possibly user-defined) color vector from heatmap_params_local$col.
  heatmap_params_local$col <- circlize::colorRamp2(limits, heatmap_params_local$col)
  heatmap_params_local$limits <- NULL  # Remove 'limits'

  heatmap_obj_local <- ComplexHeatmap::Heatmap(
    effect_matrix,
    show_row_names = TRUE,
    col = heatmap_params_local$col,
    name = heatmap_params_local$name,
    cluster_rows = heatmap_params_local$cluster_rows,
    cluster_columns = heatmap_params_local$cluster_columns,
    column_title = final_title,
    column_title_gp = grid::gpar(fontsize = 13, fontface = "bold"),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid::grid.text(sprintf("%.2f", effect_matrix[i, j]), x, y,
                      gp = grid::gpar(fontsize = 10, col = "black"))
    }
  )

  ComplexHeatmap::draw(heatmap_obj_local)

  return(invisible(list(plot = heatmap_obj_local, data = effect_values)))
}
