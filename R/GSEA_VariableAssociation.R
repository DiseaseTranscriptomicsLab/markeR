#' GSEA Variable Association
#'
#' This function assesses the association between gene expression (or another molecular score)
#' and metadata variables using differential expression (DE) analysis and Gene Set Enrichment Analysis (GSEA).
#' It generates all possible contrasts for categorical variables and uses linear modeling for continuous variables.
#'
#' @param data A matrix or data frame containing gene expression data, where rows represent genes and columns represent samples.
#' @param metadata A data frame containing sample metadata with at least one column corresponding to the variables of interest.
#' @param cols A character vector specifying the metadata columns (variables) to analyse.
#' @param stat Optional. The statistic to use for ranking genes before GSEA. If `NULL`, it is automatically determined based on the gene set:
#'   - `"B"` for gene sets with **no known direction** (vectors).
#'   - `"t"` for **unidirectional** or **bidirectional** gene sets (data frames).
#'   - If provided, this argument overrides the automatic selection.
#' @param modeA string specifying the level of detail for contrasts. Options are:
#' - `"simple"`: Performs the minimal number of pairwise comparisons between individual group levels (e.g., A - B, A - C). Default.
#' - `"medium"`: Includes comparisons between one group and the union of all other groups (e.g., A - (B + C + D)), enabling broader contrasts beyond simple pairs.
#' - `"extensive"`: Allows for all possible algebraic combinations of group levels (e.g., (A + B) - (C + D)), supporting flexible and complex contrast definitions.
#' @param gene_set A named list defining the gene sets for GSEA. **(Required)**
#'   - If using **unidirectional** gene sets, provide a list where each element is a vector of gene names representing a signature.
#'   - If using **bidirectional** gene sets, provide a list where each element is a data frame:
#'    - The **first column** should contain gene names.
#'    - The **second column** should indicate the expected direction of enrichment (`1` for upregulated, `-1` for downregulated).
#' @param signif_color A string specifying the color for the low end of the adjusted p-value gradient until the value chosen for significance (\code{sig_threshold}). Default is `"red"`.
#' @param nonsignif_color A string specifying the color for the middle of the adjusted p-value gradient. Default is `"white"`. Lower limit correspond to the value of \code{sig_threshold}.
#' @param sig_threshold A numeric value specifying the threshold for significance visualization in the plot. Default: `0.05`.
#' @param saturation_value A numeric value specifying the lower limit of the adjusted p-value gradient, below which the color will correspond to \code{signif_color}. Default is the results' minimum, unless that
#' value is above the sig_threshold; in that case, it is 0.001.
#' @param widthlabels An integer controlling the maximum width of contrast labels before text wrapping. Default: `18`.
#' @param labsize An integer controlling the axis text size in the plot. Default: `10`.
#' @param titlesize An integer specifying the plot title size. Default: `14`.
#' @param pointSize Numeric. The size of points in the lollipop plot (default is 5).
#'
#' @return A list with two elements:
#'   - `data`: A data frame containing the GSEA results, including normalized enrichment scores (NES), adjusted p-values, and contrasts.
#'   - `plot`: A ggplot2 object visualizing the GSEA results as a lollipop plot.
#'
#' @examples
#' # Example usage with random data
#' set.seed(42)  # For reproducibility
#'
#' # Create random gene expression data
#' data <- matrix(rnorm(1000), ncol = 10)
#'
#' # Assign gene identifiers as row names (e.g., Gene1, Gene2, ...)
#' rownames(data) <- paste0("Gene", 1:nrow(data))
#'
#' # Create metadata (e.g., group variable)
#' metadata <- data.frame(group = rep(c("A", "B"), each = 5))
#'
#' # Define a gene set
#' gene_set <- list(SampleSet = c("Gene1", "Gene2", "Gene3"))
#'
#' # Call the GSEA_VariableAssociation function
#' results <- GSEA_VariableAssociation(data, metadata, cols = "group", gene_set = gene_set)
#'
#' # View results
#' print(results$data)
#' print(results$plot)
#'
#' @export
GSEA_VariableAssociation <- function(data, metadata, cols, stat=NULL, mode=c("simple","medium","extensive"), gene_set,nonsignif_color = "grey", signif_color = "red", saturation_value=NULL,sig_threshold = 0.05, widthlabels=18, labsize=10, titlesize=14, pointSize=5) {
  mode <- match.arg(mode)
  metadata <- metadata[, cols %in% colnames(metadata), drop = FALSE]

  # Identify variable types
  variable_types <- identify_variable_type(metadata, cols)

  # Store results and plots
  all_results <- list()  # Store GSEA results
  all_plots <- list()    # Store plots

  cont_vec <- c()

  for (var in cols) {

    # Check if the variable exists in metadata
    if (!(var %in% colnames(metadata))) {
      warning(paste0("Variable ", var, " not in metadata"))
      next  # Skip to the next iteration if the variable is not in metadata
    }

    # Handle based on variable type
    if (variable_types[var] == "Numeric") {

      # Use a model matrix for continuous variables
      design <- model.matrix(as.formula(paste("~1+", var)), data = metadata)

      DEGs_var <- calculateDE(data = data, metadata = metadata, modelmat =  design, contrasts = c(var))
      #cont_vec <- c(cont_vec,c(paste0("intercept_",var),var))
      cont_vec <- c(cont_vec, var)

    } else {
      # For categorical variables, generate contrasts
      uniquevalues_var <- unique(metadata[, var])
      uniquevalues_var <- gsub(" ", "", uniquevalues_var)

      contrasts <- generate_all_contrasts(uniquevalues_var, mode = mode)

      # Calculate differential expression results for each contrast
      DEGs_var <- calculateDE(data = data, metadata = metadata, variables = var, contrasts = contrasts)
      cont_vec <- c(cont_vec,contrasts)
    }
    set.seed("20032025")
    # Perform GSEA for each case
    GSEA_results <- runGSEA(DEGs_var, gene_sets = gene_set, stat=stat)

    # Collect GSEA results in a list
    all_results[[var]] <- GSEA_results
  }

  # Combine results into a data frame
  combined_results <- do.call(rbind, lapply(all_results, function(x) do.call(rbind, x)))
  combined_results$Contrast <- cont_vec

  # correct adjusted p value to correct for multiple testing for the contrasts?
  combined_results$padj <- p.adjust(combined_results$padj, method = "BH")


  combined_results_toreturn <- combined_results


  # Ensure contrast ordering
  combined_results$Contrast <- sapply(combined_results$Contrast, function(x) wrap_title(x, widthlabels))
  combined_results$Contrast <- factor(combined_results$Contrast, levels = combined_results$Contrast[order(combined_results$NES)])

#
#   plot <- ggplot2::ggplot(combined_results, ggplot2::aes(x = NES, y = Contrast,fill = padj)) +
#     ggplot2::geom_segment(ggplot2::aes(yend = Contrast, xend = 0), size = .5) +
#     ggplot2::geom_point( shape = 21, stroke = 1.2, color="black", size=4) +
#     ggplot2::scale_fill_gradient2(low = low_color, mid = mid_color, high = high_color, midpoint = sig_threshold, limits = padj_limit, na.value=high_color) +
#     ggplot2::labs(title = unique(combined_results$pathway), x = "Normalized Enrichment Score (NES)", y = "Contrast", color = "Adj. p-value", fill = "Adj. p-value") +
#     ggplot2::theme_minimal() +
#     ggplot2::theme(
#       plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size=titlesize),
#       legend.position = "right",
#       axis.text =  ggplot2::element_text(size=labsize)
#     )


  if(is.null(saturation_value)){
    if (min(combined_results$padj)>sig_threshold){
      limit_pval <- 0.001
    } else{
      limit_pval <- min(combined_results$padj)
    }

  } else {
    limit_pval <- saturation_value
  }


  plot <- ggplot2::ggplot(combined_results, ggplot2::aes(x = NES, y = Contrast, fill = -log10(padj))) +
    ggplot2::geom_segment(ggplot2::aes(
      yend = Contrast,
      xend = 0,
      linetype = ifelse(stat_used == "B" & NES < 0, "dashed", "solid"),
      color = ifelse(stat_used == "B" & NES < 0, "grey", "black")
    ), size = .5) +
    ggplot2::geom_point(ggplot2::aes(
      stroke = 1.2,
      color = ifelse(stat_used == "B" & NES < 0, "grey", "black")
    ), shape = 21, size = pointSize) +
    #ggplot2::scale_fill_gradient2(low = signif_color, mid = nonsignif_color, high = nonsignif_color, midpoint = sig_threshold, limits = padj_limit, na.value = high_color) +
    ggplot2::scale_fill_gradient2(low = nonsignif_color,
                                  mid = nonsignif_color,
                                  high = signif_color,
                                  midpoint = -log10(sig_threshold),
                                  limits=c(0,-log10(limit_pval)),
                                  na.value = signif_color)+
  ggplot2::scale_linetype_identity() +
    ggplot2::scale_color_identity() +
    ggplot2::labs(
      title = unique(combined_results$pathway),
      subtitle = ifelse(any(combined_results$stat_used == "B"), "Altered", "Enriched/Depleted"),
      x = "Normalized Enrichment Score (NES)",
      y = "Contrast",
      color = "-log10(Adj. p-value)",
      fill = "-log10(Adj. p-value)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = titlesize),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, face = "italic", size = titlesize - 2),
      legend.position = "right",
      axis.text = ggplot2::element_text(size = labsize)
    )



  return(list(plot=plot,
              data=combined_results_toreturn))  # Return list if not in grid mode

}
#' Generate All Possible Unique Contrasts Between Groups
#'
#' This function creates statistical contrasts between levels of a categorical variable.
#' Users can choose the level of complexity:
#' - `"simple"`: Pairwise comparisons (e.g., A - B).
#' - `"medium"`: Pairwise comparisons plus comparisons against the mean of other groups.
#' - `"extensive"`: All possible groupwise contrasts, ensuring balance in the number of terms on each side.
#'
#' @param levels A character vector of unique group levels.
#' @param mode A string specifying the level of detail for contrasts.
#' Options are `"simple"` (pairwise only), `"medium"` (pairwise + vs. mean of others), or `"extensive"` (all possible balanced groupwise contrasts). Default is `"extensive"`.
#'
#' @return A character vector of unique contrast expressions.
#' @examples
#' \dontrun{
#' levels <- c("A", "B", "C", "D")
#' generate_all_contrasts(levels, mode = "simple")    # Pairwise only
#' generate_all_contrasts(levels, mode = "medium")    # Pairwise + mean comparisons
#' generate_all_contrasts(levels, mode = "extensive") # All balanced contrasts
#' }
#' @keywords internal
generate_all_contrasts <- function(levels, mode = "simple") {
  levels <- unique(levels)  # Ensure unique levels
  n <- length(levels)       # Total number of levels

  if (!mode %in% c("simple", "medium", "extensive")) {
    stop("Invalid mode. Choose 'simple', 'medium', or 'extensive'.")
  }

  contrasts <- c()  # Store contrasts

  # Helper function to enforce a consistent order
  normalize_contrast <- function(left, right) {
    left_sorted <- sort(left)
    right_sorted <- sort(right)

    left_expr <- if (length(left_sorted) > 1) paste0("(", paste(left_sorted, collapse = " + "), ")/", length(left_sorted)) else left_sorted
    right_expr <- if (length(right_sorted) > 1) paste0("(", paste(right_sorted, collapse = " + "), ")/", length(right_sorted)) else right_sorted

    if (paste(left_sorted, collapse = " ") < paste(right_sorted, collapse = " ")) {
      return(paste(left_expr, "-", right_expr))
    } else {
      return(paste(right_expr, "-", left_expr))
    }
  }

  # 1. Pairwise comparisons (A - B, but not B - A)
  pairwise_contrasts <- combn(levels, 2, function(x) normalize_contrast(x[1], x[2]), simplify = TRUE)

  if (mode == "simple") {
    return(unname(unique(pairwise_contrasts)))
  }

  # 2. Comparisons against the mean of other groups (medium & extensive)
  mean_contrasts <- unlist(lapply(levels, function(x) {
    others <- setdiff(levels, x)
    normalize_contrast(x, others)
  }))

  if (mode == "medium") {
    return(unname(unique(c(pairwise_contrasts, mean_contrasts))))
  }

  # 3. Groupwise comparisons (extensive mode)
  group_contrasts <- c()
  for (i in 1:(n-1)) {
    left_groups <- combn(levels, i, simplify = FALSE)  # Subsets for the first group
    for (left in left_groups) {
      right <- setdiff(levels, left)  # Remaining elements for the second group
      if (length(right) > 1 && length(left) > 1) {
        group_contrasts <- c(group_contrasts, normalize_contrast(left, right))
      }
    }
  }

  # Combine all contrasts and return unique values
  contrasts <- unique(c(pairwise_contrasts, mean_contrasts, group_contrasts))
  return(unname(contrasts))
}
