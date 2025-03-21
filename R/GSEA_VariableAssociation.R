#' GSEA Variable Association
#'
#' This function assesses the association between gene expression (or another molecular score)
#' and metadata variables using differential expression (DE) analysis and Gene Set Enrichment Analysis (GSEA).
#' It generates all possible contrasts for categorical variables and uses linear modeling for continuous variables.
#'
#' @param data A matrix or data frame containing gene expression data, where rows represent genes and columns represent samples.
#' @param metadata A data frame containing sample metadata with at least one column corresponding to the variables of interest.
#' @param cols A character vector specifying the metadata columns (variables) to analyze.
#' @param stat Optional. The statistic to use for ranking genes before GSEA. If `NULL`, it is automatically determined based on the gene set:
#'   - `"B"` for gene sets with **no known direction** (vectors).
#'   - `"t"` for **unidirectional** or **bidirectional** gene sets (data frames).
#'   - If provided, this argument overrides the automatic selection.
#' @param mode A character string specifying the contrast generation method for categorical variables. Options: `"simple"`, `"medium"`, `"extensive"`. Default is `"simple"`.
#' @param gene_set A named list defining the gene sets for GSEA. **(Required)**
#'   - If using **unidirectional** gene sets, provide a list where each element is a vector of gene names representing a signature.
#'   - If using **bidirectional** gene sets, provide a list where each element is a data frame:
#'    - The **first column** should contain gene names.
#'    - The **second column** should indicate the expected direction of enrichment (`1` for upregulated, `-1` for downregulated).
#' @param padj_limit A numeric vector of length 2 specifying the color scale limits for adjusted p-values in the plot. Default: `c(0, 0.1)`.
#' @param low_color The color representing low adjusted p-values in the plot. Default: `"blue"`.
#' @param mid_color The color representing the chosen value for `sig_threshold`. Default: `"white"`.
#' @param high_color The color representing high adjusted p-values in the plot. Default: `"red"`.
#' @param sig_threshold A numeric value specifying the threshold for significance visualization in the plot. Default: `0.05`.
#' @param widthlabels An integer controlling the maximum width of contrast labels before text wrapping. Default: `18`.
#' @param labsize An integer controlling the axis text size in the plot. Default: `10`.
#' @param titlesize An integer specifying the plot title size. Default: `14`.
#'
#' @return A list with two elements:
#'   - `data`: A data frame containing the GSEA results, including normalized enrichment scores (NES), adjusted p-values, and contrasts.
#'   - `plot`: A ggplot2 object visualizing the GSEA results as a lollipop plot.
#'
#' @examples
#' # Example usage with random data
#' data <- matrix(rnorm(1000), ncol = 10)
#' metadata <- data.frame(group = rep(c("A", "B"), each = 5))
#' gene_set <- list(SampleSet = c("Gene1", "Gene2", "Gene3"))
#' results <- GSEA_VariableAssociation(data, metadata, cols = "group", gene_set = gene_set)
#'
#' # View results
#' print(results$data)
#' print(results$plot)
#'
#' @export
GSEA_VariableAssociation <- function(data, metadata, cols, stat=NULL, mode=c("simple","medium","extensive"), gene_set, padj_limit = c(0, 0.1), low_color = "blue", mid_color = "white", high_color = "red", sig_threshold = 0.05, widthlabels=18, labsize=10, titlesize=14) {
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
      # Use a linear model expression for continuous variables
      DEGs_var <- calculateDE(data = data, metadata = metadata, lmexpression = paste0("~", var))
      cont_vec <- c(cont_vec,c(paste0("intercept_",var),var))
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

  plot_list <- list()


  # Ensure contrast ordering
  combined_results$Contrast <- sapply(combined_results$Contrast, function(x) wrap_title(x, widthlabels))
  combined_results$Contrast <- factor(combined_results$Contrast, levels = combined_results$Contrast[order(combined_results$NES)])


  plot <- ggplot2::ggplot(combined_results, ggplot2::aes(x = NES, y = Contrast,fill = padj)) +
    ggplot2::geom_segment(ggplot2::aes(yend = Contrast, xend = 0), size = .5) +
    ggplot2::geom_point( shape = 21, stroke = 1.2, color="black", size=4) +
    ggplot2::scale_fill_gradient2(low = low_color, mid = mid_color, high = high_color, midpoint = sig_threshold, limits = padj_limit, na.value=high_color) +
    ggplot2::labs(title = unique(combined_results$pathway), x = "Normalized Enrichment Score (NES)", y = "Contrast", color = "Adj. p-value", fill = "Adj. p-value") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size=titlesize),
      legend.position = "right",
      axis.text =  ggplot2::element_text(size=labsize)
    )



  return(list(plot=plot,
              data=combined_results_toreturn))  # Return list if not in grid mode

}


#' Generate All Possible Contrasts Between Groups
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
#' @return A character vector of contrast expressions (without names).
#' @examples
#' levels <- c("A", "B", "C", "D")
#' generate_all_contrasts(levels, mode = "simple")    # Pairwise only
#' generate_all_contrasts(levels, mode = "medium")    # Pairwise + mean comparisons
#' generate_all_contrasts(levels, mode = "extensive") # All balanced contrasts
#' @export
generate_all_contrasts <- function(levels, mode = "extensive") {
  levels <- unique(levels)  # Ensure unique levels
  n <- length(levels)       # Total number of levels

  if (!mode %in% c("simple", "medium", "extensive")) {
    stop("Invalid mode. Choose 'simple', 'medium', or 'extensive'.")
  }

  contrasts <- c()  # Store contrasts

  # 1. Pairwise comparisons (included in all modes, both A - B and B - A)
  pairwise_contrasts <- combn(levels, 2, function(x) c(paste(x[1], "-", x[2]), paste(x[2], "-", x[1])), simplify = TRUE)

  if (mode == "simple") {
    return(unname(as.vector(pairwise_contrasts)))
  }

  # 2. Comparisons against the mean of other groups (included in medium & extensive)
  mean_contrasts <- unlist(lapply(levels, function(x) {
    others <- setdiff(levels, x)
    if (length(others) > 1) {
      return(c(
        paste(x, "- (", paste(others, collapse = " + "), ")/", length(others)),
        paste("( ", paste(others, collapse = " + "), ")/", length(others), "-", x)
      ))
    } else {
      return(c(paste(x, "-", others), paste(others, "-", x)))  # Simple case: keep both orders
    }
  }))

  if (mode == "medium") {
    return(unname(c(as.vector(pairwise_contrasts), mean_contrasts)))
  }

  # 3. Groupwise comparisons (only included in extensive mode)
  group_contrasts <- c()
  for (i in 1:(n-1)) {
    left_groups <- combn(levels, i, simplify = FALSE)  # Subsets for the first group
    for (left in left_groups) {
      right <- setdiff(levels, left)  # Remaining elements for the second group
      if (length(right) > 1 && length(left) > 1) {
        left_weighted <- paste0("(", paste(left, collapse = " + "), ")/", length(left))
        right_weighted <- paste0("(", paste(right, collapse = " + "), ")/", length(right))
        contrast1 <- paste(left_weighted, "-", right_weighted)
        contrast2 <- paste(right_weighted, "-", left_weighted)
        group_contrasts <- c(group_contrasts, contrast1, contrast2)  # Keep both orders
      }
    }
  }

  # Combine all contrasts and return as an unnamed vector
  contrasts <- unique(c(as.vector(pairwise_contrasts), mean_contrasts, group_contrasts))  # Remove duplicates
  return(unname(contrasts))
}
