#' Generate Heatmaps for Cohen\'s d Effect Sizes using ggplot2
#'
#' This function computes Cohen\'s d effect sizes and corresponding p-values for multiple gene signatures and produces individual heatmaps. Each heatmap displays cell text showing the Cohen\'s d value along with its p-value. The heatmaps are then arranged in a grid layout.
#'
#' @param cohenlist A named list where each element corresponds to a gene signature. Output of `CohenD_allConditions`. Each signature element is a list with three components:
#' \describe{
#'   \item{CohenD}{A data frame where rows are methods and columns are group contrasts (formatted as \"Group1:Group2\"),
#'   containing the computed Cohen\'s d effect sizes.}
#'   \item{PValue}{A data frame with the same structure as \code{CohenD} containing the corresponding p-values.}
#'   \item{padj}{A data frame with the same structure as \code{PValue} containing the corresponding p-values corrected using the BH method, for all signatures and contrasts,
#'   and by method.}
#' }
#' @param nrow Optional. An integer specifying the number of rows in the heatmap grid. If \code{NULL}, the number of rows
#'   is computed automatically.
#' @param ncol Optional. An integer specifying the number of columns in the heatmap grid. If \code{NULL}, the number of columns
#'   is computed automatically.
#' @param limits Optional. A numeric vector of length 2 specifying the color scale limits (e.g., \code{c(min, max)}). If \code{NULL},
#'   the limits are determined from the data.
#' @param widthTitle An integer specifying the width used for wrapping gene set signature names in the heatmap titles. Default is 22.
#' @param titlesize An integer specifying the text size for each of the heatmap titles. Default is 12.
#' @param ColorValues A character vector specifying the colors for the gradient fill in the heatmaps. Default is \code{c("#F9F4AE", "#B44141")}.
#' @param title Title for the grid of plots.
#' @return A list with two elements:
#' \describe{
#'   \item{plt}{A combined heatmap arranged in a grid using \code{ggpubr::ggarrange}.}
#'   \item{data}{A list containing the Cohen\'s d effect sizes and p-values for each gene signature, as computed by \code{CohenD_allConditions}.}
#' }
#'
#' @details
#' The function first calculates Cohen\'s d effect sizes and corresponding p-values for each gene signature using \code{CohenD_allConditions} (assumed to be defined elsewhere in the package). The resulting matrices are converted to a long format so that each cell in the heatmap can display the Cohen\'s d value and its associated p-value (formatted as \code{Cohen\'s d (p-value)}).
#'
#' The heatmaps are then adjusted to display axis text and ticks only for the left-most column and bottom row, and combined into a grid layout. If neither \code{nrow} nor \code{ncol} are specified, the layout is automatically determined to best approximate a square grid.
#'
#' @examples
#' \dontrun{
#'   # Assuming gene_data, sample_metadata, and gene_sets are defined:
#'   result <- Heatmap_CohenD_ggplot(
#'     data = gene_data,
#'     metadata = sample_metadata,
#'     gene_sets = gene_sets,
#'     variable = "Condition",
#'     nrow = 2,
#'     ncol = 3,
#'     limits = c(-1, 1),
#'     widthTitle = 30,
#'     titlesize = 14,
#'     ColorValues = c("#F9F4AE", "#B44141")
#'   )
#'   print(result$plt)
#' }
#'
#' @seealso \code{\link{CohenD_allConditions}}, \code{\link{wrap_title}}
#'
#' @importFrom ggplot2 ggplot geom_tile geom_text labs scale_fill_gradientn theme_minimal element_text element_blank element_line margin
#' @importFrom ggpubr ggarrange
#'
#' @export
Heatmap_CohenD <- function(cohenlist, nrow = NULL, ncol = NULL, limits = NULL, widthTitle = 22, titlesize = 12, ColorValues = NULL,title=NULL ) {

  heatmaps <- list()

  if (is.null(ColorValues)) ColorValues <- c("#F9F4AE", "#B44141")

  for (signature in names(cohenlist)) {
    cohen_d_mat <- t(as.matrix(cohenlist[[signature]]$CohenD))
    p_value_mat <- t(as.matrix(cohenlist[[signature]]$padj))

    # Convert to long format manually
    long_data <- data.frame(
      Var1 = rep(rownames(cohen_d_mat), times = ncol(cohen_d_mat)),
      Var2 = rep(colnames(cohen_d_mat), each = nrow(cohen_d_mat)),
      CohenD = abs(as.vector(cohen_d_mat)), # absolute value
      PValue = as.vector(p_value_mat),
      stringsAsFactors = FALSE
    )

    # Generate text labels
    long_data$label <- paste0(sprintf("%.2f", long_data$CohenD), "\n(", format.pval(long_data$PValue, digits = 1), ")")

    # Wrap the signature title using an internal helper function
    signature_title <- wrap_title(signature, widthTitle)

    # Create heatmap using ggplot2
    p <- ggplot2::ggplot(long_data, ggplot2::aes(x = Var2, y = Var1, fill = CohenD)) +
      ggplot2::geom_tile() +
      ggplot2::geom_text(aes(label = label), color = "black", size = 3) +
      ggplot2::scale_fill_gradientn(colors = ColorValues, limits = limits) +
      ggplot2::labs(title = signature_title, x = NULL, y = NULL, fill = "|Cohen\'s D|") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     plot.title = ggplot2::element_text(hjust = 0.5, size = titlesize) )

    heatmaps[[signature]] <- p
  }

  # Determine grid layout if not provided
  num_signatures <- length(heatmaps)
  if (is.null(nrow) & is.null(ncol)) {
    ncol <- ceiling(sqrt(num_signatures))
    nrow <- ceiling(num_signatures / ncol)
  } else if (is.null(nrow)) {
    nrow <- ceiling(num_signatures / ncol)
  } else if (is.null(ncol)) {
    ncol <- ceiling(num_signatures / nrow)
  }

  # Adjust axis labels based on grid position
  for (i in seq_along(heatmaps)) {
    row_idx <- ceiling(i / ncol)  # Current row number
    col_idx <- (i - 1) %% ncol + 1  # Current column number

    p <- heatmaps[[i]] +
      theme(
        axis.text.y = if (col_idx == 1) ggplot2::element_text() else ggplot2::element_blank(),
        axis.ticks.y = if (col_idx == 1) ggplot2::element_line() else ggplot2::element_blank(),
        plot.margin = ggplot2::margin(4, 0, 0, 0)  # Adjust plot margins
      )

    heatmaps[[i]] <- p
  }

  # Compute dynamic column widths: first column is wider
  widths <- c(1.5, rep(1, ncol - 1))

  # Combine heatmaps into a single grid plot using ggpubr
  plt <- ggpubr::ggarrange(
    plotlist = heatmaps,
    ncol = ncol,
    nrow = nrow,
    common.legend = TRUE,
    legend = "right",
    align = "h",
    widths = widths
  )

  plt <- ggpubr::annotate_figure(plt, top = grid::textGrob(title, gp = grid::gpar(cex = 1.3, fontsize = titlesize+2)))


  return(list(plt = plt, data = cohenlist))
}





#' Compute Cohen\'s d for All Gene Signatures Across Conditions
#'
#' Computes Cohen\'s d effect sizes and corresponding p-values for all gene signatures using scores calculated
#' by various methods. The function first computes gene signature scores using \code{CalculateScores} with the "all"
#' option, flattens the results, and then computes pairwise comparisons for a specified grouping variable.
#'
#' @param data A data frame of gene expression data, with genes as rows and samples as columns.
#' @param metadata A data frame containing sample metadata. The first column should contain sample identifiers matching
#'   the column names of \code{data}.
#' @param gene_sets A named list of gene sets. For unidirectional gene sets, each element is a vector of gene names;
#'   for bidirectional gene sets, each element is a data frame where the first column contains gene names and the second
#'   column indicates the expected direction (1 for upregulated, -1 for downregulated).
#' @param variable A string specifying the grouping variable in \code{metadata} used to compare scores between conditions.
#' @param mode A string specifying the level of detail for contrasts.
#' Options are:
#' - `"simple"`: Pairwise comparisons (e.g., A - B).
#' - `"medium"`: Pairwise comparisons plus comparisons against the mean of other groups.
#' - `"extensive"`: All possible groupwise contrasts, ensuring balance in the number of terms on each side.
#'
#' @return A named list where each element corresponds to a gene signature. Each signature element is a list with three components:
#' \describe{
#'   \item{CohenD}{A data frame where rows are methods and columns are group contrasts (formatted as \"Group1:Group2\"),
#'   containing the computed Cohen\'s d effect sizes.}
#'   \item{PValue}{A data frame with the same structure as \code{CohenD} containing the corresponding p-values.}
#'   \item{padj}{A data frame with the same structure as \code{PValue} containing the corresponding p-values corrected using the BH method, for all signatures and contrasts,
#'   and by method.}
#' }
#'
#' @examples
#' \dontrun{
#'   # Assume gene_data is your gene expression data frame, sample_metadata is your metadata, and
#'   # gene_sets is a named list of gene sets.
#'   results <- CohenD_allConditions(data = gene_data, metadata = sample_metadata,
#'                                    gene_sets = gene_sets, variable = \"Condition\")
#'   # Access Cohen\'s d for a specific signature:
#'   results$Signature_A$CohenD
#' }
#'
#' @export
CohenD_allConditions <- function(data, metadata, gene_sets, variable, mode = c("simple","medium","extensive")) {

  # Step 1: Check if variable exists in metadata
  if (!variable %in% colnames(metadata)) {
    stop(paste("Error: Variable", variable, "not found in metadata."))
  }

  # Step 2: Compute scores for all methods and signatures
  listScores <- CalculateScores(data = data, metadata = metadata, gene_sets = gene_sets, method = "all")

  # Step 3: Flatten into a data frame with method & signature columns
  dfScores <- flatten_results(listScores)

  # Step 4: Initialize result list
  result_list <- list()

  # Step 5: Loop through each signature
  for (signature in unique(dfScores$signature)) {

    # Filter for the specific signature
    df_subset <- dfScores[dfScores$signature == signature, ]

    # Initialize storage for this signature
    cohen_d_results <- list()
    p_value_results <- list()

    # Step 6: Loop through each method
    for (method in unique(df_subset$method)) {

      # Subset data for this method
      df_method <- df_subset[df_subset$method == method, ]

      # Compute Cohen\'s d and p-values
      cohen_results <- compute_cohen_d(df_method, variable, quantitative_var = "score", mode=mode)

      # Convert to named vectors (column names = comparisons)
      # cohen_d_results[[method]] <- setNames(cohen_results$CohenD, paste0(cohen_results$Group1, " - ", cohen_results$Group2))
      # p_value_results[[method]] <- setNames(cohen_results$PValue, paste0(cohen_results$Group1, " - ", cohen_results$Group2))
      cohen_d_results[[method]] <- setNames(cohen_results$CohenD, cohen_results$contrast)
      p_value_results[[method]] <- setNames(cohen_results$PValue, cohen_results$contrast)
    }

    # Convert lists to data frames
    cohen_d_df <- as.data.frame(do.call(rbind, cohen_d_results))
    p_value_df <- as.data.frame(do.call(rbind, p_value_results))

    # Store results in the named list
    result_list[[signature]] <- list(CohenD = cohen_d_df, PValue = p_value_df)
  }

  # Step 7: Correct for multiple testing across all signatures and contrasts (but NOT methods)
  # Initialize storage for adjusted p-values
  all_pvalues <- list()

  # Step 1: Extract all p-values grouped by method
  for (signature in names(result_list)) {
    for (method in rownames(result_list[[signature]]$PValue)) {
      if (!method %in% names(all_pvalues)) {
        all_pvalues[[method]] <- c()
      }
      all_pvalues[[method]] <- c(all_pvalues[[method]], as.vector(result_list[[signature]]$PValue[method, ]))
    }
  }

  # Step 2: Apply BH correction within each method
  all_padj <- lapply(all_pvalues, function(pvals) p.adjust(pvals, method = "BH"))

  # Step 3: Store corrected p-values back into result_list
  index_tracker <- list()  # Track index position for each method
  for (signature in names(result_list)) {
    padj_matrix <- matrix(nrow = nrow(result_list[[signature]]$PValue),
                          ncol = ncol(result_list[[signature]]$PValue),
                          dimnames = dimnames(result_list[[signature]]$PValue))

    for (method in rownames(result_list[[signature]]$PValue)) {
      # Initialize index tracker for the method
      if (is.null(index_tracker[[method]])) index_tracker[[method]] <- 1

      # Assign adjusted p-values in order
      num_vals <- length(result_list[[signature]]$PValue[method, ])
      padj_matrix[method, ] <- all_padj[[method]][index_tracker[[method]]:(index_tracker[[method]] + num_vals - 1)]

      # Update index tracker
      index_tracker[[method]] <- index_tracker[[method]] + num_vals
    }

    # Store in result_list under padj
    result_list[[signature]]$padj <- as.data.frame(padj_matrix)
  }

  # Step 8: Return the list
  return(result_list)
}



#' Compute Cohen\'s d Effect Size
#'
#' Computes the absolute Cohen\'s d effect size between two numeric vectors. This function returns
#' the absolute value of the difference in means divided by the pooled standard deviation.
#'
#' @param x A numeric vector representing the values for group 1.
#' @param y A numeric vector representing the values for group 2.
#'
#' @return A numeric value representing Cohen\'s d. Returns NA if either group has fewer than two observations
#'   or if the pooled standard deviation is zero.
#'
#' @keywords internal
cohen_d <- function(x, y) {
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
  #d <- abs(d) # dor this application we don't need the direction of the classification, just if it works
  return(d)
}


#' Compute Pairwise Cohen\'s d and P-Values
#'
#' Computes Cohen\'s d effect sizes and corresponding p-values for all pairwise comparisons of a grouping variable
#' in a data frame.
#'
#' @param dfScore A data frame containing at least one numeric column and a grouping variable. Output from flatten_results.
#' @param variable A string specifying the name of the categorical grouping column in \code{dfScore}.
#' @param quantitative_var A string specifying the name of the numeric column (default is \code{"score"}).
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{Group1}{The first group in the pair.}
#'   \item{Group2}{The second group in the pair.}
#'   \item{CohenD}{The computed Cohen\'s d effect size for the comparison.}
#'   \item{PValue}{The p-value from a t-test comparing the two groups.}
#' }
#'
#' @keywords internal
compute_cohen_d <- function(dfScore, variable, quantitative_var="score", mode = c("simple","medium","extensive")) {

  # Get unique group values
  unique_groups <- unique(dfScore[[variable]])

  # Ensure there are at least two groups
  if (length(unique_groups) < 2) {
    stop("The grouping column must have at least two unique values.")
  }


  # Store results
  results <- data.frame(Group1 = character(), Group2 = character(),
                        CohenD = numeric(), PValue = numeric(), stringsAsFactors = FALSE)

  # Compute Cohen\'s d and p-value for all unique pairs
  #combs <- combn(unique_groups, 2, simplify = FALSE)
  combs <- generate_all_contrasts(unique(dfScore[[variable]]), mode = mode)
  combs <- remove_division(combs)

  for (pair in combs) {

    dfScore_subset <- create_contrast_column(dfScore, variable, pair) # subsets metadata and adds new column "cohentest" with the two parts of the contrast
    group1 <- levels(dfScore_subset$cohentest)[1]
    group2 <- levels(dfScore_subset$cohentest)[2]

    # group1 <- unique(dfScore_subset$cohentest)[1]
    # group2 <- unique(dfScore_subset$cohentest)[2]

    x <- dfScore_subset[dfScore_subset[["cohentest"]] == group1, quantitative_var, drop = TRUE]
    y <- dfScore_subset[dfScore_subset[["cohentest"]] == group2, quantitative_var, drop = TRUE]

    d <- cohen_d(x, y)
    set.seed("03042025")
    p_val <- t.test(x, y, var.equal = TRUE)$p.value

    results <- rbind(results, data.frame(Group1 = group1, Group2 = group2,
                                         CohenD = d, PValue = p_val, contrast = pair, stringsAsFactors = FALSE))
  }

  return(results)
}


#' Flatten a Nested List of Results into a Data Frame
#'
#' Converts a nested list (where the first level is a method, the second level is a gene signature,
#' and the third level is a data frame) into a single data frame. Additional columns for method and signature
#' are added to the data frame.
#'
#' @param nested_list A nested list with structure: \code{list(method = list(signature = data.frame(...)))}.
#'
#' @return A data frame combining all the nested data frames, with added columns \code{method} and \code{signature}.
#'
#'
#' @keywords internal
flatten_results <- function(nested_list) {
  # Create an empty list to store results
  result_list <- list()

  # Iterate over methods (first level)
  for (method in names(nested_list)) {
    # Iterate over gene signatures (second level)
    for (signature in names(nested_list[[method]])) {
      # Extract the data frame (third level)
      df <- nested_list[[method]][[signature]]

      # Add columns for method and signature
      df$method <- method
      df$signature <- signature

      # Store the modified data frame in the list
      result_list[[paste(method, signature, sep = "_")]] <- df
    }
  }

  # Combine all data frames into one
  final_df <- do.call(rbind, result_list)

  # Reset row names
  rownames(final_df) <- NULL

  return(final_df)
}



