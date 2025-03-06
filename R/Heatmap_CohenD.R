#' Generate Heatmaps for Cohen's d Effect Sizes
#'
#' Creates heatmaps for Cohen's d effect sizes (and corresponding p-values, shown in cell text) for each gene signature,
#' and combines them according to the specified orientation.
#'
#' @param data A data frame of gene expression data with genes as rows and samples as columns.
#'   The row names should contain gene names and the column names sample identifiers.
#' @param metadata A data frame of sample metadata. The first column must contain sample identifiers matching \code{data}.
#' @param gene_sets A named list of gene sets. For unidirectional gene sets, each element is a vector of gene names;
#'   for bidirectional gene sets, each element is a data frame with the first column as gene names and the second column as the expected direction (1 for upregulated, -1 for downregulated).
#' @param variable A string specifying the grouping variable in \code{metadata} used for computing Cohen's d comparisons.
#' @param orientation A string indicating the layout of the combined heatmaps. Options are \code{\"grid\"}, \code{\"horizontal\"}, or \code{\"vertical\"}.
#'   Defaults to \code{\"grid\"}.
#' @param limits An optional numeric vector specifying the color scale limits (e.g., \code{c(min, max)}). If not provided, the limits are computed from the data.
#' @param cluster_columns Logical; whether to cluster columns in the heatmaps. Defaults to \code{TRUE}.
#' @param cluster_rows Logical; whether to cluster rows in the heatmaps. Defaults to \code{TRUE}.
#' @param subtitlewidth An integer specifying the width used in wrapping the gene set signature names for the heatmap subtitle.
#'
#' @return Depending on the orientation, a combined heatmap object is returned:
#' \describe{
#'   \item{grid}{A grid-arranged grob (produced by \code{grid.arrange}) of heatmaps.}
#'   \item{horizontal or vertical}{A combined Heatmap object (if not converted to a grob).}
#' }
#'
#' @details
#' This function first computes Cohen's d effect sizes and p-values for all gene signatures using \code{CohenD_allConditions}.
#' It then generates individual heatmaps (with cell text showing \"Cohen's d (p-value)\") for each signature.
#' Finally, based on the specified \code{orientation}, the heatmaps are combined into a single plot using either
#' \code{grid.arrange} (for a grid layout) or horizontal/vertical stacking operators.
#'
#' @examples
#' \dontrun{
#'   # Assume gene_data, sample_metadata, and gene_sets are defined
#'   Heatmap_CohenD(data = gene_data, metadata = sample_metadata,
#'                  gene_sets = gene_sets, variable = \"Condition\", orientation = \"grid\")
#' }
#'
#' @importFrom circlize colorRamp2
#' @importFrom gridExtra grid.arrange
#' @importFrom grid grid.text gpar
#'
#' @export
Heatmap_CohenD <- function(data, metadata, gene_sets, variable, orientation=c("grid","horizontal","vertical"), limits=NULL, cluster_columns = TRUE, cluster_rows= TRUE, widthTitle=20) {

  orientation <- match.arg(orientation)

  cohenlist <- CohenD_allConditions(data=data, metadata=metadata, gene_sets=gene_sets, variable=variable)

  # Create a list to store the heatmap objects
  heatmaps <- list()
  # Generate a heatmap for each signature
  for (signature in names(cohenlist)) {

    # Ensure Cohen's d and p-value are matrices
    cohen_d_mat <- as.matrix(cohenlist[[signature]]$CohenD)
    p_value_mat <- as.matrix(cohenlist[[signature]]$PValue)

    # Format text labels
    text_labels <- matrix(
      paste0(
        sprintf("%.2f", cohen_d_mat),  # Round Cohen's d to 2 decimals
        " \n(",
        format.pval(p_value_mat, digits = 1),  # Format p-value
        ")"
      ),
      nrow = nrow(cohen_d_mat),
      ncol = ncol(cohen_d_mat)
    )

    # Use user-provided limits or default ones
    if (is.null(limits)) limits <- c(min(cohen_d_mat), max(cohen_d_mat))

    # Create heatmap
    col_fun <- circlize::colorRamp2(limits, c("#F9F4AE" ,"#B44141"))
    ht <- ComplexHeatmap::Heatmap(
      t(cohen_d_mat),
      name = "Cohen's d",
      column_title = wrap_title(signature, subtitlewidth),
      show_row_names = TRUE,
      show_column_names = TRUE,
      cluster_rows = cluster_rows,
      cluster_columns = cluster_columns,
      cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(t(text_labels)[i, j], x, y, gp = gpar(fontsize = 10, col = "black"))
      },
      col = col_fun
    )



    if (orientation=="grid"){
      # Add the heatmap to the list; but in grob so we can combine
      heatmaps[[signature]] <- invisible(grid.grabExpr(ComplexHeatmap::draw(ht)))
    }  else {
      heatmaps[[signature]] <- ht
    }

  }

  # combine all heatmaps with a common color scale in only one row

  # The user can choose to do it in only one row or only one column
  if (orientation=="horizontal"){

    for (i in 1:length(heatmaps)) {
      heatmap_row <- heatmaps   # Extract remaining heatmaps
      Heatmap_Final <- heatmaps[[1]]  # Start with the first heatmap in the row
      while (length(heatmap_row) > 0) {
        Heatmap_Final <- Heatmap_Final + heatmap_row[[1]]  # Combine horizontally
        heatmap_row <- heatmap_row[-1]  # Remove the used heatmap
      }
    }

  } else if (orientation == "vertical"){

    for (i in 1:length(heatmaps)) {
      heatmap_row <- heatmaps   # Extract remaining heatmaps
      Heatmap_Final <- heatmaps[[1]]  # Start with the first heatmap in the row
      while (length(heatmap_row) > 0) {
        Heatmap_Final <- Heatmap_Final %v% heatmap_row[[1]]  # Combine horizontally
        heatmap_row <- heatmap_row[-1]  # Remove the used heatmap
      }
    }


  } else if (orientation == "grid"){

    # Determine grid layout
    num_signatures <- length(cohenlist)
    if (is.null(ncol) && is.null(nrow)) {
      ncol <- ceiling(sqrt(num_signatures))
      nrow <- ceiling(num_signatures / ncol)
    } else if (is.null(ncol)) {
      ncol <- ceiling(num_signatures / nrow)
    } else if (is.null(nrow)) {
      nrow <- ceiling(num_signatures / ncol)
    }


    Heatmap_Final <- gridExtra::grid.arrange(grobs = heatmaps,
                                  ncol = ncol,
                                  nrow = nrow)

  }

  return(Heatmap_Final)
}




#' Compute Cohen's d for All Gene Signatures Across Conditions
#'
#' Computes Cohen's d effect sizes and corresponding p-values for all gene signatures using scores calculated
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
#'
#' @return A named list where each element corresponds to a gene signature. Each signature element is a list with two components:\n
#' \describe{\n
#'   \item{CohenD}{A data frame where rows are methods and columns are group contrasts (formatted as \"Group1:Group2\"),
#'   containing the computed Cohen's d effect sizes.}\n
#'   \item{PValue}{A data frame with the same structure as \code{CohenD} containing the corresponding p-values.}\n
#' }
#'
#' @examples
#' \dontrun{\n
#'   # Assume gene_data is your gene expression data frame, sample_metadata is your metadata, and\n
#'   # gene_sets is a named list of gene sets.\n
#'   results <- CohenD_allConditions(data = gene_data, metadata = sample_metadata,\n
#'                                    gene_sets = gene_sets, variable = \"Condition\")\n
#'   # Access Cohen's d for a specific signature:\n
#'   results$Signature_A$CohenD\n
#' }
#'
#' @export
CohenD_allConditions <- function(data, metadata, gene_sets, variable) {

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

      # Compute Cohen's d and p-values
      cohen_results <- compute_cohen_d(df_method, variable)

      # Convert to named vectors (column names = comparisons)
      cohen_d_results[[method]] <- setNames(cohen_results$CohenD, paste0(cohen_results$Group1, ":", cohen_results$Group2))
      p_value_results[[method]] <- setNames(cohen_results$PValue, paste0(cohen_results$Group1, ":", cohen_results$Group2))
    }

    # Convert lists to data frames
    cohen_d_df <- as.data.frame(do.call(rbind, cohen_d_results))
    p_value_df <- as.data.frame(do.call(rbind, p_value_results))

    # Store results in the named list
    result_list[[signature]] <- list(CohenD = cohen_d_df, PValue = p_value_df)
  }

  # Step 7: Return the list
  return(result_list)
}



#' Compute Cohen's d Effect Size
#'
#' Computes the absolute Cohen's d effect size between two numeric vectors. This function returns
#' the absolute value of the difference in means divided by the pooled standard deviation.
#'
#' @param x A numeric vector representing the values for group 1.
#' @param y A numeric vector representing the values for group 2.
#'
#' @return A numeric value representing Cohen's d. Returns NA if either group has fewer than two observations
#'   or if the pooled standard deviation is zero.
#'
#' @examples
#' x <- c(1, 2, 3, 4, 5)
#' y <- c(6, 7, 8, 9, 10)
#' cohen_d(x, y)
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
  d <- abs(d) # dor this application we don't need the direction of the classification, just if it works
  return(d)
}


#' Compute Pairwise Cohen's d and P-Values
#'
#' Computes Cohen's d effect sizes and corresponding p-values for all pairwise comparisons of a grouping variable
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
#'   \item{CohenD}{The computed Cohen's d effect size for the comparison.}
#'   \item{PValue}{The p-value from a t-test comparing the two groups.}
#' }
#'
#' @examples
#' df <- data.frame(Group = rep(c("A", "B"), each = 10), score = rnorm(20))
#' compute_cohen_d(df, "Group")
#'
#' @keywords internal
compute_cohen_d <- function(dfScore, variable, quantitative_var="score") {

  # Get unique group values
  unique_groups <- unique(dfScore[[variable]])

  # Ensure there are at least two groups
  if (length(unique_groups) < 2) {
    stop("The grouping column must have at least two unique values.")
  }



  # Store results
  results <- data.frame(Group1 = character(), Group2 = character(),
                        CohenD = numeric(), PValue = numeric(), stringsAsFactors = FALSE)

  # Compute Cohen's d and p-value for all unique pairs
  combs <- combn(unique_groups, 2, simplify = FALSE)
  for (pair in combs) {
    x <- dfScore[dfScore[[variable]] == pair[1], quantitative_var, drop = TRUE]
    y <- dfScore[dfScore[[variable]] == pair[2], quantitative_var, drop = TRUE]

    d <- cohen_d(x, y)
    p_val <- t.test(x, y, var.equal = TRUE)$p.value

    results <- rbind(results, data.frame(Group1 = pair[1], Group2 = pair[2],
                                         CohenD = d, PValue = p_val, stringsAsFactors = FALSE))
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
#' @examples
#' # Suppose nested_list is structured as follows:
#' nested_list <- list(\n  ssGSEA = list(Signature_A = data.frame(sample = c(\"S1\", \"S2\"), score = c(0.5, 0.7))),\n  logmedian = list(Signature_A = data.frame(sample = c(\"S1\", \"S2\"), score = c(0.3, 0.4)))\n)\n  flatten_results(nested_list)
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



