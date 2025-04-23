#' Compute ROC Curves and AUC Values for Gene Signature Scores
#'
#' This function calculates Receiver Operating Characteristic (ROC) curves and
#' Area Under the Curve (AUC) values for gene signature scores across different
#' contrasts of a given categorical variable.
#'
#' @param data A matrix or data frame of gene expression data (genes as rows, samples as columns).
#' @param metadata A data frame containing sample metadata, including the grouping variable.
#' @param gene_sets A named list of gene sets, where each entry is a character vector of gene names.
#' @param method A character string specifying the score calculation method. Options: `"logmedian"`, `"ssGSEA"`, `"ranking"`, or `"all"`.
#' @param variable A character string specifying the categorical variable for group comparisons.#'
#' @param mode A string specifying the level of detail for contrasts.
#' Options are:
#' - `"simple"`: Pairwise comparisons (e.g., A - B).
#' - `"medium"`: Pairwise comparisons plus comparisons against the mean of other groups.
#' - `"extensive"`: All possible groupwise contrasts, ensuring balance in the number of terms on each side.
#'
#' @return A named list containing ROC curve data and AUC values for each method, signature, and contrast.
#'
#' @importFrom pROC roc auc
#' @importFrom stats as.formula
#' @importFrom utils combn
#'
#' @keywords internal
ROCAUC_Scores_Calculate <- function(data, metadata, gene_sets, method = c("logmedian", "ssGSEA", "ranking", "all"), variable, mode = c("simple","medium","extensive")) {

  method <- match.arg(method)
  mode <- match.arg(mode)

  if (method != "all") {
    listScores <- list(CalculateScores(data = data, metadata = metadata, gene_sets = gene_sets, method = method))
    names(listScores) <- method
  } else {
    listScores <- CalculateScores(data = data, metadata = metadata, gene_sets = gene_sets, method = method)
  }

  # Initialize result list
  result_list <- list()

  for (method_name in names(listScores)) {
    method_scores <- listScores[[method_name]]

    for (signature in names(method_scores)) {
      df <- method_scores[[signature]]

      if (!variable %in% colnames(df)) {
        stop(paste("Variable", variable, "not found in metadata!"))
      }

      # Generate all possible contrasts
      contrasts <- generate_all_contrasts(unique(df[[variable]]), mode = mode )
      contrasts <- remove_division(contrasts)

      for (contrast in contrasts) {

        df_subset <- create_contrast_column(df, variable, contrast) # column cohentest
        colnames(df_subset)[colnames(df_subset)=="cohentest"] <- "group_roc"

        # Compute ROC curve
        roc_curve <- pROC::roc(group_roc ~ score , data = df_subset, quiet=T)
        auc_value <- pROC::auc(roc_curve)

        # Flip AUC if needed (AUC must be >0.5 for meaningful interpretation)
        auc_value <- ifelse(auc_value < 0.5, 1 - auc_value, auc_value)

        # Store results
        result_list[[method_name]][[signature]][[contrast]] <- list(ROC = roc_curve, AUC = auc_value)
      }
    }
  }

  return(result_list)
}





#' Plot ROC Curves for Gene Signature Scores
#'
#' This function generates ROC curve plots for different gene signatures across multiple scoring methods.
#'
#' @param data A matrix or data frame of gene expression data.
#' @param metadata A data frame containing sample metadata.
#' @param gene_sets A named list of gene sets.
#' @param method A character string specifying the scoring method(s) (`"logmedian"`, `"ssGSEA"`, `"ranking"`, or `"all"`).
#' @param variable A character string specifying the categorical variable for group comparisons.
#' @param colors A named vector specifying colors for each method. Only one color is allowed, if method != "all".
#' Default colors are `c(logmedian = "#3E5587", ssGSEA = "#B65285", ranking = "#B68C52")`.
#' @param grid Logical; if `TRUE`, arranges plots in a grid.
#' @param spacing_annotation numeric value specifying the spacing between labels of AUC values. Default is 0.3.
#' @param modeA string specifying the level of detail for contrasts. Options are:
#' - `"simple"`: Performs the minimal number of pairwise comparisons between individual group levels (e.g., A - B, A - C). Default.
#' - `"medium"`: Includes comparisons between one group and the union of all other groups (e.g., A - (B + C + D)), enabling broader contrasts beyond simple pairs.
#' - `"extensive"`: Allows for all possible algebraic combinations of group levels (e.g., (A + B) - (C + D)), supporting flexible and complex contrast definitions.
#' @param ncol Optional numeric value specifying the number of columns in the grid layout for the combined plots.
#'   If `NULL`, there will be as many columns as contrasts. If this number is 1, then a near-square grid is computed.
#' @param nrow Optional numeric value specifying the number of rows in the grid layout.
#' If `NULL`, there will be as many columns as gene sets. If this number is 1, then a near-square grid is computed.
#' @param widthTitle Optional integer specifying the maximum width of the title before inserting line breaks.
#'   Titles break at `_`, `-`, or `:` where possible, or at the exact width if no such character is found.
#'   Default is `18`.
#' @param titlesize An integer specifying the text size for each of the heatmap titles. Default is 12.
#' @param title Title for the grid of plots.
#' @return A `ggplot2` or `ggarrange` object containing the ROC curve plots.
#'
#' @importFrom ggplot2 ggplot geom_line aes labs theme scale_color_manual
#' @importFrom ggpubr annotate_figure ggarrange
#'
#'@examples
#' # Example data
#' data <- as.data.frame(abs(matrix(rnorm(1000), ncol = 10)))
#' rownames(data) <- paste0("Gene", 1:100)  # Name columns as Gene1, Gene2, ..., Gene10
#' colnames(data) <- paste0("Sample", 1:10)  # Name rows as Sample1, Sample2, ..., Sample100
#'
#' # Metadata with sample ID and condition
#' metadata <- data.frame(
#'   SampleID = colnames(data),  # Sample ID matches the colnames of the data
#'   Condition = rep(c("A", "B"), each = 5)  # Two conditions (A and B)
#' )
#'
#' # Example gene set
#' gene_sets <- list(Signature1 = c("Gene1", "Gene2", "Gene3"))  # Example gene set
#'
#' # Call ROC_Scores function
#' ROC_Scores(data, metadata, gene_sets, method = "ssGSEA", variable = "Condition")
#'
#' @export
#'
ROC_Scores <- function(data, metadata, gene_sets, method = c("logmedian","ssGSEA","ranking","all"), variable,
                       colors = c(logmedian = "#3E5587", ssGSEA = "#B65285", ranking = "#B68C52"), grid = TRUE, spacing_annotation=0.3, ncol=NULL, nrow=NULL, mode=c("simple","medium","extensive"), widthTitle = 18, title=NULL, titlesize=12) {

  data_ROCAUC <- ROCAUC_Scores_Calculate(data = data, metadata = metadata, gene_sets = gene_sets, method = method, variable = variable, mode = mode)

  plot_list <- list()

  for (signature in names(data_ROCAUC[[1]])) {  # Iterate over signatures
    for (contrast in names(data_ROCAUC[[1]][[signature]])) {  # Iterate over contrasts

      # Initialize an empty data frame to store all methods
      combined_df <- data.frame()
      auc_values <- list()

      for (method_name in names(data_ROCAUC)) {  # Iterate over methods

        if (length(names(data_ROCAUC)) == 1){
          if (is.na(colors[method_name])) names(colors) <- method_name #if the user changed the color to only one, not named
        }


        roc_data <- data_ROCAUC[[method_name]][[signature]][[contrast]]

        # Create a data frame with FPR, TPR, and Method
        temp_df <- data.frame(
          FPR = rev(1 - roc_data$ROC$specificities),
          TPR = rev(roc_data$ROC$sensitivities),
          Method = method_name
        )

        # Combine into one large data frame
        combined_df <- rbind(combined_df, temp_df)

        # Calculate AUC for this method and contrast
        auc_value <- roc_data$AUC
        auc_values[[method_name]] <- auc_value

      }

      # Create the ROC plot with all methods on the same plot
      p <- ggplot2::ggplot(combined_df, ggplot2::aes(x = FPR, y = TPR, color = Method)) +
        ggplot2::geom_line(size = 1) +  # Plot all ROC curves on the same plot
        ggplot2::scale_color_manual(values = colors) +  # Ensure correct color mapping for each method
        ggplot2::labs(title = wrap_title(signature,widthTitle), subtitle= wrap_title(contrast,widthTitle), x = "False Positive Rate", y = "True Positive Rate") +
        ggplot2::theme_classic() +
        ggplot2::theme(legend.position = "none") + # Remove the default legend
        ggplot2::geom_abline(linetype = "dashed", color = "gray") +
        ggplot2::theme( plot.title = ggplot2::element_text(hjust = 0.5 ),
                        plot.subtitle = ggplot2::element_text(hjust = 0.5 ))

      # Add AUC text labels to the bottom-right corner
      auc_texts <- data.frame(Method = names(auc_values),
                              AUC = unlist(auc_values),
                              x = rep(1, length(auc_values)),  # Place all text at x = 1 (right edge)
                              y = seq(0.05, spacing_annotation, length.out = length(auc_values)))  # Adjust the vertical positions

      p <- p + ggplot2::geom_label(data = auc_texts,
                                  ggplot2::aes(x = x, y = y, label = paste0("AUC ", Method, " = ", round(AUC, 2), ""), color = Method),
                                  size = 3,
                                  vjust = 0,  # Adjust vertical position
                                  hjust = 1,  # Adjust horizontal position to align to the right
                                  inherit.aes = FALSE,
                                  fill = "white")  # Prevent inheritance of global aes from the main plot


      # Store plot
      plot_list[[paste(signature, contrast)]] <- p
    }
  }



  if (grid) {

    # Get number of signatures and contrasts
    n_signatures <- length(names(data_ROCAUC[[1]]))
    n_contrasts <- length(names(data_ROCAUC[[1]][[names(data_ROCAUC[[1]])[1]]]))

    # Case 1: If both nrow and ncol are provided, use them as is
    if (!is.null(nrow) && !is.null(ncol)) {
      # Use the provided nrow and ncol values as is
      combined_plot <- ggpubr::ggarrange(plotlist = plot_list, ncol = ncol, nrow = nrow, align = "h")

      # Case 2: If only one of nrow or ncol is provided, adjust the other for a square-like grid
    } else if (!is.null(nrow) || !is.null(ncol)) {
      if (is.null(nrow)) {
        nrow <- ceiling(length(plot_list) / ncol)  # Calculate nrow to make the grid more square-like
      }
      if (is.null(ncol)) {
        ncol <- ceiling(length(plot_list) / nrow)  # Calculate ncol to make the grid more square-like
      }
      combined_plot <- ggpubr::ggarrange(plotlist = plot_list, ncol = ncol, nrow = nrow, align = "h")

      # Case 3: If neither nrow nor ncol are provided, use the number of signatures and contrasts
    } else {
      nrow <- n_signatures
      ncol <- n_contrasts

      # If either nrow or ncol is 1, adjust the grid to make it more square
      if (nrow == 1 || ncol == 1) {
        nrow <- ceiling(sqrt(length(plot_list)))
        ncol <- ceiling(length(plot_list) / nrow)
      }
      combined_plot <- ggpubr::ggarrange(plotlist = plot_list, ncol = ncol, nrow = nrow, align = "h")
    }
    # Add the title to the grid if provided
    if (!is.null(title)) {
      combined_plot <- ggpubr::annotate_figure(combined_plot, top = grid::textGrob(title, gp = grid::gpar(cex = 1.3, fontsize = titlesize + 2)))
    }

     return(combined_plot)

  } else {
     return(plot_list)
  }
}

#' Generate Heatmaps for AUC Scores using ggplot2
#'
#' This function computes AUC scores for multiple gene signatures and scoring methods, and generates a heatmap
#' for each gene signature. The heatmap displays the AUC scores, with the contrasts as rows and methods as columns.
#' The heatmaps are then arranged in a grid layout.
#'
#' @param data A data frame of gene expression data with genes as rows and samples as columns.
#'   Row names should contain gene names and column names sample identifiers.
#' @param metadata A data frame of sample metadata. The first column must contain sample
#'   identifiers matching those in \code{data}.
#' @param gene_sets A named list of gene sets.
#' @param method A character string specifying the scoring method(s) (`"logmedian"`, `"ssGSEA"`, `"ranking"`, or `"all"`).
#' @param mode A string specifying the level of detail for contrasts.
#' Options are:
#' - `"simple"`: Pairwise comparisons (e.g., A - B).
#' - `"medium"`: Pairwise comparisons plus comparisons against the mean of other groups.
#' - `"extensive"`: All possible groupwise contrasts, ensuring balance in the number of terms on each side.
#' @param variable A string specifying the grouping variable in \code{metadata} used for computing AUC comparisons.
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
#'   \item{data}{A list containing the AUC scores for each gene signature, as computed by \code{ROCAUC_Scores_Calculate}.}
#' }
#'
#' @details
#' The function first calculates AUC scores for each gene signature using \code{ROCAUC_Scores_Calculate}. The resulting matrices are converted to a long format so that each cell in the heatmap can display the AUC value. A title for each heatmap is dynamically created.
#'
#' The heatmaps are then adjusted to display axis text and ticks only for the left-most column and bottom row, and combined into a grid layout. If neither \code{nrow} nor \code{ncol} are specified, the layout is automatically determined to best approximate a square grid.
#'
#' @examples
#' \dontrun{
#'   result <- AUC_Scores(
#'     data = gene_data,
#'     metadata = sample_metadata,
#'     gene_sets = gene_sets,
#'     method = "ssGSEA",
#'     variable = "Condition",
#'     nrow = 2,
#'     ncol = 3,
#'     limits = c(0, 1),
#'     widthTitle = 30,
#'     titlesize = 14,
#'     ColorValues = c("#F9F4AE", "#B44141")
#'   )
#'   print(result$plt)
#' }
#'
#' @importFrom ggplot2 ggplot geom_tile geom_text labs scale_fill_gradientn theme_minimal element_text element_blank element_line margin
#' @importFrom ggpubr ggarrange
#'
#' @export
AUC_Scores <- function(data, metadata, gene_sets, method = c("logmedian", "ssGSEA", "ranking", "all"), mode = c("simple","medium","extensive"), variable, nrow = NULL, ncol = NULL, limits = NULL, widthTitle = 22, titlesize = 12, ColorValues = c("#F9F4AE", "#B44141"), title = NULL) {

  # Calculate the AUC scores using the given method
  auc_list <- ROCAUC_Scores_Calculate(data = data, metadata = metadata, gene_sets = gene_sets, method = method, variable = variable, mode = mode)

  heatmaps <- list()

  # Loop over the gene signatures
  for (signature_name in names(auc_list[[1]])) {

    # Initialize a matrix to store AUC values for this signature
    auc_matrix <- matrix(nrow = length(auc_list[[1]][[signature_name]]), ncol = length(auc_list))
    rownames(auc_matrix) <- names(auc_list[[1]][[signature_name]])  # Contrasts
    colnames(auc_matrix) <- names(auc_list)  # Methods

    # Fill the matrix with AUC values
    for (method_name in names(auc_list)) {
      for (contrast_name in names(auc_list[[method_name]][[signature_name]])) {
        auc_matrix[contrast_name, method_name] <- auc_list[[method_name]][[signature_name]][[contrast_name]]$AUC
      }
    }

    # Convert the AUC matrix into long format for ggplot
    long_data <- as.data.frame(as.table(auc_matrix))
    colnames(long_data) <- c("Contrast", "Method", "AUC")

    # Generate text labels showing the AUC value
    long_data$label <- sprintf("%.2f", long_data$AUC)

    # Wrap the signature title using an internal helper function
    signature_title <- wrap_title(signature_name, widthTitle)

    limits <- if (is.null(limits)) c(0.5 , 1) else limits

    # Create heatmap using ggplot2
    p <- ggplot2::ggplot(long_data, ggplot2::aes(x = Method, y = Contrast, fill = AUC)) +
      ggplot2::geom_tile() +
      ggplot2::geom_text(aes(label = label), color = "black", size = 3) +
      ggplot2::scale_fill_gradientn(colors = ColorValues, limits = limits) +
      ggplot2::labs(title = signature_title, x = NULL, y = NULL, fill = "AUC") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        plot.title = ggplot2::element_text(hjust = 0.5, size = titlesize)
      )

    heatmaps[[signature_name]] <- p
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

  # Add the title to the grid if provided
  if (!is.null(title)) {
    plt <- ggpubr::annotate_figure(plt, top = grid::textGrob(title, gp = grid::gpar(cex = 1.3, fontsize = titlesize + 2)))
  }

  return(plt)
}

