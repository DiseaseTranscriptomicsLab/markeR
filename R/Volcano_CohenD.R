#' Volcano Plot of Cohen's d Effect Sizes and Adjusted p-values
#'
#' This function computes Cohen's d effect sizes and adjusted p-values for multiple gene signatures across defined contrasts,
#' and generates a volcano plot (Cohen's d vs -log10(padj)) using `ggplot2`. Each point represents a method-signature pair,
#' faceted by contrast.
#'
#' @param cohenlist A named list where each element corresponds to a gene signature. Output of `CohenD_allConditions`. Each signature element is a list with three components:
#' \describe{
#'   \item{CohenD}{A data frame where rows are methods and columns are group contrasts (formatted as \"Group1:Group2\"),
#'   containing the computed Cohen\'s d effect sizes.}
#'   \item{PValue}{A data frame with the same structure as \code{CohenD} containing the corresponding p-values.}
#'   \item{padj}{A data frame with the same structure as \code{PValue} containing the corresponding p-values corrected using the BH method, for all signatures and contrasts,
#'   and by method.}
#' }
#' @param titlesize Integer. Size of the facet strip titles. Default is 12.
#' @param ColorValues Character vector of colors used to distinguish signatures. If NULL, colors are automatically generated.
#' @param title Optional title for the overall plot.
#' @param widthlegend Integer. Width used to wrap long signature names. Default is 22.
#' @param PointSize Numeric. Size of the points in the plot. Default is 3.
#' @param sig_threshold Numeric. Adjusted p-value threshold for significance. Default is 0.05.
#' @param cohend_threshold Numeric. Effect size threshold. Default is 0.5.
#' @param colorPalette Character. Name of RColorBrewer palette to use if `ColorValues` is not provided. Default is "Paired".
#' @param ncol Optional numeric value specifying the number of columns in the ggplot facet. If `NULL`, a near-square grid is computed.
#' @param nrow Optional numeric value specifying the number of rows in the grid layout. If `NULL`, a near-square grid is computed.
#'
#' @return A `ggplot` object showing a faceted volcano plot of Cohen's d effect sizes across signatures and methods for each contrast.
#'
#' @seealso \code{\link{CohenD_allConditions}}
#'
#' @importFrom ggplot2 ggplot geom_point geom_vline geom_hline facet_wrap labs scale_color_manual scale_shape_manual theme_bw ggtitle theme element_text element_rect
#' @importFrom RColorBrewer brewer.pal
#' @export
Volcano_CohenD <- function(cohenlist,
                           titlesize = 12,
                           ColorValues = NULL,
                           title = NULL,
                           widthlegend = 22,
                           PointSize = 3,
                           sig_threshold = 0.05,
                           cohend_threshold = 0.5,
                           colorPalette = "Paired",
                           nrow=NULL, ncol=NULL) {

  # Convert nested list to long format
  rows <- list()
  for (signature in names(cohenlist)) {
    sig_data <- cohenlist[[signature]]
    cohend_mat <- sig_data$CohenD
    padj_mat <- sig_data$padj

    for (method in rownames(cohend_mat)) {
      for (contrast in colnames(cohend_mat)) {
        rows[[length(rows) + 1]] <- data.frame(
          signature = signature,
          contrast = contrast,
          method = method,
          cohend = cohend_mat[method, contrast],
          padj = padj_mat[method, contrast],
          stringsAsFactors = FALSE
        )
      }
    }
  }

  final_df <- do.call(rbind, rows)

  # Wrap long signature names
  final_df$signature <- sapply(final_df$signature, function(x) wrap_title(x, widthlegend))

  # Handle colors
  if (is.null(ColorValues)) {
    ColorValues <- colorRampPalette(RColorBrewer::brewer.pal(12, colorPalette))(length(unique(final_df$signature)))
  }

  # Generate plot
  plt <- ggplot2::ggplot(final_df, ggplot2::aes(x = cohend, y = -log10(padj), shape = method)) +
    ggplot2::geom_point(colour = "black", size = PointSize) +
    ggplot2::geom_point(ggplot2::aes(colour = signature), size = PointSize - 1.5) +
    ggplot2::facet_wrap(. ~ contrast, scales = "free") +
    ggplot2::geom_hline(yintercept = -log10(sig_threshold), linetype = "dashed", color = "black", size = 0.5) +
    ggplot2::geom_vline(xintercept = cohend_threshold, linetype = "dashed", color = "black", size = 0.5) +
    ggplot2::scale_color_manual(values = ColorValues) +
    ggplot2::scale_shape_manual(values = 15:(15 + length(unique(final_df$method)) - 1)) +
    ggplot2::labs(
      x = "Cohen's d",
      y = "-log10(Adj. p-value)",
      color = "Signature",
      shape = "Method"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "right",
      strip.text = element_text(size = titlesize, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      strip.background = ggplot2::element_rect(fill = "white")
    ) +
    ggplot2::ggtitle(if (!is.null(title)) title else "Cohen's d Volcano Plot")



  return(list(plt=plt))
}
