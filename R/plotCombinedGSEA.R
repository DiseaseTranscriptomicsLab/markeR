#' Plot Combined GSEA Results
#'
#' This function creates a scatter plot visualizing multiple GSEA (Gene Set Enrichment Analysis) results
#' across different contrasts. Each point represents a pathway, where:
#' - The x-axis corresponds to the Normalized Enrichment Score (NES).
#' - The y-axis corresponds to the significance level (-log10 adjusted p-value).
#' - The color represents different pathways.
#' - The shape represents different contrasts.
#' - A dashed horizontal line marks the chosen significance threshold.
#'
#' @param GSEA_results A named list of data frames, where each data frame contains GSEA results for a contrast.
#' Each data frame should have the columns: `NES` (Normalized Enrichment Score), `padj` (adjusted p-value),
#' and `pathway` (pathway name). Output from \code{runGSEA}.
#' @param sig_threshold Numeric, default = 0.05. Adjusted p-value threshold for significance.
#' A dashed horizontal line is drawn at this threshold.
#' @param PointSize Numeric, default = 4. Size of the plotted points.
#' @param widthlegend Numeric, default = 16. Controls the width of pathway labels in the legend.
#'
#' @return A ggplot2 object displaying the combined GSEA results.
#'
#' @examples
#' # Example GSEA results (mock data)
#' GSEA_results <- list(
#' "Contrast1" = data.frame(
#' NES = rnorm(3),
#'   padj = runif(3),
#'   pathway = paste("Pathway", 1:3),
#'   stat_used = c("t", "B", "B")
#' ),
#' "Contrast2" = data.frame(
#'   NES = rnorm(3),
#'   padj = runif(3),
#'   pathway = paste("Pathway", 4:6),
#'   stat_used = c("t", "B", "B")
#' )
#' )
#'
#'
#' # Generate the plot
#' plotCombinedGSEA(GSEA_results, sig_threshold = 0.05, PointSize = 4)
#'
#' @import ggplot2
#' @import RColorBrewer
#' @export
plotCombinedGSEA <- function(GSEA_results, sig_threshold = 0.05, PointSize = 4, widthlegend=16) {

  # Combine all contrasts into a single data frame
  combined_data <- do.call(rbind, lapply(names(GSEA_results), function(contrast) {
    res <- GSEA_results[[contrast]]
    res$contrast <- contrast  # Add the contrast column
    return(res)
  }))

  # Add a significance column based on padj threshold
  combined_data$significance <- ifelse(combined_data$padj < sig_threshold, "Significant", "Not Significant")

  # Compute -log10(padj) for the y-axis
  combined_data$logpadj <- -log10(combined_data$padj)

  # Create a color palette for pathways
  # RColorBrewer has palettes for discrete color scales

  pathway_colors <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(length(unique(combined_data$pathway)))

  combined_data$pathway <- sapply(combined_data$pathway, function(x) wrap_title(x, widthlegend))

  # Create the plot
  plot <- ggplot2::ggplot(combined_data, ggplot2::aes(x = NES, y = logpadj,shape = contrast)) +
    ggplot2::geom_point(colour="black", size = PointSize) +
    ggplot2::geom_point(ggplot2::aes(colour = factor(pathway)) , size = PointSize-2.5) +
    ggplot2::geom_hline(yintercept = -log10(sig_threshold), linetype = "dashed", color = "black", size = .5)  +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = .5)  +
    ggplot2::scale_color_manual(values = pathway_colors) +  # Color pathways distinctly (fill)
    ggplot2::scale_shape_manual(values = c(15:(15 -1 + length((unique(combined_data$contrast)))))) +  # Different shapes for significant and non-significant points
    ggplot2::labs(x = "Normalized Enrichment Score (NES)",
         y = "-log10(Adj. p-value)",
         color = "Gene Set",
         shape = "Contrast") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "right",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    ) +
    #ggplot2::ggtitle("Combined GSEA Results by Contrast") +   # Title for the plot
    ggplot2::facet_grid(.~stat_used,
                        labeller = ggplot2::labeller(stat_used = c("t" = "Enriched/Depleted", "B" = "Altered")),   scales = "free", switch = "y" ) +
    theme(strip.background =element_rect(fill="white"))

  return(plot)
}
