#' Create a Lollipop Plot for GSEA Results
#'
#' This function generates a lollipop plot to visualize Gene Set Enrichment Analysis (GSEA) results.
#' Pathways are shown on the y-axis, while the Normalized Enrichment Score (NES) is shown on the x-axis.
#' The color of the lollipops represents the adjusted p-values (`padj`), with a custom color gradient.
#' It supports multiple contrasts and can combine individual plots into a grid layout.
#'
#' @param GSEA_results A named list of data frames, each containing the GSEA results for a specific contrast. Output from \code{runGSEA}.
#' Each data frame must include the following columns:
#' \describe{
#'   \item{pathway}{A character vector of pathway names.}
#'   \item{NES}{A numeric vector of Normalized Enrichment Scores for the pathways.}
#'   \item{padj}{A numeric vector of adjusted p-values for the pathways.}
#' }
#' @param padj_limit A numeric vector of length 2 that defines the limits for the adjusted p-value (`padj`) color scale.
#' The first value is the minimum limit, and the second value is the maximum limit. Default is `c(0, 0.1)`.
#' If a pathway’s `padj` exceeds the maximum value in this range, it will be assigned the `high_color`.
#' @param low_color A string specifying the color for the low end of the adjusted p-value gradient. Default is `"blue"`.
#' @param mid_color A string specifying the color for the middle of the adjusted p-value gradient. Default is `"white"`. Will correspond to the value of \code{sig_threshold}.
#' @param high_color A string specifying the color for the high end of the adjusted p-value gradient. Default is `"red"`.
#' @param sig_threshold A numeric value that sets the midpoint for the color scale. Typically used for the significance threshold. Default is `0.05`.
#' @param grid A logical value indicating whether to arrange individual plots into a grid layout. If `TRUE`, the function combines all plots into a grid. Default is `FALSE`.
#' @param nrow A numeric value specifying the number of rows to arrange the plots into if `grid = TRUE`. If `NULL`, the function calculates this automatically. Default is `NULL`.
#' @param ncol A numeric value specifying the number of columns to arrange the plots into if `grid = TRUE`. If `NULL`, the function calculates this automatically. Default is `NULL`.
#' @param widthlabels A numeric value specifying the maximum width for pathway names. If a pathway name exceeds this width, it will be wrapped to fit. Default is `18`.
#' @param title A character string for the title of the combined plot (used only when `grid = TRUE`). Default is `NULL`.
#' @param titlesize A numeric value specifying the font size for the title (used only when `grid = TRUE`). Default is `12`.
#'
#' @return If `grid = FALSE`, a list of `ggplot` objects is returned, each corresponding to a contrast. If `grid = TRUE`, a single `ggplot` object is returned, representing the combined grid of plots.
#'
#' @details
#' The function creates a lollipop plot for each contrast in the `GSEA_results` list. Each plot includes:
#' \itemize{
#'   \item A vertical segment for each pathway, where the x-coordinate represents the NES and the y-coordinate represents the pathway.
#'   \item A colored point at the end of each segment, where the color represents the adjusted p-value (`padj`), mapped using a custom color gradient.
#' }
#'
#' If a pathway's `padj` value exceeds the maximum value in `padj_limit`, the corresponding pathway is colored using the `high_color`.
#' Additionally, missing values (`NA`) for `padj` are assigned the `high_color` by setting `na.value = high_color`.
#' Pathway names are wrapped using the `wrap_title` function to fit within the specified width (`widthlabels`).
#'
#' @examples
#' # Example GSEA results (mock data)
#' GSEA_results <- list(
#'   contrast1 = data.frame(
#'     pathway = c("Pathway A", "Pathway B", "Pathway C"),
#'     NES = c(1.5, -2.3, 1.2),
#'     padj = c(0.02, 0.08, 0.01)
#'   ),
#'   contrast2 = data.frame(
#'     pathway = c("Pathway D", "Pathway E", "Pathway F"),
#'     NES = c(-1.8, 2.5, -1.3),
#'     padj = c(0.04, 0.06, 0.10)
#'   )
#' )
#'
#' # Generate individual plots without grid
#' plot_list <- plotNESlollipop(GSEA_results)
#'
#' # Generate combined grid of plots with custom title
#' combined_plot <- plotNESlollipop(GSEA_results, grid = TRUE, title = "GSEA Results Overview", titlesize = 14)
#'
#' @import ggplot2
#' @importFrom ggpubr annotate_figure ggarrange
#' @import grid
#' @export
plotNESlollipop <- function(GSEA_results, padj_limit = c(0, 0.1), low_color = "blue", mid_color = "white", high_color = "red", sig_threshold = 0.05, grid = FALSE, nrow = NULL, ncol = NULL, widthlabels=18, title=NULL, titlesize=12) {


  plot_list <- list()

  for (contrast in names(GSEA_results)) {
    res <- GSEA_results[[contrast]]


    # Ensure contrast ordering
    res$pathway <- sapply(res$pathway, function(x) wrap_title(x, widthlabels))
    res$pathway <- factor(res$pathway, levels = res$pathway[order(res$NES)])


    plot <- ggplot2::ggplot(res, ggplot2::aes(x = NES, y = pathway,fill = padj)) +
      ggplot2::geom_segment(ggplot2::aes(yend = pathway, xend = 0), size = .5) +
      ggplot2::geom_point( shape = 21, stroke = 1.2, color="black", size=4) +
      ggplot2::scale_fill_gradient2(low = low_color, mid = mid_color, high = high_color, midpoint = sig_threshold, limits = padj_limit, na.value=high_color) +
      ggplot2::labs(title = contrast, x = "Normalized Enrichment Score (NES)", y = "Gene Set", color = "Adj. p-value", fill = "Adj. p-value") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        legend.position = "right"
      )

    plot_list[[contrast]] <- plot
  }

  # Arrange plots in grid if requested
  if (grid && length(plot_list) > 1) {
    n <- length(plot_list)

    if (is.null(ncol) && is.null(nrow)) {
      ncol <- ceiling(sqrt(n))
      nrow <- ceiling(n / ncol)
    } else if (is.null(ncol)) {
      ncol <- ceiling(n / nrow)
    } else if (is.null(nrow)) {
      nrow <- ceiling(n / ncol)
    }

    combined_plot <- ggpubr::ggarrange(plotlist = plot_list, ncol = ncol, nrow = nrow, align = "hv", common.legend = T, legend = "right")
    if (!is.null(title)) combined_plot <- ggpubr::annotate_figure(combined_plot,
                                                                  top = grid::textGrob(title, gp = grid::gpar(cex = 1.3, fontsize = titlesize)))

    return(combined_plot)
  }

  return(plot_list)  # Return list if not in grid mode
}
