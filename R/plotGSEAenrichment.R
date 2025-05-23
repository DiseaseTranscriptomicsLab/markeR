#' Plot GSEA Enrichment Results
#'
#' This function generates enrichment plots for gene sets using the `fgsea::plotEnrichment()` function.
#' It supports both individual plots (returned as a list) and a grid layout using `ggpubr::ggarrange()`.
#'
#' @param GSEA_results A named list of data frames containing GSEA results for each contrast.
#' Each data frame should have a column named `pathway` specifying the gene set, and columns `NES` and `padj` for results.
#' Output from \code{runGSEA}.
#' @param DEGList A named list of data frames containing differentially expressed genes (DEGs) for each contrast.
#' Each data frame must include a column named `t` with t-statistics for ranking genes.
#' Output from \code{calculateDE}.
#' @param gene_sets A named list of gene sets, where each entry is either:
#'   - A vector of gene names (unidirectional gene set)
#'   - A data frame with two columns: gene names and direction (+1 for enriched and -1 for depleted).
#' @param widthTitle Integer. The maximum width (in characters) for wrapping plot titles. Default is 24.
#' @param grid Logical. If `TRUE`, plots are arranged in a grid using `ggpubr::ggarrange()`. Default is `FALSE`.
#' @param nrow Integer. Number of rows for the grid layout (used only if `grid = TRUE`). If `NULL`, it is auto-calculated.
#' @param ncol Integer. Number of columns for the grid layout (used only if `grid = TRUE`). If `NULL`, it is auto-calculated.
#' @param titlesize Integer. Font size for plot titles. Default is 12.
#'
#' @return
#' If `grid = FALSE`, returns a named list of ggplot objects (each plot corresponding to a contrast-signature pair).
#' If `grid = TRUE`, returns a single ggplot object with all enrichment plots arranged in a grid.
#'
#' @examples
#' # Example GSEA results (mock data, missing columns if running by runGSEA)
#'
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
#' plot <- plotCombinedGSEA(GSEA_results, sig_threshold = 0.05, PointSize = 4)
#' print(plot)
#'
#' @import ggplot2 ggpubr fgsea
#' @export
plotGSEAenrichment <- function(GSEA_results, DEGList, gene_sets, widthTitle = 24, grid = FALSE, nrow=NULL, ncol=NULL, titlesize=12) {
  plot_list <- list()

  if (length(unique(GSEA_results[[1]]$pathway))==1 & length(names(GSEA_results))==1) {
    adjustpval <- FALSE
  } else {
    adjustpval <- TRUE
  }

  for (contrast in names(GSEA_results)) {
    deg_df <- DEGList[[contrast]]


    for (signature in unique(GSEA_results[[contrast]]$pathway)) {
      if (!(signature %in% names(gene_sets))) next  # Skip missing sets



      gs <- gene_sets[[signature]]



      # Retrieve NES and adjusted p-value for this contrast-signature pair
      gsea_res <- GSEA_results[[contrast]]
      gsea_row <- gsea_res[gsea_res$pathway == signature, ]

      stat_used <- gsea_row$stat_used

      # order ranks by stat used
      ranks <- setNames(deg_df[,stat_used, drop=T], rownames(deg_df))
      ranks <- sort(ranks, decreasing = TRUE)

      nes_value <- round(gsea_row$NES, 2)
      padj_value <- signif(gsea_row$padj, 3)
      if(stat_used=="B"){
        if (adjustpval){
          subtitle_text <- paste0("Altered Gene Set\nNES: ", nes_value, " | adj. p-value: ", padj_value)
        } else {
          subtitle_text <- paste0("Altered Gene Set\nNES: ", nes_value, " | p-value: ", padj_value)
        }

      } else {

        if (adjustpval){
          subtitle_text <- paste0("Enriched/Depleted Gene Set\nNES: ", nes_value, " | adj. p-value: ", padj_value)
        } else {
          subtitle_text <- paste0("Enriched/Depleted Gene Set\nNES: ", nes_value, " | p-value: ", padj_value)
        }

      }


      # Handle bidirectional gene sets
      if (is.data.frame(gs)) {
        gs_genes <- as.character(gs[[1]])
        directions <- as.numeric(gs[[2]])
        ranks_adjusted <- ranks
        idx <- which(names(ranks_adjusted) %in% gs_genes)
        lookup <- setNames(directions, gs_genes)
        ranks_adjusted[idx] <- ranks_adjusted[idx] * lookup[names(ranks_adjusted)[idx]]

        plot <- fgsea::plotEnrichment(gs_genes, sort(ranks_adjusted, decreasing = TRUE))

      } else if (is.vector(gs)) {
        # Unidirectional gene set
        gs_genes <- as.character(gs)
        plot <- fgsea::plotEnrichment(gs_genes, sort(ranks, decreasing = TRUE))

      } else {
        warning("Gene set '", signature, "' is not in a recognized format. Skipping.")
        next
      }

      plot <- plot +
        ggtitle(
          paste0(
            wrap_title(paste0("Enrichment: ", signature), width = widthTitle),
            wrap_title(paste0(" \n(", contrast, ")"), width = widthTitle)
          )
        ) +
        labs(subtitle = subtitle_text) +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = titlesize),
          plot.subtitle = element_text(hjust = 0.5, size = titlesize-1, color = "gray40")
        )

      plot_list[[paste(contrast, signature, sep = "_")]] <- plot
    }
  }

  # If grid = TRUE, arrange plots in a grid using ggarrange
  if (grid && length(plot_list) > 0) {

    n <- length(plot_list)

    # Determine grid layout.
    if (is.null(ncol) && is.null(nrow)) {
      ncol <- ceiling(sqrt(n))
      nrow <- ceiling(n / ncol)
    } else if (is.null(ncol)) {
      ncol <- ceiling(n / nrow)
    } else if (is.null(nrow)) {
      nrow <- ceiling(n / ncol)
    }

    return(ggarrange(plotlist = plot_list, ncol = ncol, nrow = nrow, align = "hv"))
  }

  return(plot_list)  # Return a list of individual plots if grid = FALSE
}
