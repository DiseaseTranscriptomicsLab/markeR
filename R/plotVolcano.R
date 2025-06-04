#' Volcano Plots from Differential Expression Results
#'
#' This function creates a composite volcano plot grid from a list of differential expression results., or a single
#' volcano if no genes to highlight are provided and no more than one contrast is used.
#' For each contrast (provided in \code{DEResultsList}) and gene signature (from the \code{genes} argument),
#' a volcano plot is generated using the specified x and y statistics. By default, if \code{invert = FALSE}
#' and more than one gene signature is provided (i.e. the names in \code{genes} are not "ALL" or "genes"),
#' the plots are arranged with gene signatures in rows and contrasts in columns. When \code{invert = TRUE},
#' the arrangement is reversed (signatures in columns and contrasts in rows). If only one gene signature is provided,
#' an automatic grid is computed.
#'
#' @param DEResultsList A named list of data frames containing differential expression results for each contrast.
#'   Each data frame should have row names corresponding to gene names and include columns for the x and y statistics.
#'   Output from \code{calculateDE}.
#' @param genes Optional. A list of gene signatures to highlight. Each element may be a data frame (in which case its first column is extracted)
#'   or a vector of gene names. If \code{NULL}, no genes will be highlighted.
#' @param N Optional. An integer specifying the number of top (and bottom) genes to annotate with text labels.
#' @param x Character. The column name in the differential expression results to use for the x-axis (default is \code{"logFC"}).
#' @param y Character. The column name to use for the y-axis (default is \code{"-log10(adj.P.Val)"}). When using this default,
#'   threshold values for \code{threshold_y} should be provided in non-log scale (e.g., 0.05).
#' @param pointSize Numeric. The size of points in the volcano plots (default is 2).
#' @param color Character. The color used to highlight interesting genes based on thresholds (default is \code{"#6489B4"}).
#' @param highlightcolor Character. The color used to highlight genes belonging to the specified gene signatures (default is \code{"#05254A"}), if direction is not known or not specified.
#' @param highlightcolor_upreg Character. The color used to highlight upregulated genes belonging to the specified gene signatures (default is \code{"#038C65"}).
#' @param highlightcolor_downreg  Character. The color used to highlight downregulated genes belonging to the specified gene signatures (default is \code{"#8C0303"}).
#' @param nointerestcolor Character. The color for non-interesting genes (default is \code{"#B7B7B7"}).
#' @param threshold_y Numeric. A threshold value for the y-axis statistic. If \code{y} is \code{"-log10(adj.P.Val)"},
#'   the value should be provided as a non-log value (e.g., 0.05) and will be transformed internally.
#' @param threshold_x Numeric. A threshold value for the x-axis statistic.
#' @param xlab Optional. A label for the x-axis; if \code{NULL}, the value of \code{x} is used.
#' @param ylab Optional. A label for the y-axis; if \code{NULL}, the value of \code{y} is used.
#' @param ncol Optional. The number of columns for arranging plots in the grid. Only applicable if \code{genes} is \code{NULL}.
#' @param nrow Optional. The number of rows for arranging plots in the grid.
#' @param title Optional. A main title for the entire composite plot.
#' @param labsize Numeric. The font size for label annotations (default is 10). The title size will be this value + 4.
#' @param widthlabs Numeric. The width parameter to pass to the \code{wrap_title()} function for wrapping long labels (default is 20).
#' @param invert Logical. If \code{FALSE} (default), the grid is arranged with gene signatures in rows and contrasts in columns.
#'   If \code{TRUE}, the arrangement is inverted (gene signatures in columns and contrasts in rows).
#'
#' @return A composite plot (a ggplot object) arranged as a grid of volcano plots with annotated labels.
#'
#' @details
#' This function generates a volcano plot for each combination of gene signature (from \code{genes}) and contrast
#' (from \code{DEResultsList}). It uses the specified \code{x} and \code{y} statistics to plot points via \code{ggplot2}.
#' Non-interesting genes are plotted using \code{nointerestcolor}, while genes in the specified gene signature (if not "ALL")
#' are highlighted using \code{highlightcolor}. Optionally, the top and bottom \code{N} genes can be annotated with text labels
#' (using \code{ggrepel::geom_text_repel}). Threshold lines for the x and/or y axes are added if \code{threshold_x} or \code{threshold_y}
#' are provided. The individual plots are arranged into a grid using \code{ggpubr::ggarrange} and annotated with labels using
#' \code{ggpubr::annotate_figure} and \code{grid::textGrob}. The custom \code{wrap_title()} function is used to wrap long labels.
#'
#'
#'Additionally, the function allows:
#'
#' \itemize{
#'   \item Plotting of differentially expressed genes based on provided statistics (e.g., \code{x = "logFC"} and \code{y = "-log10(adj.P.Val)"}).
#'   \item Coloring of non-interesting genes and highlighting genes belonging to specific gene signatures.
#'   \item Annotation of the top \code{N} genes with text labels (using \code{ggrepel::geom_text_repel}).
#'   \item Addition of threshold lines for the x and/or y axes.
#' }
#'
#' @examples
#' \dontrun{
#' # Create example data:
#' data <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' rownames(data) <- paste0("gene", 1:100)
#' colnames(data) <- paste0("sample", 1:10)
#' metadata <- data.frame(sample = colnames(data), X = rep(c("A", "B"), each = 5))
#'
#' # Example differential expression results (for two contrasts):
#' de_results <- list(
#'   Contrast1 = data.frame(logFC = rnorm(100), `-log10(adj.P.Val)` = runif(100, 0, 5),
#'                          row.names = paste0("gene", 1:100)),
#'   Contrast2 = data.frame(logFC = rnorm(100), `-log10(adj.P.Val)` = runif(100, 0, 5),
#'                          row.names = paste0("gene", 1:100))
#' )
#'
#' # Basic volcano plot grid with default settings:
#' plotVolcano(de_results, genes = NULL, N = NULL,
#'             x = "logFC", y = "-log10(adj.P.Val)", pointSize = 2,
#'             color = "#6489B4", highlightcolor = "#05254A", nointerestcolor = "#B7B7B7",
#'             threshold_y = NULL, threshold_x = NULL,
#'             xlab = NULL, ylab = NULL, ncol = NULL, nrow = NULL, title = "Volcano Plot Grid",
#'             labsize = 10, widthlabs = 20, invert = FALSE)
#' }
#'
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggpubr ggarrange annotate_figure
#' @importFrom grid textGrob
#' @importFrom gridExtra arrangeGrob
#' @export
#'
plotVolcano <- function(DEResultsList, genes = NULL, N = NULL,
                        x = "logFC", y = "-log10(adj.P.Val)", pointSize = 2,
                        color = "#6489B4", highlightcolor = "#05254A", highlightcolor_upreg = "#038C65", highlightcolor_downreg = "#8C0303",nointerestcolor = "#B7B7B7",
                        threshold_y = NULL, threshold_x = NULL,
                        xlab = NULL, ylab = NULL, ncol = NULL, nrow = NULL, title = NULL,
                        labsize = 10, widthlabs = 20, invert = FALSE) {

  ## Helper: extract genes from list elements (if data.frame, extract first column)
  extract_direction <- function(lst) {
    lapply(lst, function(x) {
      if (is.data.frame(x)) {
        return(x[[2]])
      } else if (is.vector(x)) {
        vec <- rep("test", length(x))
        return(vec)
      } else {
        warning("Unexpected list element type: ", class(x))
        return(NULL)
      }
    })
  }

  extract_genes <- function(lst) {
    lapply(lst, function(x) {
      if (is.data.frame(x)) {
        return(x[[1]])
      } else if (is.vector(x)) {
        return(x)
      } else {
        warning("Unexpected list element type: ", class(x))
        return(NULL)
      }
    })
  }

  ## Process genes argument
  if (is.null(genes)) {
    genes <- list(ALL = row.names(DEResultsList[[1]]))
    direction <- NULL
  } else if (!is.list(genes)) {
    genes <- list(genes = genes)
    direction <- NULL
  } else if (is.list(genes)) {
    direction <- extract_direction(genes)
    genes <- extract_genes(genes)
  }

  plotList_signatures <- list()

  ## Loop over each signature (names in 'genes')
  for (sig in names(genes)) {
    genessig <- genes[[sig]]
    if (!is.null(direction)) directionsig <- direction[[sig]]

    plotList_contrasts <- list()

    ## Loop over each contrast (names in DEResultsList)
    for (contrast in names(DEResultsList)) {
      fit <- DEResultsList[[contrast]]

      ## Base plot with provided x and y aesthetics
      p <- ggplot2::ggplot(fit, ggplot2::aes_string(x = x, y = y))
      if (is.null(xlab)) xlab <- x
      if (is.null(ylab)) ylab <- y
      p <- p +
        ggplot2::theme_bw() +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                       plot.subtitle = ggplot2::element_text(hjust = 0.5),
                       axis.text.y = ggplot2::element_text(angle = 90, hjust = 0.5)) +
        ggplot2::labs(x = xlab, y = ylab, title=contrast)

      p <- p + ggplot2::geom_point(alpha = 0.4, color = nointerestcolor, size = pointSize)

      ## Highlight genes if signature is not "ALL"
      if (sig != "ALL") {
      if (is.null(direction)){
        p <- p + ggplot2::geom_point(data = fit[row.names(fit) %in% genessig, ],
                                     color = highlightcolor, size = pointSize)
      } else {

        upreg_genes <- genessig[directionsig==1]
        downreg_genes <- genessig[directionsig==-1]
        other_genes <- genessig[!directionsig %in% c(1,-1)]
        #  fit_subset <- fit[row.names(fit) %in% genessig,]
        #  fit_subset$genes <- row.names(fit_subset)
        #  fit_subset$Direction <- ifelse(fit_subset$genes %in% upreg_genes, "Upregulated", ifelse(fit_subset$genes %in% downreg_genes, "Downregulated", "No Information"))
        #
        # p <- p + ggplot2::geom_point(data = fit_subset, aes(color = Direction), size = pointSize) +
        #   ggplot2::scale_color_manual(values=c("Downregulated" = highlightcolor_downreg,
        #                                 "Upregulated" = highlightcolor_upreg,
        #                                 "No Information" = highlightcolor)) +
        #   ggplot2::theme(legend.position=NULL)

        # Working version, but upreg genes are masked by downreg genes
        # if (length(upreg_genes)>0)  p <- p + ggplot2::geom_point(data = fit[row.names(fit) %in% upreg_genes, ], color = highlightcolor_upreg, size = pointSize, alpha=0.8)
        # if (length(downreg_genes)>0)  p <- p + ggplot2::geom_point(data = fit[row.names(fit) %in% downreg_genes, ], color = highlightcolor_downreg, size = pointSize, alpha=0.8)
        # if (length(other_genes)>0)  p <- p + ggplot2::geom_point(data = fit[row.names(fit) %in% other_genes, ], color = highlightcolor, size = pointSize, alpha=0.8)

        # Approach where genes are plotted in a random order

        combined <- list()

        # Only include non-empty subsets
        if (length(upreg_genes) > 0) {
          subset_up <- fit[row.names(fit) %in% upreg_genes, ]
          subset_up$Direction <- "Upregulated"
          combined[[length(combined) + 1]] <- subset_up
        }

        if (length(downreg_genes) > 0) {
          subset_down <- fit[row.names(fit) %in% downreg_genes, ]
          subset_down$Direction <- "Downregulated"
          combined[[length(combined) + 1]] <- subset_down
        }

        if (length(other_genes) > 0) {
          subset_other <- fit[row.names(fit) %in% other_genes, ]
          subset_other$Direction <- "No Information"
          combined[[length(combined) + 1]] <- subset_other
        }

        # Combine non-empty subsets and shuffle
        if (length(combined) > 0) {
          plot_data <- do.call(rbind, combined)
          set.seed(1234)  # for reproducibility
          plot_data <- plot_data[sample(nrow(plot_data)), ]

          # Add to plot
          p <- p + ggplot2::geom_point(data = plot_data,
                                       aes(color = Direction),
                                       size = pointSize,
                                       alpha = 0.8) +
            ggplot2::scale_color_manual(values = c(
              "Upregulated" = highlightcolor_upreg,
              "Downregulated" = highlightcolor_downreg,
              "No Information" = highlightcolor
            )) +
            ggplot2::theme(legend.position = "none")
        }


      }

      }

      ## Annotate top N genes if requested
      if (!is.null(N)) {
        genes_stat <- fit[order(fit[, x], decreasing = TRUE), ]
        annotationgenes <- row.names(genes_stat)[c(1:N, (nrow(genes_stat) - N + 1):nrow(genes_stat))]
        p <- p +
          ggplot2::geom_point(data = fit[annotationgenes, ], color = highlightcolor, size = pointSize) +
          ggrepel::geom_text_repel(data = fit[annotationgenes, ],
                                   label = row.names(fit[annotationgenes, ]),
                                   nudge_x = 0, nudge_y = 1, max.overlaps = N, force = 10)
      }

      ## Add threshold lines if specified
      # color  interesting genes based on thresholds
      interesting_genes_df <- fit
      #
      # if (is.null(threshold_y) && y=="-log10(adj.P.Val)") {
      #   y <- "adj.P.Val"
      #   threshold_y <- -log10(0.05)
      #   threshold_y_subset <- 0.05
      #   interesting_genes_df <- interesting_genes_df[interesting_genes_df[y] <= threshold_y_subset,]
      #   # p <- p +
      #   #   ggplot2::geom_point(data= fit[fit[y] <= threshold_y_subset,], color=nointerestcolor, size=pointSize,alpha=0.4) +
      #   #   geom_hline(yintercept = threshold_y, linetype="dashed", size=1)
      #
      # } else

      if (!is.null(threshold_y) && y=="-log10(adj.P.Val)"){ # user specifies new value for threshold; should be without log10
        y <- "adj.P.Val"
        threshold_y_subset <- threshold_y
        threshold_y <- -log10(threshold_y)
        interesting_genes_df <- interesting_genes_df[interesting_genes_df[y] <= threshold_y_subset,]
        p <- p +
          geom_hline(yintercept = threshold_y, linetype="dashed", size=1)

      } else if (!is.null(threshold_y)){
        interesting_genes_df <- interesting_genes_df[interesting_genes_df[y] >= threshold_x,]
        p <- p +
          geom_hline(yintercept = threshold_y, linetype="dashed", size=1)
      }
      if (!is.null(threshold_x)){
        interesting_genes_df <- interesting_genes_df[abs(interesting_genes_df[x]) >= threshold_x,]
        p <- p +
          geom_vline(xintercept = threshold_x, linetype="dashed", size=1) +
          geom_vline(xintercept = -threshold_x, linetype="dashed", size=1)
      }

      if (!is.null(threshold_x) || !is.null(threshold_y)){
        p <- p +
          ggplot2::geom_point(data= interesting_genes_df, color=color, size=pointSize,alpha=0.4)
      }


      plotList_contrasts[[contrast]] <- p
    } # end contrast loop

    ## Arrange contrast plots for current signature
    if (any(names(genes) %in% c("ALL", "genes"))) {
      n <- length(plotList_contrasts)
      if (is.null(ncol) && is.null(nrow)) {
        ncol <- ceiling(sqrt(n))
        nrow <- ceiling(n / ncol)
      } else if (is.null(ncol)) {
        ncol <- ceiling(n / nrow)
      } else if (is.null(nrow)) {
        nrow <- ceiling(n / ncol)
      }
      plotList_signatures[[1]] <- ggpubr::ggarrange(plotlist = plotList_contrasts,
                                                    ncol = ncol, nrow = nrow, align = "h")
    } else {

      ## Remove the individual plot title now (we'll add common contrast labels later)
      for (contrast in names(plotList_contrasts)){
        plotList_contrasts[[contrast]] <- plotList_contrasts[[contrast]] + ggplot2::labs(title = NULL)
      }

      if (invert) {
        ## Signatures in columns: arrange contrast plots vertically
        arranged <- ggpubr::ggarrange(plotlist = plotList_contrasts,
                                      ncol = 1, nrow = length(plotList_contrasts),  align = "v")
        ## Annotate signature on top with wrapped title
        arranged <- ggpubr::annotate_figure(arranged,
                                            top = grid::textGrob(wrap_title(sig, width = widthlabs),
                                                                 gp = grid::gpar(cex = 1.3, fontsize = labsize)))
      } else {
        ## Signatures in rows: arrange contrast plots horizontally
        arranged <- ggpubr::ggarrange(plotlist = plotList_contrasts,
                                      ncol = length(plotList_contrasts), nrow = 1,  align = "h")
        ## Annotate signature on left with wrapped title and rotated text
        arranged <- ggpubr::annotate_figure(arranged,
                                            left = grid::textGrob(wrap_title(sig, width = widthlabs),
                                                                  rot = 90, vjust = 0.5,
                                                                  gp = grid::gpar(cex = 1.3, fontsize = labsize)))
      }
      plotList_signatures[[sig]] <- arranged
    }

  } # end signature loop

  ## Combine signature-level plots
  if (length(plotList_signatures) == 1) {
    combined_plot <- plotList_signatures[[1]]
    if (!is.null(title))
      combined_plot <- ggpubr::annotate_figure(combined_plot,
                                               top = grid::textGrob(wrap_title(title, width = widthlabs),
                                                                    gp = grid::gpar(cex = 1.3, fontsize = labsize + 4,  fontface = "bold")))
  } else {
    if (invert) {
      ## Signatures as columns
      combined_plot <- ggpubr::ggarrange(plotlist = plotList_signatures,
                                         ncol = length(plotList_signatures), nrow = 1,  align = "v")
      ## Add a left annotation for contrast names using wrap_title()
      contrast_names <- names(DEResultsList)
      contrast_grobs <- lapply(contrast_names, function(cn) {
        grid::textGrob(wrap_title(cn, width = widthlabs), rot = 90, vjust = 0.5, gp = grid::gpar(cex = 1.5, fontsize = labsize))
      })
      left_annotation <- gridExtra::arrangeGrob(grobs = contrast_grobs, ncol = 1)
      combined_plot <- ggpubr::annotate_figure(combined_plot, left = left_annotation)
    } else {
      ## Signatures as rows
      combined_plot <- ggpubr::ggarrange(plotlist = plotList_signatures,
                                         ncol = 1, nrow = length(plotList_signatures), align = "h")
      ## Add a top annotation for contrast names using wrap_title()
      contrast_names <- names(DEResultsList)
      contrast_grobs <- lapply(contrast_names, function(cn) {
        grid::textGrob(wrap_title(cn, width = widthlabs), gp = grid::gpar(cex = 1.3, fontsize = labsize))
      })
      top_annotation <- gridExtra::arrangeGrob(grobs = contrast_grobs, ncol = length(contrast_grobs))
      combined_plot <- ggpubr::annotate_figure(combined_plot, top = top_annotation)
    }

    if (!is.null(title))
      combined_plot <- ggpubr::annotate_figure(combined_plot,
                                               top = grid::textGrob(wrap_title(title, width = widthlabs),
                                                                    gp = grid::gpar(cex = 1.3, fontsize = labsize + 4,  fontface = "bold")))
  }

  return(combined_plot)
}
