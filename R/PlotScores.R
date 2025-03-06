#' Plot Gene Signature Scores as Violin Plots
#'
#' This function generates violin plots (with overlaid jittered points, median summary crossbars,
#' and optional group connection lines) for each sample based on one or more predefined gene sets
#' (signatures). Four methods are available:
#'
#'   - **ssGSEA**: Uses the single-sample Gene Set Enrichment Analysis (ssGSEA) method to compute an enrichment score
#'     for each signature in each sample using the `gsva()` function from the `GSVA` package.
#'   - **logmedian**: Computes the score as the sum of the normalized (log2-median-centered) expression values of the
#'     signature genes divided by the number of genes in the signature.
#'   - **ranking**: Computes gene signature scores for each sample by ranking the expression of signature genes
#'     in the dataset and normalizing the score based on the total number of genes.
#'   - **all**: Computes gene signature scores using all three methods (`ssGSEA`, `logmedian`, and `ranking`).
#'     Returns a heatmap summarizing Cohen's D for all metric combinations of the variables of interest.
#'
#' @param data A data frame of normalized (non-transformed) counts where each row is a gene and each column is a sample.
#'   Row names should contain gene names, and column names should contain sample identifiers. **(Required)**
#' @param metadata A data frame describing the attributes of each sample, where each row corresponds to a sample and each column to an attribute.
#'   The first column should contain sample identifiers (i.e., the column names of `data`). **(Required if method = "all")**
#' @param gene_sets Gene set input. **(Required)**
#'   - **Unidirectional gene sets**: Provide a named list where each element is a vector of gene names representing a gene signature.
#'   - **Bidirectional gene sets**: Provide a named list where each element is a data frame with two columns:
#'       - The **first column** contains gene names.
#'       - The **second column** indicates the expected direction of enrichment (1 for upregulated genes, -1 for downregulated genes).
#' @param method A character string indicating the scoring method to use. Options are `"ssGSEA"`,
#'   `"logmedian"`, `"ranking"`, or `"all"` (to compute scores using all methods). Defaults to `"logmedian"`.
#' @param ColorVariable A character string indicating the column name used to color the points.
#'   If `NULL` (default), the "Paired" brewer palette is applied. Only applicable if `method != "all"`.
#' @param GroupingVariable A character string indicating the column name in `metadata` used for grouping
#'   on the x-axis. **(Required)** If `method == "all"`, this variable is used for generating comparisons.
#' @param ColorValues An optional named vector mapping unique values of `ColorVariable` to specific colors.
#'   If not provided, the "Paired" brewer palette is used. Only applicable if `method != "all"`.
#' @param ConnectGroups Logical, indicating whether to connect groups using lines across the x-axis.
#'   If `TRUE`, a line connecting median values across groups is drawn, colored by `ColorVariable`.
#'   Default is `FALSE`. Only applicable if `method != "all"`.
#' @param ncol Optional numeric value specifying the number of columns in the grid layout for the combined plots.
#'   If `NULL`, a near-square grid is computed.
#' @param nrow Optional numeric value specifying the number of rows in the grid layout. If `NULL`, it is computed based on `ncol`.
#' @param title A string specifying the main title of the grid of plots.
#' @param titlesize Numeric; font size of the main title of the grid of plots (default = `14`).
#' @param widthTitle Optional integer specifying the maximum width of the title before inserting line breaks.
#'   Titles break at `_`, `-`, or `:` where possible, or at the exact width if no such character is found.
#'   Default is `10`.
#' @param limits Optional numeric vector of length 2 specifying the y-axis limits (if `method != "all"`) or the limits of the color scale (if `method == "all"`).
#' If `NULL`, the limits adjusts automatically.
#' @param legend_nrow Optional numeric value specifying the number of rows in the legend. If `NULL`, determined by ggplot2.
#'   Only applicable if `method != "all"`.
#' @param pointSize Optional numeric value specifying the point size. Default is `2`.
#'   Only applicable if `method != "all"`.
#' @param xlab Optional character string specifying the x-axis label. Default is the name of `GroupingVariable`.
#'   Only applicable if `method != "all"`.
#' @param cond_cohend Optional named list specifying two groups for which Cohen's d effect size should be calculated.
#'   The list should contain exactly two named elements (e.g., `list("GroupA" = c("Condition1"), "GroupB" = c("Condition2", "Condition3"))`).
#'   If not provided, Cohen's d is not computed. Currently works only for two groups. Only applicable if `method != "all"`.
#' @param labsize Numeric; font size of the plot's axis labels (default = `14`). Only applicable if `method != "all"`.
#' @param cluster_columns Logical; whether to cluster columns in the heatmap (only used if `method == "all"`). Default is `TRUE`.
#' @param cluster_rows Logical; whether to cluster rows in the heatmap (only used if `method == "all"`). Default is `TRUE`.
#' @param orientation A character string indicating heatmap orientation ("grid", "horizontal", or "vertical"). Only used if `method == "all"`.
#'
#' @return
#' - If `method != "all"`: A combined ggplot object (using `ggpubr::ggarrange` and `ggpubr::annotate_figure`) displaying violin plots for each gene signature.
#' - If `method == "all"`: A heatmap summarizing Cohen's D for all metric combinations of the `GroupingVariable`.
#'
#' @details
#' For each gene signature in `ResultsList`, the function creates a violin plot using `ggplot2`.
#' The x-axis is determined by the grouping variable (e.g., `Condition`), while the y-axis shows the signature score.
#' Jittered points are overlaid and optionally colored by `ColorVariable`. A median summary is added as a crossbar.
#' If `ConnectGroups = TRUE`, median values across groups are connected with lines, colored by `ColorVariable`.
#'
#' If `method == "all"`, the function instead returns a heatmap displaying Cohen's D effect sizes for all variable combinations.
#'
#' @importFrom ggpubr ggarrange
#' @importFrom ggpubr annotate_figure
#' @import ggplot2
#' @importFrom grid textGrob
#' @importFrom grid gpar
#' @importFrom rstatix cohens_d
#'
#' @export
PlotScores <- function(data, metadata, gene_sets,
                       method = c("ssGSEA", "logmedian", "ranking"),
                       ColorVariable = NULL, GroupingVariable,
                       ColorValues = NULL, ConnectGroups = FALSE, ncol = NULL, nrow = NULL, title=NULL,
                       widthTitle = 10, titlesize=12,limits = NULL, legend_nrow = NULL, pointSize=2, xlab=NULL, labsize=10,cond_cohend=NULL,
                       cluster_columns=T, cluster_rows=T, orientation=c("grid","horizontal","vertical")) {

  if (method == "all"){

    # if user wants "all" methods, a heatmap of cohen D's is returned, for all combination of variables in GroupingVariable
    Heatmap_Final <- Heatmap_CohenD(data = data,
                   metadata = metadata,
                   gene_sets = gene_sets,
                   variable=GroupingVariable,
                   orientation = orientation,
                   limits = c(0,2),
                   cluster_columns = cluster_columns,
                   cluster_rows = cluster_rows,
                   widthTitle = 20)
    return(Heatmap_Final)

  } else {

    ResultsList <- CalculateScores(data=data,
                                   metadata=metadata,
                                   gene_sets=gene_sets,
                                   method = method)

    # Initialize an empty list to store individual ggplot objects.
    plot_list <- list()

    # Loop over each gene signature in the ResultsList.
    for (signature in names(ResultsList)) {
      # Extract the data frame for the current signature.
      df <- ResultsList[[signature]]

      # using factors so we can retrieve the first condition for cohens d if none is specified
      df[,GroupingVariable] <- factor(df[,GroupingVariable], levels = unique(df[,GroupingVariable]))

      # Wrap the signature name using the helper function
      wrapped_title <- wrap_title(signature, width = widthTitle)

      # Create a base ggplot object with the specified grouping on the x-axis and score on the y-axis.
      p <- ggplot2::ggplot(df, ggplot2::aes_string(x = GroupingVariable, y = "score"))

      # Add jittered points, optionally colored by ColorVariable.Default: Brewer Pallette "Paired"
      if (!is.null(ColorVariable)) {
        p <- p + ggplot2::geom_jitter(ggplot2::aes_string(color = ColorVariable), size = pointSize, alpha = 0.5)
      } else {
        p <- p + ggplot2::geom_jitter(size = pointSize, alpha = 0.5) + ggplot2::scale_color_brewer(palette = "Paired")
      }

      # Overlay violin plots.
      p <- p + ggplot2::geom_violin(alpha = 0.5, scale = "width")

      # Add median summary crossbar.
      p <- p + ggplot2::stat_summary(fun = median, fun.min = median, fun.max = median,
                                     geom = "crossbar", width = 0.25,
                                     position = ggplot2::position_dodge(width = 0.13))

      # add stats
      # Compute Cohen's d

      if (!is.null(cond_cohend)){

        if (sum(unlist(cond_cohend) %in% unique(df[,GroupingVariable])) != length(unique(df[,GroupingVariable]))) warning("Warning: Not all conditions of GroupingVariable  were specified for Cohen's d calculation")

        df$cohen <- ifelse(df[,GroupingVariable] %in% cond_cohend[[1]], names(cond_cohend)[1],  names(cond_cohend)[2]  )

        cohen_d_results <- rstatix::cohens_d(df, formula = score ~ cohen)

        subtitle <- wrap_title(paste0("Cohen's d = ", round(cohen_d_results$effsize,3)), width = widthTitle)
      } else {

        # cond_cohend <- list(A=unique(df[,GroupingVariable])[1], # if no variable is defined, will be the first that appears in the ggplot
        #                     B=unique(df[,GroupingVariable])[-1])
        subtitle <- NULL
      }




      # If ConnectGroups is TRUE, add a line connecting medians across groups
      if (ConnectGroups && !is.null(ColorVariable)) {
        p <- p + ggplot2::stat_summary(ggplot2::aes_string(group = ColorVariable, color = ColorVariable),
                                       fun.y = median, geom = "line", size = 1.5, alpha=0.75)
      }

      # Customize the plot appearance.
      p <- p + ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                       plot.title = ggplot2::element_text(hjust = 0.5, size = 8),
                       plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 8, face = "italic")) +
        ggplot2::labs(title = wrapped_title, subtitle=subtitle,color = "", x = "", y = "")

      # If ColorValues is provided, use a manual color scale; otherwise, if ColorVariable is provided,
      # use a default brewer palette.
      if (!is.null(ColorValues)) {
        p <- p + ggplot2::scale_color_manual(values = ColorValues)
      } else if (!is.null(ColorVariable)) {
        p <- p + ggplot2::scale_color_brewer(palette = "Paired")
      }

      # If limits is specified, crop the plot without adjusting the data (violins).
      if (!is.null(limits)) {
        p <- p + ggplot2::coord_cartesian(ylim = limits)
      }

      # Adjust legend rows if legend_nrow is specified
      if (!is.null(legend_nrow)) {
        p <- p + ggplot2::guides(color = ggplot2::guide_legend(nrow = legend_nrow))
      }

      # Store the plot in the list.
      plot_list[[signature]] <- p
    }


    n <- length(plot_list)

    # Determine grid layout
    if (is.null(ncol) && is.null(nrow)) {

      ncol <- ceiling(sqrt(n))
      nrow <- ceiling(n / ncol)

    } else if (is.null(ncol)){

      ncol <- ceiling(n / nrow)

    } else if (is.null(nrow)){

      nrow <- ceiling(n / ncol)

    }


    # Combine plots
    combined_plot <- ggpubr::ggarrange(plotlist = plot_list, ncol = ncol, nrow = nrow, common.legend = TRUE, align = "h")

    # Annotate with axis labels
    # Change labels
    if(is.null(xlab)){
      xlab <- GroupingVariable
    }

    if (!is.null(title)) title <- wrap_title(title, width = widthTitle)

    # create label for y axis
    if (method == "ssGSEA"){
      ylab <- "ssGSEA Enrichment Score"
    } else if (method == "logmedian"){
      ylab <- "Normalized Signature Score"
    } else if (method == "ranking"){
      ylab <- "Signature Genes' Ranking"
    }

    combined_plot <- ggpubr::annotate_figure(combined_plot,
                                             left = grid::textGrob(ylab,
                                                                   rot = 90, vjust = 1, gp = grid::gpar(cex = 1.3, fontsize = labsize)),
                                             bottom = grid::textGrob(xlab, gp = grid::gpar(cex = 1.3, fontsize = labsize)),
                                             top = grid::textGrob(title, gp = grid::gpar(cex = 1.3, fontsize = titlesize)))

    return(combined_plot)

  }
}
