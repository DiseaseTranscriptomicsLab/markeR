#' Plot Gene Signature Scores as Violin Plots
#'
#' This function generates violin plots (with overlaid jittered points, median summary crossbars,
#' and optional group connection lines) for each gene signature provided in \code{ResultsList}.
#' Each element in \code{ResultsList} should be a data frame with at least the following columns:
#' \code{sample} (sample identifier) and \code{score} (the calculated gene signature score).
#' Any additional metadata columns will be added to the plots if available.
#'
#' @param ResultsList A named list where each element is a data frame containing the calculated scores for a gene signature.
#'   Each data frame must have a column named \code{sample} (matching the sample identifiers from the expression data)
#'   and a column named \code{score}. Additional columns (e.g., metadata such as \code{Condition}, \code{Genotype}, etc.)
#'   can be included, by running \code{CalculateScores} accordingly. **(Required)**
#' @param ColorVariable A character string indicating the column name used to color the points.
#'   If \code{NULL} (default), the "Paired" brewer pallette is applied.
#' @param GroupingVariable A character string indicating the column name in each data frame that will be used for grouping
#'   on the x-axis. **(Required)**
#' @param method A character string indicating the scoring method used when running \code{CalculateScores}. Options are \code{"ssGSEA"},
#'   \code{"logmedian"} or \code{"ranking"}. **This argument is mandatory.**
#' @param ColorValues An optional named vector that maps the unique values of \code{ColorVariable} to specific colors.
#'   If not provided, a default brewer palette ("Paired") will be used.
#' @param ConnectGroups A logical value indicating whether to connect groups using lines across the x-axis.
#'   If \code{TRUE}, a line connecting the median values across groups will be drawn, colored by \code{ColorVariable}.
#'   Default is \code{FALSE}. **(Optional)**
#' @param ncol An optional numeric value specifying the number of columns in the grid layout for the combined plots.
#'   If not provided, the function will compute a near-square grid (e.g., for 20 signatures, a 4x5 or 5x4 grid).
#' @param nrow An optional numeric value specifying the number of rows in the grid layout. If not provided,
#'   it is computed based on the number of signatures and \code{ncol}.
#' @param title A string specifying the main title of the grid of plots.
#' @param titlesize Numeric; font size of the main title of the grid of plots (default = `14`).
#' @param widthTitle An optional integer specifying the maximum width of the title before inserting line breaks.
#'   Titles will break at `_`, `-`, or `:` where possible, or at the exact width if no such character is found.
#'   Default is 10.
#' @param y_limits An optional numeric vector of length 2 specifying the limits of the y-axis. If \code{NULL}, the y-axis will adjust automatically.
#'   **(Optional)**
#' @param legend_nrow An optional numeric value specifying the number of rows in the legend. If \code{NULL}, the default number of rows is determined by ggplot2.
#'   **(Optional)**
#' @param pointSize An optional numeric value specifying the size of the points. If \code{NULL}, the default number is 2.
#'   **(Optional)**
#' @param xlab Optional parameter, specifying the name of the x label. Default is the name of the Grouping Variable
#'   **(Optional)**
#' @param cond_cohend An optional named list specifying the two groups for which Cohen's d effect size should be calculated.
#'   The list should contain exactly two named elements (e.g., \code{list("GroupA" = c("Condition1"), "GroupB" = c("Condition2", "Condition3"))}).
#'   If not provided, Cohen's d is not computed. Currently only working for two groups. **(Optional)**
#' @param labsize Numeric; font size of the plot's axis labels (default = `14`).
#'
#' @return A combined ggplot object (created using \code{ggpubr::ggarrange} and \code{ggpubr::annotate_figure})
#'   that displays a grid of violin plots. Each plot corresponds to one gene signature.
#'
#' @details
#' For each gene signature in \code{ResultsList}, the function creates a violin plot using \code{ggplot2}.
#' The x-axis is determined by the grouping variable (e.g., \code{Condition}), while the y-axis shows the signature score.
#' Jittered points are overlaid and optionally colored by \code{ColorVariable}. A median summary is added as a crossbar.
#' If \code{ConnectGroups = TRUE}, median values across groups are connected with lines, colored by \code{ColorVariable}.
#'
#' The individual plots are arranged into a grid using \code{ggpubr::ggarrange}, with a near-square layout if
#' \code{ncol} and \code{nrow} are not explicitly provided.
#'
#' @importFrom ggpubr ggarrange
#' @importFrom ggpubr annotate_figure
#' @import ggplot2
#' @importFrom grid textGrob
#' @importFrom grid gpar
#' @importFrom rstatix cohens_d
#' @export
PlotScores <- function(ResultsList, ColorVariable = NULL, GroupingVariable, method = c("ssGSEA", "logmedian", "ranking"),
                       ColorValues = NULL, ConnectGroups = FALSE, ncol = NULL, nrow = NULL, title = NULL,
                       widthTitle = 10, titlesize = 12, y_limits = NULL, legend_nrow = NULL, pointSize = 2,
                       xlab = NULL, labsize = 10, cond_cohend = NULL, pvalcalc = FALSE) {

  # Initialize an empty list to store individual ggplot objects.
  plot_list <- list()

  # Loop over each gene signature in the ResultsList.
  for (signature in names(ResultsList)) {
    # Extract the data frame for the current signature.
    df <- ResultsList[[signature]]

    # Using factors so we can retrieve the first condition for Cohen's d if none is specified.
    df[, GroupingVariable] <- factor(df[, GroupingVariable], levels = unique(df[, GroupingVariable]))

    # Wrap the signature name using the helper function.
    wrapped_title <- wrap_title(signature, width = widthTitle)

    # Create a base ggplot object with the specified grouping on the x-axis and score on the y-axis.
    p <- ggplot2::ggplot(df, ggplot2::aes_string(x = GroupingVariable, y = "score"))

    # Add jittered points, optionally colored by ColorVariable. Default: Brewer Palette "Paired"
    if (!is.null(ColorVariable)) {
      p <- p + ggplot2::geom_jitter(ggplot2::aes_string(color = ColorVariable), size = pointSize, alpha = 0.5)
    } else {
      p <- p + ggplot2::geom_jitter(size = pointSize, alpha = 0.5) +
        ggplot2::scale_color_brewer(palette = "Paired")
    }

    # Overlay violin plots.
    p <- p + ggplot2::geom_violin(alpha = 0.5, scale = "width")

    # Add median summary crossbar.
    p <- p + ggplot2::stat_summary(fun = median, fun.min = median, fun.max = median,
                                   geom = "crossbar", width = 0.25,
                                   position = ggplot2::position_dodge(width = 0.13))

    # Add stats: Compute Cohen's d (and optionally the p-value).
    if (!is.null(cond_cohend)) {

      if (sum(unlist(cond_cohend) %in% unique(df[, GroupingVariable])) != length(unique(df[, GroupingVariable]))) {
        warning("Warning: Not all conditions of GroupingVariable were specified for Cohen's d calculation")
      }

      # Create a new grouping variable for effect size calculations.
      df$cohen <- ifelse(df[, GroupingVariable] %in% cond_cohend[[1]],
                         names(cond_cohend)[1], names(cond_cohend)[2])

      # Compute Cohen's d.
      cohen_d_results <- rstatix::cohens_d(df, formula = score ~ cohen)

      # If pvalcalc is TRUE, compute the p-value via a t-test.
      if (pvalcalc) {
        ttest_results <- rstatix::t_test(df, formula = score ~ cohen)
        p_val <- ttest_results$p[1]
        line1 <- wrap_title(paste0("Cohen's d = ", round(cohen_d_results$effsize, 3)), width = widthTitle)
        line2 <- wrap_title(paste0("p = ", round(p_val, 3)), width = widthTitle)
        subtitle <- paste(line1, line2, sep = "\n")
      } else {
        subtitle <- wrap_title(paste0("Cohen's d = ", round(cohen_d_results$effsize, 3)), width = widthTitle)
      }
    } else {
      subtitle <- NULL
    }

    # If ConnectGroups is TRUE, add a line connecting medians across groups.
    if (ConnectGroups && !is.null(ColorVariable)) {
      p <- p + ggplot2::stat_summary(ggplot2::aes_string(group = ColorVariable, color = ColorVariable),
                                     fun.y = median, geom = "line", size = 1.5, alpha = 0.75)
    }

    # Customize the plot appearance.
    p <- p + ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     plot.title = ggplot2::element_text(hjust = 0.5, size = 8),
                     plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 8, face = "italic")) +
      ggplot2::labs(title = wrapped_title, subtitle = subtitle, color = "", x = "", y = "")

    # If ColorValues is provided, use a manual color scale; otherwise, if ColorVariable is provided,
    # use a default brewer palette.
    if (!is.null(ColorValues)) {
      p <- p + ggplot2::scale_color_manual(values = ColorValues)
    } else if (!is.null(ColorVariable)) {
      p <- p + ggplot2::scale_color_brewer(palette = "Paired")
    }

    # If y_limits is specified, crop the plot without adjusting the data (violins).
    if (!is.null(y_limits)) {
      p <- p + ggplot2::coord_cartesian(ylim = y_limits)
    }

    # Adjust legend rows if legend_nrow is specified.
    if (!is.null(legend_nrow)) {
      p <- p + ggplot2::guides(color = ggplot2::guide_legend(nrow = legend_nrow))
    }

    # Store the plot in the list.
    plot_list[[signature]] <- p
  }

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

  # Combine plots.
  combined_plot <- ggpubr::ggarrange(plotlist = plot_list, ncol = ncol, nrow = nrow, common.legend = TRUE, align = "h")

  # Annotate with axis labels.
  if (is.null(xlab)) {
    xlab <- GroupingVariable
  }

  if (!is.null(title)) title <- wrap_title(title, width = widthTitle)

  # Set y-axis label based on method.
  if (method == "ssGSEA") {
    ylab <- "ssGSEA Enrichment Score"
  } else if (method == "logmedian") {
    ylab <- "Normalized Signature Score"
  } else if (method == "ranking") {
    ylab <- "Signature Genes' Ranking"
  }

  combined_plot <- ggpubr::annotate_figure(combined_plot,
                                           left = grid::textGrob(ylab,
                                                                 rot = 90, vjust = 1, gp = grid::gpar(cex = 1.3, fontsize = labsize)),
                                           bottom = grid::textGrob(xlab, gp = grid::gpar(cex = 1.3, fontsize = labsize)),
                                           top = grid::textGrob(title, gp = grid::gpar(cex = 1.3, fontsize = titlesize)))

  return(combined_plot)
}

