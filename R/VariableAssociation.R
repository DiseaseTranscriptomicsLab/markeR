#' Identify Variable Types
#'
#' Determines the type of each variable in a given data frame.
#' Variables are classified as "Numeric", "Categorical Bin" (binary categorical),
#' or "Categorical Multi" (multi-level categorical). Warnings are issued if
#' categorical variables have more than 10 unique values.
#'
#' @param df A data frame containing the variables to classify.
#' @param cols Optional. A character vector of column names to consider. If NULL, all columns in `df` are used.
#'
#' @return A named character vector where names correspond to column names
#' and values indicate the variable type: "Numeric", "Categorical Bin", or "Categorical Multi".
#'
#' @examples
#' df <- data.frame(
#'   age = c(25, 30, 35, 40),
#'   gender = c("Male", "Female", "Female", "Male"),
#'   score = c(80, 85, 90, 95)
#' )
#' identify_variable_type(df)
#'
#' @export
identify_variable_type <- function(df, cols = NULL) {

  # Define only cols of interest
  if (!is.null(cols)) df <- df[, cols, drop = FALSE]

  variable_types <- sapply(names(df), function(col_name) {

    col <- df[[col_name]]
    unique_vals <- length(unique(col))

    if (is.numeric(col) | is.integer(col)) {
      if (unique_vals > 10) {
        return("Numeric")
      } else if (unique_vals == 2) {
        return("Categorical Bin")
      } else {
        return("Categorical Multi")
      }
    } else if (is.character(col) || is.factor(col)) {
      if (unique_vals == 2) {
        return("Categorical Bin")
      } else if (unique_vals > 10) {
        warning(paste0("Warning: Number of unique values in '", col_name, "' is too high (>10). Consider removing this variable from the analysis."))
        return("Categorical Multi")
      } else {
        return("Categorical Multi")
      }
    }
    return("Unknown")
  }, USE.NAMES = TRUE)

  return(variable_types)
}




#' Compute Statistical Tests for Variable Associations with a Target Variable
#'
#' Performs statistical tests to assess the relationship between predictor variables
#' and a target variable, selecting appropriate methods based on variable types.
#' Returns a list of data frames containing metric values and p-values.
#'
#' ## **Variable Classification**
#' - **Numeric**: Continuous numeric or integer variables with more than 10 unique values.
#' - **Categorical Bin**: Binary categorical variables (factors, characters, or integers with exactly 2 unique values).
#' - **Categorical Multi**: Categorical variables with more than 2 unique values (up to 10 levels recommended).
#'   A warning is issued for categorical variables with more than 10 unique values.
#'
#' ## **Statistical Tests Applied**
#' - **Numeric Predictors**: Pearson, Spearman, or Kendall correlation.
#' - **Categorical Bin Predictors**: T-test or Wilcoxon rank-sum test.
#' - **Categorical Multi Predictors**: ANOVA (default) or Kruskal-Wallis test.
#'   If ANOVA is used, Tukey’s HSD post-hoc test is performed for multiple comparisons.
#'
#' @param df A data frame containing the target variable and predictors.
#' @param target_var A string specifying the dependent variable.
#' @param cols Optional. A character vector of predictor variables.
#'   If `NULL`, all variables except `target_var` are used.
#' @param numeric The correlation method for numeric predictors.
#'   Options: `"pearson"` (default), `"spearman"`, `"kendall"`.
#' @param categorical_bin The statistical test for binary categorical variables.
#'   Options: `"t.test"` (default) or `"wilcoxon"`.
#' @param categorical_multi The statistical test for multi-level categorical variables.
#'   Options: `"anova"` (default) or `"kruskal-wallis"`.
#'
#' @return A named list (one entry per variable being analysed) where each element is a data frame with:
#'   - **Metric**: The test statistic (correlation coefficient, t-statistic, ANOVA F-value, etc.).
#'   - **p-value**: The significance value of the test.
#'   - For **Categorical Multi**, multiple rows are included for pairwise comparisons (Tukey HSD results).
#'
#' @details
#' The function automatically detects variable types and applies the appropriate test.
#' If a categorical variable has more than 10 unique levels, a warning is issued.
#' If an invalid statistical test is requested, the function stops with an error message.
#'
#' @examples
#' df <- data.frame(
#'   score = c(80, 85, 90, 95, 100),
#'   age = c(25, 30, 35, 40, 45),
#'   gender = c("Male", "Female", "Male", "Female", "Male"),
#'   group = factor(c("A", "B", "A", "B", "C"))
#' )
#'
#' results <- compute_stat_tests(df, target_var = "score")
#' print(results)
#'
#' @importFrom stats cor.test t.test wilcox.test aov TukeyHSD
#'
#' @export
compute_stat_tests <- function(df, target_var, cols = NULL,
                               numeric = "pearson",
                               categorical_bin = "t.test",
                               categorical_multi = "anova") {

  # Ensure only one method is selected per variable type
  if (length(numeric) > 1 | length(categorical_bin) > 1 | length(categorical_multi) > 1) {
    stop("Error: Please select only one method per variable type.")
  }

  results <- list()
  variable_types <- identify_variable_type(df, cols = cols)

  for (var in names(variable_types)) {
    if (var == target_var) next  # Skip the target variable itself

    test_result <- NULL
    test_df <- NULL
    method_used <- NULL

    if (variable_types[var] == "Numeric") {
      test_result <- stats::cor.test(df[[target_var]], df[[var]], method = numeric)
      test_df <- data.frame(metric = test_result$estimate, p_value = test_result$p.value)
      row.names(test_df) <- numeric
      method_used <- "Correlation"

    } else if (variable_types[var] == "Categorical Bin") {
      if (categorical_bin == "t.test") {
        test_result <- stats::t.test(df[[target_var]] ~ df[[var]])
      } else if (categorical_bin == "wilcoxon") {
        test_result <- stats::wilcox.test(df[[target_var]] ~ df[[var]])
      }

      test_df <- data.frame(metric = test_result$statistic, p_value = test_result$p.value)
      row.names(test_df) <- categorical_bin
      method_used <- "Binary Comparison"

    } else if (variable_types[var] == "Categorical Multi") {
      if (categorical_multi == "anova") {
        test_result <- stats::aov(df[[target_var]] ~ df[[var]])
        anova_p <- summary(test_result)[[1]][["Pr(>F)"]][1]
        anova_m <- summary(test_result)[[1]]["F value"][1,1]
        test_df <- data.frame(metric = anova_m,
                              p_value = anova_p)
        row.names(test_df) <- "ANOVA"

        # Tukey's HSD for multiple comparisons
        tukey_result <- stats::TukeyHSD(test_result)
        tukey_df <- as.data.frame(tukey_result[[1]])[,c(1,4)] # diff & p-value
        colnames(tukey_df) <- c("metric", "p_value")
        test_df <- rbind(test_df, tukey_df)

      } else if (categorical_multi == "kruskal-wallis") {
        test_result <- kruskal.test(df[[target_var]] ~ df[[var]])
        test_df <- data.frame(metric = test_result$statistic, p_value = test_result$p.value)
        row.names(test_df) <- "Kruskal-Wallis"
      }

      method_used <- "Multi-Category Comparison"
    }


    # scientific notation
    test_df$metric <- formatC(test_df$metric, format = "e", digits = 2)
    # correct for multiple testing per variable
    test_df$p_value <- p.adjust(test_df$p_value, method = "BH")
    test_df$p_value <- formatC(test_df$p_value, format = "e", digits = 3)


    results[[var]] <-  test_df



  }

  return(results)
}



#' Visualize Statistical Associations Between Variables and a Target Score
#'
#' Creates visual representations of statistical relationships between predictor variables
#' and a target variable (`target_var`). The function generates scatter plots with regression
#' lines for numeric variables and density plots for categorical variables, incorporating
#' statistical test results from `compute_stat_tests()`.
#'
#' ## **Plot Types**
#' - **Numeric Predictors**: Scatter plots with regression lines, annotated with correlation metrics.
#' - **Categorical Predictors**: Density plots colored by factor levels, annotated with test statistics.
#'
#' ## **Variable Classification & Statistical Methods**
#' - **Numeric**: Continuous numeric or integer variables with more than 10 unique values.
#'   - Test: Pearson, Spearman, or Kendall correlation.
#' - **Categorical Bin**: Binary categorical variables (factors, characters, or integers with exactly 2 unique values).
#'   - Test: T-test or Wilcoxon rank-sum test.
#' - **Categorical Multi**: Categorical variables with more than 2 unique values (up to 10 levels recommended).
#'   - Test: ANOVA (default) or Kruskal-Wallis.
#'   - If ANOVA is used, Tukey’s HSD post-hoc test is performed.
#'   - A warning is issued if a categorical variable has more than 10 unique values.
#'
#' ## **Color Customization**
#' - **Continuous Variables**: `continuous_color` specifies the color of scatter plot points.
#' - **Categorical Variables**:
#'   - User can provide a named list (`discrete_colors`) with custom colors for factor levels.
#'   - If `discrete_colors` is `NULL`, colors are chosen from an `RColorBrewer` palette (`color_palette`).
#'
#' @param df A data frame containing the target variable and predictors.
#' @param cols A character vector of predictor variables to include in the plots.
#' @param target_var The dependent variable to be plotted.
#' @param targetvar_lab A string specifying the label for the target variable. Default is `"Score"`.
#' @param discrete_colors Optional. A named list specifying custom colors for categorical variables.
#'   Each element should be a named vector where names correspond to factor levels.
#' @param continuous_color The color for numeric variables in scatter plots. Default is `"#8C6D03"`.
#' @param color_palette A color palette from RColorBrewer for categorical variables. Default is `"Set2"`.
#' @param sizeannot The font size for p-value annotations in plots. Default is `3`.
#' @param ncol Number of columns in the arranged plot grid. If `NULL`, layout is auto-determined.
#' @param nrow Number of rows in the arranged plot grid. If `NULL`, layout is auto-determined.
#' @param numeric The correlation method for numeric predictors.
#'   Options: `"pearson"` (default), `"spearman"`, `"kendall"`.
#' @param categorical_bin The statistical test for binary categorical variables.
#'   Options: `"t.test"` (default) or `"wilcoxon"`.
#' @param categorical_multi The statistical test for multi-level categorical variables.
#'   Options: `"anova"` (default) or `"kruskal-wallis"`.
#' @param title A string specifying the main title of the grid of plots.
#' @param titlesize Numeric; font size of the main title of the grid of plots (default = `14`).
#' @param signif_color A string specifying the color for the low end of the adjusted p-value gradient until the value chosen for significance (\code{sig_threshold}). Default is `"blue"`.
#' @param nonsignif_color A string specifying the color for the middle of the adjusted p-value gradient. Default is `"white"`. Lower limit correspond to the value of \code{sig_threshold}.
#' @param sig_threshold A numeric value specifying the threshold for significance visualization in the plot. Default: `0.05`.
#' @param saturation_value A numeric value specifying the lower limit of the adjusted p-value gradient, below which the color will correspond to \code{signif_color}. Default is the results' minimum, unless that
#' value is above the sig_threshold; in that case, it is 0.001.
#' @param widthlabels An integer controlling the maximum width of contrast labels before text wrapping. Default: `18`.
#' @param pointSize Numeric. The size of points in the lollipop plot (default is 5).
#' @param widths Numerical vector of relative columns widths. Should be the same length as the number of variables provided. Default is that each plot has the same width.
#' @param heights Numerical vector of relative columns heights  Should be the same length as the number of variables provided. Default is that each plot has the same height
#'
#' @return A named list with two elements:
#' \itemize{
#'   \item \code{plot}: A `ggarrange` object displaying the arranged statistical plots in a grid format.
#'   \item \code{data}: A named list containing the statistical test results for each variable, where each entry is a data frame with metric values and p-values.
#' }
#'
#' @details
#' - Each plot is annotated with the corresponding statistical test result (e.g., correlation coefficient,
#'   ANOVA F-value, t-test statistic) and p-values.
#' - If a categorical variable has more than 10 unique levels, a warning is issued.
#' - If `discrete_colors` is provided, it overrides `color_palette` for the specified variables.
#' - If `ncol` and `nrow` are `NULL`, the function automatically determines an optimal grid layout.
#'
#' @examples
#' df <- data.frame(
#'   score = c(80, 85, 90, 95),
#'   age = c(25, 30, 35, 40),
#'   gender = c("Male", "Female", "Female", "Male")
#' )
#' plot_stat_tests(df, cols = c("age", "gender"), target_var = "score")
#'
#' @import ggplot2 patchwork
#' @importFrom ggpubr ggarrange annotate_figure
#' @importFrom RColorBrewer brewer.pal
#'
#' @export
VariableAssociation <- function(df, cols, target_var, targetvar_lab="Score",
                                discrete_colors = NULL, continuous_color = "#8C6D03",
                                color_palette = "Set2",
                                sizeannot=3, ncol=NULL, nrow=NULL,
                                numeric = "pearson",
                                categorical_bin = "t.test",
                                categorical_multi = "anova", title=NULL, titlesize=14,
                                nonsignif_color = "grey", signif_color = "red", saturation_value=NULL,sig_threshold = 0.05, widthlabels=18, pointSize=5,
                                widths=1,heights=1) {

  results <- compute_stat_tests(df, target_var, cols = cols,
                                numeric =numeric,
                                categorical_bin = categorical_bin,
                                categorical_multi = categorical_multi)

  plot_list <- list()  # Store individual plots
  variable_types <- identify_variable_type(df, cols = names(results))

  for (var in names(results)) {
    p <- NULL
    test_results <- results[[var]]
    test_results$metric <- as.numeric(test_results$metric)
    test_results$p_value <- as.numeric(test_results$p_value)
    test_results$contrast <- row.names(test_results)
    test_results$contrast <- sapply(test_results$contrast, function(x) wrap_title(x, widthlabels))

    # Adding p-value in parenthesis next to the metric
    if (variable_types[var] == "Numeric") {
      # Create a label for the metric with p-value in parentheses (e.g., "pearson: 4.56e-02 (p = 0.03)")
      metric_label <- paste0(row.names(test_results), ": ", test_results$metric,
                             " (p = ", test_results$p_value, ")")

      p <- ggplot2::ggplot(df, ggplot2::aes_string(x = var, y = target_var)) +
        ggplot2::geom_point(alpha = 0.6, size=4, color = continuous_color) +  # Use continuous_color
        ggplot2::geom_smooth(method = "lm", col = "black", se = FALSE, size=2) +
        # Place the label in the top-left corner using -Inf/Inf coordinates
        ggplot2::annotate("text",
                          x = -Inf, y = Inf,
                          label = metric_label,
                          hjust = -0.1, vjust = 1.1,
                          size = sizeannot, color = "black") +
        ggplot2::coord_cartesian(clip = "off") +  # Allow text outside the plot area
        ggplot2::theme_classic() +
        ggplot2::ggtitle(var) +
        ggplot2::labs(y=targetvar_lab) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

    } else if (variable_types[var] %in% c("Categorical Bin")) {

      # Combine multiple test results into one multiline label if needed
      metric_labels <- paste(row.names(test_results), ": ", test_results$metric,
                             " (p = ", test_results$p_value, ")", collapse = "\n")

      # Check if user provided a custom named list for discrete colors
      if (!is.null(discrete_colors) && var %in% names(discrete_colors)) {
        # Use user-specified colors for the specific variable
        colors <- discrete_colors[[var]]
      } else {
        num_levels <- length(unique(df[[var]]))
        colors <- colorRampPalette(RColorBrewer::brewer.pal(8, color_palette))(num_levels)

      }

      p <-ggplot2:: ggplot(df, ggplot2::aes_string(x = target_var, fill = var)) +
        ggplot2::geom_density(alpha = 0.6) +
        ggplot2::scale_fill_manual(values = colors) +  # Apply custom discrete colors
        # Place the label in the top-left corner
        ggplot2::annotate("text",
                          x = -Inf, y = Inf,
                          label = metric_labels,
                          hjust = -0.1, vjust = 1.1,
                          size = sizeannot, color = "black" ) +
        ggplot2::coord_cartesian(clip = "off") +
        ggplot2::theme_classic() +
        ggplot2::ggtitle(var) +
        ggplot2::labs(x = targetvar_lab, y = "Density", fill="") +
        ggplot2::theme(legend.position = "right",
                       plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

    } else if (variable_types[var] %in% c("Categorical Multi")){

      # Check if user provided a custom named list for discrete colors
      if (!is.null(discrete_colors) && var %in% names(discrete_colors)) {
        # Use user-specified colors for the specific variable
        colors <- discrete_colors[[var]]
      } else {
        num_levels <- length(unique(df[[var]]))
        colors <- colorRampPalette(RColorBrewer::brewer.pal(8, color_palette))(num_levels)

      }

      density <- ggplot2::ggplot(df, ggplot2::aes_string(y = target_var, fill = var)) +
        ggplot2::geom_density(alpha = 0.6) +
        ggplot2::scale_fill_manual(values = colors) +
        #ggplot2::coord_cartesian(clip = "off") +
        ggplot2::theme_classic() +
        ggplot2::labs(y = targetvar_lab, x = "", fill = "") +
        ggplot2::theme(
          legend.position = "right",
          axis.title.x = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank(),
          axis.line.x = ggplot2::element_blank()
        )



      if(is.null(saturation_value)){
        if (min(test_results$p_value)>sig_threshold){
          limit_pval <- 0.001
        } else{
          limit_pval <- min(test_results$p_value)
        }

      } else {
        limit_pval <- saturation_value
      }


     lolli <-  ggplot2::ggplot(test_results, ggplot2::aes(x = metric, y = contrast, fill = -log10(p_value))) +
        ggplot2::geom_segment(ggplot2::aes(
          yend = contrast,
          xend = 0
        ), size = .5) +
        ggplot2::geom_point(ggplot2::aes(
          stroke = 1.2
        ), shape = 21, size = pointSize) +
        ggplot2::scale_linetype_identity() +
        ggplot2::scale_color_identity() +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5, face = "bold" ),
          legend.position = "right"
        ) +
        ggplot2::scale_fill_gradient2(low = nonsignif_color,
                                      mid = nonsignif_color,
                                      high = signif_color,
                                      midpoint = -log10(sig_threshold),
                                      limits=c(0,-log10(limit_pval)),
                                      na.value = signif_color)+
       ggplot2::ggtitle(var) +
       ylab("Contrasts") + xlab("Metrics")


     lolli_adj <- lolli + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
     density_adj <- density + theme(plot.margin = unit(c(0,0,0,0), "cm"))

     p <- (lolli_adj | density_adj) +
       patchwork::plot_layout(guides = "collect", widths = c(0.8, 0.2)) &
       ggplot2::theme(legend.position = "right")

    }

    if (!is.null(p)) {
      plot_list[[var]] <- p
    }


  }

  n <- length(plot_list)
  if (is.null(ncol) && is.null(nrow)) {
    ncol <- ceiling(sqrt(n))
    nrow <- ceiling(n / ncol)
  } else if (is.null(ncol)) {
    ncol <- ceiling(n / nrow)
  } else if (is.null(nrow)) {
    nrow <- ceiling(n / ncol)
  }

  # Arrange the individual plots in a grid.
  plt <- ggpubr::ggarrange(plotlist = plot_list,
                   ncol = ncol,
                   nrow = nrow,
                   widths=widths,
                   heights=heights)

  if (!is.null(title)) plt <- ggpubr::annotate_figure(plt, top = grid::textGrob(title, gp = grid::gpar(cex = 1.3, fontsize = titlesize)))

  print(plt)

  invisible(list(plot=plt,
                 data=results))
}

