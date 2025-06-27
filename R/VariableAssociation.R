#' Identify Variable Types
#'
#' Determines the type of each variable in a given data frame.
#' Variables are classified as "Numeric", "Categorical Bin" (binary categorical),
#' or "Categorical Multi" (multi-level categorical). Warnings are issued if
#' categorical variables have more than 10 unique values.
#'
#' @param df A data frame containing the variables to classify.
#' @param cols A character vector of column names to consider.
#'
#' @return A named character vector where names correspond to column names
#' and values indicate the variable type: "Numeric", "Categorical Bin", or "Categorical Multi".
#'
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   age = c(25, 30, 35, 40),
#'   gender = c("Male", "Female", "Female", "Male"),
#'   score = c(80, 85, 90, 95)
#' )
#' identify_variable_type(df)
#'}
#' @keywords internal
identify_variable_type <- function(df, cols = NULL) {

  # Define only cols of interest
  #if (!is.null(cols)) df <- df[, cols, drop = FALSE]

  if (is.null(cols)) return("Unknown")

  if (!is.null(cols)) df <- df[, cols, drop = FALSE]

  variable_types <- sapply(names(df), function(col_name) {

    col <- df[[col_name]]
    unique_vals <- length(unique(col))

    if (is.numeric(col) | is.integer(col)) {
 #     if (unique_vals > 10) {
        return("Numeric")
      # } else if (unique_vals == 2) {
      #   return("Categorical Bin")
      # } else {
      #   return("Categorical Multi")
      # }
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
#' \dontrun{
#' df <- data.frame(
#'   score = c(80, 85, 90, 95, 100),
#'   age = c(25, 30, 35, 40, 45),
#'   gender = c("Male", "Female", "Male", "Female", "Male"),
#'   group = factor(c("A", "B", "A", "B", "C"))
#' )
#'
#' results <- compute_stat_tests(df, target_var = "score")
#' print(results)
#'}
#' @importFrom stats cor.test t.test wilcox.test aov TukeyHSD
#'
#' @keywords internal
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






#' Variable Association Analysis
#'
#' This unified function evaluates associations between gene expression and sample metadata
#' using multiple methods: score-based (logmedian, ssGSEA, ranking) or GSEA-based association.
#' The function returns statistical results and visualizations summarizing effect sizes and significance.
#'
#' @param method Character string specifying the method to use. One of:
#'   - `"logmedian"`
#'   - `"ssGSEA"`
#'   - `"ranking"`
#'   - `"GSEA"`
#'
#' @section Shared Arguments (All Methods):
#' @param data A data frame with gene expression data (genes as rows, samples as columns).
#' @param metadata A data frame containing sample metadata; the first column should be the sampleID.
#' @param cols Character vector of metadata column names to analyze.
#' @param gene_set A named list of gene sets:
#'   - For score-based methods: list of gene vectors.
#'   - For GSEA: list of vectors (unidirectional) or data frames (bidirectional).
#' @param mode Contrast mode: `"simple"` (default), `"medium"`, or `"extensive"`.
#' @param signif_color Color used for significant associations (default: `"red"`).
#' @param nonsignif_color Color used for non-significant associations (default: `"grey"`).
#' @param sig_threshold Numeric significance cutoff (default: `0.05`).
#' @param saturation_value Lower limit for p-value coloring (default: auto).
#' @param widthlabels Integer for contrast label width before wrapping (default: `18`).
#' @param labsize Axis text size (default: `10`).
#' @param titlesize Plot title size (default: `14`).
#' @param pointSize Size of plot points (default: `5`).
#' @param printplt Logical. If `TRUE`, plots are printed. Default: `TRUE`.
#' @param discrete_colors (Score-based only) Optional named list mapping factor levels to colors.
#' @param continuous_color (Score-based only) Color for continuous variable points (default: `"#8C6D03"`).
#' @param color_palette (Score-based only) ColorBrewer palette name for categorical variables (default: `"Set2"`).
#' @param stat (GSEA only) Optional. Statistic for ranking genes (`"B"` or `"t"`). Auto-detected if `NULL`.
#' @param ignore_NAs (GSEA only) Logical. If `TRUE`, rows with NA metadata are removed. Default: `FALSE`.
#'
#' @return A list with method-specific results and ggplot2-based visualizations:
#'
#' **For score-based methods (`logmedian`, `ssGSEA`, `ranking`):**
#'
#' - `Overall`: Data frame of effect sizes (Cohen's f) and p-values for each metadata variable.
#' - `Contrasts`: Data frame of Cohen’s d values and adjusted p-values for pairwise comparisons (based on `mode`).
#' - `plot`: A combined visualization including:
#'     - Lollipop plots of Cohen’s f,
#'     - Distribution plots by variable (density or scatter),
#'     - Lollipop plots of Cohen’s d for contrasts.
#' - `plot_contrasts`: Lollipop plots of Cohen’s d effect sizes, colored by adjusted p-values (BH).
#' - `plot_overall`: Lollipop plot of Cohen’s f, colored by p-values.
#' - `plot_distributions`: List of distribution plots of scores by variable.
#'
#' **For GSEA-based method (`GSEA`):**
#'
#' - `data`: A data frame with GSEA results, including normalized enrichment scores (NES), adjusted p-values, and contrasts.
#' - `plot`: A ggplot2 lollipop plot of GSEA enrichment across contrasts.
#'
#' @export
VariableAssociation <- function(method = c("ssGSEA", "logmedian", "ranking", "GSEA"),
                                data, metadata, cols, gene_set, mode = c("simple", "medium", "extensive"),
                                stat = NULL, ignore_NAs = FALSE,
                                signif_color = "red", nonsignif_color = "grey", sig_threshold = 0.05,
                                saturation_value = NULL, widthlabels = 18, labsize = 10, titlesize = 14,
                                pointSize = 5,
                                discrete_colors = NULL, continuous_color = "#8C6D03", color_palette = "Set2",
                                printplt = TRUE) {
  method <- match.arg(method)
  mode <- match.arg(mode)

  if (method == "GSEA") {
    result <- GSEA_VariableAssociation(
      data = data,
      metadata = metadata,
      cols = cols,
      stat = stat,
      mode = mode,
      gene_set = gene_set,
      signif_color = signif_color,
      nonsignif_color = nonsignif_color,
      sig_threshold = sig_threshold,
      saturation_value = saturation_value,
      widthlabels = widthlabels,
      labsize = labsize,
      titlesize = titlesize,
      pointSize = pointSize,
      ignore_NAs = ignore_NAs
    )

  } else if (method %in% c("ssGSEA", "logmedian", "ranking")) {
    result <- Score_VariableAssociation(
      data = data,
      metadata = metadata,
      cols = cols,
      method = method,
      gene_set = gene_set,
      mode = mode,
      signif_color = signif_color,
      nonsignif_color = nonsignif_color,
      sig_threshold = sig_threshold,
      saturation_value = saturation_value,
      widthlabels = widthlabels,
      labsize = labsize,
      titlesize = titlesize,
      pointSize = pointSize,
      discrete_colors = discrete_colors,
      continuous_color = continuous_color,
      color_palette = color_palette,
      printplt = printplt
    )
  }

  return(result)
}
