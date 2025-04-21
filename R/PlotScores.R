#' Plot gene signature scores using various methods.
#'
#' Computes and visualizes gene signature scores using one or more methods, returning plots
#' such as scatter plots, violin plots, heatmaps, or volcano plots depending on inputs.
#'
#' @param data A data frame of normalized (non-transformed) counts where each row is a gene and each column is a sample.
#'   Row names should contain gene names, and column names should contain sample identifiers. **(Required)**
#' @param metadata A data frame with sample-level attributes. Each row corresponds to a sample, with the first column containing sample IDs
#'   that match `colnames(data)`. **Required if `method = "all"` or if metadata-derived groupings or colors are used.**
#' @param gene_sets A named list of gene sets to score.
#'   - For **unidirectional gene sets**: a list of character vectors.
#'   - For **bidirectional gene sets**: a list of data frames with two columns: gene names and direction (1 = up, -1 = down). **(Required)**
#' @param method Scoring method to use. One of `"ssGSEA"`, `"logmedian"`, `"ranking"`, or `"all"` (default = `"logmedian"`).
#'   - `"all"` triggers a full analysis returning both heatmap and volcano plots.
#'   - Other values return single-score plots depending on `Variable` type.
#' @param ColorVariable Name of a metadata column to color points by. Used in **single-method mode** (`"ssGSEA"`, etc.). Ignored in `"all"` mode.
#' @param Variable Metadata column to define groups or numeric comparisons.
#'   - **Required if `method = "all"`** (used to compute and compare effect sizes).
#'   - If `NULL` and `method != "all"`, density plots of each signature scores across samples are shown (no grouping or comparison).
#' @param ColorValues Optional. A named vector or list of colors used to control the coloring of plot elements across different methods and variable types.
#' Behavior depends on the combination of `method` and `Variable`:
#'
#'   - **If `method != "all"`**:
#'     - If `Variable` is `NULL`:
#'       - Used in density plots; a **single color** will be applied (default: `"#ECBD78"` if `ColorValues` is not specified).
#'     - If `Variable` is **categorical**:
#'       - A **named vector** should map each level of `Variable` (or `ColorVariable`) to a specific color.
#'       - Overrides the palette specified by `colorPalette`.
#'     - If `Variable` is **numeric**:
#'       - A **single color** is applied to all points in the scatter plot (default: `"#5264B6"`).
#'
#'   - **If `method == "all"`**:
#'     - `ColorValues` can be a **named list** with two elements:
#'       - `heatmap`: A vector of **two colors** used as a diverging scale for the heatmap of effect sizes (default: `c("#F9F4AE", "#B44141")`).
#'       - `volcano`: A named vector of colors used for labeling or grouping gene signatures (e.g., in the volcano plot).
#'     - If not provided, defaults will be used for both components.
#'
#'   In all cases, `ColorValues` takes precedence over the default `colorPalette` setting if specified.
#'
#' @param ConnectGroups Logical. If `TRUE`, connects points by sample ID across conditions (used for categorical variables and `method != "all"`).
#' @param ncol Number of columns for facet layout (used in both heatmaps and score plots).
#' @param nrow Number of rows for facet layout (used in both heatmaps and score plots).
#' @param title Plot title (optional).
#' @param widthTitle Width allocated for title (affects alignment).
#' @param titlesize Font size for plot title.
#' @param limits Y-axis limits (numeric vector of length 2).
#' @param legend_nrow Number of rows for plot legend (used in single-method plots).
#' @param pointSize Numeric. Size of points in **score plots** (violin or scatter),
#'   used when plotting individual sample scores for both categorical and numeric variables,
#'   including when `method = "all"`.
#' @param xlab Label for x-axis (optional; defaults to `Variable`).
#' @param labsize Font size for axis and facet labels.
#' @param compute_cohen Logical. Whether to compute Cohen's effect sizes in **score plots** (`method != "all"`).
#'   - **Only applies when `method != "all"`**; ignored otherwise.
#'   - If the variable is **categorical** and `cond_cohend` is specified, computes **Cohen's d** for the specified comparison.
#'   - If the variable is **categorical** and `cond_cohend` is not specified, computes:
#'       - **Cohen's d** if there are exactly two groups.
#'       - **Cohen's f** if there are more than two groups.
#'   - If the variable is **numeric**, computes **Cohen's f** regardless of `cond_cohend`.
#' @param cond_cohend Optional. List of length 2 with the two groups being used to compute effect size. The values in each entry should be levels of `Variable (used with `compute_cohen = TRUE`).
#' @param cohen_threshold Effect size threshold shown as a **guide line** in volcano plots. Used only when `method = "all"`.
#' @param pvalcalc Logical. If `TRUE`, computes p-values between groups.
#' @param mode A string specifying the level of detail for contrasts, if `method = "all"`.
#' Options are:
#' - `"simple"`: Pairwise comparisons (e.g., A - B).
#' - `"medium"`: Pairwise comparisons plus comparisons against the mean of other groups.
#' - `"extensive"`: All possible groupwise contrasts, ensuring balance in the number of terms on each side.
#' @param widthlegend Width of the legend in **volcano plots** (used only if `method = "all"`) and violin score plots.
#' @param sig_threshold P-value cutoff shown as a **guide line** in volcano plots. Only applies when `method = "all"`.
#' @param colorPalette Name of an RColorBrewer palette used to assign colors in plots. Applies to all methods. Default is "Set3".
#'   If `ColorValues` is provided, it overrides this palette.
#'   - If `Variable` is `NULL` and `method != "all"` (i.e., for density plots), a default color `"#ECBD78"` is used.
#'   - If `method = "all"` (i.e., for heatmaps and volcano plots), a default diverging color scale is used: `c("#F9F4AE", "#B44141")`, unless `ColorValues` is manually specified.
#' @param cor Correlation method for numeric variables. One of `"pearson"` (default), `"spearman"`, or `"kendall"`.
#'   Only applies when the variable is numeric and `method != "all"`.
#'
#' @return Depending on `method`:
#'   - If `method = "all"`: returns a list with `heatmap` and `volcano` ggplot objects.
#'   - If `method` is a single method: returns a single ggplot object (scatter or violin plot depending on variable type).
#'
#' @details
#' **Behavior Based on `method`:**
#'
#' - `"all"`:
#'   - Requires `metadata` and `Variable`.
#'   - Computes scores using all available methods and returns:
#'     - A heatmap of Cohen’s effect sizes.
#'     - A volcano plot showing effect size vs p-value across gene signatures.
#'   - Uses additional parameters:
#'     - `mode`: defines how contrasts between groups are constructed.
#'     - `sig_threshold` and `cohend_threshold`: add **guide dashed lines** to the volcano plot (do not affect point coloring).
#'     - `widthlegend`: controls width of the volcano plot legend.
#'     - `pointSize`: controls dot size for signature points in the volcano plot.
#'   - `ColorValues` can be a **named list**:
#'     - `heatmap`: two-color gradient for effect sizes (default: `c("#F9F4AE", "#B44141")`).
#'     - `signatures`: named vector of colors for gene signatures in the volcano plot (default is color palette "Set3").
#'
#' - `"ssGSEA"`, `"logmedian"`, or `"ranking"`:
#'   - The type of `Variable` determines the plot:
#'     - If **categorical**: produces violin plots with optional group comparisons.
#'     - If **numeric**: produces scatter plots with correlation.
#'     - If `Variable` is `NULL`: produces density plots for each signature across all samples.
#'   - Additional arguments:
#'     - `ColorVariable` and `ColorValues`: control coloring of points or violins.
#'     - `colorPalette`: default palette (overridden by `ColorValues` if present).
#'     - `ConnectGroups`: links samples by ID (for categorical `Variable` only).
#'     - `cor`: specifies correlation method for numeric `Variable`.
#'     - `pvalcalc`: enables group-wise p-value calculations (categorical only).
#'     - `compute_cohen`: calculates effect sizes when applicable.
#'     - `cond_cohend`: focuses Cohen’s d calculation on a specific comparison.
#'
#' **Behavior Based on `Variable` Type:**
#'
#' - **If `Variable` is numeric**:
#'   - Outputs scatter plots (in single-method mode).
#'   - Computes correlation (`cor`).
#'   - Ignores `compute_cohen`, `cond_cohend`, and `pvalcalc`.
#'   - Color is uniform (default: `"#5264B6"`) unless overridden via `ColorValues`.
#'   - Cohen’s f effect size estimation (`compute_cohen = TRUE`) and significance if `pvalcalc` is `TRUE`.
#'
#' - **If `Variable` is categorical**:
#'   - Outputs violin plots (in single-method mode).
#'   - Supports:
#'     - p-value comparisons (`pvalcalc = TRUE`),
#'     - optional connection lines (`ConnectGroups = TRUE`),
#'     - Cohen’s effect size estimation (`compute_cohen = TRUE`) and significance (`pvalcalc` is `TRUE`):
#'       - If `cond_cohend` is specified, computes **Cohen's d** for that comparison.
#'       - If not specified:
#'         - Computes **Cohen’s d** if 2 groups.
#'         - Computes **Cohen’s f** if >2 groups.
#'   - Colors are matched to factor levels using `ColorValues` or `colorPalette`.
#'
#' - **If `Variable` is NULL and `method != "all"`**:
#'   - Produces density plots of signature scores.
#'   - Uses a **single fill color** (`"#ECBD78"` by default or from `ColorValues`).
#'
#'
#' @export
PlotScores <- function(data, metadata, gene_sets,
                       method = c("ssGSEA", "logmedian", "ranking", "all"),
                       ColorVariable = NULL, Variable = NULL,
                       ColorValues = NULL, ConnectGroups = FALSE, ncol = NULL, nrow = NULL, title = NULL,
                       widthTitle = 20, titlesize = 12, limits = NULL, legend_nrow = NULL, pointSize = 4,
                       xlab = NULL, labsize = 10, compute_cohen=TRUE, cond_cohend = NULL, pvalcalc = FALSE, mode = c("simple","medium","extensive"),
                       widthlegend=22, sig_threshold=0.05, cohen_threshold=0.5, colorPalette="Set3", cor=c("pearson","spearman","kendall")) {

  method <- match.arg(method)

  type <- identify_variable_type(metadata, Variable)#[Variable]

  if (method == "all") { # returns heatmap

    if (type =="Numeric"){

      cohenlist <- CohenF_allConditions(data = data, metadata = metadata, gene_sets = gene_sets, variable = Variable )

    } else {

      cohenlist <- CohenD_allConditions(data = data, metadata = metadata, gene_sets = gene_sets, variable = Variable, mode = mode)

    }

    # if user wants "all" methods, a heatmap of Cohen's d's is returned, for all combination of variables in GroupingVariable
    Heatmap_Final <- Heatmap_Cohen(cohenlist = cohenlist,
                                   nrow = nrow,
                                   ncol = ncol,
                                   limits = limits,
                                   widthTitle = widthTitle,
                                   titlesize = titlesize,
                                   ColorValues = ColorValues,
                                   title = title )

    Volcano_Cohen <- Volcano_Cohen(cohenlist = cohenlist,
                                   titlesize = 12,
                                   ColorValues = ColorValues,
                                   title = title,
                                   widthlegend = widthlegend,
                                   pointSize = pointSize,
                                   sig_threshold = sig_threshold,
                                   cohen_threshold = cohen_threshold,
                                   colorPalette =colorPalette,
                                   ncol = ncol,
                                   nrow = nrow)

    return(list(heatmap=Heatmap_Final$plt,
                volcano=Volcano_Cohen$plt))

  } else {



    if (type!="Numeric"){

      return(

        PlotScores_Categorical(data=data, metadata=metadata, gene_sets=gene_sets,
                               method = method,
                               ColorVariable = ColorVariable, GroupingVariable = Variable,
                               ColorValues = ColorValues, ConnectGroups = ConnectGroups, ncol = ncol, nrow = nrow, title = title,
                               widthTitle = widthTitle, titlesize = titlesize, limits = limits, legend_nrow = legend_nrow, pointSize = pointSize,
                               xlab = xlab, labsize = labsize, compute_cohen=compute_cohen, cond_cohend = cond_cohend, pvalcalc = pvalcalc, mode = mode,
                               widthlegend=widthlegend, cohen_threshold=cohen_threshold, colorPalette=colorPalette)

      )

    } else {

      return(

        PlotScores_Numeric(data=data,
                           metadata=metadata,
                           gene_sets=gene_sets,
                           method = method,
                           Variable = Variable,
                           ColorValues = ColorValues,
                           ncol = ncol,
                           nrow = nrow,
                           title = title,
                           widthTitle = widthTitle,
                           titlesize = titlesize,
                           limits = limits,
                           pointSize = pointSize,
                           xlab = xlab,
                           labsize = labsize,
                           compute_cohen = compute_cohen,
                           pvalcalc = pvalcalc,
                           colorPalette = colorPalette,
                           cor=cor)

      )

    }


  }


}
#' Plot Gene Set Scores by Group or Continuous Variable
#'
#' This function computes and visualizes gene set enrichment scores using various methods, optionally comparing across groups or numeric variables.
#' It supports categorical and numeric comparisons, statistical testing, Cohen's d effect sizes, and visualizations such as heatmaps and volcano plots.
#'
#' Four methods are available:
#'
#'   - **ssGSEA**: Uses the single-sample Gene Set Enrichment Analysis (ssGSEA) method to compute an enrichment score
#'     for each signature in each sample using an adaptation of the `gsva()` function from the `GSVA` package.
#'   - **logmedian**: Computes the score as the sum of the normalized (log2-median-centered) expression values of the
#'     signature genes divided by the number of genes in the signature.
#'   - **ranking**: Computes gene signature scores for each sample by ranking the expression of signature genes
#'     in the dataset and normalizing the score based on the total number of genes.
#'   - **all**: Computes gene signature scores using all three methods (`ssGSEA`, `logmedian`, and `ranking`).
#'     Returns a heatmap summarizing Cohen's D for all metric combinations of the variables of interest.
#'
#' Depending on the method and the type of variable (categorical, numeric, or `NULL`), the function produces different plots:
#' - **If `method = "all"`** and the variable is **categorical**, a heatmap of Cohen's d or F statistics and a volcano plot showing contrasts between all groups of that variable are produced.
#' - **If `method = "all"`** and the variable is **numeric**, a heatmap of Cohen's F and a volcano plot are produced.
#' - **If `method != "all"`** and the variable is **categorical**, a violin plot for each signature is generated.
#' - **If `method != "all"`** and the variable is `NULL`, a density plot of the score distribution is displayed.
#' - **If `method != "all"`** and the variable is **numeric**, a scatter plot is created to show the relationship between the scores and the numeric variable.
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
#' @param method A character string indicating the scoring method to use. Options are `"ssGSEA"`, `"logmedian"` or `"ranking"`. Defaults to `"logmedian"`.
#' @param ColorVariable Optional. Name of the metadata column to use for point color in plots.
#' @param GroupingVariable Optional. Name of the metadata column to use for group comparison.
#' @param ColorValues Optional. Named vector of colors to use for each group in `ColorVariable` or `GroupingVariable`.
#' @param ConnectGroups Logical. If TRUE, connects points of the same sample across conditions.
#' @param ncol Number of columns in the facet layout of the plot.
#' @param nrow Number of rows in the facet layout of the plot.
#' @param title Optional. Main title of the plot.
#' @param widthTitle Numeric. Width of the title area (for alignment purposes).
#' @param titlesize Numeric. Font size of the title text.
#' @param limits Optional numeric vector of length 2 specifying y-axis limits.
#' @param legend_nrow Optional. Number of rows in the plot legend.
#' @param pointSize Numeric. Size of the points in the plots.
#' @param xlab Optional. Label for the x-axis.
#' @param labsize Numeric. Font size for axis and facet labels.
#' @param compute_cohen Logical. If TRUE, computes Cohen's d effect sizes between groups.
#' @param cond_cohend Optional. Specify a condition or comparison subset for calculating Cohen’s d.
#' @param pvalcalc Logical. If TRUE, computes p-values for group comparisons.
#' @param mode Character string indicating comparison complexity. Options: `"simple"`, `"medium"`, `"extensive"`.
#' @param widthlegend Numeric. Width of the legend area in volcano plots.
#' @param cohen_threshold Numeric. Cohen's d threshold to highlight effect size in volcano plots (default = 0.6).
#' @param colorPalette Character. Name of RColorBrewer palette for coloring (default = `"Set3"`).
#'
#' @details
#' - **If `method = "all"`** and the variable is **categorical**, the function returns a heatmap of Cohen's d or F statistics and a volcano plot showing contrasts between all groups of that variable.
#' - **If `method = "all"`** and the variable is **numeric**, a heatmap of Cohen's F and a volcano plot will be produced.
#' - **If `method != "all"`** and the variable is **categorical**, a violin plot for each signature will be displayed.
#' - **If `method != "all"`** and the variable is `NULL`, a density plot of the score distribution will be displayed.
#' - **If `method != "all"`** and the variable is **numeric**, a scatter plot will be generated to show the relationship between the scores and the numeric variable.
#'
#' @import ggplot2
#' @importFrom grid textGrob
#' @importFrom grid gpar
#' @importFrom rstatix t_test
#' @export
PlotScores_Categorical <- function(data, metadata, gene_sets,
                                   method = c("ssGSEA", "logmedian", "ranking"),
                                   ColorVariable = NULL, GroupingVariable = NULL,
                                   ColorValues = NULL, ConnectGroups = FALSE, ncol = NULL, nrow = NULL, title = NULL,
                                   widthTitle = 10, titlesize = 12, limits = NULL, legend_nrow = NULL, pointSize = 2,
                                   xlab = NULL, labsize = 10, compute_cohen=TRUE, cond_cohend = NULL, pvalcalc = FALSE, mode = c("simple","medium","extensive"),
                                   widthlegend=22, cohen_threshold=0.6, colorPalette="Set3") {

  method <- match.arg(method)

  ResultsList <- CalculateScores(data = data,
                                 metadata = metadata,
                                 gene_sets = gene_sets,
                                 method = method)

  # if grouping variable is NULL, then the function displays a density / distribution of scores
  if (is.null(GroupingVariable) | is.null(metadata)) {

    plot_list <- list()

    for (signature in names(ResultsList)) {

      df <- ResultsList[[signature]]
      # Wrap the signature name using the helper function
      wrapped_title <- wrap_title(signature, width = widthTitle)

      ColorValues <- if (is.null(ColorValues)) "#ECBD78" else ColorValues

      p <- ggplot2::ggplot(df, ggplot2::aes(x = score)) +
        ggplot2::geom_density(fill = ColorValues, alpha = 0.5) +
        ggplot2::labs(title = "Density Plot of Score", x = xlab, y = "Density")

      # Customize the plot appearance.
      p <- p + ggplot2::theme_classic() +
        ggplot2::labs(title = wrapped_title, color = "", x = "", y = "") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = labsize - .5),
                       axis.text.y = ggplot2::element_text(  size = labsize - .5),
                       plot.title = ggplot2::element_text(hjust = 0.5, size = titlesize-1),
                       plot.subtitle = ggplot2::element_text(hjust = 0.5, size = titlesize - 1.5, face = "italic"))

      # If limits is specified, crop the plot without adjusting the data (violins).
      if (!is.null(limits)) {
        p <- p + ggplot2::coord_cartesian(xlim = limits)
      }

      plot_list[[signature]] <- p

    }

    n <- length(plot_list)

    # Determine grid layout
    if (is.null(ncol) && is.null(nrow)) {
      ncol <- ceiling(sqrt(n))
      nrow <- ceiling(n / ncol)
    } else if (is.null(ncol)) {
      ncol <- ceiling(n / nrow)
    } else if (is.null(nrow)) {
      nrow <- ceiling(n / ncol)
    }

    # create label for y axis
    if (method == "ssGSEA") {
      xlab <- "ssGSEA Enrichment Score"
    } else if (method == "logmedian") {
      xlab <- "Normalized Signature Score"
    } else if (method == "ranking") {
      xlab <- "Signature Genes' Ranking"
    }

    combined_plot <- ggpubr::ggarrange(plotlist = plot_list, ncol = ncol, nrow = nrow, common.legend = TRUE, align = "h")
    combined_plot <- ggpubr::annotate_figure(combined_plot,
                                             left = grid::textGrob("Density",
                                                                   rot = 90, vjust = 1, gp = grid::gpar(cex = 1.3, fontsize = labsize)),
                                             bottom = grid::textGrob(xlab, gp = grid::gpar(cex = 1.3, fontsize = labsize)),
                                             top = grid::textGrob(title, gp = grid::gpar(cex = 1.3, fontsize = titlesize + 2)))
    return(combined_plot)
  }

  if (!(GroupingVariable %in% colnames(metadata)))
    stop(paste0(GroupingVariable, " not in metadata columns. Please check metadata."))

  # Initialize an empty list to store individual ggplot objects.
  plot_list <- list()

  # Loop over each gene signature in the ResultsList.
  for (signature in names(ResultsList)) {
    # Extract the data frame for the current signature.
    df <- ResultsList[[signature]]

    # Using factors so we can retrieve the first condition for Cohen's d if none is specified.
    df[, GroupingVariable] <- factor(df[, GroupingVariable],
                                     levels = sort(unique(as.character(df[, GroupingVariable]))))

    # Wrap the signature name using the helper function.
    wrapped_title <- wrap_title(signature, width = widthTitle)

    # Create a base ggplot object with the specified grouping on the x-axis and score on the y-axis.
    p <- ggplot2::ggplot(df, ggplot2::aes_string(x = GroupingVariable, y = "score"))

    # Add jittered points, optionally colored by ColorVariable.
    if (!is.null(ColorVariable)) {
      p <- p + ggplot2::geom_jitter(ggplot2::aes_string(color = ColorVariable), size = pointSize, alpha = 0.5)
    } else {
      p <- p + ggplot2::geom_jitter(size = pointSize, alpha = 0.5) +
        ggplot2::scale_color_brewer(palette = colorPalette)
    }

    # Overlay violin plots.
    p <- p + ggplot2::geom_violin(alpha = 0.5, scale = "width")

    # Add median summary crossbar.
    p <- p + ggplot2::stat_summary(fun = median, fun.min = median, fun.max = median,
                                   geom = "crossbar", width = 0.25,
                                   position = ggplot2::position_dodge(width = 0.13))

    # Add stats: Compute Cohen's d (and optionally p‑value)
    if(compute_cohen){
      if (!is.null(cond_cohend)){
        # can be of the following form:
        # cond_cohend <- list(A=c("Senescent"),
        #                     B=c("Proliferative","Quiescent"))

        if (sum(unlist(cond_cohend) %in% unique(df[, GroupingVariable])) != length(unique(df[, GroupingVariable])))
          warning("Warning: Not all conditions of GroupingVariable were specified for Cohen's d calculation")

        x <- df[df[[GroupingVariable]] %in% cond_cohend[[1]], "score", drop = TRUE]
        y <- df[df[[GroupingVariable]]  %in% cond_cohend[[2]], "score", drop = TRUE]

        cohen_d_results <- cohen_d(x, y)

        # df$cohen <- ifelse(df[, GroupingVariable] %in% cond_cohend[[1]], names(cond_cohend)[1], names(cond_cohend)[2])
        # cohen_d_results <- rstatix::cohens_d(df, formula = score ~ cohen)

        if (pvalcalc) {
          df$cohen <- ifelse(df[, GroupingVariable] %in% cond_cohend[[1]], names(cond_cohend)[1], names(cond_cohend)[2])
          ttest_results <- rstatix::t_test(df, formula = score ~ cohen)
          p_val <- ttest_results$p[1]
          line1 <- wrap_title(paste0("Cohen's d = ", round(cohen_d_results, 3)), width = widthTitle)
          line2 <- wrap_title(paste0("p = ", round(p_val, 3)), width = widthTitle)
          subtitle <- paste(line1, line2, sep = "\n")
        } else {
          subtitle <- wrap_title(paste0("Cohen's d = ", round(cohen_d_results, 3)), width = widthTitle)
        }


      } else {

        if(length(unique(df[, GroupingVariable])) < 2){

          warning("Not enough conditions available to report Cohen's d.")

        } else if(length(unique(df[, GroupingVariable])) == 2) {

          # Calculate Cohen's d based on ordering of the x axis
          group1 <- levels(df[, GroupingVariable])[1]
          group2 <- levels(df[, GroupingVariable])[2]

          x <- df[df[[GroupingVariable]] == group1, "score", drop = TRUE]
          y <- df[df[[GroupingVariable]] == group2, "score", drop = TRUE]

          cohen_d_results <- cohen_d(x, y)

          if (pvalcalc) {
            ttest_results <- rstatix::t_test(df, formula = score ~ GroupingVariable)
            p_val <- ttest_results$p[1]
            line1 <- wrap_title(paste0("Cohen's d = ", round(cohen_d_results$effsize, 3)), width = widthTitle)
            line2 <- wrap_title(paste0("p = ", round(p_val, 3)), width = widthTitle)
            subtitle <- paste(line1, line2, sep = "\n")
          } else {
            subtitle <- wrap_title(paste0("Cohen's d = ", round(cohen_d_results$effsize, 3)), width = widthTitle)
          }



        } else if(length(unique(df[, GroupingVariable])) > 2){

          # Calculate Cohen's f
          type <- identify_variable_type(df, GroupingVariable)[GroupingVariable]
          #Without scaling, the coefficient represents the change in score per unit increase in the variable (if numeric, the unit of the variable. Makes sense to not scale...)
          model <- lm(score ~ get(GroupingVariable), data = df)
          results_var <- compute_cohens_f_pval(model, type)


          if (pvalcalc) {
            line1 <- wrap_title(paste0("Cohen's f = ", round(results_var["Cohen_f"], 3)), width = widthTitle)
            line2 <- wrap_title(paste0("p = ", round(results_var["P_Value"], 3)), width = widthTitle)
            subtitle <- paste(line1, line2, sep = "\n")
          } else {
            subtitle <- wrap_title(paste0("Cohen's f = ", round(results_var["Cohen_f"], 3)), width = widthTitle)
          }

        }


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
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = labsize),
                     axis.text.y = ggplot2::element_text(  size = labsize),
                     plot.title = ggplot2::element_text(hjust = 0.5, size = titlesize-1),
                     plot.subtitle = ggplot2::element_text(hjust = 0.5, size = titlesize - 1.5, face = "italic")) +
      ggplot2::labs(title = wrapped_title, subtitle = subtitle, color = "", x = "", y = "")

    # If ColorValues is provided, use a manual color scale; otherwise, if ColorVariable is provided,
    # use a default brewer palette.
    if (!is.null(ColorValues)) {
      p <- p + ggplot2::scale_color_manual(values = ColorValues)
    } else if (!is.null(ColorVariable)) {
      p <- p + ggplot2::scale_color_brewer(palette = colorPalette)
    }

    # If limits is specified, crop the plot without adjusting the data (violins).
    if (!is.null(limits)) {
      p <- p + ggplot2::coord_cartesian(ylim = limits)
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

  # Create label for y axis based on method.
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


#' Plot Gene Signature Scores with Continuous Variables
#'
#' This function visualizes gene signature scores using scatter plots and regression lines across a continuous metadata variable.
#' Signature scores are computed per sample using one of three methods: \code{"ssGSEA"}, \code{"logmedian"}, or \code{"ranking"}.
#' Optionally, the effect size (Cohen's f) and p-value for the association between the signature score and the continuous variable can be computed and displayed.
#'
#' @param data A data frame of normalized (non-transformed) gene expression counts. Rows are genes, columns are samples. Row names should be gene names, and column names should match sample identifiers in \code{metadata}.
#' @param metadata A data frame where each row corresponds to a sample and contains sample-level attributes (e.g., clinical or experimental metadata). Must include a column matching the sample IDs in \code{data}.
#' @param gene_sets A list of gene sets (signatures). Each element is either a character vector of gene names or a data frame with gene names and enrichment direction (1 for upregulated, -1 for downregulated).
#' @param method Scoring method to use. One of \code{"ssGSEA"}, \code{"logmedian"}, or \code{"ranking"}. Default is \code{"logmedian"}.
#' @param Variable Name of the continuous variable in \code{metadata} to use on the x-axis for scoring association.
#' @param ColorValues (Optional) A named vector defining the color for the plotted points. If NULL, defaults to a preset color.
#' @param ncol,nrow Number of columns and rows in the facet grid layout. If NULL, computed automatically.
#' @param title Optional string for the overall title of the plot grid.
#' @param widthTitle Maximum character width for titles before inserting line breaks. Default is 10.
#' @param titlesize Numeric value for the font size of plot titles. Default is 12.
#' @param limits Optional numeric vector of length 2 to define y-axis limits.
#' @param pointSize Size of the plotted points. Default is 2.
#' @param xlab Optional label for the x-axis. If NULL, defaults to the name of \code{Variable}.
#' @param labsize Font size for axis labels. Default is 10.
#' @param compute_cohen Logical. If TRUE (default), computes Cohen's f effect size for the association between signature score and the continuous variable.
#' @param pvalcalc Logical. If TRUE, includes the p-value in the plot subtitle. Default is FALSE.
#' @param colorPalette Name of the RColorBrewer palette for coloring. Default is "Set3". Currently unused but kept for consistency.
#' @param cor Character string indicating the correlation method to be used in `ggpubr::stat_cor()`. Options are "pearson" (default), "kendall", or "spearman".

#' @return A ggplot2 object or a multi-plot figure showing scatter plots for each gene signature, with linear regression lines and optional statistical annotations.
#'
#' @details
#' For each gene signature, the function:
#' \itemize{
#'   \item Computes a signature score per sample using the selected method.
#'   \item Plots the score against a continuous metadata variable (\code{Variable}).
#'   \item Adds a regression line and optionally computes and displays Cohen’s f effect size and p-value.
#'   \item Returns a faceted grid of ggplots, arranged by \code{ncol} and \code{nrow}.
#' }
#'
#' This version of the function is specifically tailored for use with continuous variables.
#'
#' @import ggplot2
#' @importFrom ggpubr ggarrange annotate_figure stat_cor
#' @importFrom grid textGrob gpar
#' @export
PlotScores_Numeric <- function(data,
                               metadata,
                               gene_sets,
                               method = c("ssGSEA", "logmedian", "ranking"),
                               Variable = NULL,
                               ColorValues = NULL,
                               ncol = NULL,
                               nrow = NULL,
                               title = NULL,
                               widthTitle = 10,
                               titlesize = 12,
                               limits = NULL,
                               pointSize = 2,
                               xlab = NULL,
                               labsize = 10,
                               compute_cohen = TRUE,
                               pvalcalc = FALSE,
                               colorPalette = "Set3",
                               cor = c("pearson","spearman","kendall")) {

  method <- match.arg(method)


  ResultsList <- CalculateScores(data = data,
                                 metadata = metadata,
                                 gene_sets = gene_sets,
                                 method = method)


  # if grouping variable is NULL, then the function displays a density / distribution of scores
  if (is.null(Variable) | is.null(metadata)) {

    plot_list <- list()

    for (signature in names(ResultsList)) {

      df <- ResultsList[[signature]]
      # Wrap the signature name using the helper function
      wrapped_title <- wrap_title(signature, width = widthTitle)

      ColorValues <- if (is.null(ColorValues)) "#ECBD78" else ColorValues

      p <- ggplot2::ggplot(df, ggplot2::aes(x = score)) +
        ggplot2::geom_density(fill = ColorValues, alpha = 0.5) +
        ggplot2::labs(title = "Density Plot of Score", x = xlab, y = "Density")

      # Customize the plot appearance.
      p <- p + ggplot2::theme_classic() +
        ggplot2::labs(title = wrapped_title, color = "", x = "", y = "") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = labsize - .5),
                       axis.text.y = ggplot2::element_text(  size = labsize - .5),
                       plot.title = ggplot2::element_text(hjust = 0.5, size = titlesize-1),
                       plot.subtitle = ggplot2::element_text(hjust = 0.5, size = titlesize - 1.5, face = "italic"))



      # If limits is specified, crop the plot without adjusting the data (violins).
      if (!is.null(limits)) {
        p <- p + ggplot2::coord_cartesian(xlim = limits)
      }

      plot_list[[signature]] <- p

    }

    n <- length(plot_list)

    # Determine grid layout
    if (is.null(ncol) && is.null(nrow)) {
      ncol <- ceiling(sqrt(n))
      nrow <- ceiling(n / ncol)
    } else if (is.null(ncol)) {
      ncol <- ceiling(n / nrow)
    } else if (is.null(nrow)) {
      nrow <- ceiling(n / ncol)
    }

    # create label for y axis
    if (method == "ssGSEA") {
      xlab <- "ssGSEA Enrichment Score"
    } else if (method == "logmedian") {
      xlab <- "Normalized Signature Score"
    } else if (method == "ranking") {
      xlab <- "Signature Genes' Ranking"
    }

    combined_plot <- ggpubr::ggarrange(plotlist = plot_list, ncol = ncol, nrow = nrow, common.legend = TRUE, align = "h")
    combined_plot <- ggpubr::annotate_figure(combined_plot,
                                             left = grid::textGrob("Density",
                                                                   rot = 90, vjust = 1, gp = grid::gpar(cex = 1.3, fontsize = labsize)),
                                             bottom = grid::textGrob(xlab, gp = grid::gpar(cex = 1.3, fontsize = labsize)),
                                             top = grid::textGrob(title, gp = grid::gpar(cex = 1.3, fontsize = titlesize + 2)))
    return(combined_plot)
  }

  if (!(Variable %in% colnames(metadata)))
    stop(paste0(Variable, " not in metadata columns. Please check metadata."))

  # Initialize an empty list to store individual ggplot objects.
  plot_list <- list()

  # Loop over each gene signature in the ResultsList.
  for (signature in names(ResultsList)) {
    # Extract the data frame for the current signature.
    df <- ResultsList[[signature]]
    #
    #       # Using factors so we can retrieve the first condition for Cohen's d if none is specified.
    #       df[, Variable] <- factor(df[, Variable],
    #                                        levels = sort(unique(as.character(df[, Variable]))))

    # Wrap the signature name using the helper function.
    wrapped_title <- wrap_title(signature, width = widthTitle)

    # Create a base ggplot object with the specified grouping on the x-axis and score on the y-axis.
    p <- ggplot2::ggplot(df, ggplot2::aes_string(x = Variable, y = "score"))

    #add points
    # If ColorValues is provided, use a manual color scale;
    # use a default brewer palette.
    ColorValues <- if (is.null(ColorValues)) "#5264B6" else ColorValues
    p <- p + ggplot2::geom_point(size = pointSize, alpha = 0.5, color=ColorValues[1])

    # Add  line
    p <- p + ggplot2::geom_smooth(method = "lm", col = "black", se = FALSE, size=2) + ggpubr::stat_cor(method=cor) # cor in "pearson" (default), "kendall", or "spearman".

    # Add stats: Compute Cohen's f (and optionally p‑value)
    if(compute_cohen){

      # Calculate Cohen's f
      type <- identify_variable_type(df, Variable)[Variable]
      #Without scaling, the coefficient represents the change in score per unit increase in the variable (if numeric, the unit of the variable. Makes sense to not scale...)
      model <- lm(score ~ get(Variable), data = df)
      results_var <- compute_cohens_f_pval(model, type)

      if (pvalcalc) {
        line1 <- wrap_title(paste0("Cohen's f = ", round(results_var["Cohen_f"], 3)), width = widthTitle)
        line2 <- wrap_title(paste0("p = ", round(results_var["P_Value"], 3)), width = widthTitle)
        subtitle <- paste(line1, line2, sep = "\n")

      } else {
        subtitle <- wrap_title(paste0("Cohen's f = ", round(results_var["Cohen_f"], 3)), width = widthTitle)
      }

    } else {
      subtitle <- NULL
    }


    # Customize the plot appearance.
    p <- p + ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = labsize),
                     axis.text.y = ggplot2::element_text(  size = labsize),
                     plot.title = ggplot2::element_text(hjust = 0.5, size = titlesize-1),
                     plot.subtitle = ggplot2::element_text(hjust = 0.5, size = titlesize - 1.5, face = "italic")) +
      ggplot2::labs(title = wrapped_title, subtitle = subtitle, color = "", x = "", y = "")



    # If limits is specified, crop the plot without adjusting the data (violins).
    if (!is.null(limits)) {
      p <- p + ggplot2::coord_cartesian(ylim = limits)
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
    xlab <- Variable
  }

  if (!is.null(title)) title <- wrap_title(title, width = widthTitle)

  # Create label for y axis based on method.
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


