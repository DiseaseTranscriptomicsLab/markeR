#' Compute Cohen's f and p-value for a given model and predictor
#'
#' This function calculates the Cohen's f effect size and the corresponding p-value for a given linear model
#' or ANOVA model based on the predictor variable type (numeric or categorical).
#'
#' @param model A linear model (\code{lm}) or ANOVA model (\code{aov}) fitted to the data.
#' @param predictor A character string specifying the name of the predictor variable in the model.
#' @param type A string indicating whether the predictor is numeric or categorical. Options are "Numeric" or "Categorical".
#'
#' @details
#' Cohen's f effect size is computed from the eta-squared (\eqn{\eta^2}) value.
#' For numeric predictors (continuous variables), the p-value is obtained from the t-test in \code{summary(lm(...))}.
#' For categorical predictors (binary or multi-level), the p-value is obtained from the F-test in \code{anova(lm(...))}.
#'
#' @return A named vector with two elements:
#' \itemize{
#'   \item \code{Cohen_f}: The Cohen's f effect size value.
#'   \item \code{P_Value}: The p-value from the statistical test.
#' }
#'
#' @examples
#' model <- lm(Sepal.Length ~ Petal.Width, data = iris)
#' compute_cohens_f_pval(model, "Petal.Width", type = "Numeric")
#'
#' @importFrom effectsize eta_squared
#' @importFrom stats anova lm
#'
#' @export
compute_cohens_f_pval <- function(model, predictor, type) {
  eta_sq_value <- effectsize::eta_squared(model, partial = FALSE)$Eta2
  cohen_f <- sqrt(eta_sq_value / (1 - eta_sq_value))  # Convert η² to Cohen's f

  # Get p-value
  if (type=="Numeric") {
    # Numeric predictor: get p-value from regression summary
    #For numeric predictors (continuous variables): The t-test from summary(lm(...)) is appropriate because it directly tests the significance of the slope.
    p_value <- summary(model)$coefficients[2, 4]
  } else {
    # Categorical predictor: get p-value from ANOVA
    #For categorical predictors (binary or multi-level): The F-test from anova(lm(...)) is preferred because it tests for overall group differences.
    p_value <- anova(model)$`Pr(>F)`[1]
  }

  return(c(Cohen_f = cohen_f, P_Value = p_value))
}


#' Remove Division Notation in Contrast Labels
#'
#' This function removes division notation (e.g., `/2`, `/3`) after closing parentheses in contrast labels.
#'
#' @param contrasts A character vector containing contrast labels.
#' @return A character vector with division notation removed.
#' @export
remove_division <- function(contrasts) {
  gsub("\\)/[0-9]+", ")", contrasts)  # Removes "/2", "/3", etc. after parentheses
}

#' Create Contrast Column in Metadata
#'
#' This function extracts and processes contrast groups from a given contrast string,
#' then assigns contrast labels to a metadata subset based on the variable of interest.
#'
#' @param metadata A data frame containing sample metadata.
#' @param variable_name A character string specifying the column name in `metadata` that represents the variable of interest.
#' @param contrast A character string representing the contrast in the form "(A + B) - (C + D)" (e.g.).
#' @return A subset of `metadata` with an added `cohentest` column, indicating group membership based on the contrast.
#' @export
create_contrast_column <- function(metadata, variable_name, contrast) {
  # Extract left and right sides of the contrast
  sides <- strsplit(contrast, " - ")[[1]]

  left_side <- gsub("[()]", "", sides[1])  # Remove parentheses
  right_side <- gsub("[()]", "", sides[2]) # Remove parentheses

  # Split into individual levels
  left_levels <- unlist(strsplit(left_side, " \\+ "))
  right_levels <- unlist(strsplit(right_side, " \\+ "))

  # Keep only relevant rows (samples in the contrast)
  relevant_metadata <- metadata[metadata[[variable_name]] %in% c(left_levels, right_levels), ]

  # Create new column with labels based on the contrast
  relevant_metadata$cohentest <- ifelse(relevant_metadata[[variable_name]] %in% right_levels,
                                        right_side,
                                        left_side)

  return(relevant_metadata)  # Return only the relevant subset
}

#' Score Variable Association
#'
#' This function evaluates the association between gene expression scores and metadata variables.
#' It uses linear modeling to get Cohen's F, and contrast-based comparisons for categorical variables to compute Cohen's D.
#' The function generates plots summarizing the results.
#'
#' @param data A data frame or matrix containing gene expression data.
#' @param metadata A data frame containing sample metadata with at least one column corresponding to the variables of interest.
#' @param cols A character vector specifying metadata columns to analyze.
#' @param method A character string specifying the scoring method (`"logmedian"`, `"ssGSEA"`, or `"ranking"`).
#' @param gene_set A named list containing one gene set for scoring.
#' @param mode A character string specifying the contrast generation method (`"simple"`, `"medium"`, `"extensive"`). Four methods are available:
#'
#'   - **ssGSEA**: Uses the single-sample Gene Set Enrichment Analysis (ssGSEA) method to compute an enrichment score
#'     for each signature in each sample using an adaptation of the `gsva()` function from the `GSVA` package.
#'   - **logmedian**: Computes the score as the sum of the normalized (log2-median-centered) expression values of the
#'     signature genes divided by the number of genes in the signature.
#'   - **ranking**: Computes gene signature scores for each sample by ranking the expression of signature genes
#'     in the dataset and normalizing the score based on the total number of genes.
#'
#' @param nonsignif_color A string specifying the color for non-significant results. Default: `"grey"`.
#' @param signif_color A string specifying the color for significant results. Default: `"red"`.
#' @param saturation_value A numeric value for color saturation threshold. Default: `NULL` (auto-determined).
#' @param sig_threshold A numeric value specifying the significance threshold. Default: `0.05`.
#' @param widthlabels An integer controlling contrast label wrapping. Default: `18`.
#' @param labsize An integer controlling axis text size. Default: `10`.
#' @param title A string specifying the plot title. Default: `NULL`.
#' @param titlesize An integer specifying the title size. Default: `14`.
#' @param pointSize A numeric value for point size in plots. Default: `5`.
#' @param discrete_colors A named list mapping categorical variable levels to colors.
#' Each element should be a named vector where names correspond to factor levels. Default: `NULL`.
#' @param continuous_color A string specifying the color for continuous variables. Default: `"#8C6D03"`.
#' @param color_palette A string specifying the color palette for discrete variables. Default: `"Set2"`.
#'
#' @return A list with:
#'   - `Overall`: Data frame with Cohen's F and p-values for overall associations.
#'   - `Contrasts`: Data frame with Cohen's D and p-values for contrasts.
#'   - `plot`: Combined visualization of results.
#'   - `plot_contrasts`: Plot of contrast results.
#'   - `plot_overall`: Plot of overall associations.
#'   - `plot_distributions`: List of variable distribution plots.
#'
#' @examples
#' data <- matrix(rnorm(1000), ncol = 10)
#' metadata <- data.frame(group = rep(c("A", "B"), each = 5))
#' gene_set <- list(SampleSet = c("Gene1", "Gene2", "Gene3"))
#' results <- Score_VariableAssociation(data, metadata, cols = "group", gene_set = gene_set)
#' print(results$plot)
#'
#' @importFrom ggpubr ggarrange
#' @importFrom ggpubr annotate_figure
#' @import ggplot2
#' @importFrom grid textGrob
#' @importFrom grid gpar
#'
#' @export
#'
#'
Score_VariableAssociation <- function(data, metadata, cols, method=c("logmedian","ssGSEA","ranking"), gene_set,
                                      mode=c("simple","medium","extensive"),
                                      nonsignif_color = "grey", signif_color = "red", saturation_value=NULL,sig_threshold = 0.05,
                                      widthlabels=18, labsize=10, title=NULL, titlesize=14, pointSize=5, discrete_colors=NULL,
                                      continuous_color = "#8C6D03", color_palette = "Set2"){

  # calculate scores for a given metric
  df_ranking <- CalculateScores(data = data, metadata = metadata, method = method, gene_sets = gene_set)[[1]] # if more than one gene set is provided, only results for the first one will be returned


  # Initialise results data frame
  df_results_overall <- data.frame(Variable = NULL,
                                   cohen_f = NULL,
                                   p_value = NULL)

  df_results_contrast <- data.frame(Contrast = NULL,
                                    Group1 = NULL,
                                    Group2 = NULL,
                                    CohenD = NULL,
                                    PValue = NULL)

  for (var in cols) {

    # Check if the variable exists in metadata
    if (!(var %in% colnames(metadata))) {
      warning(paste0("Variable ", var, " not in metadata"))
      next  # Skip to the next iteration if the variable is not in metadata
    }

    ########### OVERALL MODE ############

    # for each variable, fit a linear model between the score and that variable
    # Define variable type (numeric or non numeric)
    type <- identify_variable_type(metadata, var)[var]
    #Without scaling, the coefficient represents the change in score per unit increase in the variable (if numeric, the unit of the variable. Makes sense to not scale...)
    model <- lm(score ~ get(var), data = df_ranking)
    results_var <- compute_cohens_f_pval(model, predictor, type)

    df_results_overall <- rbind(
      df_results_overall,
      data.frame(
        Variable = var,
        Cohen_f = results_var["Cohen_f"],
        P_Value = results_var["P_Value"],
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    )

    ########### CONTRASTS MODE ############

    if (type =="Numeric"){
      next
    }
    contrasts <- generate_all_contrasts(unique(df_ranking[[var]]), mode = mode)
    contrasts <- remove_division(contrasts)

    for (cont in contrasts){

      metadata_cohentest <- create_contrast_column(df_ranking, var, cont) # subsets metadata and adds new column "cohentest" with the two parts of the contrast
      result_cohen_contrast <- compute_cohen_d(metadata_cohentest, "cohentest", quantitative_var="score")
      df_results_contrast <-  rbind(df_results_contrast,
                                    data.frame(
                                      Variable = var,
                                      Contrast = paste0("(",result_cohen_contrast$Group1,") - (",result_cohen_contrast$Group2,")"),
                                      Group1 = result_cohen_contrast$Group1,
                                      Group2 = result_cohen_contrast$Group2,
                                      CohenD = result_cohen_contrast$CohenD,
                                      PValue = result_cohen_contrast$PValue),
                                    stringsAsFactors = FALSE,
                                    row.names = NULL)



    }


  }

  df_results_contrast$padj <- p.adjust(df_results_contrast$PValue, method = "BH")


  ########### CONTRAST MODE PLOT ############

  # Ensure contrast ordering
  df_results_contrast$Contrast <- sapply(df_results_contrast$Contrast, function(x) wrap_title(x, widthlabels))
  df_results_contrast$Contrast <- factor(df_results_contrast$Contrast, levels = df_results_contrast$Contrast[order(df_results_contrast$CohenD)])


  if(is.null(saturation_value)){
    if (min(df_results_contrast$padj)>sig_threshold){
      limit_pval <- 0.001
    } else{
      limit_pval <- min(df_results_contrast$padj)
    }

  } else {
    limit_pval <- saturation_value
  }


  plot_contrasts <- ggplot2::ggplot(df_results_contrast, ggplot2::aes(x = CohenD, y = Contrast, fill = -log10(padj))) +
    ggplot2::geom_segment(ggplot2::aes(
      yend = Contrast,
      xend = 0,
      linetype =  "solid",
      color =  "black"
    ), size = .5) +
    ggplot2::geom_point(ggplot2::aes(
      stroke = 1.2,
      color =   "black"
    ), shape = 21, size = pointSize) +
    ggplot2::scale_fill_gradient2(low = nonsignif_color,
                                  mid = nonsignif_color,
                                  high = signif_color,
                                  midpoint = -log10(sig_threshold),
                                  limits=c(0,-log10(limit_pval)),
                                  na.value = signif_color)+
    ggplot2::scale_linetype_identity() +
    ggplot2::scale_color_identity() +
    ggplot2::labs(
      x = "Cohen's D",
      y = "Contrast",
      color = "-log10(Adj. p-value)",
      fill = "-log10(Adj. p-value)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = titlesize),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, face = "italic", size = titlesize - 2),
      legend.position = "top",
      axis.text = ggplot2::element_text(size = labsize)
    ) +
    ggplot2::facet_grid(Variable ~.,   scales = "free", switch = "y", space = "free" ) +
    theme(strip.background =element_rect(fill="white"))


  ########### OVERALL MODE PLOT ############

  plot_overall <- ggplot2::ggplot(df_results_overall, ggplot2::aes(x = Cohen_f, y = Variable, fill = -log10(P_Value))) +
    ggplot2::geom_segment(ggplot2::aes(
      yend = Variable,
      xend = 0,
      linetype =  "solid",
      color =  "black"
    ), size = .5) +
    ggplot2::geom_point(ggplot2::aes(
      stroke = 1.2,
      color =   "black"
    ), shape = 21, size = pointSize) +
    ggplot2::scale_fill_gradient2(low = nonsignif_color,
                                  mid = nonsignif_color,
                                  high = signif_color,
                                  midpoint = -log10(sig_threshold),
                                  limits=c(0,-log10(limit_pval)),
                                  na.value = signif_color)+
    ggplot2::scale_linetype_identity() +
    ggplot2::scale_color_identity() +
    ggplot2::labs(
      x = "Cohen's F",
      y = "Variable",
      color = "-log10(p-value)",
      fill = "-log10(p-value)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = titlesize),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, face = "italic", size = titlesize - 2),
      legend.position = "top",
      axis.text = ggplot2::element_text(size = labsize)
    )



  ########### DISTRIBUTION VARS PLOT ############


  variable_types <- identify_variable_type(metadata, cols = cols)
  list_plts_var_distribution <- list()
  for (var in cols){

    # Check if the variable exists in metadata
    if (!(var %in% colnames(metadata))) {
      warning(paste0("Variable ", var, " not in metadata"))
      next  # Skip to the next iteration if the variable is not in metadata
    }

    # Adding p-value in parenthesis next to the metric
    if (variable_types[var] == "Numeric") {

      list_plts_var_distribution[[var]] <- ggplot2::ggplot(df_ranking, ggplot2::aes_string(x = var, y = "score")) +
        ggplot2::geom_point(alpha = 0.6, size=4, color = continuous_color) +  # Use continuous_color
        ggplot2::geom_smooth(method = "lm", col = "black", se = FALSE, size=2) +
        ggplot2::coord_cartesian(clip = "off") +  # Allow text outside the plot area
        ggplot2::theme_classic() +
        ggplot2::ggtitle(var) +
        ggplot2::labs(y=paste0("Score (",method,")")) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

    } else {

      # Check if user provided a custom named list for discrete colors
      if (!is.null(discrete_colors) && var %in% names(discrete_colors)) {
        # Use user-specified colors for the specific variable
        colors <- discrete_colors[[var]]
      } else {
        num_levels <- length(unique(df_ranking[[var]]))
        colors <- colorRampPalette(RColorBrewer::brewer.pal(8, color_palette))(num_levels)

      }

      list_plts_var_distribution[[var]] <- ggplot2:: ggplot(df_ranking, ggplot2::aes_string(x = "score", fill = var)) +
        ggplot2::geom_density(alpha = 0.6) +
        ggplot2::scale_fill_manual(values = colors) +  # Apply custom discrete colors
        ggplot2::coord_cartesian(clip = "off") +
        ggplot2::theme_classic() +
        ggplot2::ggtitle(var) +
        ggplot2::labs(x = paste0("Score (",method,")"), y = "Density", fill="") +
        ggplot2::theme(legend.position = "right",
                       plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
    }
  }

  plot_distributions <- ggarrange(plotlist = list_plts_var_distribution, ncol=1)

  ########### ALL ############

  row1 <- ggpubr::ggarrange(plot_distributions,plot_contrasts, nrow=1, widths=c(0.4,0.6))
  plotfinal <- ggpubr::ggarrange(plot_overall, row1, ncol=1, heights=c(0.3,0.7))
  plotfinal <- ggpubr::annotate_figure(plotfinal, top = grid::textGrob(names(gene_set)[1], gp = grid::gpar(cex = 1.3, fontsize = titlesize, face="bold")))


  print(plotfinal)

  invisible(list(Overall=df_results_overall,
                 Contrasts = df_results_contrast,
                 plot=plotfinal,
                 plot_contrasts=plot_contrasts,
                 plot_overall=plot_overall,
                 plot_distributions=list_plts_var_distribution))


}
