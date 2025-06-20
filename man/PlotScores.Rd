% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PlotScores.R
\name{PlotScores}
\alias{PlotScores}
\title{Plot gene signature scores using various methods.}
\usage{
PlotScores(
  data,
  metadata,
  gene_sets,
  method = c("ssGSEA", "logmedian", "ranking", "all"),
  ColorVariable = NULL,
  Variable = NULL,
  ColorValues = NULL,
  ConnectGroups = FALSE,
  ncol = NULL,
  nrow = NULL,
  title = NULL,
  widthTitle = 20,
  titlesize = 12,
  limits = NULL,
  legend_nrow = NULL,
  pointSize = 4,
  xlab = NULL,
  labsize = 10,
  compute_cohen = TRUE,
  cond_cohend = NULL,
  pvalcalc = FALSE,
  mode = c("simple", "medium", "extensive"),
  widthlegend = 22,
  sig_threshold = 0.05,
  cohen_threshold = 0.5,
  colorPalette = "Set3",
  cor = c("pearson", "spearman", "kendall")
)
}
\arguments{
\item{data}{A data frame of Normalised (non-transformed) counts where each row is a gene and each column is a sample.
Row names should contain gene names, and column names should contain sample identifiers. \strong{(Required)}}

\item{metadata}{A data frame with sample-level attributes. Each row corresponds to a sample, with the first column containing sample IDs
that match \code{colnames(data)}. \strong{Required if \code{method = "all"} or if metadata-derived groupings or colors are used.}}

\item{gene_sets}{A named list of gene sets to score.
\itemize{
\item For \strong{unidirectional gene sets}: a list of character vectors.
\item For \strong{bidirectional gene sets}: a list of data frames with two columns: gene names and direction (1 = up, -1 = down). \strong{(Required)}
}}

\item{method}{Scoring method to use. One of \code{"ssGSEA"}, \code{"logmedian"}, \code{"ranking"}, or \code{"all"} (default = \code{"logmedian"}).
\itemize{
\item \code{"all"} triggers a full analysis returning both heatmap and volcano plots.
\item Other values return single-score plots depending on \code{Variable} type.
}}

\item{ColorVariable}{Name of a metadata column to color points by. Used in \strong{single-method mode} (\code{"ssGSEA"}, etc.). Ignored in \code{"all"} mode.}

\item{Variable}{Metadata column to define groups or numeric comparisons.
\itemize{
\item \strong{Required if \code{method = "all"}} (used to compute and compare effect sizes).
\item If \code{NULL} and \code{method != "all"}, density plots of each signature scores across samples are shown (no grouping or comparison).
}}

\item{ColorValues}{Optional. A named vector or list of colors used to control the coloring of plot elements across different methods and variable types.
Behavior depends on the combination of \code{method} and \code{Variable}:
\itemize{
\item \strong{If \code{method != "all"}}:
\itemize{
\item If \code{Variable} is \code{NULL}:
\itemize{
\item Used in density plots; a \strong{single color} will be applied (default: \code{"#ECBD78"} if \code{ColorValues} is not specified).
}
\item If \code{Variable} is \strong{categorical}:
\itemize{
\item A \strong{named vector} should map each level of \code{Variable} (or \code{ColorVariable}) to a specific color.
\item Overrides the palette specified by \code{colorPalette}.
}
\item If \code{Variable} is \strong{numeric}:
\itemize{
\item A \strong{single color} is applied to all points in the scatter plot (default: \code{"#5264B6"}).
}
}
\item \strong{If \code{method == "all"}}:
\itemize{
\item \code{ColorValues} can be a \strong{named list} with two elements:
\itemize{
\item \code{heatmap}: A vector of \strong{two colors} used as a diverging scale for the heatmap of effect sizes (default: \code{c("#F9F4AE", "#B44141")}).
\item \code{volcano}: A named vector of colors used for labeling or grouping gene signatures (e.g., in the volcano plot).
}
\item If not provided, defaults will be used for both components.
}
}

In all cases, \code{ColorValues} takes precedence over the default \code{colorPalette} setting if specified.}

\item{ConnectGroups}{Logical. If \code{TRUE}, connects points by sample ID across conditions (used for categorical variables and \code{method != "all"}).}

\item{ncol}{Number of columns for facet layout (used in both heatmaps and score plots).}

\item{nrow}{Number of rows for facet layout (used in both heatmaps and score plots).}

\item{title}{Plot title (optional).}

\item{widthTitle}{Width allocated for title (affects alignment).}

\item{titlesize}{Font size for plot title.}

\item{limits}{Y-axis limits (numeric vector of length 2).}

\item{legend_nrow}{Number of rows for plot legend (used in single-method plots).}

\item{pointSize}{Numeric. Size of points in \strong{score plots} (violin or scatter),
used when plotting individual sample scores for both categorical and numeric variables,
including when \code{method = "all"}.}

\item{xlab}{Label for x-axis (optional; defaults to \code{Variable}).}

\item{labsize}{Font size for axis and facet labels.}

\item{compute_cohen}{Logical. Whether to compute Cohen's effect sizes in \strong{score plots} (\code{method != "all"}).
\itemize{
\item \strong{Only applies when \code{method != "all"}}; ignored otherwise.
\item If the variable is \strong{categorical} and \code{cond_cohend} is specified, computes \strong{Cohen's d} for the specified comparison.
\item If the variable is \strong{categorical} and \code{cond_cohend} is not specified, computes:
\itemize{
\item \strong{Cohen's d} if there are exactly two groups.
\item \strong{Cohen's f} if there are more than two groups.
}
\item If the variable is \strong{numeric}, computes \strong{Cohen's f} regardless of \code{cond_cohend}.
}}

\item{cond_cohend}{Optional. List of length 2 with the two groups being used to compute effect size. The values in each entry should be levels of \verb{Variable (used with }compute_cohen = TRUE`).}

\item{pvalcalc}{Logical. If \code{TRUE}, computes p-values between groups.}

\item{mode}{A string specifying the contrast mode when \code{method = "all"}.
Determines the complexity and breadth of comparisons performed between group levels.
Options are:
\itemize{
\item \code{"simple"}: Performs the minimal number of pairwise comparisons between individual group levels (e.g., A - B, A - C). Default.
\item \code{"medium"}: Includes comparisons between one group and the union of all other groups (e.g., A - (B + C + D)), enabling broader contrasts beyond simple pairs.
\item \code{"extensive"}: Allows for all possible algebraic combinations of group levels (e.g., (A + B) - (C + D)), supporting flexible and complex contrast definitions.
}}

\item{widthlegend}{Width of the legend in \strong{volcano plots} (used only if \code{method = "all"}) and violin score plots.}

\item{sig_threshold}{P-value cutoff shown as a \strong{guide line} in volcano plots. Only applies when \code{method = "all"}.}

\item{cohen_threshold}{Effect size threshold shown as a \strong{guide line} in volcano plots. Used only when \code{method = "all"}.}

\item{colorPalette}{Name of an RColorBrewer palette used to assign colors in plots. Applies to all methods. Default is "Set3".
If \code{ColorValues} is provided, it overrides this palette.
\itemize{
\item If \code{Variable} is \code{NULL} and \code{method != "all"} (i.e., for density plots), a default color \code{"#ECBD78"} is used.
\item If \code{method = "all"} (i.e., for heatmaps and volcano plots), a default diverging color scale is used: \code{c("#F9F4AE", "#B44141")}, unless \code{ColorValues} is manually specified.
}}

\item{cor}{Correlation method for numeric variables. One of \code{"pearson"} (default), \code{"spearman"}, or \code{"kendall"}.
Only applies when the variable is numeric and \code{method != "all"}.}
}
\value{
Depending on \code{method}:
\itemize{
\item If \code{method = "all"}: returns a list with \code{heatmap} and \code{volcano} ggplot objects.
\item If \code{method} is a single method: returns a single ggplot object (scatter or violin plot depending on variable type).
}
}
\description{
Computes and visualizes gene signature scores using one or more methods, returning plots
such as scatter plots, violin plots, heatmaps, or volcano plots depending on inputs.
}
\details{
\strong{Behavior Based on \code{method}:}
\itemize{
\item \code{"all"}:
\itemize{
\item Requires \code{metadata} and \code{Variable}.
\item Computes scores using all available methods and returns:
\itemize{
\item A heatmap of Cohen’s effect sizes.
\item A volcano plot showing effect size vs p-value across gene signatures.
}
\item Uses additional parameters:
\itemize{
\item \code{mode}: defines how contrasts between groups are constructed.
\item \code{sig_threshold} and \code{cohend_threshold}: add \strong{guide dashed lines} to the volcano plot (do not affect point coloring).
\item \code{widthlegend}: controls width of the volcano plot legend.
\item \code{pointSize}: controls dot size for signature points in the volcano plot.
}
\item \code{ColorValues} can be a \strong{named list}:
\itemize{
\item \code{heatmap}: two-color gradient for effect sizes (default: \code{c("#F9F4AE", "#B44141")}).
\item \code{signatures}: named vector of colors for gene signatures in the volcano plot (default is color palette "Set3").
}
}
\item \code{"ssGSEA"}, \code{"logmedian"}, or \code{"ranking"}:
\itemize{
\item The type of \code{Variable} determines the plot:
\itemize{
\item If \strong{categorical}: produces violin plots with optional group comparisons.
\item If \strong{numeric}: produces scatter plots with correlation.
\item If \code{Variable} is \code{NULL}: produces density plots for each signature across all samples.
}
\item Additional arguments:
\itemize{
\item \code{ColorVariable} and \code{ColorValues}: control coloring of points or violins.
\item \code{colorPalette}: default palette (overridden by \code{ColorValues} if present).
\item \code{ConnectGroups}: links samples by ID (for categorical \code{Variable} only).
\item \code{cor}: specifies correlation method for numeric \code{Variable}.
\item \code{pvalcalc}: enables group-wise p-value calculations (categorical only).
\item \code{compute_cohen}: calculates effect sizes when applicable.
\item \code{cond_cohend}: focuses Cohen’s d calculation on a specific comparison.
}
}
}

\strong{Behavior Based on \code{Variable} Type:}
\itemize{
\item \strong{If \code{Variable} is numeric}:
\itemize{
\item Outputs scatter plots (in single-method mode).
\item Computes correlation (\code{cor}).
\item Ignores \code{compute_cohen}, \code{cond_cohend}, and \code{pvalcalc}.
\item Color is uniform (default: \code{"#5264B6"}) unless overridden via \code{ColorValues}.
\item Cohen’s f effect size estimation (\code{compute_cohen = TRUE}) and significance if \code{pvalcalc} is \code{TRUE}.
}
\item \strong{If \code{Variable} is categorical}:
\itemize{
\item Outputs violin plots (in single-method mode).
\item Supports:
\itemize{
\item p-value comparisons (\code{pvalcalc = TRUE}),
\item optional connection lines (\code{ConnectGroups = TRUE}),
\item Cohen’s effect size estimation (\code{compute_cohen = TRUE}) and significance (\code{pvalcalc} is \code{TRUE}):
\itemize{
\item If \code{cond_cohend} is specified, computes \strong{Cohen's d} for that comparison.
\item If not specified:
\itemize{
\item Computes \strong{Cohen’s d} if 2 groups.
\item Computes \strong{Cohen’s f} if >2 groups.
}
}
}
\item Colors are matched to factor levels using \code{ColorValues} or \code{colorPalette}.
}
\item \strong{If \code{Variable} is NULL and \code{method != "all"}}:
\itemize{
\item Produces density plots of signature scores.
\item Uses a \strong{single fill color} (\code{"#ECBD78"} by default or from \code{ColorValues}).
}
}
}
