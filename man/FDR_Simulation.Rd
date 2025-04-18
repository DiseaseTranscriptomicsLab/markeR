% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FDR_Simulation.R
\name{FDR_Simulation}
\alias{FDR_Simulation}
\title{FDR Simulation Plot}
\usage{
FDR_Simulation(
  data,
  metadata,
  original_signatures,
  Variable,
  gene_list = NULL,
  number_of_sims = 10,
  title = NULL,
  widthTitle = 30,
  titlesize = 12,
  pointSize = 2,
  labsize = 10,
  mode = c("none", "simple", "medium", "extensive"),
  ColorValues = NULL,
  ncol = NULL,
  nrow = NULL
)
}
\arguments{
\item{data}{A data frame or matrix of gene expression values (genes as rows, samples as columns).}

\item{metadata}{A data frame containing metadata for the samples (columns of \code{data}).}

\item{original_signatures}{A named list of gene signatures. Each element can be either:
\itemize{
\item A vector of gene names (unidirectional), or
\item A data frame with columns \code{"Gene"} and \code{"Signal"} for bidirectional signatures.
}}

\item{Variable}{A column in \code{metadata} indicating the variable of interest for grouping or regression.
This can be categorical or numeric.}

\item{gene_list}{A character vector of gene names from which simulated signatures are generated by sampling.
Default is all genes in \code{data}.}

\item{number_of_sims}{Integer. Number of simulated gene signatures to generate per original signature.}

\item{widthTitle}{Integer. Max width for wrapping the title text (default: 30).}

\item{titlesize}{Numeric. Font size for the title text (default: 12).}

\item{pointSize}{Numeric. Size of the points representing simulations (default: 2).}

\item{labsize}{Numeric. Font size for axis labels (default: 10).}

\item{mode}{Character. Type of contrast calculation if \code{Variable} is categorical. Options are:
\itemize{
\item \code{"simple"}: basic pairwise contrasts,
\item \code{"medium"}: pairwise + comparisons vs group mean,
\item \code{"extensive"}: all possible balanced contrasts,
\item \code{"none"}: Comparing all levels of \code{Variable} (default)
}}

\item{ColorValues}{Named vector of colors for plot points, typically \code{Original} and \code{Simulated}. If \code{NULL}, default colors are used.}

\item{ncol, }{nrow Integers. Number of columns/rows for arranging signature plots in a grid layout. If \code{NULL}, layout is auto-calculated.}

\item{title_for_plot}{Character. Title prefix used in the overall plot.}
}
\value{
Invisibly returns the combined \code{ggplot} object showing observed vs simulated effect sizes.
One violin plot is generated per signature and contrast. Observed values are highlighted and compared
to the simulated distribution. Significance (adjusted p-value ≤ 0.05) is indicated by point shape.
}
\description{
This function simulates false discovery rates (FDR) by generating simulated gene signatures and comparing
the observed effect size values (Cohen's \emph{d} or \emph{f}) of the original signatures to those from simulated signatures.
The effect size is computed using three scoring methods (\code{ssGSEA}, \code{logmedian}, and \code{ranking}), and the results
are visualized as violin plots with overlaid observed values.
}
\details{
The function supports both categorical and numeric variables:
\itemize{
\item For \strong{categorical variables}, Cohen's \emph{d} is used and contrasts are defined by the \code{mode} parameter, if \code{mode!=none}.
\item For \strong{numeric variables}, Cohen's \emph{f} is used to quantify associations through linear modeling.
}

For each original gene signature, a number of simulated signatures are created by sampling genes from \code{gene_list}.
Each simulated signature is scored using three methods, and its effect size is computed relative to the variable of interest.
The resulting distributions are shown as violins, overlaid with the observed value from the original signature.
A red dashed line marks the 95th percentile of the simulated distribution per method.

The function internally uses \code{CohenD_allConditions()} and \code{CohenF_allConditions()} depending on variable type.
}
\examples{
\dontrun{
FDR_Simulation(
  data = expression_data,
  metadata = sample_metadata,
  original_signatures = my_signatures,
  Variable = "condition",
  number_of_sims = 100,
  title_for_plot = "Simulation for FDR"
)
}

}
