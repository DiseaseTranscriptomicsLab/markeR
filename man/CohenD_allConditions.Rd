% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Heatmap_Cohen.R
\name{CohenD_allConditions}
\alias{CohenD_allConditions}
\title{Compute Cohen\'s d for All Gene Signatures Across Conditions}
\usage{
CohenD_allConditions(
  data,
  metadata,
  gene_sets,
  variable,
  mode = c("simple", "medium", "extensive")
)
}
\arguments{
\item{data}{A data frame of gene expression data, with genes as rows and samples as columns.}

\item{metadata}{A data frame containing sample metadata. The first column should contain sample identifiers matching
the column names of \code{data}.}

\item{gene_sets}{A named list of gene sets. For unidirectional gene sets, each element is a vector of gene names;
for bidirectional gene sets, each element is a data frame where the first column contains gene names and the second
column indicates the expected direction (1 for upregulated, -1 for downregulated).}

\item{variable}{A string specifying the grouping variable in \code{metadata} used to compare scores between conditions.}

\item{mode}{A string specifying the level of detail for contrasts.
Options are:
\itemize{
\item \code{"simple"}: Pairwise comparisons (e.g., A - B).
\item \code{"medium"}: Pairwise comparisons plus comparisons against the mean of other groups.
\item \code{"extensive"}: All possible groupwise contrasts, ensuring balance in the number of terms on each side.
}}
}
\value{
A named list where each element corresponds to a gene signature. Each signature element is a list with three components:
\describe{
\item{CohenD}{A data frame where rows are methods and columns are group contrasts (formatted as \"Group1:Group2\"),
containing the computed Cohen\'s d effect sizes.}
\item{PValue}{A data frame with the same structure as \code{CohenD} containing the corresponding p-values.}
\item{padj}{A data frame with the same structure as \code{PValue} containing the corresponding p-values corrected using the BH method, for all signatures and contrasts,
and by method.}
}
}
\description{
Computes Cohen\'s d effect sizes and corresponding p-values for all gene signatures using scores calculated
by various methods. The function first computes gene signature scores using \code{CalculateScores} with the "all"
option, flattens the results, and then computes pairwise comparisons for a specified grouping variable.
}
\examples{
\dontrun{
  # Assume gene_data is your gene expression data frame, sample_metadata is your metadata, and
  # gene_sets is a named list of gene sets.
  results <- CohenD_allConditions(data = gene_data, metadata = sample_metadata,
                                  gene_sets = gene_sets, variable = "Condition")
  # Access Cohen's d for a specific signature:
  results$Signature_A$CohenD
}

}
\keyword{internal}
