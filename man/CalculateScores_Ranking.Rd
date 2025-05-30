% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CalculateScores_Ranking.R
\name{CalculateScores_Ranking}
\alias{CalculateScores_Ranking}
\title{Calculate Gene Signature Scores using Ranking Approach}
\usage{
CalculateScores_Ranking(data, metadata = NULL, gene_sets)
}
\arguments{
\item{data}{A data frame where rows represent genes, columns represent samples, and values
correspond to gene expression levels. \strong{(Required)}}

\item{metadata}{A data frame containing sample metadata. The first column must contain sample
names. \strong{(Optional)}}

\item{gene_sets}{A named list of gene sets. \strong{(Required)}
\itemize{
\item \strong{Unidirectional gene sets}: Provide a named list where each element is a vector of gene names.
\item \strong{Bidirectional gene sets}: Provide a named list where each element is a data frame.
\itemize{
\item The \strong{first column} should contain gene names.
\item The \strong{second column} should indicate the expected direction of enrichment (1 for upregulated genes, -1 for downregulated genes).
}
}}
}
\value{
A named list of data frames, where each data frame contains:
\itemize{
\item \code{sample}: Sample name.
\item \code{score}: Normalized ranking score for the given gene signature.
\item Additional metadata columns (if \code{metadata} is provided).
}
}
\description{
Computes gene signature scores for each sample by ranking the expression of signature genes
in the dataset and normalizing the score based on the total number of genes.
}
\details{
\itemize{
\item The function first validates inputs and extracts relevant genes from the dataset.
\item For \strong{unidirectional} signatures, it computes rankings based on gene expression levels.
\item For \strong{bidirectional} signatures, it computes separate rankings for upregulated and downregulated genes,
then calculates a final score by subtracting downregulated rankings from upregulated rankings.
\item The final scores are normalized by dividing by the total number of genes.
\item This metric is not suitable to compare absolute values between different gene sets, i.e. should be used only for
relative comparisons between samples when using the same gene set.
}
}
\examples{
\dontrun{
# Example dataset with 5 genes (rows) and 3 samples (columns)
set.seed(123)
data <- as.data.frame(matrix(runif(15, 1, 100), nrow = 5, ncol = 3))
rownames(data) <- paste0("Gene_", 1:5)
colnames(data) <- paste0("Sample_", 1:3)

# Unidirectional gene set example
gene_sets <- list(Signature1 = c("Gene_1", "Gene_3", "Gene_5"))

# Compute scores
scores <- CalculateScores_Ranking(data, gene_sets = gene_sets)
print(scores)
}
}
\keyword{internal}
