#' Calculate Gene Signature Scores using Log-Median Approach
#'
#' Computes a log2-median-centered score for each sample based on signature gene expression.
#'
#' @param data A data frame of normalized counts where each row is a gene and each column is a sample.
#' @param metadata A data frame containing sample metadata (optional).
#' @param gene_sets A named list of gene signatures.
#'
#' @return A list of data frames containing log-median scores for each signature.
#' @keywords internal
CalculateScores_logmedian <- function(data, metadata = NULL, gene_sets) {

  ResultsList <- list()

  for (sig in names(gene_sets)) {
    signature <- gene_sets[[sig]]
    data_subset <- na.omit(subset(log2(data + 1), row.names(log2(data + 1)) %in% signature))
    data_subset <- data_subset - apply(data_subset, 1, median)  # Center gene in its log2 median

    dfScore <- colSums(data_subset) / nrow(data_subset)  # Normalize by signature size
    dfScore <- data.frame(sample = names(dfScore), score = dfScore)

    if (!is.null(metadata)) dfScore <- merge(dfScore, metadata, by = "sample")

    row.names(dfScore) <- NULL
    ResultsList[[sig]] <- dfScore
  }

  return(ResultsList)
}
