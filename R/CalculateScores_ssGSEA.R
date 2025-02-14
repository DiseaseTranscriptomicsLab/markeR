#' Calculate Gene Signature Scores using ssGSEA
#'
#' Computes an enrichment score for each gene signature in each sample using the single-sample Gene Set Enrichment Analysis (ssGSEA).
#'
#' @param data A data frame of normalized (non-transformed) counts where each row is a gene and each column is a sample.
#' @param metadata A data frame containing sample metadata (optional).
#' @param gene_sets A named list of gene signatures.
#'
#' @return A list of data frames containing ssGSEA scores for each signature.
#' @keywords internal
CalculateScores_ssGSEA <- function(data, metadata = NULL, gene_sets) {

  ResultsList <- list()

  for (sig in names(gene_sets)) {
    siglist <- list(gene_sets[[sig]])
    names(siglist) <- c(sig)

    mtx <- as.matrix(log2(data))
    ssgsea_results <- GSVA::gsva(expr = mtx, gset.idx.list = siglist, method = "ssgsea", kcdf = "Gaussian",verbose=FALSE)

    # Format results
    ssgsea_results <- as.data.frame(ssgsea_results)
    ssgsea_results$signature <- row.names(ssgsea_results)
    ssgsea_results <- reshape2::melt(ssgsea_results)
    colnames(ssgsea_results) <- c("signature", "sample", "score")
    ssgsea_results$signature <- NULL

    # Merge with metadata
    if (!is.null(metadata)) ssgsea_results <- merge(ssgsea_results, metadata, by = "sample")

    row.names(ssgsea_results) <- NULL
    ResultsList[[sig]] <- ssgsea_results
  }

  return(ResultsList)
}
