#' Calculate Gene Signature Scores using Log-Median Approach
#'
#' Computes log2-median-centered scores for each sample based on gene signature expression.
#'
#' @param data A data frame of normalized counts where each row is a gene and each column is a sample.
#' @param metadata A data frame containing sample metadata (optional). If provided, the resulting scores will be merged with metadata.
#' @param gene_sets A named list representing gene sets. **(Required)**
#'
#' - **Unidirectional gene sets:** Each element should be a vector of gene names representing a signature.
#'   - The names of the list elements serve as labels for the signatures.
#' - **Bidirectional gene sets:** Each element should be a data frame with:
#'   - The **first column** containing gene names.
#'   - The **second column** specifying the expected direction of enrichment:
#'     - `1` for upregulated genes.
#'     - `-1` for downregulated genes.
#'
#' @return A list of data frames containing log-median scores for each signature. If `metadata` is provided, it is merged with the scores.
#' @examples
#' \dontrun{
#' data <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' colnames(data) <- paste0("Sample_", 1:10)
#' rownames(data) <- paste0("Gene_", 1:100)
#' gene_sets <- list(
#'   Signature_A = sample(rownames(data), 10),
#'   Signature_B = data.frame(Gene = sample(rownames(data), 10), Direction = sample(c(1, -1), 10, replace = TRUE))
#' )
#' scores <- CalculateScores_logmedian(data, gene_sets = gene_sets)
#' }
#' @export
CalculateScores_logmedian <- function(data, metadata = NULL, gene_sets) {
  ResultsList <- list()

  for (sig in names(gene_sets)) {
    signature <- gene_sets[[sig]]

    if (is.data.frame(signature)) {  # If a data frame, check enrichment direction

      # number of values in the enrichment column
      nb_factors <- length(unique(signature[, 2]))

      if (nb_factors > 1) {  # Both up (1) and down (-1) genes
        message(paste0("Considering bidirectional gene signature mode for signature ", sig))
        dfScore <- calculateScore_logmedian_bidirectional(data, signature)
      } else {  # Only one known direction (up or down)
        message(paste0("Considering unidirectional gene signature mode for signature ", sig))
        dfScore <- calculateScore_logmedian_unidirectional(data, signature)
      }
    } else {  # If vector of genes (unidirectional)
      message(paste0("Considering unidirectional gene signature mode for signature ", sig))
      dfScore <- calculateScore_logmedian_unidirectional(data, signature)
    }

    if (!is.null(metadata)) dfScore <- merge(dfScore, metadata, by = "sample")

    row.names(dfScore) <- NULL
    ResultsList[[sig]] <- dfScore
  }

  return(ResultsList)
}

#' Calculate Log-Median Scores for Bidirectional Gene Sets
#'
#' Computes gene signature scores considering both upregulated and downregulated genes separately,
#' then calculates a differential score by subtracting downregulated from upregulated scores.
#'
#' @param data A data frame of normalized counts (genes as rows, samples as columns).
#' @param signature A data frame with:
#'   - The **first column** containing gene names.
#'   - The **second column** specifying enrichment direction (`1` for upregulated, `-1` for downregulated).
#'
#' @return A named vector with log-median-centered scores per sample.
#' @keywords internal
calculateScore_logmedian_bidirectional <- function(data, signature) {
  signature[, 2] <- as.numeric(signature[, 2])

  if (!all(unique(signature[, 2]) %in% c(-1, 1))) {
    stop("Values for enrichment not supported. Please use 1 for enriched genes and -1 for depleted genes.")
  }

  signaturegenes_up <- signature[signature[, 2] == 1, ]
  signaturegenes_down <- signature[signature[, 2] == -1, ]

  # Compute log-median scores for upregulated genes
  data_subset_up <- na.omit(subset(log2(data + 1), row.names(data) %in% signaturegenes_up[, 1]))
  data_subset_up <- data_subset_up - apply(data_subset_up, 1, median)
  score_up <- colSums(data_subset_up) / nrow(data_subset_up)

  # Compute log-median scores for downregulated genes
  data_subset_down <- na.omit(subset(log2(data + 1), row.names(data) %in% signaturegenes_down[, 1]))
  data_subset_down <- data_subset_down - apply(data_subset_down, 1, median)
  score_down <- colSums(data_subset_down) / nrow(data_subset_down)

  score <- score_up - score_down
  dfScore <- data.frame(sample = names(score), score = score)

  return(dfScore)
}

#' Calculate Log-Median Scores for Unidirectional Gene Sets
#'
#' Computes log-median-centered scores for gene signatures where all genes are expected to be enriched in
#' the same direction, or when direction is not known.
#'
#' @param data A data frame of normalized counts (genes as rows, samples as columns).
#' @param signature A vector of gene names or a data frame where the first column contains gene names.
#'
#' @return A named vector with log-median-centered scores per sample.
#' @keywords internal
calculateScore_logmedian_unidirectional <- function(data, signature) {
  if (is.data.frame(signature)) signature <- as.vector(signature[, 1])

  data_subset <- na.omit(subset(log2(data + 1), row.names(data) %in% signature))

  data_subset <- data_subset - apply(data_subset, 1, median)  # Center gene in its log2 median
  dfScore <- colSums(data_subset) / nrow(data_subset)  # Normalize by signature size
  dfScore <- data.frame(sample = names(dfScore), score = dfScore)

  return(dfScore)
}
