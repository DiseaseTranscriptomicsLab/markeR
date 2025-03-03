#' Calculate Gene Signature Scores using ssGSEA
#'
#' Computes an enrichment score for each gene signature in each sample using the single-sample Gene Set Enrichment Analysis (ssGSEA).
#'
#' @param data A data frame of normalized (non-transformed) counts where each row is a gene and each column is a sample.
#' @param metadata A data frame containing sample metadata (optional).
#' @param gene_sets Gene set input. **(Required)**
#'
#' If using **unidirectional** gene sets, provide a named list where each element is a vector of gene names representing a gene signature.
#' The names of the list elements should correspond to the labels for each signature.
#'
#' If using **bidirectional** gene sets, provide a named list where each element is a data frame. The names of the list elements should
#' correspond to the labels for each signature, and each data frame should contain the following structure:
#' - The **first column** should contain gene names.
#' - The **second column** should indicate the expected direction of enrichment (1 for upregulated genes, -1 for downregulated genes).
#'
#'
#' @importFrom GSVA gsva
#' @importFrom reshape2 melt
#' @importFrom reshape2 melt
#'
#' @return A list of data frames containing ssGSEA scores for each signature.
#' @keywords internal

CalculateScores_ssGSEA <- function(data, metadata = NULL, gene_sets) {
  ResultsList <- list()

  for (sig in names(gene_sets)) {
    signature <- gene_sets[[sig]]

    if (is.data.frame(signature)) {  # If a data frame, check enrichment direction

      # number of values in the enrichment column
      nb_factors <- length(unique(signature[, 2]))

      if (nb_factors > 1) {  # Both up (1) and down (-1) genes
        message(paste0("Considering bidirectional gene signature mode for signature ", sig))
        dfScore <- CalculateScores_ssGSEA_bidirectional(data, signature)
      } else {  # Only one known direction (up or down)
        message(paste0("Considering unidirectional gene signature mode for signature ", sig))
        dfScore <- CalculateScores_ssGSEA_unidirectional(data, signature)
      }
    } else {  # If vector of genes (unidirectional)
      message(paste0("Considering unidirectional gene signature mode for signature ", sig))
      dfScore <- CalculateScores_ssGSEA_unidirectional(data, signature)
    }

    if (!is.null(metadata)) dfScore <- merge(dfScore, metadata, by = "sample")

    row.names(dfScore) <- NULL
    ResultsList[[sig]] <- dfScore
  }

  return(ResultsList)
}


CalculateScores_ssGSEA_unidirectional <- function (data, metadata = NULL, gene_sets)
{
  ResultsList <- list()
    for (sig in names(gene_sets)) {
      siglist <- list(gene_sets[[sig]]$Gene)
      names(siglist) <- c(sig)
      mtx <- log2(data + 1)
      mtx <- as.matrix(mtx)
      ssgsea_results <- GSVA::gsva(expr = mtx, gset.idx.list = siglist,
                                   method = "ssgsea", kcdf = "Gaussian", verbose = FALSE)
      ssgsea_results <- as.data.frame(ssgsea_results)
      ssgsea_results$signature <- row.names(ssgsea_results)
      ssgsea_results <- reshape2::melt(ssgsea_results)
      colnames(ssgsea_results) <- c("signature", "sample",
                                    "score")
      ssgsea_results$signature <- NULL
      if (!is.null(metadata))
        ssgsea_results <- merge(ssgsea_results, metadata,
                                by = "sample")
      row.names(ssgsea_results) <- NULL
      ResultsList[[sig]] <- ssgsea_results
    }
  return(ResultsList)
}


CalculateScores_ssGSEA_bidirectional <- function (data, metadata = NULL, gene_sets)
{
  ResultsList <- list()

    for (sig in names(gene_sets)) {
      mtx <- log2(data + 1)
      mtx <- as.matrix(mtx)


      up_gene_sets <- subset(gene_sets[[sig]], Signal== 1)
      down_gene_sets <- subset(gene_sets[[sig]], Signal== -1)

      # ssGSEA for UP genes
      up_siglist <- list(up_gene_sets$Gene)
      names(up_siglist) <- sig
      up_results <- GSVA::gsva(expr = mtx, gset.idx.list = up_siglist,
                               method = "ssgsea", kcdf = "Gaussian", verbose = FALSE)

      up_results <- as.data.frame(up_results)
      up_results$signature <- row.names(up_results)
      up_results <- reshape2::melt(up_results)
      colnames(up_results) <- c("signature", "sample", "score_up")
      up_results$signature <- NULL

      # ssGSEA for DOWN genes
      if (nrow(down_gene_sets)>1){
        down_siglist <- list(down_gene_sets$Gene)
        names(down_siglist) <- sig
        down_results <- GSVA::gsva(expr = mtx, gset.idx.list = down_siglist,
                                   method = "ssgsea", kcdf = "Gaussian", verbose = FALSE)

        down_results <- as.data.frame(down_results)
        down_results$signature <- row.names(down_results)
        down_results <- reshape2::melt(down_results)
        colnames(down_results) <- c("signature", "sample", "score_down")
        down_results$signature <- NULL

        # Merge UP & DOWN results
        merged_results <- merge(up_results, down_results, by = "sample", all.x = TRUE)

        # Calculate Final Score (UP - DOWN)
        merged_results$score <- (merged_results$score_up * (nrow(up_gene_sets)/nrow(gene_sets[[sig]]))) - (merged_results$score_down * (nrow(down_gene_sets)/nrow(gene_sets[[sig]])))
        merged_results$score_up <- NULL
        merged_results$score_down <- NULL
      }

      else{
        merged_results <- up_results
        colnames(merged_results) <- c("sample", "score")
        merged_results$score <- (merged_results$score)
      }
      # Merge with metadata if provided
      if (!is.null(metadata)) {
        merged_results <- merge(merged_results, metadata, by = "sample")
      }

      row.names(merged_results) <- NULL
      ResultsList[[sig]] <- merged_results
    }
  return(ResultsList)
}
