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

      colnames(signature) <- c("Gene","Signal")

      # number of values in the enrichment column
      nb_factors <- length(unique(signature[, 2]))

      if (nb_factors > 1) {  # Both up (1) and down (-1) genes
        message(paste0("Considering bidirectional gene signature mode for signature ", sig))
        dfScore <- CalculateScores_ssGSEA_bidirectional(data, signature)
      } else {  # Only one known direction (up or down)
        message(paste0("Considering unidirectional gene signature mode for signature ", sig))
        signature <- signature$Gene
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


#' Calculate ssGSEA Scores for Unidirectional Gene Signatures
#'
#' Computes single-sample Gene Set Enrichment Analysis (ssGSEA) scores for each sample using
#' a unidirectional gene signature.
#'
#' @param data A data frame of normalized (non-transformed) counts where rows are genes and columns are samples.
#' @param signature A vector of gene names representing a unidirectional gene signature.
#'
#' @importFrom reshape2 melt
#'
#' @return A data frame containing:
#' - `sample`: Sample name.
#' - `score`: ssGSEA enrichment score for the gene signature.
#'
#'
#' @examples
#' # Example dataset with 5 genes (rows) and 3 samples (columns)
#' set.seed(123)
#' data <- matrix(runif(15, 1, 100), nrow = 5, ncol = 3)
#' rownames(data) <- paste0("Gene_", 1:5)
#' colnames(data) <- paste0("Sample_", 1:3)
#'
#' # Define a unidirectional gene signature
#' signature <- c("Gene_1", "Gene_3", "Gene_5")
#'
#' # Compute scores
#' scores <- CalculateScores_ssGSEA_unidirectional(data, signature = signature)
#' print(scores)
#'
#' @export
CalculateScores_ssGSEA_unidirectional <- function(data, signature) {
  ResultsList <- list()

  siglist <- list(signature)
  mtx <- log2(data)
  mtx <- as.matrix(mtx)

  # ssgsea_results <- GSVA::gsva(expr = mtx, gset.idx.list = siglist,
  #                              method = "ssgsea", kcdf = "Gaussian", verbose = FALSE)
  ssgsea_results <- ssGSEA_alternative(X = mtx, gene_sets = siglist)

  ssgsea_results <- as.data.frame(ssgsea_results)
  ssgsea_results <- reshape2::melt(ssgsea_results)
  colnames(ssgsea_results) <- c("sample", "score")
  row.names(ssgsea_results) <- NULL

  return(ssgsea_results)
}


#' Calculate ssGSEA Scores for Bidirectional Gene Signatures
#'
#' Computes single-sample Gene Set Enrichment Analysis (ssGSEA) scores for each sample using
#' a bidirectional gene signature (separating upregulated and downregulated genes).
#'
#' @param data A data frame of normalized (non-transformed) counts where rows are genes and columns are samples.
#' @param signature A data frame with:
#' - The **first column** containing gene names.
#' - The **second column** (`Signal`) indicating the expected direction of enrichment (1 for upregulated genes, -1 for downregulated genes).
#'
#' @importFrom reshape2 melt
#'
#' @return A data frame containing:
#' - `sample`: Sample name.
#' - `score`: Final ssGSEA enrichment score (computed as the difference between upregulated and downregulated scores).
#'
#' @details
#' - The input gene expression matrix (`data`) is log2-transformed before applying ssGSEA.
#' - Upregulated and downregulated genes are analyzed separately.
#' - As both upregulated and downregulated genes are present, the final score is computed as:
#'   \deqn{score = (score_{up} \* \frac{|up\_genes|}{|total\_genes|}) - (score_{down} \* \frac{|down\_genes|}{|total\_genes|})}
#' - If no downregulated genes are present, only the upregulated score is used.
#' - The results are reshaped into a long-format data frame with one score per sample.
#'
#' @examples
#' # Example dataset with 5 genes (rows) and 3 samples (columns)
#' set.seed(123)
#' data <- matrix(runif(15, 1, 100), nrow = 5, ncol = 3)
#' rownames(data) <- paste0("Gene_", 1:5)
#' colnames(data) <- paste0("Sample_", 1:3)
#'
#' # Define a bidirectional gene signature
#' signature <- data.frame(Gene = c("Gene_1", "Gene_3", "Gene_5"),
#'                         Signal = c(1, -1, 1))
#'
#' # Compute scores
#' scores <- CalculateScores_ssGSEA_bidirectional(data, signature = signature)
#' print(scores)
#'
#' @export
CalculateScores_ssGSEA_bidirectional <- function(data, signature) {
  ResultsList <- list()

  mtx <- log2(data)
  mtx <- as.matrix(mtx)

  up_genes <- subset(signature, Signal == 1)$Gene
  down_genes <- subset(signature, Signal == -1)$Gene

  ################## ssGSEA for UP genes ##################

  #up_siglist <- list(up_gene_sets)
  #names(up_siglist) <- sig
  # up_results <- GSVA::gsva(expr = mtx, gset.idx.list = list(up_genes),
  #                          method = "ssgsea", kcdf = "Gaussian", verbose = FALSE)
  up_results <- ssGSEA_alternative(X = mtx, gene_sets = list(up_genes))

  up_results <- as.data.frame(up_results)
  up_results <- reshape2::melt(up_results)
  colnames(up_results) <- c("sample", "score_up")
  row.names(up_results) <- NULL

  ################## ssGSEA for DOWN genes ##################

  # down_results <- GSVA::gsva(expr = mtx, gset.idx.list = list(down_genes),
  #                            method = "ssgsea", kcdf = "Gaussian", verbose = FALSE)

  down_results <- ssGSEA_alternative(X = mtx, gene_sets = list(down_genes))

  down_results <- as.data.frame(down_results)
  down_results <- reshape2::melt(down_results)
  colnames(down_results) <- c("sample", "score_down")
  row.names(down_results) <- NULL


  ################## Merge ssGSEA results ##################

  ssGSEAresults <- merge(up_results, down_results, by = "sample")

  # Calculate Final Score
  frac_upgenes <- length(up_genes)/nrow(signature) # fraction of genes for the up part of the signature
  frac_downgenes <- length(down_genes)/nrow(signature) # fraction of genes for the down part of the signature

  # weight ssGSEA results based on the proportion of genes from each direction; calculate difference - if negative, means that the
  # signature is doing the opposite that it should do (either DOWN is positive, UP is negative; both...)
  ssGSEAresults$score <- (ssGSEAresults$score_up * frac_upgenes) - (ssGSEAresults$score_down * frac_downgenes)
  ssGSEAresults$score_up <- NULL
  ssGSEAresults$score_down <- NULL


  colnames(ssGSEAresults) <- c("sample", "score")


  return(ssGSEAresults)
}

