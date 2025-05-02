#' Alternative Implementation of Single-Sample Gene Set Enrichment Analysis (ssGSEA)
#'
#' This function computes an enrichment score for each sample using an alternative single-sample Gene Set Enrichment Analysis (ssGSEA) method.
#' It first maps gene sets to the gene indices present in the expression matrix, then ranks the genes for each sample, and finally calculates
#' a weighted enrichment score based on the cumulative differences between in-set and out-of-set gene ranks.
#' Source: https://rpubs.com/pranali018/SSGSEA
#'
#' @param X A numeric matrix of gene expression values with rows representing genes and columns representing samples. Row names should correspond to gene identifiers.
#' @param gene_sets A list of gene sets, where each element is a vector of gene identifiers. The function will match these identifiers with the row names of \code{X}.
#' @param alpha A numeric value specifying the exponent used to weight the ranking scores. Default is \code{0.25}.
#' @param scale Logical; if \code{TRUE}, the cumulative difference is normalized by the total number of genes. Default is \code{TRUE}.
#' @param norm Logical; if \code{TRUE}, the enrichment scores are further normalized by the absolute difference between the maximum and minimum scores. Default is \code{FALSE}.
#' @param single Logical; if \code{TRUE}, the function returns the sum of the cumulative differences as the enrichment score. If \code{FALSE}, the maximum absolute cumulative difference is used. Default is \code{TRUE}.
#'
#' @return A matrix of enrichment scores with rows corresponding to gene sets and columns corresponding to samples.
#'
#' @details The function performs the following steps:
#' \enumerate{
#'   \item Maps each gene set to the indices of genes in \code{X} by matching gene identifiers.
#'   \item Computes column-wise rankings for the gene expression matrix using a ranking method (via the \code{colRanking} function) with tie resolution set to \code{'average'}.
#'   \item For each sample, orders the gene ranks in decreasing order.
#'   \item For each gene set in the sample, calculates:
#'     \itemize{
#'       \item The weighted contribution (\code{rank_alpha}) for genes in the set raised to the power of \code{alpha}.
#'       \item The cumulative distribution functions (CDFs) for genes within the gene set (\code{step_cdf_pos}) and those not in the gene set (\code{step_cdf_neg}).
#'       \item The difference between these CDFs, optionally scaled by the number of genes if \code{scale = TRUE}.
#'       \item Depending on the \code{single} parameter, either the sum of the differences (if \code{TRUE}) or the maximum absolute difference (if \code{FALSE}) is used as the enrichment score for that gene set.
#'     }
#'   \item Optionally normalizes the final enrichment scores by the range of values if \code{norm = TRUE}.
#' }
#'
#' @examples
#' \dontrun{
#'   # Create a sample gene expression matrix:
#'   X <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#'   rownames(X) <- paste0("gene", 1:100)
#'
#'   # Define example gene sets:
#'   gene_sets <- list(
#'     set1 = sample(rownames(X), 10),
#'     set2 = sample(rownames(X), 15)
#'   )
#'
#'   # Compute the ssGSEA enrichment scores:
#'   es <- ssGSEA_alternative(X, gene_sets, alpha = 0.25, scale = TRUE, norm = FALSE, single = TRUE)
#'   print(es)
#' }
#'
#' @keywords internal
ssGSEA_alternative <- function(X, gene_sets, alpha = 0.25, scale = T, norm = F, single = T) {
  row_names = rownames(X)
  num_genes = nrow(X)
  gene_sets = lapply(gene_sets, function(genes) { which(row_names %in% genes) })

  # Ranks for genes
  R = colRanking(X, ties.method = 'average')

  # Calculate enrichment score (es) for each sample (column)
  es = apply(R, 2, function(R_col) {
    gene_ranks = order(R_col, decreasing = TRUE)

    # Calc es for each gene set
    es_sample = sapply(gene_sets, function(gene_set_idx) {
      # pos: match (within the gene set)
      # neg: non-match (outside the gene set)
      indicator_pos = gene_ranks %in% gene_set_idx
      indicator_neg = !indicator_pos

      rank_alpha  = (R_col[gene_ranks] * indicator_pos) ^ alpha

      step_cdf_pos = cumsum(rank_alpha)    / sum(rank_alpha)
      step_cdf_neg = cumsum(indicator_neg) / sum(indicator_neg)

      step_cdf_diff = step_cdf_pos - step_cdf_neg

      # Normalize by gene number
      if (scale) step_cdf_diff = step_cdf_diff / num_genes

      # Use ssGSEA or not
      if (single) {
        sum(step_cdf_diff)
      } else {
        step_cdf_diff[which.max(abs(step_cdf_diff))]
      }
    })
    unlist(es_sample)
  })

  if (length(gene_sets) == 1) es = matrix(es, nrow = 1)

  # Normalize by absolute diff between max and min
  if (norm) es = es / diff(range(es))

  # Prepare output
  rownames(es) = names(gene_sets)
  colnames(es) = colnames(X)
  return(es)
}


#' Compute Independent Column-wise Ranks of Matrix Elements
#'
#' This function computes the rank of each element in every column of a numeric matrix independently.
#' For each column, the smallest element receives a rank of 1, the second smallest a rank of 2, and so on.
#'
#' @param x A numeric matrix.
#' @param ties.method A character string specifying the method used for tie-breaking.
#' Options include \code{"average"}, \code{"first"}, \code{"random"}, \code{"max"}, or \code{"min"}.
#' The default is \code{"average"}.
#'
#' @return A numeric matrix of the same dimensions as \code{x} where each column contains the ranks of the corresponding column's elements.
#'
#' @keywords internal
colRanking <- function(x, ties.method = "average") {
  apply(x, 2, rank, ties.method = ties.method)
}
