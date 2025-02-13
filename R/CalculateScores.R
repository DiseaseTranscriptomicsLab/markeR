#' Calculate Gene Signature Scores using Score-Based Approaches
#'
#' This function calculates a gene signature score for each sample based on one or more predefined gene sets
#' (signatures).
#'
#' \describe{
#' This function calculates a gene signature score for each sample based on one or more predefined gene sets
#' (signatures). Two methods are available:
#'
#'   \item{\code{ssGSEA}}{
#'     Uses the single-sample Gene Set Enrichment Analysis (ssGSEA) method to compute an enrichment score
#'     for each signature in each sample. This method uses the \code{gsva()} function from the \code{GSVA} package
#'     to compute an enrichment score, representing the absolute enrichment of each gene set in each sample.
#'   }
#'   \item{\code{logmedian}}{
#'     Computes, for each sample, the score as the sum of the normalized (log2-median-centered)
#'     expression values of the signature genes divided by the number of genes in the signature.
#'   }
#' }
#'
#' @param data A data frame of normalised (non-transformed) counts where each row is a gene and each column is a sample. The row names should contain gene names
#'   and the column names should contain sample identifiers. **(Required)**
#' @param metadata A data frame describing the attributes of each sample. Each row corresponds to a sample and each column to an attribute.
#'   The row names of \code{metadata} must match the sample identifiers (i.e., the column names of \code{data}).
#'   Defaults to \code{NULL} if no metadata is provided.
#' @param gene_sets A named list where each element is a vector of gene names corresponding to a gene signature.
#'   The name of each list element should be the label for that signature. **(Required)**
#' @param method A character string indicating the scoring method to use. Options are \code{"ssGSEA"}
#'   or \code{"logmedian"}. Defaults to \code{"logmedian"}.
#'
#' @return A list containing the calculated scores for each gene signature. Each element of the list
#' corresponds to one signature (named accordingly) and is a data frame with the following attributes:
#' \describe{
#'   \item{sample}{The sample identifier (matching the column names of the input data).}
#'   \item{score}{The calculated gene signature score for the corresponding sample.}
#'   \item{(metadata)}{Any additional columns from the \code{metadata} data frame provided by the user, if available.}
#' }
#'
#'
#' @examples
#' \dontrun{
#'   # Assume 'gene_data' is your gene expression data frame and 'sample_metadata'
#'   # is your metadata. Define a list of gene signatures as follows:
#'   gene_sets <- list(
#'     "Signature_A" = c("Gene1", "Gene5", "Gene10", "Gene20"),
#'     "Signature_B" = c("Gene2", "Gene6", "Gene15", "Gene30")
#'   )
#'
#'   # Using the ssGSEA method:
#'   scores_ssgsea <- calculate_signature_score(data = gene_data,
#'                                              metadata = sample_metadata,
#'                                              gene_sets = gene_sets,
#'                                              method = "ssGSEA")
#'
#'   # Using the logmedian method (default):
#'   scores_logmedian <- calculate_signature_score(data = gene_data,
#'                                                 gene_sets = gene_sets)
#' }
#'
#' @export
CalculateScores <- function(data, metadata, gene_sets, method = c("ssGSEA", "logmedian")) {
  method <- match.arg(method)  # Validate method input

  if (!is.data.frame(data)) stop("Error: data must be a data-frame")
  if (!is.null(metadata) && !is.data.frame(metadata)) stop("Error: metadata must be a data-frame")
  if (!is.list(gene_sets)) stop("Error: gene_sets must be a list")

  # Call appropriate function based on method
  if (method == "ssGSEA") {
    return(CalculateScores_ssGSEA(data, metadata, gene_sets))
  } else if (method == "logmedian") {
    return(CalculateScores_logmedian(data, metadata, gene_sets))
  }
}
