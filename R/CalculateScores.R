#' Calculate Gene Signature Scores using Score-Based Approaches
#'
#' This function calculates a gene signature score for each sample based on one
#' or more predefined gene sets (signatures).
#'
#' \describe{
#' This function calculates a gene signature score for each sample based on one or
#' more predefined gene sets (signatures). Four methods are available:
#'
#'   \item{\code{ssGSEA}}{
#'     Uses the single-sample Gene Set Enrichment Analysis (ssGSEA) method
#'     to compute an enrichment score for each signature in each sample.
#'     This method uses an adaptation from the the \code{gsva()} function from
#'     the \code{GSVA} package to compute an enrichment score, representing the
#'     absolute enrichment of each gene set in each sample.
#'   }
#'   \item{\code{logmedian}}{
#'     Computes, for each sample, the score as the sum of the normalized
#'     (log2-median-centered) expression values of the signature genes divided
#'     by the number of genes in the signature.
#'   }
#'   \item{\code{ranking}}{
#'     Computes gene signature scores for each sample by ranking the expression
#'     of signature genes in the dataset and normalizing the score based on the
#'     total number of genes.
#'   }
#'   \item{\code{all}}{
#'     Computes gene signature scores using all three methods
#'     (\code{ssGSEA}, \code{logmedian}, and \code{ranking}).
#'     The function returns a list containing the results of each method.
#'   }
#' }
#'
#' @param data A data frame of normalized (non-transformed) counts where each
#'   row is a gene and each column is a sample. The row names should contain
#'   gene names, and the column names should contain sample identifiers.
#'   **(Required)**
#' @param metadata A data frame describing the attributes of each sample. Each
#'   row corresponds to a sample and each column to an attribute. The first
#'   column of \code{metadata} should be the sample identifiers (i.e., the
#'   column names of \code{data}). Defaults to \code{NULL} if no metadata is
#'   provided.
#' @param gene_sets Gene set input. **(Required)**
#'
#'   If using **unidirectional** gene sets, provide a named list where each
#'   element is a vector of gene names representing a gene signature. The names
#'   of the list elements should correspond to the labels for each signature.
#'
#'   If using **bidirectional** gene sets, provide a named list where each
#'   element is a data frame. The names of the list elements should correspond
#'   to the labels for each signature, and each data frame should contain the
#'   following structure:
#' - The **first column** should contain gene names.
#' - The **second column** should indicate the expected direction of enrichment
#' (1 for upregulated genes, -1 for downregulated genes).
#'
#' @param method A character string indicating the scoring method to use.
#'   Options are \code{"ssGSEA"}, \code{"logmedian"}, \code{"ranking"}, or
#'   \code{"all"} (to compute scores using all methods). Defaults to
#'   \code{"logmedian"}.
#'
#' @return If a single method is chosen, a data frame containing the calculated
#'   scores for each gene signature, including metadata if provided. If
#'   \code{method = "all"}, a list is returned where each element corresponds to
#'   a scoring method and contains the respective data frame of scores.
#'
#' \describe{
#'   \item{sample}{The sample identifier (matching the column names of the input data).}
#'   \item{score}{The calculated gene signature score for the corresponding sample.}
#'   \item{(metadata)}{Any additional columns from the \code{metadata} data frame
#'   provided by the user, if available.}
#' }
#'
#' @examples
#' # Simulate positive gene expression data (genes as rows, samples as columns)
#' set.seed(42)
#' expr <- as.data.frame(matrix(rexp(60, rate = 0.2), nrow = 6, ncol = 10))
#' rownames(expr) <- paste0("Gene", 1:6)
#' colnames(expr) <- paste0("Sample", 1:10)
#'
#' # Simulate metadata for samples
#' metadata <- data.frame(
#'   sample = colnames(expr),
#'   Group = rep(c("A", "B"), each = 5)
#' )
#'
#' # Define two simple gene sets
#' gene_sets <- list(
#'   Signature1 = c("Gene1", "Gene2", "Gene3"),
#'   Signature2 = c("Gene4", "Gene5", "Gene6")
#' )
#'
#' # Calculate logmedian scores
#' scores_logmedian <- CalculateScores(
#'   data = expr,
#'   metadata = metadata,
#'   gene_sets = gene_sets,
#'   method = "logmedian"
#' )
#' head(scores_logmedian)
#'
#' # Calculate all score types
#' scores_all <- CalculateScores(
#'   data = expr,
#'   metadata = metadata,
#'   gene_sets = gene_sets,
#'   method = "all"
#' )
#' lapply(scores_all, head)
#'
#' @export
CalculateScores <- function(data, metadata, gene_sets,
                            method = c("ssGSEA", "logmedian","ranking", "all")) {
  data <- as.data.frame(data) # Ensure data is a data frame
  method <- match.arg(method)  # Validate method input

  # if (!is.data.frame(data)) stop("Error: data must be a data-frame")
  # if (!is.null(metadata) && !is.data.frame(metadata)) stop(
  #   "Error: metadata must be a data-frame")
  # if (!is.list(gene_sets)) stop("Error: gene_sets must be a list")
  
  if (!is.data.frame(data) || (!is.null(metadata) && !is.data.frame(metadata)) || !is.list(gene_sets)) {
    stop("Invalid input: 'data' must be a data-frame, 'metadata' (if provided) must be a data-frame, and 'gene_sets' must be a list.")
  }
  

  # Change first column name to default name "sample", for merging purposes
  if (!is.null(metadata)) colnames(metadata)[1] <- "sample"


  # Call appropriate function based on method
  if (method == "ssGSEA") {
    return(CalculateScores_ssGSEA(data, metadata, gene_sets))
  } else if (method == "logmedian") {
    return(CalculateScores_logmedian(data, metadata, gene_sets))
  } else if (method == "ranking"){
    return(CalculateScores_Ranking(data, metadata, gene_sets))
  } else if (method == "all"){
    return(list(ssGSEA=CalculateScores_ssGSEA(data, metadata, gene_sets),
                logmedian=CalculateScores_logmedian(data, metadata, gene_sets),
                ranking=CalculateScores_Ranking(data, metadata, gene_sets)))

  }

}
