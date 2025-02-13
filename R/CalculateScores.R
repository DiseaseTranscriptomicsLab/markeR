#' Calculate Gene Signature Scores using Score-Based Approaches
#'
#' This function calculates a gene signature score for each sample based on one or more predefined gene sets
#' (signatures). Two methods are available:
#'
#' \describe{
#'   \item{\code{ssGSEA}}{
#'     Uses the single-sample Gene Set Enrichment Analysis (ssGSEA) method from the
#'     \code{GSVA} package to compute an enrichment score for each signature in each sample.
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
#' corresponds to one signature (named accordingly) and is a data frame with the following columns:
#' \describe{
#'   \item{sample}{The sample identifier (matching the column names of the input data).}
#'   \item{score}{The calculated gene signature score for the corresponding sample.}
#'   \item{<metadata>}{Any additional columns from the \code{metadata} data frame provided by the user, if available.}
#' }
#'
#' @details
#' \strong{ssGSEA:} This method uses the \code{gsva()} function from the \code{GSVA} package to compute an enrichment score,
#' representing the absolute enrichment of each gene set in each sample.
#'
#' \strong{logmedian:} This method log2-transforms the expression values, centers each gene at its median,
#' and calculates the average of these centered values for the genes in each signature.
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

  # Compares the provided method value against the allowed options (c("ssGSEA", "logmedian")).
  # If the user does not explicitly provide a value, it will choose the default (the first element in the allowed values, if not overridden).
  method <- match.arg(method)

  # Ensure data is in the proper format
  if (!is.list(data)) stop("Error: data must be a data-frame")
  if (!is.list(metadata)) stop("Error: metadata must be a data-frame")
  if (!is.list(gene_sets)) stop("Error: gene_sets must be a list")

  # Create new variable with sample name, will be used for merging purposes
  metadata$sample <- row.names(metadata)

  # Initialise results list
  ResultsList <- list()

  if(method == "ssGSEA"){

    # For each signature in gene_sets, compute ssGSEA scores.
    for (sig in names(gene_sets)){

      # create a list per signature, so that the results are one data frame per signature
      siglist <- list(gene_sets[[sig]])
      names(siglist) <- c(sig)

      mtx <- as.matrix(log2(data))
      ssgsea_results <- gsva(expr=mtx,
                             gset.idx.list=siglist,
                             method="ssgsea",
                             kcdf="Gaussian") #suitable when input expression values are continuous, such as RNA-seq log-CPMs, log-RPKMs or log-TPMs.

      # Formatting results to a data-frame
      ssgsea_results <- as.data.frame(ssgsea_results)
      ssgsea_results$signature <- row.names(ssgsea_results)
      ssgsea_results <- melt(ssgsea_results)
      colnames(ssgsea_results) <- c("signature","sample","score")
      ssgsea_results$signature <- NULL # remove, resulting list will have one entry per signature

      # merge metadata with results, if metadata is available
      if(!is.null(metadata)) ssgsea_results <- merge(ssgsea_results,metadata, by="sample")
      # remove redundant information
      row.names(ssgsea_results) <- NULL
      # create new entry for each signature
      ResultsList[[sig]] <- ssgsea_results

    }



  } else if(method == "logmedian"){

    for (sig in names(gene_sets)){

      signature <- gene_sets[[sig]]
      data_subset <- na.omit(subset(log2(data+1), row.names(log2(data+1)) %in% signature)) # log(data + 1) to not have Inf values
      data_subset <- data_subset-apply(data_subset, 1, median) # center gene in its log2 median

      dfScore <- colSums(data_subset) # sum based on its enrichment (if depleted, -1 of the expression; otherwise, sum; only works because expression in zero is zero - log2(x+1))
      dfScore <- dfScore/nrow(data_subset) # vector
      dfScore <- data.frame(sample=names(dfScore), # create data frame
                            score=dfScore)

      if(!is.null(metadata)){
        # merge metadata with results, if metadata is available
        dfScore <- merge(dfScore, metadata, by="sample")
      }

      # remove redundant information
      row.names(dfScore) <- NULL
      # create new entry for each signature
      ResultsList[[sig]] <- dfScore

    }


  }

  return(ResultsList)
}
