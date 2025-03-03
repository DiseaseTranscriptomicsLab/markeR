#' Calculate Gene Signature Scores using Ranking Approach
#'
#' Computes gene signature scores for each sample by ranking the expression of signature genes
#' in the dataset and normalizing the score based on the total number of genes.
#'
#' @param data A data frame where rows represent genes, columns represent samples, and values
#'   correspond to gene expression levels. **(Required)**
#' @param metadata A data frame containing sample metadata. The first column must contain sample
#'   names. **(Optional)**
#' @param gene_sets A named list of gene sets. **(Required)**
#'
#' - **Unidirectional gene sets**: Provide a named list where each element is a vector of gene names.
#' - **Bidirectional gene sets**: Provide a named list where each element is a data frame.
#'   - The **first column** should contain gene names.
#'   - The **second column** should indicate the expected direction of enrichment (1 for upregulated genes, -1 for downregulated genes).
#'
#' @return A named list of data frames, where each data frame contains:
#' - `sample`: Sample name.
#' - `score`: Normalized ranking score for the given gene signature.
#' - Additional metadata columns (if `metadata` is provided).
#'
#' @details
#' - The function first validates inputs and extracts relevant genes from the dataset.
#' - For **unidirectional** signatures, it computes rankings based on gene expression levels.
#' - For **bidirectional** signatures, it computes separate rankings for upregulated and downregulated genes,
#'   then calculates a final score by subtracting downregulated rankings from upregulated rankings.
#' - The final scores are normalized by dividing by the total number of genes.
#' - This metric is not suitable to compare absolute values between different gene sets, i.e. should be used only for
#'   relative comparisons between samples when using the same gene set.
#'
#' @examples
#' # Example dataset with 5 genes (rows) and 3 samples (columns)
#' set.seed(123)
#' data <- matrix(runif(15, 1, 100), nrow = 5, ncol = 3)
#' rownames(data) <- paste0("Gene_", 1:5)
#' colnames(data) <- paste0("Sample_", 1:3)
#'
#' # Unidirectional gene set example
#' gene_sets <- list(Signature1 = c("Gene_1", "Gene_3", "Gene_5"))
#'
#' # Compute scores
#' scores <- CalculateScores_Ranking(data, gene_sets = gene_sets)
#' print(scores)
#'
#' @export
CalculateScores_Ranking <- function(data, metadata = NULL, gene_sets) {

  ResultsList <- list()

  if (!is.data.frame(data)) stop("Error: data must be a data frame")
  if (!is.null(metadata) && !is.data.frame(metadata)) stop("Error: metadata must be a data frame")
  if (!is.list(gene_sets)) stop("Error: gene_sets must be a list")

  # Change first column name to default name "sample" for merging purposes
  if (!is.null(metadata)) colnames(metadata)[1] <- "sample"

  # Define universe of genes
  universe_genes <- row.names(data)

  for (sig in names(gene_sets)) {

    signature <- gene_sets[[sig]]

    if (is.data.frame(signature)) {  # If a data frame, check enrichment direction

      nb_factors <- length(unique(signature[, 2]))

      if (nb_factors > 1) {  # Both up (1) and down (-1) genes

        message(paste0("Considering bidirectional gene signature mode for signature ", sig))

        signaturegenes_up <- signature[signature[,2] == 1, 1]
        signaturegenes_down <- signature[signature[,2] == -1, 1]

        # Apply getRanking function to each sample (column)
        rankings_up <- sapply(colnames(data), function(sample) getRanking(data, sample, signaturegenes_up))
        rankings_down <- sapply(colnames(data), function(sample) getRanking(data, sample, signaturegenes_down))

        ranking_final <- (rankings_up - rankings_down) / length(universe_genes)
        ranking_final <- data.frame(sample = colnames(data), score = ranking_final)

      } else {  # Only one known direction (up or down)

        message(paste0("Considering unidirectional gene signature mode for signature ", sig))

        signaturegenes <- signature[, 1]

        # Apply getRanking function to each sample (column)
        rankings <- sapply(colnames(data), function(sample) getRanking(data, sample, signaturegenes))

        ranking_final <- rankings / length(universe_genes)
        ranking_final <- data.frame(sample = colnames(data), score = ranking_final)

      }
    } else {  # If vector of genes (unidirectional)
      message(paste0("Considering unidirectional gene signature mode for signature ", sig))

      signaturegenes <- signature

      # Apply getRanking function to each sample (column)
      rankings <- sapply(colnames(data), function(sample) getRanking(data, sample, signaturegenes))

      ranking_final <- rankings / length(universe_genes)
      ranking_final <- data.frame(sample = colnames(data), score = ranking_final)

    }

    if (!is.null(metadata)) ranking_final <- merge(ranking_final, metadata, by = "sample")

    row.names(ranking_final) <- NULL

    ResultsList[[sig]] <- ranking_final
  }

  return(ResultsList)
}

#' Get Gene Expression Ranking
#'
#' Computes the rank sum of a given gene set within a sample based on its expression level.
#'
#' @param data A data frame where rows represent genes, columns represent samples, and values correspond to expression levels.
#' @param sample A character string specifying the sample name (column in `data`).
#' @param geneset A vector of gene names to be ranked.
#'
#' @return The sum of the ranks of the genes found in the sample.
#' @details
#' - The function orders gene expression levels from lowest to highest.
#' - It then determines the rank of each gene in `geneset` and returns the sum of these ranks.
#' - If some genes are missing, they are omitted from the ranking calculation.
#'
#' @examples
#' # Example dataset with 5 genes and 3 samples
#' set.seed(123)
#' data <- matrix(runif(15, 1, 100), nrow = 5, ncol = 3)
#' rownames(data) <- paste0("Gene_", 1:5)
#' colnames(data) <- paste0("Sample_", 1:3)
#'
#' # Define gene set
#' geneset <- c("Gene_1", "Gene_3", "Gene_5")
#'
#' # Compute ranking for Sample_1
#' rank_score <- getRanking(data, "Sample_1", geneset)
#' print(rank_score)
#'
#' @export
getRanking <- function(data, sample, geneset) {

  expressiongene <- data[, sample]  # Isolate one sample and get the expression of all genes
  names(expressiongene) <- row.names(data)  # Name vector
  expressiongene <- expressiongene[order(expressiongene, decreasing = FALSE)]  # Order from least to most expressed
  ranking <- match(geneset, names(expressiongene))  # Find gene positions in ordered list
  ranking <- as.vector(na.omit(ranking))  # Remove missing genes

  return(sum(ranking))  # Return sum of ranks
}
