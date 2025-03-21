#' Run Gene Set Enrichment Analysis (GSEA) for Multiple Contrasts
#'
#' This function performs GSEA using `fgsea` for each contrast in a list of differential expression results.
#' It automatically determines the appropriate ranking statistic based on the gene set format unless specified by the user.
#'
#' @param DEGList A named list where each element represents a contrast and contains a data frame of differential expression results.
#'   - Each data frame must include at least the `"t"` statistic and the `"B"` statistic for each gene.
#'   - Row names should correspond to gene identifiers.
#'
#' @param gene_sets A named list where each element represents a gene set. Each gene set can be:
#'   - A **vector of gene names** (for unidirectional gene sets).
#'   - A **data frame** with two columns:
#'     - Column 1: Gene names.
#'     - Column 2: Expected direction (`1` for upregulated genes, `-1` for downregulated genes).
#'
#' @param stat Optional. The statistic to use for ranking genes before GSEA. If `NULL`, it is automatically determined based on the gene set:
#'   - `"B"` for gene sets with **no known direction** (vectors).
#'   - `"t"` for **unidirectional** or **bidirectional** gene sets (data frames).
#'   - If provided, this argument overrides the automatic selection.
#'
#' @return A named list where each element corresponds to a contrast. Each contrast contains a **single data frame** with GSEA results for all gene sets.
#'   The result includes the standard `fgsea` output plus two additional columns:
#'   - `pathway`: The name of the gene set.
#'   - `stat_used`: The statistic used for ranking genes in that analysis (`"t"` or `"B"`).
#'
#' @examples
#' # Example input data
#' DEGList <- list(
#'   Contrast1 = data.frame(t = rnorm(100), B = rnorm(100), row.names = paste0("Gene", 1:100)),
#'   Contrast2 = data.frame(t = rnorm(100), B = rnorm(100), row.names = paste0("Gene", 1:100))
#' )
#'
#' gene_sets <- list(
#'   UnidirectionalSet = c("Gene1", "Gene5", "Gene20"),
#'   BidirectionalSet = data.frame(Gene = c("Gene2", "Gene10", "Gene15"), Direction = c(1, -1, 1))
#' )
#'
#' results <- runGSEA(DEGList, gene_sets)
#'
#' @importFrom fgsea fgsea
#' @export
runGSEA <- function(DEGList, gene_sets, stat = NULL) {

  # Initialize storage for results across contrasts
  results_by_contrast <- list()

  # Loop over each contrast in DEGList
  for (contrast in names(DEGList)) {
    deg_df <- DEGList[[contrast]]

    # Initialize list to store results for all gene sets in this contrast
    contrast_results <- list()

    for (set_name in names(gene_sets)) {
      gs <- gene_sets[[set_name]]

      # Automatically determine which statistic to use unless specified
      current_stat <- if (is.null(stat)) {
        if (is.data.frame(gs)) "t" else "B"
      } else {
        stat
      }

      # Create the ranking vector for GSEA
      ranks <- setNames(deg_df[[current_stat]], rownames(deg_df))

      if (is.data.frame(gs)) {
        # Bidirectional gene set: modify rankings based on expected direction
        gs_genes <- as.character(gs[[1]])
        directions <- as.numeric(gs[[2]])

        # Ensure that all genes in gs_genes are present in ranks
        matched_genes <- gs_genes[gs_genes %in% names(ranks)]
        matched_directions <- directions[gs_genes %in% names(ranks)]

        # Multiply ranks by direction
        ranks_adjusted <- ranks
        ranks_adjusted[matched_genes] <- ranks_adjusted[matched_genes] * matched_directions

        # Sort after modifying ranks
        ranks_sorted <- sort(ranks_adjusted, decreasing = TRUE)

        set.seed("20032025")
        fgseaRes <- fgsea::fgsea(pathways = list(temp = matched_genes),
                                 stats = ranks_sorted)

      } else if (is.vector(gs)) {
        # Unidirectional gene set: use the ranking vector directly
        gs_genes <- as.character(gs)

        # Sort ranks before running fgsea
        ranks_sorted <- sort(ranks, decreasing = TRUE)
        set.seed("20032025")
        fgseaRes <- fgsea::fgsea(pathways = list(temp = gs_genes),
                                 stats = ranks_sorted)

      } else {
        warning("Gene set '", set_name, "' is not in a recognized format. Skipping.")
        next
      }

      # Add gene set name and ranking method used
      fgseaRes$pathway <- set_name
      fgseaRes$stat_used <- current_stat

      contrast_results[[set_name]] <- fgseaRes
    }

    # Combine all gene sets into a single data frame for this contrast
    results_by_contrast[[contrast]] <- do.call(rbind, contrast_results)
  }

  return(results_by_contrast)
}
