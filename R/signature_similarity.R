#' Plot Signature Similarity via Jaccard Index
#'
#' Visualises the similarity between user-defined gene signatures and either
#' other user-defined signatures or MSigDB gene sets, using the Jaccard index.
#' Produces a heatmap of pairwise Jaccard indices between reference and comparison
#' signatures.
#'
#' @param signatures A named list of character vectors representing reference gene signatures.
#' @param other_user_signatures Optional. A named list of character vectors representing other user-defined signatures to compare against. If `NULL`, MSigDB gene sets (if provided) will be used.
#' @param collection Optional. MSigDB collection name (e.g., `"H"` for hallmark, `"C2"` for curated gene sets).
#' @param subcollection Optional. Subcategory within an MSigDB collection (e.g., `"CP:KEGG"`).
#' @param limits Numeric vector of length 2. Limits for the Jaccard index color scale. Default is `c(0, 1)`.
#' @param title_size Integer specifying the font size for the plot title. Default is `12`.
#' @param color_values Character vector of colors used in the Jaccard fill gradient. Default is `c("#F9F4AE", "#B44141")`.
#' @param title Optional. Custom title for the plot. If `NULL`, the title defaults to `"Signature Overlap"`.
#' @param num_sigs_toplot Optional. Integer. Maximum number of comparison signatures (including user and MSigDB) to display.
#' @param jaccard_threshold Numeric. Minimum Jaccard index required for an MSigDB gene set to be included in the plot. Default is `0`.
#' @param msig_subset Optional. Character vector of MSigDB gene set names to subset from the specified collection. Useful to restrict analysis to a specific set of pathways. If supplied, other filters will apply only to this subset.
#'
#' @return A `ggplot` object showing a heatmap of Jaccard similarities between reference and comparison signatures.
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom purrr map_dfr
#' @importFrom tibble tibble
#' @importFrom msigdbr msigdbr
#' @importFrom scales squish
#'
#' @examples
#' sig1 <- list(A = c("TP53", "BRCA1", "EGFR"))
#' sig2 <- list(B = c("TP53", "MYC", "EGFR"), C = c("GATA3", "STAT3"))
#'
#' signature_similarity(
#'   signatures = sig1,
#'   other_user_signatures = sig2
#' )
#'
#' signature_similarity(
#'   signatures = sig1,
#'   collection = "H",
#'   jaccard_threshold = 0.1,
#'   num_sigs_toplot = 5
#' )
#'
#' @export


signature_similarity <- function(
  signatures,
  other_user_signatures = NULL,
  collection = NULL,
  subcollection = NULL,
  limits = c(0, 1),
  title_size = 12,
  color_values = c("#F9F4AE", "#B44141"),
  title = NULL,
  num_sigs_toplot = NULL,
  jaccard_threshold = 0,
  msig_subset = NULL
) {

  # Convert gene names to uppercase for consistency
  signatures <- lapply(signatures, toupper)
  if (!is.null(other_user_signatures)) {
    other_user_signatures <- lapply(other_user_signatures, toupper)
  }

  # Initialize MSigDB sets
  gsets_filtered <- list()

  # If MSigDB collection is specified
  if (!is.null(collection)) {
    gs <- msigdbr::msigdbr(
      species = "Homo sapiens",
      category = collection,
      subcategory = subcollection
    )

    gsets <- base::split(gs$gene_symbol, gs$gs_name)
    gsets <- lapply(gsets, toupper)

    # Warn if msig_subset has pathways not in MSigDB
    if (!is.null(msig_subset)) {
      missing <- base::setdiff(msig_subset, base::names(gsets))
      if (base::length(missing) > 0) {
        base::message("The following pathways from msig_subset were not found in the MSigDB collection:\n", base::paste(missing, collapse = ", "))
      }
      # Keep only valid subset pathways
      gsets <- gsets[base::names(gsets) %in% msig_subset]
    }

    # Compute Jaccard indices
    jaccard_msigdb <- purrr::map_dfr(
      base::names(signatures),
      function(ref_name) {
        purrr::map_dfr(
          base::names(gsets),
          function(comp_name) {
            sig1 <- signatures[[ref_name]]
            sig2 <- gsets[[comp_name]]
            jaccard <- base::length(base::intersect(sig1, sig2)) / base::length(base::union(sig1, sig2))
            tibble::tibble(
              Reference_Signature = ref_name,
              Compared_Signature = comp_name,
              Jaccard = jaccard
            )
          }
        )
      }
    )

    # Rank by max Jaccard
    msigdb_ranks <- jaccard_msigdb %>%
      dplyr::group_by(Compared_Signature) %>%
      dplyr::summarise(Max_Jaccard = max(Jaccard)) %>%
      dplyr::arrange(dplyr::desc(Max_Jaccard))

    # Filter by threshold
    msigdb_ranks <- msigdb_ranks %>%
      dplyr::filter(Max_Jaccard >= jaccard_threshold)

    if (nrow(msigdb_ranks) == 0) {
      base::message("No MSigDB pathways passed the Jaccard threshold of ", jaccard_threshold, ". Only user-defined signatures will be plotted.")
    } else {
      # Determine number to plot
      num_custom <- if (!is.null(other_user_signatures)) base::length(other_user_signatures) else 0
      num_remaining <- if (!is.null(num_sigs_toplot)) base::max(num_sigs_toplot - num_custom, 0) else nrow(msigdb_ranks)
      top_msigdb_names <- utils::head(msigdb_ranks$Compared_Signature, num_remaining)
      gsets_filtered <- gsets[top_msigdb_names]
    }
  }

  # Combine user-defined and MSigDB-filtered
  if (is.null(other_user_signatures)) {
    other_user_signatures <- gsets_filtered
  } else {
    other_user_signatures <- base::c(other_user_signatures, gsets_filtered)
  }

  # Compute Jaccard for plot
  jaccard_df <- purrr::map_dfr(
    base::names(signatures),
    function(ref_name) {
      purrr::map_dfr(
        base::names(other_user_signatures),
        function(comp_name) {
          sig1 <- signatures[[ref_name]]
          sig2 <- other_user_signatures[[comp_name]]
          jaccard <- base::length(base::intersect(sig1, sig2)) / base::length(base::union(sig1, sig2))
          tibble::tibble(
            Reference_Signature = ref_name,
            Compared_Signature = comp_name,
            Jaccard = jaccard
          )
        }
      )
    }
  )

  # Fallback title
  plot_title <- if (is.null(title)) "Signature Overlap" else title

  # Plot
  ggplot2::ggplot(jaccard_df, ggplot2::aes(x = Reference_Signature, y = Compared_Signature, fill = Jaccard)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = base::sprintf("%.2f", Jaccard)), color = "black") +
    ggplot2::scale_fill_gradientn(colors = color_values, limits = c(0, 1), oob = scales::squish) +
    ggplot2::labs(
      x = "",
      y = "Compared Signature",
      fill = "Jaccard Index",
      title = plot_title
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(hjust = 0.5, size = title_size)
    )
}
