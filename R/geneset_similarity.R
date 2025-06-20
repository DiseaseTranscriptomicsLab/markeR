#' Plot Signature Similarity via Jaccard Index or Fisher's Odds Ratio
#'
#' Visualizes similarity between user-defined gene signatures and either
#' other user-defined signatures or MSigDB gene sets, using either the Jaccard index
#' or Fisher's Odds Ratio. Produces a heatmap of pairwise similarity metrics.
#'
#' @param signatures A named list of character vectors representing reference gene signatures.
#' @param other_user_signatures Optional. A named list of character vectors representing other user-defined signatures to compare against.
#' @param collection Optional. MSigDB collection name (e.g., `"H"` for hallmark, `"C2"` for curated gene sets). Use msigdbr::msigdbr_collections() for the available options.
#' @param subcollection Optional. Subcategory within an MSigDB collection (e.g., `"CP:REACTOME"`). Use msigdbr::msigdbr_collections() for the available options.
#' @param metric Character. Either "jaccard" or "odds_ratio".
#' @param universe Character vector. Background gene universe. Required for odds ratio.
#' @param or_threshold Numeric. Minimum Odds Ratio required for a gene set to be included in the plot. Default is 1.
#' @param pval_threshold Numeric. Maximum adjusted p-value to show a label. Default is 0.05.
#' @param limits Numeric vector of length 2. Limits for color scale.
#' @param title_size Integer specifying the font size for the plot title. Default is `12`.
#' @param color_values Character vector of colors used for the fill gradient. Default is `c("#F9F4AE", "#B44141")`.
#' @param title Optional. Custom title for the plot. If `NULL`, the title defaults to `"Signature Overlap"`.
#' @param num_sigs_toplot Optional. Integer. Maximum number of comparison signatures (including user and MSigDB) to display.
#' @param jaccard_threshold Numeric. Minimum Jaccard index required for a gene set to be included in the plot. Default is `0`.
#' @param msig_subset Optional. Character vector of MSigDB gene set names to subset from the specified collection. Useful to restrict analysis to a specific set of pathways.
#'                    If supplied, other filters will apply only to this subset. Use "collection = "all" to mix gene sets from different collections.
#' @param width_text Integer. Character wrap width for labels.
#' @param na_color Character. Color for NA values in the heatmap. Default is `"grey90"`.
#'
#' @return A ggplot heatmap object.
#'
#' @import ggplot2
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
#' @export
geneset_similarity <- function(
    signatures,
    other_user_signatures = NULL,
    collection = NULL,
    subcollection = NULL,
    metric = c("jaccard","odds_ratio"),
    universe = NULL,
    or_threshold = 1,
    pval_threshold = 0.05,
    limits = NULL,
    title_size = 12,
    color_values = c("#F9F4AE", "#B44141"),
    title = NULL,
    num_sigs_toplot = NULL,
    jaccard_threshold = 0,
    msig_subset = NULL,
    width_text = 20,
    na_color = "grey90"
) {
  if (is.null(signatures) || length(signatures) == 0) {
    stop("You must provide at least one signature.")
  }
  if (!is.list(signatures) || !all(sapply(signatures, is.character))) {
    stop("Signatures must be a named list of character vectors.")
  }
  if (!is.null(other_user_signatures) && (!is.list(other_user_signatures) || !all(sapply(other_user_signatures, is.character)))) {
    stop("Other user signatures must be a named list of character vectors.")
  }
  if (!is.null(collection) && !is.character(collection)) {
    stop("Collection must be a character string or NULL.")
  }
  if (!is.null(subcollection) && !is.character(subcollection)) {
    stop("Subcollection must be a character string or NULL.")
  }
  if (!is.null(universe) && !is.character(universe)) {
    stop("Universe must be a character vector or NULL.")
  }

  if (!is.numeric(or_threshold) || or_threshold < 0) {
    stop("or_threshold must be a non-negative numeric value.")
  }

  if (!is.numeric(pval_threshold) || pval_threshold < 0 || pval_threshold > 1) {
    stop("pval_threshold must be a numeric value between 0 and 1.")
  }

  if (!is.null(limits) && (!is.numeric(limits) || length(limits) != 2)) {
    stop("limits must be a numeric vector of length 2.")
  }

  if (!is.numeric(title_size) || title_size <= 0) {
    stop("title_size must be a positive numeric value.")
  }

  if (!is.character(color_values) || length(color_values) < 2) {
    stop("color_values must be a character vector with two colors.")
  }

  if (!is.null(title) && !is.character(title)) {
    stop("title must be a character string or NULL.")
  }

  if (!is.null(num_sigs_toplot) && (!is.numeric(num_sigs_toplot) || num_sigs_toplot <= 0)) {
    stop("num_sigs_toplot must be a positive numeric value or NULL.")
  }

  if (!is.numeric(jaccard_threshold) || jaccard_threshold < 0 || jaccard_threshold > 1) {
    stop("jaccard_threshold must be a numeric value between 0 and 1.")
  }
  if (!is.null(msig_subset) && (!is.character(msig_subset) || length(msig_subset) == 0)) {
    stop("msig_subset must be a character vector or NULL.")
  }

  if (!is.character(metric) || length(metric) != 1) {
    stop("metric must be a single character string.")
  }
  metric <- tolower(metric)
  if (is.null(metric) || metric == "") {
    stop("You must specify a metric: 'jaccard' or 'odds_ratio'.")
  } else if (!metric %in% c("jaccard", "odds_ratio")) {
    stop("Invalid metric specified. Use 'jaccard' or 'odds_ratio'.")
  }

  signatures <- lapply(signatures, toupper)
  if (!is.null(other_user_signatures)) {
    other_user_signatures <- lapply(other_user_signatures, toupper)
  }
  if (!is.null(universe)) {
    universe <- toupper(universe)
  }

  if (!is.null(collection)) {

    if (collection=="all"){
      gs <- msigdbr::msigdbr(
        species = "Homo sapiens",
        collection = NULL,
        subcollection = NULL
      )
    } else {
      gs <- msigdbr::msigdbr(
        species = "Homo sapiens",
        collection = collection,
        subcollection = subcollection
      )
    }


    gsets <- split(toupper(gs$gene_symbol), gs$gs_name)

    if (!is.null(msig_subset)) {
      gsets <- gsets[names(gsets) %in% msig_subset]
    }

  }

  if (is.null(other_user_signatures)) {
    other_user_signatures <- gsets
  } else {
    other_user_signatures <- c(other_user_signatures, gsets)
  }

  similarity_list <- list()

  for (ref_name in names(signatures)) {
    sig1 <- signatures[[ref_name]]

    for (comp_name in names(other_user_signatures)) {
      sig2 <- other_user_signatures[[comp_name]]

      if (metric == "jaccard") {
        score <- length(intersect(sig1, sig2)) / length(union(sig1, sig2))
        label <- sprintf("%.2f", score)
        pval <- NA
      } else {
        if (is.null(universe)) {
          stop("You must provide a gene universe for odds_ratio.")
        }

        a <- length(intersect(sig1, sig2))
        b <- length(setdiff(sig1, sig2))
        c <- length(setdiff(sig2, sig1))
        d <- length(setdiff(universe, union(sig1, sig2)))

        cont_tbl <- matrix(c(a, b, c, d), nrow = 2)
        ft <- fisher.test(cont_tbl)

        score <- log10(ft$estimate)
        if (!is.na(ft$p.value) && ft$p.value <= pval_threshold && ft$estimate >= or_threshold) {
          label <- sprintf("%.1f", score)
        } else {
          label <- ""
        }
        pval <- ft$p.value
      }

      row <- data.frame(
        Reference_Signature = ref_name,
        Compared_Signature = comp_name,
        Score = score,
        Label = label,
        Pval = pval,
        stringsAsFactors = FALSE
      )

      similarity_list[[length(similarity_list) + 1]] <- row
    }
  }

  # Combine all rows into one data frame
  similarity_df <- do.call(rbind, similarity_list)


  if (metric == "odds_ratio" ) {
    similarity_df <- similarity_df %>%
      dplyr::group_by(Compared_Signature) %>%
      dplyr::filter(any(10^Score >= or_threshold, na.rm = TRUE)) %>%
      dplyr::ungroup()
  }


  if (metric == "odds_ratio") {
    similarity_df <- similarity_df %>%
      dplyr::mutate(
        Label = ifelse(Pval <= pval_threshold, sprintf("%.1f", Score), "")
      )
  }

  if (metric == "jaccard" && jaccard_threshold > 0) {
    similarity_df <- similarity_df %>%
      dplyr::group_by(Compared_Signature) %>%
      dplyr::filter(any(Score >= jaccard_threshold, na.rm = TRUE)) %>%
      dplyr::ungroup()
  }


  similarity_df$Reference_Signature <- sapply(similarity_df$Reference_Signature, function(x) wrap_title(x, width_text))
  similarity_df$Compared_Signature <- sapply(similarity_df$Compared_Signature, function(x) wrap_title(x, width_text))

  if (is.null(limits)) {
    if (metric == "jaccard") {
      limits <- c(0, 1)
    } else {
      # For odds ratio, we set limits based on the data
      # but ensure they are at least 0 to avoid negative log10 values
    limits <- c(0, max(similarity_df$Score, na.rm = TRUE))
  }
  }

  plt <- ggplot(similarity_df, aes(x = Reference_Signature, y = Compared_Signature, fill = Score)) +
    geom_tile(color = "white") +
    geom_text(aes(label = Label), color = "black") +
    scale_fill_gradientn(colors = color_values, limits = limits, oob = scales::squish, na.value = na_color) +
    labs(
      x = "",
      y = "Compared Signature",
      fill = ifelse(metric == "jaccard", "Jaccard Index", "log10(OR)"),
      title = ifelse(is.null(title), paste("Signature Overlap (", metric, ")"), title)
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, size = title_size)
    )

  return(plt)
}
