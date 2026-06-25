#' Plot Signature Similarity via Jaccard Index or Fisher's Odds Ratio
#'
#' Visualizes similarity between user-defined gene signatures and either other
#' user-defined signatures or MSigDB gene sets, using either the Jaccard index
#' or Fisher's Odds Ratio. Produces a heatmap of pairwise similarity metrics.
#'
#' @param signatures A named list of character vectors representing reference
#'   gene signatures.
#' @param other_user_signatures Optional. A named list of character vectors
#'   representing other user-defined signatures to compare against.
#' @param collection Optional. MSigDB collection name (e.g., `"H"` for hallmark,
#'   `"C2"` for curated gene sets). Use msigdbr::msigdbr_collections() for the
#'   available options.
#' @param subcollection Optional. Subcategory within an MSigDB collection (e.g.,
#'   `"CP:REACTOME"`). Use msigdbr::msigdbr_collections() for the available
#'   options.
#' @param metric Character. Either "jaccard" or "odds_ratio".
#' @param universe Character vector. Background gene universe. Required for odds
#'   ratio.
#' @param or_threshold (only if method == "odds_ratio" only) Numeric. Minimum
#'   Odds Ratio required for a gene set to be included in the plot. Default is
#'   1.
#' @param pval_threshold (only if method == "odds_ratio" only) Numeric. Maximum
#'   adjusted p-value required for a gene set to be included in the plot.
#'   Default is 0.05.
#' @param limits Numeric vector of length 2. Limits for color scale. If `NULL`,
#'   is automatically set to c(0,1) for Jaccard or the range of OR for odds
#'   ratio.
#' @param title_size Integer specifying the font size for the plot title.
#'   Default is `12`.
#' @param color Character. The color for the maximum of the scale. Default is
#'   `red.`
#'   - If `method = "jaccard"`, the scale goes from `neutral_color` to `color`.
#'   - If `method = "odds_ratio"` and any OR >= 1, the scale ends at `color`.
#'   - If `method = "odds_ratio"` and all OR <= 1, `color` is not used; instead, the scale
#'   runs from `cold_color` (minimum) to `neutral_color` (OR = 1, if present;
#'   otherwise `neutral_color` is the maximum).
#' @param neutral_color Character. The neutral reference color. Default is
#'   `white`.
#'   - If `method = "jaccard"`, this is the minimum of the scale.
#'   - If `method = "odds_ratio"` and any OR >= 1, this corresponds to OR = 1 if such values exist; otherwise it is the minimum of the scale.
#'   - If `method = "odds_ratio"` and all OR <= 1, this corresponds to OR = 1 if such values exist; otherwise it is the maximum of the scale (with `cold_color` as the minimum).
#' @param cold_color Character. The color for values below OR = 1 (only used
#'   when `method = "odds_ratio"`). Default is `blue`.
#'   - If `method = "odds_ratio"` and any OR < 1, the scale runs from `cold_color`
#'   (minimum) to `neutral_color` (OR = 1 if present; otherwise `neutral_color`
#'   is the maximum).
#'   - Ignored if `method = "jaccard"` or if all OR >= 1.
#' @param title Optional. Custom title for the plot. If `NULL`, the title
#'   defaults to `"Signature Overlap"`.
#' @param jaccard_threshold (only if method == "jaccard" only) Numeric. Minimum
#'   Jaccard index required for a gene set to be included in the plot. Default
#'   is `0`.
#' @param organism Character. Species name passed to [msigdbr::msigdbr()].
#'   Defaults to `"Homo sapiens"`. Use `msigdbr::msigdbr_species()` to see all
#'   supported species (e.g. `"Mus musculus"`, `"Rattus norvegicus"`).
#'   Only relevant when `collection` is specified.
#' @param db_species Character or NULL. Controls which MSigDB database to query.
#'   `NULL` (default) omits the argument, using the msigdbr default (human gene
#'   sets with ortholog mapping). Pass `"MM"` to retrieve gene sets natively
#'   curated in mouse, or `"HS"` to explicitly request human-curated sets.
#'   Only relevant in msigdbr >= 2024 and when `collection` is specified.
#' @param msig_subset Optional. Character vector of MSigDB gene set names to
#'   subset from the specified collection. Useful to restrict analysis to a
#'   specific set of pathways. If supplied, other filters will apply only to
#'   this subset. Use "collection = "all" to mix gene sets from different
#'   collections.
#' @param width_text Integer. Character wrap width for labels.
#' @param na_color Character. Color for NA values in the heatmap. Default is
#'   `"grey90"`.
#'
#' @return Invisibly returns a list containing:
#'   \describe{
#'     \item{\code{plot}}{The \pkg{ggplot2} object of the similarity heatmap.}
#'     \item{\code{data}}{The data frame object containing the similarity
#'     scores per pair of gene sets.}
#'   }
#'
#' @import ggplot2
#' @importFrom tibble tibble
#' @importFrom msigdbr msigdbr
#' @importFrom scales squish rescale
#'
#' @examples
#' # Create two simple gene signatures
#' sig1 <- c("TP53", "BRCA1", "MYC", "EGFR", "CDK2")
#' sig2 <- c("ATXN2", "FUS", "MTOR", "CASP3")
#' signatures <- list(SignatureA = sig1, SignatureB = sig2)
#'
#' # Compare the signatures using the Jaccard index
#' plt <- geneset_similarity(
#'   signatures = signatures,
#'   metric = "jaccard",
#'   collection = "H",
#'   jaccard_threshold = 0.01
#' )
#'
#' # Print the plot (will show a small heatmap)
#' print(plt)
#'
#'
#' # Odds ratio example (requires universe)
#' gene_universe <- unique(c(
#'   sig1, sig2,
#'   msigdbr::msigdbr(species = "Homo sapiens", collection = "C2")$gene_symbol
#' ))
#'
#' plt_or <- geneset_similarity(
#'   signatures = signatures,
#'   metric = "odds_ratio",
#'   universe = gene_universe,
#'   collection = "H"
#' )
#' print(plt_or)
#'
#' @export
geneset_similarity <- function(
    signatures,
    other_user_signatures = NULL,
    collection = NULL,
    subcollection = NULL,
    organism = "Homo sapiens",
    db_species = NULL,
    metric = c("jaccard","odds_ratio"),
    universe = NULL,
    or_threshold = 1,
    pval_threshold = 0.05,
    limits = NULL,
    title_size = 12,
    color = "#B44141",         # color for the maximum of the scale
    neutral_color = "white",   # neutral reference color
    cold_color = "#4173B4",       # color for OR < 1 when applicable
    title = NULL,
    jaccard_threshold = 0,
    msig_subset = NULL,
    width_text = 20,
    na_color = "grey90"
) {
  if (is.null(signatures) || length(signatures) == 0) {
    stop("You must provide at least one signature.")
  }
  if (!is.list(signatures) || !all(vapply(signatures, is.character, logical(1)))) {
    stop("Signatures must be a named list of character vectors.")
  }
  if (!is.null(other_user_signatures) && (!is.list(other_user_signatures) ||
                                          !all(vapply(other_user_signatures, is.character, logical(1))))) {
    stop("Other user signatures must be a named list of character vectors.")
  }
  if (!is.null(collection) && !is.character(collection)) {
    stop("Collection must be a character string or NULL.")
  }
  if (!is.character(organism) || length(organism) != 1 || !nzchar(organism)) {
    stop("organism must be a non-empty character string (e.g. 'Homo sapiens', 'Mus musculus').")
  }
  if (!is.null(db_species) && (!is.character(db_species) || length(db_species) != 1 || !nzchar(db_species))) {
    stop("db_species must be a single character string (e.g. 'HS', 'MM') or NULL.")
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
  
  if (!is.null(limits) && any(limits < 0 | !is.finite(limits))) {
    warning("Limits contain negative, or non-finite values. Ensure limits are positive finite numbers.")
  }

  if (!is.numeric(title_size) || title_size <= 0) {
    stop("title_size must be a positive numeric value.")
  }
 

  if (!is.null(title) && !is.character(title)) {
    stop("title must be a character string or NULL.")
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

    .msigdbr_call <- function(collection_arg, subcollection_arg) {
      args <- list(species = organism, collection = collection_arg,
                   subcollection = subcollection_arg)
      if (!is.null(db_species)) args$db_species <- db_species
      do.call(msigdbr::msigdbr, args)
    }
    if (collection == "all") {
      gs <- .msigdbr_call(NULL, NULL)
    } else {
      gs <- .msigdbr_call(collection, subcollection)
    }


    gsets <- split(toupper(gs$gene_symbol), gs$gs_name)

    if (!is.null(msig_subset)) {
      gsets <- gsets[names(gsets) %in% msig_subset]
    }

  } else {

    gsets <- c()
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
        ft <- stats::fisher.test(cont_tbl)

        score <- log10(ft$estimate)
        # if (!is.na(ft$p.value) && ft$p.value <= pval_threshold && ft$estimate >= or_threshold) {
        #   label <- sprintf("%.1f", 10^score) # to show non log values in heatmap
        #   #label <- sprintf("%.1f", score)
        # } else {
        #   label <- ""
        # }
        pval <- ft$p.value
      }

      row <- data.frame(
        Reference_Signature = ref_name,
        Compared_Signature = comp_name,
        Score = score,
        #Label = label,
        Pval = pval,
        stringsAsFactors = FALSE
      )

      similarity_list[[length(similarity_list) + 1]] <- row
    }
  }

  # Combine all rows into one data frame
  similarity_df <- do.call(rbind, similarity_list)


  if (metric == "odds_ratio") {
    # Filter groups where any 10^Score >= threshold
    # keep_rows <- by(similarity_df, similarity_df$Compared_Signature, function(group) {
    #   any(10^group$Score >= or_threshold, na.rm = TRUE)
    # })

    keep_rows <- by(similarity_df, similarity_df$Compared_Signature, function(group) {
      any(10^group$Score >= or_threshold & group$Pval <= pval_threshold, na.rm = TRUE)
    })
    
    
    kept_signatures <- names(keep_rows[keep_rows])
    similarity_df <- similarity_df[similarity_df$Compared_Signature %in% kept_signatures, , drop = FALSE]

    # Add Label column
    # similarity_df$Label <- ifelse(
    #   similarity_df$Pval <= pval_threshold,
    #   sprintf("%.1f", similarity_df$Score),
    #   ""
    # )
  }

  if (metric == "jaccard" && jaccard_threshold > 0) {
    # Filter groups where any Score >= threshold
    keep_rows <- by(similarity_df, similarity_df$Compared_Signature, function(group) {
      any(group$Score >= jaccard_threshold, na.rm = TRUE)
    })

    kept_signatures <- names(keep_rows[keep_rows])
    similarity_df <- similarity_df[similarity_df$Compared_Signature %in% kept_signatures, , drop = FALSE]
  }

  data <- similarity_df
 
  
  
  if (nrow(similarity_df) == 0) {
    stop("No signatures passed the filtering criteria.")  
  }
  
  similarity_df$Reference_Signature <- vapply(similarity_df$Reference_Signature,
                                              function(x) wrap_title(x, width_text),
                                              character(1))
  similarity_df$Compared_Signature <- vapply(similarity_df$Compared_Signature,
                                             function(x) wrap_title(x, width_text),
                                             character(1))
# 
#   if (is.null(limits)) {
#     if (metric == "jaccard") {
#       limits <- c(0, 1)
#     } else {
#       # For odds ratio, we set limits based on the data 
#     limits <- c(min(similarity_df$Score[is.finite(similarity_df$Score)], na.rm = TRUE), max(similarity_df$Score, na.rm = TRUE))
#      
#   }
#   }
#   
#   
# 
#   plt <- ggplot(similarity_df, aes(x = .data$Reference_Signature,
#                                    y = .data$Compared_Signature, fill = .data$Score)) +
#     geom_tile(color = "white") +
#     #geom_text(aes(label = .data$Label), color = "black") +
#     scale_fill_gradientn(colors = color_values, limits = limits,
#                          oob = scales::squish, na.value = na_color) +
#     labs(
#       x = "",
#       y = "Compared Signature",
#       fill = ifelse(metric == "jaccard", "Jaccard Index", "log10(OR)"),
#       title = ifelse(is.null(title), paste("Signature Overlap (", metric, ")"), title)
#     ) +
#     theme_minimal() +
#     theme(
#       axis.text.x = element_text(angle = 45, hjust = 1),
#       plot.title = element_text(hjust = 0.5, size = title_size)
#     )
# 
#   plt
# 
#   invisible(list(plot=plt,
#                  data=data))
    
  
  # ----------------------------
  # Safe handling of limits (user provides OR)
  # ----------------------------
  if (is.null(limits)) {
    if (metric == "jaccard") {
      limits <- c(0, 1)
    } else { # odds_ratio
       
      # Extract finite scores
      finite_scores <- similarity_df$Score[is.finite(similarity_df$Score)]
      
      # Identify if there are any OR < 1 (logOR < 0)
      has_below1 <- any(finite_scores < 0)
      
      # Compute padding for -Inf (original OR = 0)
      if (has_below1) {
        # Place -Inf one log unit below the minimum finite score < 0
        min_below <- min(finite_scores[finite_scores < 0])
        pad_value <- min_below - 1
      } else {
        # All OR >= 1 → map -Inf slightly above the maximum score, then invert sign
        max_score <- max(finite_scores)
        pad_value <- -(max_score + 1)
      }
      
      # Replace -Inf in Score with computed padding
      similarity_df$Score[is.infinite(similarity_df$Score) & similarity_df$Score < 0] <- pad_value
      
      # Convert back from log
      OR_values <- 10^similarity_df$Score  
      
      # Replace any zero or negative OR with a small number
      OR_values[OR_values <= 0] <- 1e-6
      #similarity_df$Score[is.infinite(similarity_df$Score) & similarity_df$Score < 0] <- log10(1e-6)
      
      # Compute limits
      limits <- c(min(OR_values, na.rm = TRUE), max(OR_values, na.rm = TRUE))
    }
  }

  
  # Convert OR limits to log space (Score already log10 OR)
  log_limits <- if (metric == "odds_ratio") log10(limits) else limits
  if (min(log_limits) == max(log_limits)) {
    log_limits <- log_limits + c(-0.01, 0.01)  # small padding
  } 
  zero <- 0  # neutral color at OR = 1 → log10(1) = 0
  
  # ----------------------------
  # Define fill colors
  # ----------------------------
  if (metric == "jaccard") {
    fill_colors <- c(neutral_color, color)
    fill_values <- c(log_limits[1], log_limits[2])
  } else if (metric == "odds_ratio") {
    min_lim <- log_limits[1]
    max_lim <- log_limits[2]
    
    if (min_lim >= zero) {
      fill_colors <- c(neutral_color, color)
      fill_values <- c(min_lim, max_lim)
    } else if (max_lim <= zero) {
      fill_colors <- c(cold_color, neutral_color)
      fill_values <- c(min_lim, max_lim)
    } else {
      fill_colors <- c(cold_color, neutral_color, color)
      fill_values <- c(min_lim, zero, max_lim)
    }
    
    # ----------------------------
    # Safe legend breaks in OR space
    # ----------------------------
    valid_OR <- limits[limits > 0 & is.finite(limits)]
    if (length(valid_OR) == 0) valid_OR <- 1  # fallback
    
    log_breaks <- 10^seq(floor(log10(min(valid_OR))), ceiling(log10(max(valid_OR))))
    log_breaks <- log_breaks[log_breaks >= min(valid_OR) & log_breaks <= max(valid_OR)]
  }
  
  # ----------------------------
  # Build plot
  # ----------------------------
  plt <- ggplot(similarity_df, aes(
    x = .data$Reference_Signature,
    y = .data$Compared_Signature,
    fill = .data$Score)) +
    geom_tile(color = "white") +
    scale_fill_gradientn(
      colors = fill_colors,
      values = scales::rescale(fill_values),
      limits = log_limits,
      oob = scales::squish,
      na.value = na_color,
      trans = "identity",  # already logged
      breaks = if (metric == "odds_ratio") log10(log_breaks) else waiver(),
      labels = if (metric == "odds_ratio") log_breaks else waiver()
    ) +
    labs(
      x = "",
      y = "Compared Signature",
      fill = ifelse(metric == "jaccard", "Jaccard Index", "Odds Ratio"),
      title = ifelse(is.null(title), "Signature Overlap", title)
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, size = title_size)
    )
   
  if (metric == "odds_ratio") {
    data$Score <- 10^data$Score  # convert back to OR for data output
  }
  
  invisible(list(plot = plt, data = data))
  
  
}
