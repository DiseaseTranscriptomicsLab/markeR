# =============================================================================
# Gene Sets Tab – Shiny Module
# =============================================================================

# ---- Helpers -----------------------------------------------------------------

if (!exists("%||%", mode = "function"))
  `%||%` <- function(a, b) if (!is.null(a)) a else b

.gs_genes <- function(gs_entry) {
  if (is.data.frame(gs_entry)) as.character(gs_entry[[1]]) else as.character(gs_entry)
}
.gs_to_char_list <- function(gene_sets) lapply(gene_sets, .gs_genes)

# Human-readable labels for MSigDB collection identifiers.
# Covers both HS (human) collections (H, C1–C8) and MM (mouse-native) collections (MH, M1–M8).
.collection_labels <- c(
  # Human (HS) collections
  "H"  = "Hallmarks (H)",
  "C1" = "Positional (C1)",
  "C2" = "Curated gene sets (C2)",
  "C3" = "Regulatory targets (C3)",
  "C4" = "Computational (C4)",
  "C5" = "Ontology / GO (C5)",
  "C6" = "Oncogenic signatures (C6)",
  "C7" = "Immunologic (C7)",
  "C8" = "Cell type signatures (C8)",
  # Mouse-native (MM) collections
  "MH" = "Mouse Hallmarks (MH)",
  "M1" = "Positional (M1)",
  "M2" = "Curated gene sets (M2)",
  "M3" = "Regulatory targets (M3)",
  "M5" = "Ontology / GO (M5)",
  "M7" = "Immunologic (M7)",
  "M8" = "Cell type signatures (M8)"
)

# Query msigdbr at runtime to get the real available collections.
# db_species = NULL  → human MSigDB (HS), also used for all non-mouse species via ortholog mapping
# db_species = "MM"  → mouse-native MSigDB (MH instead of H; no "all" option)
# Falls back gracefully if the installed msigdbr version does not accept db_species or uses
# different column names (gs_cat in ≤7.5.x; column name varies in newer versions).
.get_msigdb_collections <- function(db_species = NULL) {
  .fallback <- function(mm = FALSE) {
    if (mm)
      c("Mouse Hallmarks (MH)" = "MH", "Curated gene sets (C2)" = "C2",
        "Ontology / GO (C5)" = "C5", "Oncogenic signatures (C6)" = "C6",
        "Immunologic (C7)" = "C7")
    else
      c("Hallmarks (H)" = "H", "Curated gene sets (C2)" = "C2",
        "Ontology / GO (C5)" = "C5", "Oncogenic signatures (C6)" = "C6",
        "Immunologic (C7)" = "C7", "All collections" = "all")
  }
  fetch <- function(ds) {
    if (!is.null(ds)) msigdbr::msigdbr_collections(db_species = ds)
    else              msigdbr::msigdbr_collections()
  }
  colls <- tryCatch(fetch(db_species), error = function(e) {
    if (!is.null(db_species)) tryCatch(fetch(NULL), error = function(e2) NULL) else NULL
  })
  if (is.null(colls) || nrow(colls) == 0) return(.fallback(identical(db_species, "MM")))

  # Column name varies by msigdbr version:
  #   "gs_collection" in msigdbr 2024+  |  "gs_cat" in msigdbr ≤7.5.x
  cat_col <- intersect(c("gs_collection", "gs_cat", "collection", "category"), colnames(colls))
  if (length(cat_col) == 0) {
    # Last resort: pick the first column that looks like short category codes
    char_cols <- colnames(colls)[sapply(colls, is.character)]
    cat_col <- if (length(char_cols) > 0) char_cols[1] else colnames(colls)[1]
  }
  cats <- sort(unique(as.character(colls[[cat_col[1]]])))
  if (length(cats) == 0) return(.fallback(identical(db_species, "MM")))

  labels  <- ifelse(cats %in% names(.collection_labels), .collection_labels[cats], cats)
  choices <- stats::setNames(cats, labels)
  # "all" (passes NULL to msigdbr) is only valid when NOT using db_species = "MM"
  if (!identical(db_species, "MM")) c(choices, "All collections" = "all") else choices
}

# Static fallback used only for the initial UI render (before the server reactive fires).
# Does NOT call msigdbr so it is safe at package load time.
.msigdb_collections <- c(
  "Hallmarks (H)"             = "H",
  "Curated gene sets (C2)"    = "C2",
  "Ontology / GO (C5)"        = "C5",
  "Oncogenic signatures (C6)" = "C6",
  "Immunologic (C7)"          = "C7",
  "All collections"           = "all"
)

.msigdb_organisms <- c(
  "Human (Homo sapiens)"           = "Homo sapiens",
  "Mouse (Mus musculus)"           = "Mus musculus",
  "Rat (Rattus norvegicus)"        = "Rattus norvegicus",
  "Zebrafish (Danio rerio)"        = "Danio rerio",
  "Fly (Drosophila melanogaster)"  = "Drosophila melanogaster",
  "Chicken (Gallus gallus)"        = "Gallus gallus",
  "Pig (Sus scrofa)"               = "Sus scrofa",
  "Macaque (Macaca mulatta)"       = "Macaca mulatta",
  "Chimpanzee (Pan troglodytes)"   = "Pan troglodytes",
  "Dog (Canis lupus familiaris)"   = "Canis lupus familiaris"
)

# Heuristic: infer organism from gene name casing conventions.
#   Human / primate  → ALL CAPS     (TP53, ACTB)
#   Mouse / rat etc. → Title case   (Trp53, Actb)
#   Zebrafish        → all lowercase (tp53, actb)
# Defaults to "Homo sapiens" when the signal is ambiguous.
.detect_organism <- function(gene_sets) {
  all_genes <- unlist(lapply(gene_sets, function(x) {
    if (is.data.frame(x)) as.character(x[[1]]) else as.character(x)
  }))
  all_genes <- all_genes[nzchar(all_genes) & !is.na(all_genes)]
  if (length(all_genes) == 0) return("Homo sapiens")
  prop_upper <- mean(grepl("^[A-Z][A-Z0-9]", all_genes))  # all-caps: human
  prop_title <- mean(grepl("^[A-Z][a-z]",    all_genes))  # title-case: mouse/rat
  prop_lower <- mean(grepl("^[a-z]",          all_genes))  # lowercase: zebrafish
  if (prop_upper > 0.5) return("Homo sapiens")
  if (prop_title > 0.5) return("Mus musculus")
  if (prop_lower > 0.5) return("Danio rerio")
  "Homo sapiens"
}

.dyn_h <- function(n, base = 200, per_item = 28, max_h = 4000L) {
  as.integer(min(max_h, base + n * per_item))
}

# Build OR heatmap from res$data — caps Inf self-comparisons and provides text aesthetic
.build_or_heatmap <- function(data, font = 10, title = "Odds Ratio",
                               cold = "#4173B4", mid = "white", hot = "#B44141") {
  df <- data
  finite_or <- df$Score[is.finite(df$Score) & df$Score > 0]
  if (length(finite_or) == 0) finite_or <- 1
  cap <- max(finite_or, na.rm = TRUE)
  df$fill_val <- log10(pmin(df$Score, cap))
  df$fill_val[!is.finite(df$fill_val)] <- log10(1e-6)
  log_range <- range(df$fill_val, na.rm = TRUE)
  if (diff(log_range) < 1e-9) log_range <- log_range + c(-0.5, 0.5)
  pval_col <- intersect(c("Pval", "pval", "p_adj", "padj"), colnames(df))
  df$hover_text <- paste0(
    "Reference: ", df$Reference_Signature,
    "\nCompared: ", df$Compared_Signature,
    "\nOR: ", ifelse(is.infinite(df$Score), "Inf", sprintf("%.2f", df$Score)),
    if (length(pval_col) > 0) paste0("\np-adj: ", sprintf("%.3e", df[[pval_col[1]]])) else ""
  )
  zero_in_range <- (log_range[1] < 0 && log_range[2] > 0)
  grad_colors <- if (zero_in_range) c(cold, mid, hot) else c(mid, hot)
  grad_values <- scales::rescale(
    if (zero_in_range) c(log_range[1], 0, log_range[2]) else log_range)
  ggplot2::ggplot(df, ggplot2::aes(x = .data$Reference_Signature,
                                   y = .data$Compared_Signature,
                                   fill = .data$fill_val,
                                   text = .data$hover_text)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradientn(colors = grad_colors, values = grad_values,
      limits = log_range, oob = scales::squish, na.value = "grey90",
      breaks = pretty(log_range, n = 5), labels = function(x) sprintf("%.1f", 10^x)) +
    ggplot2::labs(x = "", y = "Compared Signature", fill = "Odds Ratio", title = title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = font),
                   axis.text.y = ggplot2::element_text(size = font),
                   legend.text = ggplot2::element_text(size = font),
                   legend.title = ggplot2::element_text(size = font + 1),
                   legend.position = "bottom",
                   plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(barwidth = ggplot2::unit(10, "lines"),
                    barheight = ggplot2::unit(0.8, "lines"), title.position = "top"))
}

# Build Jaccard heatmap with text aesthetic for ggplotly hover
.build_ji_heatmap <- function(data, font = 10, title = "Jaccard Index") {
  df <- data
  df$hover_text <- paste0("Signature: ", df$Reference_Signature,
                          "\nPathway: ", df$Compared_Signature,
                          "\nJaccard: ", sprintf("%.4f", df$Score))
  ggplot2::ggplot(df, ggplot2::aes(x = .data$Reference_Signature,
                                   y = .data$Compared_Signature,
                                   fill = .data$Score,
                                   text = .data$hover_text)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradientn(colors = c("white", "#4173B4", "#801F4F"),
      limits = c(0, 1), oob = scales::squish, na.value = "grey90") +
    ggplot2::labs(x = "", y = "", fill = "Jaccard Index", title = title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = font),
                   axis.text.y = ggplot2::element_text(size = font),
                   legend.text = ggplot2::element_text(size = font),
                   legend.title = ggplot2::element_text(size = font + 1),
                   legend.position = "bottom",
                   plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(barwidth = ggplot2::unit(10, "lines"),
                    barheight = ggplot2::unit(0.8, "lines"), title.position = "top"))
}

# Build correlation heatmap at render time with controllable font size.
# corrmat_data: either a matrix (no separate.by) or a named list of matrices
.build_corr_heatmap <- function(corrmat_data, method = "spearman", font = 10,
                                 cold = "#4173B4", mid = "white", hot = "#B44141",
                                 show_col_names = TRUE) {
  col_fun <- circlize::colorRamp2(c(-1, 0, 1), c(cold, mid, hot))
  legend_title <- switch(method,
    spearman = "Spearman\nCorr.",
    pearson  = "Pearson\nCorr.",
    kendall  = "Kendall\nCorr.",
    "Correlation"
  )
  make_ht <- function(mat, title = NULL) {
    ComplexHeatmap::Heatmap(
      mat,
      name        = legend_title,
      col         = col_fun,
      column_title = title,
      show_row_names    = TRUE,
      show_column_names = show_col_names,
      row_names_gp      = grid::gpar(fontsize = font),
      column_names_gp   = grid::gpar(fontsize = font),
      heatmap_legend_param = list(
        title       = legend_title,
        title_gp    = grid::gpar(fontsize = font),
        labels_gp   = grid::gpar(fontsize = max(font - 2, 6))
      )
    )
  }
  if (is.matrix(corrmat_data) || is.data.frame(corrmat_data)) {
    ht <- make_ht(corrmat_data)
    ComplexHeatmap::draw(ht, merge_legend = TRUE)
  } else if (is.list(corrmat_data)) {
    conds <- names(corrmat_data)
    ht_list <- lapply(seq_along(conds), function(i) make_ht(corrmat_data[[i]], conds[i]))
    combined <- Reduce(`+`, ht_list)
    ComplexHeatmap::draw(combined, merge_legend = TRUE)
  }
}

# Styled settings panel used inside each main-panel tab
.tab_settings <- function(...) {
  shiny::div(
    style = paste("background:#f8f9fa; border:1px solid #dee2e6; border-radius:6px;",
                  "padding:12px 14px; margin-bottom:14px;"),
    ...
  )
}

# ---- UI ---------------------------------------------------------------------

#' @importFrom bslib layout_sidebar sidebar navset_card_tab nav_panel tooltip
#' @importFrom shiny NS radioButtons sliderInput selectInput numericInput
#'   actionButton uiOutput plotOutput hr h4 h5 p tags div icon checkboxInput
#' @importFrom shinyWidgets pickerInput pickerOptions
#' @importFrom plotly plotlyOutput
#' @importFrom DT DTOutput

geneSetsUI <- function(id) {
  ns <- shiny::NS(id)

  bslib::layout_sidebar(

    # ---- Sidebar: shared gene/sample selection ------------------------------
    sidebar = bslib::sidebar(
      width = 310,
      shiny::div(
        style = "padding-bottom:6px;",
        shiny::h4("Gene Sets", style = "color:#801F4F;"),
        shiny::p(shiny::em(
          "Assess overlap with MSigDB, cross-gene-set similarity,",
          "and individual gene behaviour across all plots."),
          style = "color:#6c757d; font-size:0.85em;"),
        shiny::hr()
      ),
      shiny::tags$small(shiny::tags$b("Gene / sample selection")),
      shiny::helpText(style = "margin-bottom:6px; font-size:0.82em;",
        "Shared across all individual gene-level plots."),
      shiny::uiOutput(ns("gs_picker_ui")),
      shiny::uiOutput(ns("large_gs_warning_ui")),
      shiny::uiOutput(ns("gene_selector_ui")),
      shiny::uiOutput(ns("group_var_ui")),
      shiny::uiOutput(ns("level_warning_ui")),
      shiny::uiOutput(ns("color_var_ui"))
    ),

    # ---- Main panel ---------------------------------------------------------
    bslib::navset_card_tab(

      # ---- 1. Pairwise Overlap -----------------------------------------------
      bslib::nav_panel("Pairwise Overlap",
        .tab_settings(
          shiny::div(
            style = "display:grid; grid-template-columns:1fr 1fr; gap:10px;",
            shiny::div(
              shiny::tags$label("Metric ",
                bslib::tooltip(shiny::icon("circle-question", style="color:#6c757d;cursor:help;"),
                  shiny::tags$span(
                    shiny::tags$b("Jaccard index:"), " |A ∩ B| / |A ∪ B|. JI=0: no shared genes.",
                    shiny::tags$br(),
                    shiny::tags$b("Odds ratio:"), " OR=0: no overlap. OR=1: chance. OR>1: enriched.",
                    " Self-comparisons = Inf (shown at max colour)."
                  ), placement = "right")),
              shiny::radioButtons(ns("pw_metric"), label = NULL,
                choices  = c("Jaccard index" = "jaccard", "Odds ratio" = "odds_ratio"),
                selected = "jaccard")
            ),
            shiny::div(
              shiny::numericInput(ns("pw_fontsize"), "Axis font size (pt):",
                                  value = 10, min = 6, max = 20, step = 1),
              shiny::actionButton(ns("run_pairwise"), "Run Pairwise Overlap",
                                  class = "btn-primary btn-sm", width = "100%")
            )
          )
        ),
        shiny::uiOutput(ns("pairwise_status_ui")),
        shiny::uiOutput(ns("pairwise_plot_ui")),
        shiny::uiOutput(ns("pairwise_table_hdr_ui")),
        DT::DTOutput(ns("pairwise_table"))
      ),

      # ---- 2. Similarity vs MSigDB ------------------------------------------
      bslib::nav_panel("Similarity vs MSigDB",
        .tab_settings(
          # Row 1: metric + collection
          shiny::div(
            style = "display:grid; grid-template-columns:1fr 1fr; gap:10px; margin-bottom:8px;",
            shiny::div(
              shiny::tags$label("Similarity metric ",
                bslib::tooltip(shiny::icon("circle-question", style="color:#6c757d;cursor:help;"),
                  shiny::tags$span(
                    shiny::tags$b("Jaccard:"), " JI=0: no overlap. JI=1: identical.",
                    " Pathway kept if ANY signature exceeds threshold.",
                    shiny::tags$br(),
                    shiny::tags$b("Odds ratio:"), " OR=1: chance. OR>1: enriched.",
                    " Pathway removed only if ALL fail OR+p-value thresholds."
                  ), placement = "right")),
              shiny::radioButtons(ns("gs_metric"), label = NULL,
                choices  = c("Jaccard index" = "jaccard", "Odds ratio" = "odds_ratio"),
                selected = "jaccard")
            ),
            shiny::div(
              shiny::tags$label("Organism ",
                bslib::tooltip(shiny::icon("circle-question", style = "color:#6c757d;cursor:help;"),
                  shiny::tags$span(
                    "Selects the species whose gene symbols are used from MSigDB.",
                    shiny::tags$br(),
                    shiny::tags$b("Note:"), " Pathway names (e.g. HALLMARK_APOPTOSIS) are",
                    " identical across species — only the gene symbols within each",
                    " pathway change. Check the plot title to confirm which organism was used."
                  ), placement = "right")),
              shinyWidgets::pickerInput(ns("gs_organism"), label = NULL,
                choices  = .msigdb_organisms,
                selected = "Homo sapiens",
                options  = shinyWidgets::pickerOptions(liveSearch = TRUE)),
              shiny::uiOutput(ns("gs_db_species_ui")),
              shiny::uiOutput(ns("gs_organism_info_ui")),
              shiny::uiOutput(ns("gs_collection_ui")),
              shiny::textInput(ns("gs_subcollection"), "Subcollection (optional):",
                placeholder = "e.g. CP:REACTOME")
            )
          ),
          # Row 2: pathway filter + threshold controls
          shiny::div(
            style = "display:grid; grid-template-columns:1fr 1fr; gap:10px; margin-bottom:8px;",
            shiny::div(
              shiny::textAreaInput(ns("msig_subset"),
                label = shiny::tags$span("Restrict to pathways (optional) ",
                  bslib::tooltip(shiny::icon("circle-question", style="color:#6c757d;cursor:help;"),
                    "Comma-separated exact MSigDB pathway names.", placement = "right")),
                placeholder = "HALLMARK_APOPTOSIS, ...", rows = 2)
            ),
            shiny::div(
              shiny::conditionalPanel(paste0("input['", ns("gs_metric"), "'] == 'jaccard'"),
                shiny::sliderInput(ns("jaccard_threshold"), "Min. Jaccard index:",
                                   min = 0, max = 1, value = 0, step = 0.01),
                shiny::uiOutput(ns("jaccard_size_warn_ui"))),
              shiny::conditionalPanel(paste0("input['", ns("gs_metric"), "'] == 'odds_ratio'"),
                shiny::numericInput(ns("or_threshold"), "Min. odds ratio:", value = 2, min = 0, step = 0.5),
                shiny::sliderInput(ns("or_pval"), "Max adj. p-value:",
                                   min = 0, max = 1, value = 0.05, step = 0.01))
            )
          ),
          # Row 3: font size + run button
          shiny::div(
            style = "display:grid; grid-template-columns:1fr 1fr; gap:10px; align-items:end;",
            shiny::numericInput(ns("sim_fontsize"), "Axis font size (pt):", value = 9, min = 6, max = 20, step = 1),
            shiny::actionButton(ns("run_similarity"), "Run Similarity",
                                class = "btn-primary btn-sm", width = "100%",
                                style = "margin-top:24px;")
          )
        ),
        shiny::uiOutput(ns("similarity_status_ui")),
        shiny::uiOutput(ns("similarity_plot_ui")),
        shiny::uiOutput(ns("similarity_table_hdr_ui")),
        shiny::uiOutput(ns("similarity_table_container"))
      ),

      # ---- 3. Expression (violins + heatmap) --------------------------------
      bslib::nav_panel("Expression",
        .tab_settings(
          # Gene picker tip
          shiny::div(
            style = paste("background:#fff8e1; border-left:3px solid #EBB43E;",
                          "border-radius:4px; padding:7px 11px; margin-bottom:10px;",
                          "font-size:0.84em;"),
            shiny::icon("hand-point-left", style = "color:#b08000;"), " ",
            "Choose ", shiny::strong("which genes"), " appear in both plots using the",
            shiny::strong("gene picker in the sidebar"),
            " (select all, a subset, or type to search).",
            " For best results keep the selection under 20 genes."
          ),
          shiny::uiOutput(ns("expr_gene_count_warn_ui")),
          shiny::div(
            style = "display:grid; grid-template-columns:1fr 1fr; gap:14px;",
            shiny::div(
              shiny::tags$small(shiny::tags$b("Violin options")),
              shiny::numericInput(ns("expr_fontsize"), "Axis font size (pt):",
                                  value = 10, min = 6, max = 20, step = 1),
              shiny::div(style = "margin-top:14px;",
                shiny::actionButton(ns("run_expression"), "Plot Violins",
                                    class = "btn-primary btn-sm", width = "100%"))
            ),
            shiny::div(
              shiny::tags$small(shiny::tags$b("Heatmap options")),
              shiny::div(style = "display:flex; gap:8px; margin-top:4px;",
                shiny::checkboxInput(ns("hm_cluster_rows"), "Cluster genes",    value = TRUE),
                shiny::checkboxInput(ns("hm_cluster_cols"), "Cluster samples",  value = TRUE)
              ),
              shiny::numericInput(ns("hm_titlesize"), "Title size (pt):",
                                  value = 12, min = 8, max = 22, step = 1),
              shiny::actionButton(ns("run_heatmap"), "Plot Expression Heatmap",
                                  class = "btn-outline-secondary btn-sm", width = "100%")
            )
          )
        ),
        shiny::uiOutput(ns("violin_status_ui")),
        shiny::uiOutput(ns("expr_plot_ui")),
        shiny::uiOutput(ns("heatmap_status_ui")),
        shiny::uiOutput(ns("heatmap_plot_ui"))
      ),

      # ---- 4. Correlation ---------------------------------------------------
      bslib::nav_panel("Correlation",
        .tab_settings(
          shiny::div(
            style = "display:grid; grid-template-columns:1fr 1fr; gap:10px; align-items:end;",
            shiny::div(
              shiny::selectInput(ns("corr_method"), "Correlation method:",
                choices = c("Spearman" = "spearman", "Pearson" = "pearson", "Kendall" = "kendall"),
                selected = "spearman"),
              shiny::uiOutput(ns("corr_sep_ui"))
            ),
            shiny::div(
              shiny::checkboxInput(ns("show_col_names"), "Show gene labels on x-axis", value = TRUE),
              shiny::numericInput(ns("corr_fontsize"), "Axis font size (pt):",
                                  value = 10, min = 6, max = 20, step = 1),
              shiny::actionButton(ns("run_corr"), "Plot Correlation Heatmap",
                                  class = "btn-primary btn-sm", width = "100%")
            )
          )
        ),
        shiny::uiOutput(ns("corr_status_ui")),
        shiny::uiOutput(ns("corr_plot_ui"))
      ),

      # ---- 5. Cohen's d -----------------------------------------------------
      bslib::nav_panel("Cohen's d",
        .tab_settings(
          shiny::helpText(style = "margin-bottom:8px; font-size:0.85em;",
            "Effect size: how many standard deviations the positive class differs from",
            "the rest. Requires a binary condition variable."),
          shiny::div(
            style = "display:grid; grid-template-columns:1fr 1fr; gap:14px;",
            shiny::div(
              shiny::tags$small(shiny::tags$b("Analysis parameters")),
              shiny::uiOutput(ns("cohend_cond_var_ui")),
              shiny::uiOutput(ns("cohend_cond_class_ui")),
              shiny::numericInput(ns("cohend_max_genes"),
                label = shiny::tags$span("Max genes to plot ",
                  bslib::tooltip(shiny::icon("circle-question", style="color:#6c757d;cursor:help;"),
                    shiny::tags$span(
                      "Limits the number of genes plotted. Set to 0 to plot all selected genes.",
                      shiny::tags$br(), shiny::tags$br(),
                      shiny::tags$b("Top N by |Cohen's d|:"), " runs the analysis on ALL selected",
                      " genes, then displays only those with the largest absolute effect size.",
                      " This is usually what you want.",
                      shiny::tags$br(),
                      shiny::tags$b("First N (gene set order):"), " takes the first N genes in the",
                      " order shown in the sidebar picker — i.e. the order they appear in your gene set.",
                      shiny::tags$br(),
                      shiny::tags$b("Random N:"), " samples N genes at random from the picker selection."
                    ), placement = "right")),
                value = 30, min = 0, max = 500, step = 5),
              shiny::radioButtons(ns("cohend_gene_sel"), "Gene selection:",
                choices  = c("Top N by |Cohen's d|"   = "top_d",
                             "First N (gene set order)" = "first",
                             "Random N"                 = "random"),
                selected = "top_d")
            ),
            shiny::div(
              shiny::tags$small(shiny::tags$b("Display options")),
              shiny::numericInput(ns("cohend_fontsize"), "Axis font size (pt):",
                                  value = 10, min = 6, max = 20, step = 1),
              shiny::div(style = "margin-top:14px;",
                shiny::actionButton(ns("run_cohend"), "Plot Cohen's d",
                                    class = "btn-primary btn-sm", width = "100%"))
            )
          )
        ),
        shiny::uiOutput(ns("cohend_status_ui")),
        shiny::uiOutput(ns("cohend_plot_ui"))
      ),

      # ---- 6. Gene PCA ------------------------------------------------------
      bslib::nav_panel("Gene PCA",
        .tab_settings(
          shiny::helpText(style = "margin-bottom:6px; font-size:0.85em;",
            "PCA run on the expression matrix subsetted to the selected genes.",
            "Each point is a sample; genes are the features used to build the",
            "principal components. Points are coloured by the",
            shiny::strong("Colour by"), "variable chosen in the sidebar.",
            "Hover over points to see sample names."),
          shiny::numericInput(ns("pca_fontsize"), "Axis font size (pt):", value = 10, min = 6, max = 20, step = 1),
          shiny::actionButton(ns("run_pca"), "Plot Gene PCA",
                              class = "btn-primary btn-sm", width = "100%")
        ),
        shiny::uiOutput(ns("pca_status_ui")),
        shiny::uiOutput(ns("pca_plot_ui"))
      ),

      # ---- 7. ROC / AUC -----------------------------------------------------
      bslib::nav_panel("ROC / AUC",
        .tab_settings(
          shiny::helpText(style = "margin-bottom:8px; font-size:0.85em;",
            "AUC = 0.5: gene performs no better than chance.",
            "AUC = 1: perfect separation. ROC curves show the",
            "true positive rate vs false positive rate across all thresholds.",
            " Colours are consistent between the ROC curves and the AUC barplot."),
          shiny::div(
            style = "display:grid; grid-template-columns:1fr 1fr; gap:14px;",
            shiny::div(
              shiny::tags$small(shiny::tags$b("Analysis parameters")),
              shiny::uiOutput(ns("roc_cond_var_ui")),
              shiny::uiOutput(ns("roc_cond_class_ui")),
              shiny::numericInput(ns("roc_max_genes"),
                label = shiny::tags$span("Max genes to analyse ",
                  bslib::tooltip(shiny::icon("circle-question", style="color:#6c757d;cursor:help;"),
                    shiny::tags$span(
                      "Limits the number of genes analysed and plotted.",
                      shiny::tags$br(),
                      shiny::tags$b("First N as listed:"), " takes the first N genes in the",
                      " order shown in the sidebar gene picker (i.e. the order in your gene set).",
                      shiny::tags$br(),
                      shiny::tags$b("Random N:"), " samples N genes at random.",
                      shiny::tags$br(), "Set to 0 to use all genes. Recommended ≤ 20."
                    ), placement = "right")),
                value = 20, min = 0, max = 500, step = 5),
              shiny::radioButtons(ns("roc_gene_sel"), "Gene selection:",
                choices  = c("First N (as listed)" = "first", "Random N" = "random"),
                selected = "first", inline = TRUE)
            ),
            shiny::div(
              shiny::tags$small(shiny::tags$b("Display options")),
              shiny::numericInput(ns("roc_fontsize"), "Axis font size (pt):",
                                  value = 10, min = 6, max = 20, step = 1),
              shiny::div(style = "margin-top:14px;",
                shiny::actionButton(ns("run_rocauc"), "Plot ROC / AUC",
                                    class = "btn-primary btn-sm", width = "100%"))
            )
          )
        ),
        shiny::uiOutput(ns("rocauc_status_ui")),
        shiny::uiOutput(ns("roc_curves_ui")),
        shiny::uiOutput(ns("rocauc_plot_ui"))
      )
    )
  )
}

# ---- Server -----------------------------------------------------------------

#' @importFrom shiny moduleServer reactiveVal observeEvent req renderPlot
#'   renderUI withProgress showNotification isolate outputOptions
#' @importFrom plotly renderPlotly ggplotly
#' @importFrom DT renderDT datatable
#' @importFrom scales rescale squish
#' @importFrom ComplexHeatmap draw
geneSetsServer <- function(id, get_expr, get_meta, get_gene_sets) {

  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ---- Cached results -----------------------------------------------------
    pairwise_result   <- shiny::reactiveVal(NULL)
    similarity_result <- shiny::reactiveVal(NULL)
    expression_result <- shiny::reactiveVal(NULL)
    heatmap_result    <- shiny::reactiveVal(NULL)
    corr_result       <- shiny::reactiveVal(NULL)
    cohend_result     <- shiny::reactiveVal(NULL)
    pca_result        <- shiny::reactiveVal(NULL)
    rocauc_result     <- shiny::reactiveVal(NULL)

    # ---- Per-plot reactive heights ------------------------------------------
    pw_h         <- shiny::reactiveVal(400L)
    sim_h        <- shiny::reactiveVal(500L)
    expr_h       <- shiny::reactiveVal(500L)
    hm_h         <- shiny::reactiveVal(500L)
    corr_h       <- shiny::reactiveVal(400L)
    cohend_h     <- shiny::reactiveVal(400L)
    pca_h        <- shiny::reactiveVal(600L)
    roc_curves_h <- shiny::reactiveVal(400L)
    rocauc_h     <- shiny::reactiveVal(450L)

    # ---- Reset when gene sets uploaded/changed (all results) ----------------
    shiny::observeEvent(get_gene_sets(), {
      pairwise_result(NULL); similarity_result(NULL)
      expression_result(NULL); heatmap_result(NULL)
      corr_result(NULL); cohend_result(NULL)
      pca_result(NULL); rocauc_result(NULL)
      roc_curves_h(400L)
    }, ignoreInit = TRUE)

    # ---- Reset gene-level plots when selected gene set changes --------------
    shiny::observeEvent(input$selected_gs, {
      expression_result(NULL); heatmap_result(NULL)
      corr_result(NULL); cohend_result(NULL)
      pca_result(NULL); rocauc_result(NULL)
      roc_curves_h(400L)
    }, ignoreInit = TRUE, ignoreNULL = TRUE)

    # ---- Reset violin/heatmap when grouping variable changes ----------------
    shiny::observeEvent(input$group_var, {
      expression_result(NULL); heatmap_result(NULL)
    }, ignoreInit = TRUE, ignoreNULL = TRUE)

    # ---- Reset PCA when colour variable changes -----------------------------
    shiny::observeEvent(input$color_var, {
      pca_result(NULL)
    }, ignoreInit = TRUE, ignoreNULL = TRUE)

    # ---- Reset when selected genes change (subset change) -------------------
    shiny::observeEvent(input$selected_genes, {
      expression_result(NULL); heatmap_result(NULL)
      corr_result(NULL); cohend_result(NULL)
      pca_result(NULL); rocauc_result(NULL)
    }, ignoreInit = TRUE, ignoreNULL = TRUE)

    # ---- Reset when metric changes ------------------------------------------
    shiny::observeEvent(input$pw_metric,    { pairwise_result(NULL)   }, ignoreInit = TRUE)
    shiny::observeEvent(input$gs_metric,    { similarity_result(NULL) }, ignoreInit = TRUE)
    shiny::observeEvent(input$gs_organism,  { similarity_result(NULL) }, ignoreInit = TRUE)

    # =========================================================================
    # ---- SIDEBAR: dynamic UI ------------------------------------------------
    # =========================================================================

    output$gs_picker_ui <- shiny::renderUI({
      gs <- get_gene_sets()
      if (is.null(gs) || length(gs) == 0) return(NULL)
      shinyWidgets::pickerInput(ns("selected_gs"), "Gene set:",
        choices = names(gs), selected = names(gs)[1],
        options = shinyWidgets::pickerOptions(liveSearch = TRUE))
    })

    output$large_gs_warning_ui <- shiny::renderUI({
      gs <- get_gene_sets()
      if (is.null(gs) || is.null(input$selected_gs) || !input$selected_gs %in% names(gs)) return(NULL)
      n <- length(.gs_genes(gs[[input$selected_gs]]))
      if (n > 20)
        shiny::div(class = "alert alert-warning",
                   style = "font-size:0.82em; padding:5px 9px; margin-bottom:4px;",
                   shiny::icon("triangle-exclamation"), " ",
                   shiny::strong(n), " genes — violin plots may be crowded.")
    })

    output$gene_selector_ui <- shiny::renderUI({
      gs <- get_gene_sets()
      if (is.null(gs) || is.null(input$selected_gs) || !input$selected_gs %in% names(gs)) return(NULL)
      genes <- .gs_genes(gs[[input$selected_gs]])
      shinyWidgets::pickerInput(ns("selected_genes"), "Genes to plot:",
        choices = genes, selected = genes, multiple = TRUE,
        options = shinyWidgets::pickerOptions(actionsBox = TRUE, liveSearch = TRUE,
          selectedTextFormat = "count > 4", countSelectedText = "{0} genes selected"))
    })

    output$group_var_ui <- shiny::renderUI({
      meta <- get_meta()
      if (is.null(meta)) return(NULL)
      cols <- setdiff(colnames(meta), colnames(meta)[1])
      shiny::selectInput(ns("group_var"), "Grouping variable (x-axis / violin):",
                         choices = cols, selected = cols[1])
    })

    output$level_warning_ui <- shiny::renderUI({
      meta <- get_meta()
      if (is.null(meta) || is.null(input$group_var) || !input$group_var %in% colnames(meta)) return(NULL)
      n_lev <- length(unique(as.character(meta[[input$group_var]])))
      if (n_lev > 10)
        shiny::div(class = "alert alert-warning",
                   style = "font-size:0.82em; padding:5px 9px; margin-bottom:4px;",
                   shiny::icon("triangle-exclamation"), " ", shiny::strong(n_lev),
                   " levels — plots may be hard to read.")
    })

    output$color_var_ui <- shiny::renderUI({
      meta <- get_meta()
      if (is.null(meta)) return(NULL)
      cols <- setdiff(colnames(meta), colnames(meta)[1])
      shiny::selectInput(ns("color_var"), "Colour by (optional):",
                         choices = c("None" = "", cols), selected = "")
    })

    # In-tab gene-count warning (Expression tab) — shown when > 20 genes selected
    output$expr_gene_count_warn_ui <- shiny::renderUI({
      n <- length(input$selected_genes)
      if (n > 20)
        shiny::div(
          class = "alert alert-warning",
          style = "font-size:0.82em; padding:6px 10px; margin-bottom:8px;",
          shiny::icon("triangle-exclamation"), " ",
          shiny::strong(n), " genes currently selected.",
          " Violin plots and heatmaps with more than ~20 genes can be very slow",
          " or hard to read. Use the sidebar picker to reduce the selection.")
    })

    # =========================================================================
    # ---- PER-TAB: dynamic UI ------------------------------------------------
    # =========================================================================

    # Warn when total gene set size is large (Jaccard will be low by design)
    output$jaccard_size_warn_ui <- shiny::renderUI({
      gs <- get_gene_sets()
      if (is.null(gs) || length(gs) == 0) return(NULL)
      n_total <- sum(sapply(gs, function(x) {
        if (is.data.frame(x)) nrow(x) else length(x)
      }))
      if (n_total > 500)
        shiny::div(
          class = "alert alert-warning",
          style = "font-size:0.80em; padding:5px 8px; margin-top:4px;",
          shiny::icon("triangle-exclamation"), " ",
          shiny::strong(n_total), " genes in your signature(s).",
          " Jaccard scores will be low by design — keep Min. Jaccard at 0 to see all pathways."
        )
      else NULL
    })

    # ---- Organism / db_species / collection reactives -------------------------

    # Auto-detect organism from gene name casing when gene sets are uploaded.
    # Uses session$onFlushed so the updatePickerInput reaches the browser only
    # after all outputs have been sent (avoids silent early-fire timing failures).
    shiny::observeEvent(get_gene_sets(), {
      gs <- get_gene_sets()
      if (is.null(gs) || length(gs) == 0) return()
      detected <- .detect_organism(gs)
      session$onFlushed(function() {
        shinyWidgets::updatePickerInput(session, "gs_organism", selected = detected)
      }, once = TRUE)
    }, ignoreInit = FALSE)

    # db_species radio — shown ONLY for Mus musculus because that is the only
    # species for which msigdbr ships a native database (db_species = "MM").
    # Every other non-human species uses ortholog mapping from the HS database.
    output$gs_db_species_ui <- shiny::renderUI({
      org <- input$gs_organism
      if (is.null(org) || org != "Mus musculus") return(NULL)
      shiny::radioButtons(
        ns("gs_db_species"),
        label = shiny::tags$span("Gene set database ",
          bslib::tooltip(shiny::icon("circle-question", style = "color:#6c757d;cursor:help;"),
            shiny::tags$span(
              shiny::tags$b("Native mouse (MM):"),
              " gene sets curated directly in mouse (e.g. Mouse Hallmarks = MH).",
              shiny::tags$br(),
              shiny::tags$b("Human with orthologs (HS):"),
              " human-curated sets with symbols mapped to mouse orthologs."
            ), placement = "right")),
        choices  = c("Native mouse gene sets (MM)" = "MM",
                     "Human with ortholog mapping (HS)" = "HS"),
        selected = "MM",
        inline   = FALSE
      )
    })

    # Info banner: explains which database/mapping is being used and what to expect
    output$gs_organism_info_ui <- shiny::renderUI({
      org <- input$gs_organism
      if (is.null(org)) return(NULL)
      ds  <- effective_db_species()

      lines <- if (org == "Homo sapiens") {
        list(
          shiny::tagList(shiny::icon("circle-info", style = "color:#1d4ed8;"),
            " Searching the ", shiny::strong("human MSigDB (HS)"),
            " — gene symbols matched as-is (case-insensitive).")
        )
      } else if (org == "Mus musculus" && identical(ds, "MM")) {
        list(
          shiny::tagList(shiny::icon("circle-info", style = "color:#15803d;"),
            " Searching the ", shiny::strong("mouse-native MSigDB (MM)"),
            " — gene sets curated directly in mouse. Collections: MH, M2, M3, M5, M7, M8."),
          shiny::tagList(shiny::icon("triangle-exclamation", style = "color:#b45309;"),
            " Pathway names (e.g. HALLMARK_APOPTOSIS) are the same as in human MSigDB.",
            " Only the gene symbols inside each pathway are species-specific.")
        )
      } else if (org == "Mus musculus") {
        list(
          shiny::tagList(shiny::icon("circle-info", style = "color:#b45309;"),
            " Searching the ", shiny::strong("human MSigDB (HS)"),
            " with symbols ", shiny::strong("mapped to mouse orthologs"),
            " by msigdbr. Same pathway names as human; gene content differs."),
          shiny::tagList(shiny::icon("triangle-exclamation", style = "color:#6c757d;"),
            " Results may look similar to human because many gene names are conserved",
            " across species (e.g. COL5A1 is the same symbol in both).",
            " Differences appear for genes with species-specific symbols (e.g.",
            " mouse Trp53 ≠ human TP53).")
        )
      } else {
        list(
          shiny::tagList(shiny::icon("circle-info", style = "color:#6c757d;"),
            " Searching the ", shiny::strong("human MSigDB (HS)"),
            " with symbols mapped to ", shiny::em(org), " orthologs via msigdbr.",
            " Not all genes have a 1-to-1 ortholog."),
          shiny::tagList(shiny::icon("triangle-exclamation", style = "color:#6c757d;"),
            " Results may look similar across species because the comparison is",
            " gene-symbol based: conserved genes with identical symbols always match,",
            " regardless of organism. Only genes with species-specific symbols differ.")
        )
      }

      shiny::div(
        style = paste("font-size:0.81em; padding:6px 9px; margin-top:6px;",
                      "background:#f8f9fa; border-left:3px solid #adb5bd;",
                      "border-radius:3px; display:flex; flex-direction:column; gap:3px;"),
        lapply(lines, function(l) shiny::div(l))
      )
    })

    # Effective db_species passed to msigdbr:
    #   Mouse + MM radio → "MM"
    #   Mouse + HS radio → NULL  (HS ortholog mapping is the msigdbr default)
    #   Human            → NULL
    #   Any other species→ NULL  (only HS ortholog mapping available)
    effective_db_species <- shiny::reactive({
      org <- input$gs_organism
      if (is.null(org)) return(NULL)
      if (org == "Mus musculus") {
        ds <- input$gs_db_species          # NULL until the radio renders
        if (is.null(ds) || ds == "MM") return("MM")
        return(NULL)                        # HS → omit db_species
      }
      NULL
    })

    # Dynamic collection picker — rebuilt whenever organism or db_species changes.
    # Queries msigdbr_collections() for the real available collections so the list
    # is always correct regardless of species or msigdbr version.
    output$gs_collection_ui <- shiny::renderUI({
      db_sp   <- effective_db_species()
      choices <- .get_msigdb_collections(db_species = db_sp)
      # Sensible default: first hallmark-style entry, else first available
      default_col <- if ("H"  %in% choices) "H"  else
                     if ("MH" %in% choices) "MH" else choices[[1]]
      shinyWidgets::pickerInput(ns("gs_collection"), "MSigDB collection:",
        choices  = choices,
        selected = default_col,
        options  = shinyWidgets::pickerOptions(liveSearch = TRUE))
    })

    output$corr_sep_ui <- shiny::renderUI({
      meta <- get_meta()
      if (is.null(meta)) return(NULL)
      cols <- setdiff(colnames(meta), colnames(meta)[1])
      shiny::selectInput(ns("corr_sep"), "Separate by (optional):",
                         choices = c("None" = "", cols), selected = "")
    })

    output$cohend_cond_var_ui <- shiny::renderUI({
      meta <- get_meta()
      if (is.null(meta)) return(NULL)
      cols <- setdiff(colnames(meta), colnames(meta)[1])
      cat_cols <- cols[sapply(meta[cols], function(x) is.character(x) || is.factor(x))]
      choices <- if (length(cat_cols) > 0) cat_cols else cols
      shiny::selectInput(ns("cohend_cond_var"), "Condition variable:",
                         choices = choices, selected = choices[1])
    })

    output$cohend_cond_class_ui <- shiny::renderUI({
      meta <- get_meta()
      if (is.null(meta) || is.null(input$cohend_cond_var) ||
          !input$cohend_cond_var %in% colnames(meta)) return(NULL)
      levs <- unique(as.character(meta[[input$cohend_cond_var]]))
      shiny::selectInput(ns("cohend_cond_class"), "Positive class:",
                         choices = levs, selected = levs[1])
    })

    output$roc_cond_var_ui <- shiny::renderUI({
      meta <- get_meta()
      if (is.null(meta)) return(NULL)
      cols <- setdiff(colnames(meta), colnames(meta)[1])
      cat_cols <- cols[sapply(meta[cols], function(x) is.character(x) || is.factor(x))]
      choices <- if (length(cat_cols) > 0) cat_cols else cols
      shiny::selectInput(ns("roc_cond_var"), "Condition variable:",
                         choices = choices, selected = choices[1])
    })

    output$roc_cond_class_ui <- shiny::renderUI({
      meta <- get_meta()
      if (is.null(meta) || is.null(input$roc_cond_var) ||
          !input$roc_cond_var %in% colnames(meta)) return(NULL)
      levs <- unique(as.character(meta[[input$roc_cond_var]]))
      shiny::selectInput(ns("roc_cond_class"), "Positive class:",
                         choices = levs, selected = levs[1])
    })

    # =========================================================================
    # ---- RUN OBSERVERS ------------------------------------------------------
    # =========================================================================

    shiny::observeEvent(input$run_pairwise, {
      gs <- get_gene_sets(); shiny::req(!is.null(gs), length(gs) >= 1)
      gs_char <- .gs_to_char_list(gs)
      metric  <- shiny::isolate(input$pw_metric)
      shiny::withProgress(message = "Computing pairwise overlap...", value = 0, {
        result <- tryCatch({
          shiny::incProgress(0.4)
          args <- list(signatures = gs_char, other_user_signatures = gs_char, metric = metric)
          if (metric == "odds_ratio") {
            expr_d <- get_expr(); shiny::req(!is.null(expr_d))
            args$universe <- rownames(expr_d)
          }
          shiny::incProgress(0.6)
          do.call(geneset_similarity, args)
        }, error = function(e) {
          shiny::showNotification(paste("Pairwise failed:", conditionMessage(e)),
                                  type = "error", duration = 10); NULL
        })
      })
      if (!is.null(result)) {
        result$run_params <- list(metric = metric)
        if (!is.null(result$data))
          pw_h(.dyn_h(length(unique(result$data$Compared_Signature)), base = 280, per_item = 45))
      }
      pairwise_result(result)
    })

    shiny::observeEvent(input$run_similarity, {
      shiny::req(get_gene_sets())
      similarity_result(NULL)   # clear stale result so the old table never persists
      gs_char <- .gs_to_char_list(get_gene_sets())
      metric       <- shiny::isolate(input$gs_metric)
      organism     <- shiny::isolate(input$gs_organism)
      if (is.null(organism) || !nzchar(organism)) organism <- "Homo sapiens"
      db_species   <- shiny::isolate(effective_db_species())
      # If collection picker is still rendering (NULL), derive the expected default
      # rather than cancelling with req() — this was causing first-run failures.
      collection <- shiny::isolate(input$gs_collection)
      if (is.null(collection) || !nzchar(collection)) {
        collection <- if (!is.null(db_species) && db_species == "MM") "MH" else "H"
      }
      # Guard against stale collection codes from the previous organism.
      # The collection renderUI updates asynchronously; if the user switches
      # organism and clicks Run before it finishes, the old collection code is read.
      #   • MM database  → valid codes start with "M" (MH, M2, M5 …)
      #   • HS / ortholog → valid codes start with "C" or equal "H" or "all"
      if (!is.null(db_species) && db_species == "MM" && !grepl("^M", collection)) {
        collection <- "MH"
        shiny::showNotification(
          "Collection reset to Mouse Hallmarks (MH) for the MM database.",
          type = "warning", duration = 5)
      } else if ((is.null(db_species) || db_species != "MM") && grepl("^M[^a-z]", collection)) {
        collection <- "H"
        shiny::showNotification(
          paste0("Collection reset to Hallmarks (H) — mouse-native codes (M-prefix) ",
                 "are not valid for ", organism, "."),
          type = "warning", duration = 5)
      }
      subset_raw   <- shiny::isolate(input$msig_subset)
      msig_subset  <- if (!is.null(subset_raw) && nzchar(trimws(subset_raw)))
        trimws(unlist(strsplit(subset_raw, "[,\n]+"))) else NULL
      subco_raw    <- shiny::isolate(input$gs_subcollection)
      subcollection <- if (!is.null(subco_raw) && nzchar(trimws(subco_raw))) trimws(subco_raw) else NULL
      shiny::withProgress(
        message = paste0("Fetching ", organism, " gene sets from MSigDB",
          if (!is.null(db_species)) paste0(" (", db_species, " database)") else "",
          "…"), value = 0, {
        result <- tryCatch({
          shiny::incProgress(0.3)
          args <- list(signatures = gs_char, collection = collection,
                       subcollection = subcollection, msig_subset = msig_subset, metric = metric,
                       organism = organism, db_species = db_species)
          if (metric == "jaccard") {
            args$jaccard_threshold <- shiny::isolate(input$jaccard_threshold)
          } else {
            expr_d <- get_expr(); shiny::req(!is.null(expr_d))
            args$universe       <- rownames(expr_d)
            args$or_threshold   <- shiny::isolate(input$or_threshold)
            args$pval_threshold <- shiny::isolate(input$or_pval)
          }
          shiny::incProgress(0.7)
          res <- do.call(geneset_similarity, args)
          if (!is.null(res)) {
            res$organism   <- organism
            res$run_params <- list(
              metric        = metric,
              organism      = organism,
              db_species    = db_species,
              collection    = collection,
              subcollection = subcollection,
              msig_subset   = msig_subset,
              jaccard_thr   = if (metric == "jaccard")   shiny::isolate(input$jaccard_threshold) else NULL,
              or_thr        = if (metric == "odds_ratio") shiny::isolate(input$or_threshold)      else NULL,
              pval_thr      = if (metric == "odds_ratio") shiny::isolate(input$or_pval)           else NULL
            )
          }
          res
        }, error = function(e) {
          msg <- conditionMessage(e)
          if (grepl("No signatures passed", msg, fixed = TRUE)) {
            shiny::showNotification(
              paste0("No pathways matched for ", organism, ". ",
                     "Try a different collection, lower your threshold, or check that ",
                     "the correct organism is selected."),
              type = "warning", duration = 10)
          } else {
            shiny::showNotification(paste("Similarity failed:", msg),
                                    type = "error", duration = 10)
          }
          NULL
        })
      })
      if (!is.null(result) && !is.null(result$data))
        sim_h(.dyn_h(length(unique(result$data$Compared_Signature)), base = 280, per_item = 32))
      similarity_result(result)
    })

    shiny::observeEvent(input$run_expression, {
      expr_d <- get_expr(); meta <- get_meta(); gs <- get_gene_sets()
      shiny::req(expr_d, meta, gs, input$selected_gs, input$group_var,
                 input$selected_genes, length(input$selected_genes) > 0)
      genes_plot <- shiny::isolate(input$selected_genes)
      group_var  <- shiny::isolate(input$group_var)
      color_var  <- shiny::isolate(input$color_var)
      color_var  <- if (!is.null(color_var) && nzchar(color_var)) color_var else NULL

      color_vals <- if (!is.null(color_var) && color_var %in% colnames(meta)) {
        levs <- unique(as.character(meta[[color_var]]))
        stats::setNames(.pp_palette[seq_along(levs)], levs)
      } else NULL

      n_levels  <- length(unique(as.character(meta[[group_var]])))
      ncol_opt  <- max(1L, min(5L, as.integer(ceiling(12 / n_levels))))
      n_genes   <- length(genes_plot)
      max_lbl   <- max(nchar(unique(as.character(meta[[group_var]]))), na.rm = TRUE)
      per_row_h <- 200L + ceiling(max_lbl * 4L)

      shiny::withProgress(message = "Generating violin plots...", value = 0.2, {
        result <- tryCatch({
          violin <- IndividualGenes_Violins(
            data             = expr_d,
            metadata         = meta,
            genes            = genes_plot,
            GroupingVariable = group_var,
            ColorVariable    = color_var,
            ColorValues      = color_vals,
            ncol             = ncol_opt,
            plot             = FALSE)
          shiny::incProgress(0.8)
          list(violin = violin, gs_name = shiny::isolate(input$selected_gs),
               n_genes = n_genes, group_var = group_var,
               ncol_opt = ncol_opt, color_var = color_var)
        }, error = function(e) {
          shiny::showNotification(paste("Violin plots failed:", conditionMessage(e)),
                                  type = "error", duration = 10); NULL
        })
      })
      if (!is.null(result)) {
        n_rows  <- ceiling(n_genes / ncol_opt)
        total_h <- min(20000L, max(350L, as.integer(n_rows * per_row_h + 80L)))
        expr_h(total_h)
      }
      expression_result(result)
    })

    shiny::observeEvent(input$run_heatmap, {
      expr_d <- get_expr(); meta <- get_meta(); gs <- get_gene_sets()
      shiny::req(expr_d, meta, gs, input$selected_gs, input$group_var,
                 input$selected_genes, length(input$selected_genes) > 0)
      genes_plot <- shiny::isolate(input$selected_genes)
      group_var  <- shiny::isolate(input$group_var)
      cl_rows    <- shiny::isolate(input$hm_cluster_rows)
      cl_cols    <- shiny::isolate(input$hm_cluster_cols)
      tsz        <- shiny::isolate(input$hm_titlesize)
      if (is.null(cl_rows)) cl_rows <- TRUE
      if (is.null(cl_cols)) cl_cols <- TRUE
      if (is.null(tsz) || is.na(tsz)) tsz <- 12

      levs      <- unique(as.character(meta[[group_var]]))
      annot_col <- list(stats::setNames(.pp_palette[seq_along(levs)], levs))
      names(annot_col) <- group_var

      shiny::withProgress(message = "Generating expression heatmap...", value = 0.3, {
        result <- tryCatch({
          hm <- ExpressionHeatmap(
            data              = expr_d,
            metadata          = meta,
            genes             = genes_plot,
            annotate.by       = group_var,
            annotation_colors = annot_col,
            colorlist         = list(low = "#4173B4", mid = "white", high = "#B44141"),
            cluster_rows      = cl_rows,
            cluster_columns   = cl_cols,
            title             = paste("Expression –", shiny::isolate(input$selected_gs)),
            titlesize         = tsz
          )
          shiny::incProgress(0.7)
          list(heatmap = hm, gs_name = shiny::isolate(input$selected_gs),
               n_genes = length(genes_plot))
        }, error = function(e) {
          shiny::showNotification(paste("Expression heatmap failed:", conditionMessage(e)),
                                  type = "error", duration = 10); NULL
        })
      })
      if (!is.null(result)) hm_h(.dyn_h(result$n_genes, base = 200, per_item = 22))
      heatmap_result(result)
    })

    shiny::observeEvent(input$run_corr, {
      expr_d <- get_expr(); meta <- get_meta()
      shiny::req(expr_d, meta, input$selected_genes, length(input$selected_genes) > 0)
      genes_plot  <- shiny::isolate(input$selected_genes)
      sep         <- shiny::isolate(input$corr_sep)
      sep_by      <- if (!is.null(sep) && nzchar(sep)) sep else NULL
      method      <- shiny::isolate(input$corr_method)
      show_colnam <- shiny::isolate(input$show_col_names)
      if (is.null(show_colnam)) show_colnam <- TRUE
      shiny::withProgress(message = "Computing correlations...", value = 0.3, {
        result <- tryCatch({
          # Run CorrelationHeatmap to get the data (corrmat) — store the raw
          # correlation matrices so font size can be applied at render time.
          res <- CorrelationHeatmap(
            data = expr_d, metadata = meta, genes = genes_plot,
            separate.by = sep_by, method = method,
            colorlist = list(low = "#4173B4", mid = "white", high = "#B44141"),
            limits_colorscale = c(-1, 0, 1),
            show_column_names = show_colnam,
            detailedresults   = FALSE)
          shiny::incProgress(0.7)
          # res$data: either a corrmat (no sep) or a named list of corrmat (sep)
          list(corrmat = res$data, n_genes = length(genes_plot),
               sep = sep_by, method = method, show_col_names = show_colnam)
        }, error = function(e) {
          shiny::showNotification(paste("Correlation failed:", conditionMessage(e)),
                                  type = "error", duration = 10); NULL
        })
      })
      if (!is.null(result)) corr_h(max(350L, as.integer(result$n_genes * 30 + 80)))
      corr_result(result)
    })

    shiny::observeEvent(input$run_cohend, {
      expr_d <- get_expr(); meta <- get_meta()
      shiny::req(expr_d, meta, input$selected_genes, length(input$selected_genes) > 0,
                 input$cohend_cond_var, input$cohend_cond_class)
      genes_all  <- shiny::isolate(input$selected_genes)
      max_g      <- shiny::isolate(input$cohend_max_genes)
      gene_sel   <- shiny::isolate(input$cohend_gene_sel)
      cond_var   <- shiny::isolate(input$cohend_cond_var)
      cond_class <- shiny::isolate(input$cohend_cond_class)

      # For "top_d" we run on all genes and post-filter by |d|.
      # For "first" / "random" we pre-slice before running.
      if (isTRUE(gene_sel == "top_d")) {
        genes_compute <- genes_all           # compute on the full selection
      } else {
        genes_compute <- genes_all
        if (!is.null(max_g) && max_g > 0 && length(genes_compute) > max_g) {
          genes_compute <- if (isTRUE(gene_sel == "random"))
            sample(genes_compute, max_g)
          else
            genes_compute[seq_len(max_g)]
        }
      }

      shiny::withProgress(message = "Computing Cohen's d...", value = 0.3, {
        result <- tryCatch({
          res <- CohenD_IndividualGenes(
            data = expr_d, metadata = meta, genes = genes_compute,
            condition_var = cond_var, class = cond_class, group_var = NULL,
            params = list(colors = .pp_palette[1], limits = NULL))
          shiny::incProgress(0.6)

          # Post-filter for "top_d": keep only the N genes with largest |d|.
          if (isTRUE(gene_sel == "top_d") && !is.null(max_g) && max_g > 0) {
            d_df <- if (is.list(res) && is.data.frame(res$data)) res$data else NULL
            if (!is.null(d_df) && "CohensD" %in% colnames(d_df) && nrow(d_df) > max_g) {
              top_genes <- d_df$Gene[order(abs(d_df$CohensD), decreasing = TRUE)][seq_len(max_g)]
              res$data  <- d_df[d_df$Gene %in% top_genes, , drop = FALSE]
            }
          }
          shiny::incProgress(0.4)
          n_shown <- if (is.list(res) && is.data.frame(res$data)) nrow(res$data) else length(genes_compute)
          list(res = res, n_genes = n_shown, cond_var = cond_var, cond_class = cond_class,
               gene_sel = gene_sel)
        }, error = function(e) {
          shiny::showNotification(paste("Cohen's d failed:", conditionMessage(e)),
                                  type = "error", duration = 10); NULL
        })
      })
      if (!is.null(result)) cohend_h(.dyn_h(result$n_genes, base = 200, per_item = 30))
      cohend_result(result)
    })

    shiny::observeEvent(input$run_pca, {
      expr_d <- get_expr(); meta <- get_meta()
      shiny::req(expr_d, meta, input$selected_genes, length(input$selected_genes) > 0)
      genes_plot <- shiny::isolate(input$selected_genes)
      cv         <- shiny::isolate(input$color_var)
      color_var  <- if (!is.null(cv) && nzchar(trimws(cv))) cv else NULL
      shiny::withProgress(message = "Running PCA on genes...", value = 0.3, {
        result <- tryCatch({
          res <- plotPCA(data = expr_d, metadata = meta, genes = genes_plot,
                         scale = FALSE, center = TRUE, PCs = list(c(1, 2)),
                         ColorVariable = color_var, ColorValues = NULL,
                         pointSize = 4, legend_nrow = 1, ncol = 1, nrow = 1)
          shiny::incProgress(0.7)
          list(res = res, color_var = color_var)
        }, error = function(e) {
          shiny::showNotification(paste("PCA failed:", conditionMessage(e)),
                                  type = "error", duration = 10); NULL
        })
      })
      pca_result(result)
    })

    shiny::observeEvent(input$run_rocauc, {
      expr_d <- get_expr(); meta <- get_meta()
      shiny::req(expr_d, meta, input$selected_genes, length(input$selected_genes) > 0,
                 input$roc_cond_var, input$roc_cond_class)
      genes_plot <- shiny::isolate(input$selected_genes)
      max_g      <- shiny::isolate(input$roc_max_genes)
      gene_sel   <- shiny::isolate(input$roc_gene_sel)
      if (!is.null(max_g) && max_g > 0 && length(genes_plot) > max_g) {
        genes_plot <- if (isTRUE(gene_sel == "random"))
          sample(genes_plot, max_g)
        else
          genes_plot[seq_len(max_g)]
      }
      n_genes    <- length(genes_plot)
      roc_ncol   <- max(1L, min(4L, as.integer(ceiling(sqrt(n_genes)))))
      cond_var   <- shiny::isolate(input$roc_cond_var)
      cond_class <- shiny::isolate(input$roc_cond_class)
      # Build a stable gene→colour mapping in the original gene order so that
      # the same gene always gets the same colour in BOTH the ROC curves and AUC plots.
      stable_colors <- stats::setNames(rep_len(.pp_palette, length(genes_plot)), genes_plot)

      shiny::withProgress(message = "Computing ROC / AUC...", value = 0.3, {
        tmp_pdf <- tempfile(fileext = ".pdf")
        grDevices::pdf(tmp_pdf, width = 8, height = 6)
        result <- tryCatch({
          res <- ROCandAUCplot(
            data = expr_d, metadata = meta, genes = genes_plot,
            condition_var = cond_var, class = cond_class,
            group_var = NULL, plot_type = "all",
            roc_params = list(ncol = roc_ncol))
          shiny::incProgress(0.7)
          list(res = res, n_genes = n_genes, roc_ncol = roc_ncol,
               cond_var = cond_var, cond_class = cond_class,
               stable_colors = stable_colors)
        }, error = function(e) {
          shiny::showNotification(paste("ROC/AUC failed:", conditionMessage(e)),
                                  type = "error", duration = 10); NULL
        })
        try(grDevices::dev.off(), silent = TRUE)
        unlink(tmp_pdf)
      })

      if (!is.null(result)) {
        n_rows_roc <- ceiling(result$n_genes / result$roc_ncol)
        roc_curves_h(max(350L, as.integer(n_rows_roc * 280 + 60)))
        rocauc_h(.dyn_h(result$n_genes, base = 200, per_item = 30))
      }
      rocauc_result(result)
    })

    # =========================================================================
    # ---- RENDERS ------------------------------------------------------------
    # =========================================================================

    # ---- Pairwise -----------------------------------------------------------

    output$pairwise_status_ui <- shiny::renderUI({
      if (is.null(pairwise_result()))
        return(shiny::div(class = "alert alert-info", style = "margin:10px 0;",
                   shiny::icon("circle-info"),
                   " Click ", shiny::strong("Run Pairwise Overlap"),
                   " to compare your gene sets against each other."))
      NULL
    })

    output$pairwise_plot_ui <- shiny::renderUI({
      res <- pairwise_result()
      if (is.null(res)) return(NULL)
      shiny::tagList(shiny::h5("Pairwise Gene Set Overlap", style = "margin:12px 0 4px;"),
                     plotly::plotlyOutput(ns("pairwise_plot"), height = paste0(pw_h(), "px")))
    })

    output$pairwise_plot <- plotly::renderPlotly({
      res  <- pairwise_result(); shiny::req(!is.null(res))
      font <- input$pw_fontsize
      metric <- shiny::isolate(input$pw_metric)
      p <- if (metric == "odds_ratio" && !is.null(res$data))
             .build_or_heatmap(res$data, font = font, title = "Pairwise Overlap (Odds Ratio)")
           else if (!is.null(res$data))
             .build_ji_heatmap(res$data, font = font, title = "Pairwise Overlap (Jaccard Index)")
           else shiny::req(FALSE)
      plotly::ggplotly(p, tooltip = "text") |>
        plotly::layout(hoverlabel = list(bgcolor = "white"), height = pw_h())
    })

    # Table header rendered separately; DTOutput is STATIC in UI so binding
    # always exists in the DOM. outputOptions(suspendWhenHidden=FALSE) ensures
    # renderDT fires even when the tab is not currently active.
    output$pairwise_table_hdr_ui <- shiny::renderUI({
      res <- pairwise_result()
      if (is.null(res) || !is.data.frame(res$data)) return(NULL)
      p   <- res$run_params
      metric_label <- if (!is.null(p$metric) && p$metric == "odds_ratio") "Odds Ratio" else "Jaccard Index"
      n_rows <- nrow(res$data)
      shiny::tagList(
        shiny::hr(),
        shiny::h5("Pairwise Overlap Data"),
        shiny::div(
          style = "margin-bottom:6px; font-size:0.85em; color:#555;",
          shiny::icon("circle-info", style = "color:#0d6efd; margin-right:4px;"),
          sprintf("Showing %d pairs — metric: %s. No threshold filtering applied.",
                  n_rows, metric_label)
        )
      )
    })
    output$pairwise_table <- DT::renderDT({
      res <- pairwise_result()
      if (is.null(res) || !is.data.frame(res$data)) return(NULL)
      metric <- shiny::isolate(input$pw_metric)
      df     <- res$data
      metric_label <- if (isTRUE(metric == "odds_ratio")) "OR" else "Jaccard"
      if ("Score" %in% colnames(df)) colnames(df)[colnames(df) == "Score"] <- metric_label
      df[sapply(df, is.numeric)] <- lapply(df[sapply(df, is.numeric)], round, 4)
      DT::datatable(df, rownames = FALSE, filter = "top",
                    extensions = "Buttons",
                    options = list(pageLength = 15, scrollX = TRUE,
                                   dom = "Bfrtip", buttons = c("csv", "excel")))
    })
    shiny::outputOptions(output, "pairwise_table", suspendWhenHidden = FALSE)

    # ---- Similarity vs MSigDB -----------------------------------------------

    output$similarity_status_ui <- shiny::renderUI({
      if (is.null(similarity_result()))
        return(shiny::div(class = "alert alert-info", style = "margin:10px 0;",
                   shiny::icon("circle-info"),
                   " Configure the settings above and click ",
                   shiny::strong("Run Similarity"), "."))
      NULL
    })

    output$similarity_plot_ui <- shiny::renderUI({
      res <- similarity_result()
      if (is.null(res)) return(NULL)
      shiny::tagList(shiny::h5("Similarity Heatmap", style = "margin:12px 0 4px;"),
                     plotly::plotlyOutput(ns("similarity_plot"), height = paste0(sim_h(), "px")))
    })

    output$similarity_plot <- plotly::renderPlotly({
      res  <- similarity_result(); shiny::req(!is.null(res))
      font <- input$sim_fontsize
      metric <- shiny::isolate(input$gs_metric)
      org_label <- if (!is.null(res$organism)) res$organism else "Homo sapiens"
      p <- if (metric == "odds_ratio" && !is.null(res$data))
             .build_or_heatmap(res$data, font = font,
               title = paste0("Similarity vs MSigDB (Odds Ratio) — ", org_label))
           else if (!is.null(res$data))
             .build_ji_heatmap(res$data, font = font,
               title = paste0("Similarity vs MSigDB (Jaccard Index) — ", org_label))
           else shiny::req(FALSE)
      plotly::ggplotly(p, tooltip = "text") |>
        plotly::layout(hoverlabel = list(bgcolor = "white"), height = sim_h())
    })

    # Similarity table header — shows every active filter/threshold from the last run
    output$similarity_table_hdr_ui <- shiny::renderUI({
      res <- similarity_result()
      if (is.null(res) || !is.data.frame(res$data) || nrow(res$data) == 0) return(NULL)
      p <- res$run_params
      if (is.null(p)) return(shiny::tagList(shiny::hr(), shiny::h5("Similarity Data Table")))

      badge  <- function(txt, cls = "bg-secondary") {
        shiny::tags$span(class = paste("badge me-1", cls), style = "font-weight:500;", txt)
      }
      wbadge <- function(txt) badge(txt, "bg-warning text-dark")

      tags <- list()

      # ── Organism + database ─────────────────────────────────────────────────
      org_label <- p$organism %||% "Homo sapiens"
      db_sfx    <- if (!is.null(p$db_species) && p$db_species == "MM") " · MM native"
                   else if (org_label != "Homo sapiens")                " · HS orthologs"
                   else                                                  ""
      tags <- c(tags, list(badge(paste0(org_label, db_sfx), "bg-primary")))

      # ── Collection / subcollection ──────────────────────────────────────────
      coll_disp <- if (!is.null(p$collection) && p$collection != "all") p$collection
                   else "All collections"
      tags <- c(tags, list(badge(paste0("Collection: ", coll_disp))))
      if (!is.null(p$subcollection) && nzchar(p$subcollection))
        tags <- c(tags, list(badge(paste0("Sub-collection: ", p$subcollection))))

      # ── Pathway name filter ─────────────────────────────────────────────────
      if (!is.null(p$msig_subset) && length(p$msig_subset) > 0)
        tags <- c(tags, list(wbadge(paste0("Pathway filter: ",
                                           paste(p$msig_subset, collapse = ", ")))))

      # ── Metric ─────────────────────────────────────────────────────────────
      metric_disp <- if (!is.null(p$metric) && p$metric == "odds_ratio") "Odds Ratio" else "Jaccard Index"
      tags <- c(tags, list(badge(paste0("Metric: ", metric_disp))))

      # ── Thresholds (all that are applicable) ────────────────────────────────
      if (!is.null(p$metric) && p$metric == "jaccard") {
        thr <- p$jaccard_thr %||% 0
        tags <- c(tags, list(wbadge(
          if (thr > 0) paste0("Jaccard ≥ ", thr, " (pathway kept if ANY signature qualifies)")
          else         "Jaccard threshold: none (all pathways kept)"
        )))
      } else if (!is.null(p$metric) && p$metric == "odds_ratio") {
        or_thr   <- p$or_thr   %||% 1
        pval_thr <- p$pval_thr %||% 0.05
        tags <- c(tags, list(wbadge(paste0("OR ≥ ", or_thr))))
        tags <- c(tags, list(wbadge(paste0("adj. p-value ≤ ", pval_thr))))
        tags <- c(tags, list(badge("Pathway kept if ANY signature passes both OR + p-value", "bg-secondary")))
      }

      n_pathways <- length(unique(res$data$Compared_Signature))
      n_pairs    <- nrow(res$data)
      shiny::tagList(
        shiny::hr(),
        shiny::h5("Similarity Data Table"),
        shiny::div(
          style = "margin-bottom:8px; font-size:0.85em;",
          shiny::icon("circle-info", style = "color:#0d6efd; margin-right:4px;"),
          shiny::tags$span(style = "color:#555; margin-right:6px;",
            sprintf("%d pathway(s), %d row(s) — run with:", n_pathways, n_pairs)),
          do.call(shiny::tagList, tags)
        )
      )
    })

    # Similarity table — wrapped in renderUI so the DT widget is fully destroyed
    # and recreated each time the result changes. This is the reliable pattern:
    # req() inside renderDT causes a "silent error" that leaves the old DT in
    # place; putting DTOutput inside renderUI forces a fresh DOM binding instead.
    output$similarity_table_container <- shiny::renderUI({
      res <- similarity_result()
      if (is.null(res) || !is.data.frame(res$data) || nrow(res$data) == 0) return(NULL)
      DT::DTOutput(ns("similarity_table"))
    })

    output$similarity_table <- DT::renderDT({
      res <- similarity_result()
      shiny::req(!is.null(res), is.data.frame(res$data), nrow(res$data) > 0)
      df <- res$data

      # Infer metric from the data (Pval is NA for Jaccard, numeric for OR)
      is_or        <- "Pval" %in% colnames(df) && any(!is.na(df$Pval))
      metric_label <- if (is_or) "Odds Ratio" else "Jaccard Index"

      # Rename columns for display
      colnames(df)[colnames(df) == "Reference_Signature"] <- "Your Gene Set"
      colnames(df)[colnames(df) == "Compared_Signature"]  <- "MSigDB Pathway"
      colnames(df)[colnames(df) == "Score"]               <- metric_label
      if ("Pval" %in% colnames(df)) colnames(df)[colnames(df) == "Pval"] <- "Adj. p-value"

      # Round numerics
      df[sapply(df, is.numeric)] <- lapply(df[sapply(df, is.numeric)],
                                            function(x) round(x, 4))

      score_col <- which(colnames(df) == metric_label)
      DT::datatable(
        df,
        rownames   = FALSE,
        filter     = "top",
        extensions = "Buttons",
        options    = list(
          pageLength = 20,
          scrollX    = TRUE,
          dom        = "Bfrtip",
          buttons    = c("csv", "excel"),
          order      = if (length(score_col) == 1L) list(list(score_col - 1L, "desc")) else list()
        )
      )
    }, server = FALSE)

    # ---- Expression (violins + heatmap) ------------------------------------

    output$violin_status_ui <- shiny::renderUI({
      if (is.null(expression_result()))
        return(shiny::div(class = "alert alert-info", style = "margin:10px 0;",
                   shiny::icon("circle-info"),
                   " Select a gene set and grouping variable in the sidebar,",
                   " then click ", shiny::strong("Plot Violins"), "."))
      NULL
    })

    output$expr_plot_ui <- shiny::renderUI({
      res <- expression_result()
      if (is.null(res)) return(NULL)
      shiny::tagList(
        shiny::h5(paste("Per-gene Violin Plots –", res$gs_name), style = "margin:12px 0 4px;"),
        shiny::plotOutput(ns("violin_plot"), height = paste0(expr_h(), "px")))
    })

    output$violin_plot <- shiny::renderPlot({
      res <- expression_result(); shiny::req(!is.null(res), !is.null(res$violin))
      font <- input$expr_fontsize
      p    <- res$violin$plot
      for (i in seq_along(p$layers)) {
        layer_class <- class(p$layers[[i]]$geom)[1]
        if (layer_class %in% c("GeomJitter", "GeomPoint")) {
          mapped_aes <- names(p$layers[[i]]$mapping)
          if (!any(c("colour", "color") %in% mapped_aes))
            p$layers[[i]]$aes_params$colour <- "grey55"
        }
      }
      print(p + ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = font),
        axis.text.y = ggplot2::element_text(size = font),
        axis.title  = ggplot2::element_text(size = font + 1),
        strip.text  = ggplot2::element_text(size = font + 1, face = "bold"),
        legend.text = ggplot2::element_text(size = font)
      ))
    }, height = function() expr_h(), res = 150)

    output$heatmap_status_ui <- shiny::renderUI({
      if (is.null(heatmap_result()) && !is.null(expression_result()))
        return(shiny::div(class = "alert alert-secondary",
                   style = "margin:12px 0; font-size:0.88em;",
                   shiny::icon("circle-info"),
                   " Click ", shiny::strong("Plot Expression Heatmap"),
                   " to view a sample × gene heatmap."))
      NULL
    })

    output$heatmap_plot_ui <- shiny::renderUI({
      res <- heatmap_result()
      if (is.null(res)) return(NULL)
      shiny::tagList(
        shiny::h5(paste("Expression Heatmap –", res$gs_name), style = "margin:12px 0 4px;"),
        shiny::plotOutput(ns("heatmap_plot"), height = paste0(hm_h(), "px")))
    })

    output$heatmap_plot <- shiny::renderPlot({
      res <- heatmap_result(); shiny::req(!is.null(res), !is.null(res$heatmap))
      ComplexHeatmap::draw(res$heatmap$plot)
    }, height = function() hm_h(), res = 150)

    # ---- Correlation --------------------------------------------------------

    output$corr_status_ui <- shiny::renderUI({
      if (is.null(corr_result()))
        return(shiny::div(class = "alert alert-info", style = "margin:10px 0;",
                   shiny::icon("circle-info"),
                   " Configure settings above and click ",
                   shiny::strong("Plot Correlation Heatmap"), "."))
      NULL
    })

    output$corr_plot_ui <- shiny::renderUI({
      res <- corr_result()
      if (is.null(res)) return(NULL)
      shiny::tagList(shiny::h5("Pairwise Gene Correlation", style = "margin:12px 0 4px;"),
                     shiny::plotOutput(ns("corr_plot"), height = paste0(corr_h(), "px")))
    })

    # renderPlot reads input$corr_fontsize and input$show_col_names directly,
    # so changes to those inputs cause an instant re-draw WITHOUT re-computing
    # the correlation matrices (expensive step already done in run_corr).
    output$corr_plot <- shiny::renderPlot({
      res <- corr_result(); shiny::req(!is.null(res), !is.null(res$corrmat))
      font     <- if (!is.null(input$corr_fontsize))  input$corr_fontsize  else 10L
      show_col <- if (!is.null(input$show_col_names)) input$show_col_names else TRUE
      tryCatch(
        .build_corr_heatmap(res$corrmat, method = res$method,
                             font = font, show_col_names = show_col),
        error = function(e) shiny::showNotification(
          paste("Correlation plot failed:", e$message), type = "warning", duration = 8))
    }, height = function() corr_h(), res = 150)

    # ---- Cohen's d ----------------------------------------------------------

    output$cohend_status_ui <- shiny::renderUI({
      if (is.null(cohend_result()))
        return(shiny::div(class = "alert alert-info", style = "margin:10px 0;",
                   shiny::icon("circle-info"),
                   " Set condition variable and positive class, then click ",
                   shiny::strong("Plot Cohen's d"), "."))
      NULL
    })

    output$cohend_plot_ui <- shiny::renderUI({
      res <- cohend_result()
      if (is.null(res)) return(NULL)
      shiny::tagList(
        shiny::h5(paste0("Cohen's d per Gene — ", res$cond_class,
                         " vs rest (", res$cond_var, ")"), style = "margin:12px 0 4px;"),
        plotly::plotlyOutput(ns("cohend_plot"), height = paste0(cohend_h(), "px")))
    })

    output$cohend_plot <- plotly::renderPlotly({
      res  <- cohend_result(); shiny::req(!is.null(res))
      font <- input$cohend_fontsize
      df   <- if (is.list(res$res) && is.data.frame(res$res$data)) res$res$data else NULL
      shiny::req(!is.null(df), "CohensD" %in% colnames(df))
      df$hover <- paste0("Gene: ", df$Gene, "<br>Cohen's d: ", sprintf("%.3f", df$CohensD))
      df$Gene  <- stats::reorder(df$Gene, df$CohensD)
      levs     <- levels(df$Gene)
      col_vals <- stats::setNames(rep_len(.pp_palette, length(levs)), levs)
      p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$CohensD, y = .data$Gene,
                                             fill = .data$Gene, text = .data$hover)) +
        ggplot2::geom_col() +
        ggplot2::scale_fill_manual(values = col_vals) +
        ggplot2::labs(x = "Cohen's d", y = "Gene") +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "none",
                       axis.text = ggplot2::element_text(size = font),
                       axis.title = ggplot2::element_text(size = font + 1))
      plotly::ggplotly(p, tooltip = "text") |> plotly::layout(height = cohend_h())
    })

    # ---- Gene PCA -----------------------------------------------------------

    output$pca_status_ui <- shiny::renderUI({
      if (is.null(pca_result()))
        return(shiny::div(class = "alert alert-info", style = "margin:10px 0;",
                   shiny::icon("circle-info"),
                   " Click ", shiny::strong("Plot Gene PCA"), ".",
                   " Set 'Colour by' in the sidebar to colour points by a metadata variable."))
      NULL
    })

    output$pca_plot_ui <- shiny::renderUI({
      res <- pca_result()
      if (is.null(res)) return(NULL)
      shiny::tagList(shiny::h5("Gene PCA", style = "margin:12px 0 4px;"),
                     plotly::plotlyOutput(ns("pca_plot"), height = paste0(pca_h(), "px")))
    })

    output$pca_plot <- plotly::renderPlotly({
      res       <- pca_result(); shiny::req(!is.null(res))
      font      <- input$pca_fontsize
      raw       <- res$res
      color_var <- res$color_var
      pca_df    <- raw$data
      shiny::req(!is.null(pca_df), "PC1" %in% colnames(pca_df), "PC2" %in% colnames(pca_df))

      sample_lbl <- if ("Sample" %in% colnames(pca_df)) as.character(pca_df$Sample) else
                      as.character(seq_len(nrow(pca_df)))

      if (!is.null(color_var) && color_var %in% colnames(pca_df)) {
        levs       <- unique(as.character(pca_df[[color_var]]))
        color_vals <- stats::setNames(.pp_palette[seq_along(levs)], levs)
        pca_df[[color_var]] <- factor(pca_df[[color_var]], levels = levs)
        hover_text <- paste0("Sample: ", sample_lbl,
                             "<br>", color_var, ": ", as.character(pca_df[[color_var]]))
        p <- ggplot2::ggplot(pca_df,
               ggplot2::aes(x = .data$PC1, y = .data$PC2,
                            fill = .data[[color_var]], text = hover_text)) +
          ggplot2::geom_point(size = 4, alpha = 0.8, shape = 21, color = "black") +
          ggplot2::scale_fill_manual(values = color_vals, name = color_var)
      } else {
        p <- ggplot2::ggplot(pca_df,
               ggplot2::aes(x = .data$PC1, y = .data$PC2,
                            text = paste0("Sample: ", sample_lbl))) +
          ggplot2::geom_point(size = 4, alpha = 0.8, shape = 21,
                              color = "black", fill = .pp_palette[1])
      }

      p <- p +
        ggplot2::labs(x = "PC1", y = "PC2", title = "Gene PCA") +
        ggplot2::geom_vline(xintercept = 0, linetype = "dotted", color = "grey60") +
        ggplot2::geom_hline(yintercept = 0, linetype = "dotted", color = "grey60") +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text = ggplot2::element_text(size = font),
                       axis.title = ggplot2::element_text(size = font + 1),
                       legend.text = ggplot2::element_text(size = font),
                       legend.position = "bottom",
                       plot.title = ggplot2::element_text(hjust = 0.5))

      plotly::ggplotly(p, tooltip = "text") |>
        plotly::layout(height = pca_h(), legend = list(orientation = "h"))
    })

    # ---- ROC / AUC ----------------------------------------------------------

    output$rocauc_status_ui <- shiny::renderUI({
      if (is.null(rocauc_result()))
        return(shiny::div(class = "alert alert-info", style = "margin:10px 0;",
                   shiny::icon("circle-info"),
                   " Set condition variable and positive class, then click ",
                   shiny::strong("Plot ROC / AUC"), ".",
                   " ROC curves and an AUC barplot will both be shown."))
      NULL
    })

    output$roc_curves_ui <- shiny::renderUI({
      res <- rocauc_result()
      if (is.null(res)) return(NULL)
      shiny::tagList(
        shiny::h5(paste0("ROC Curves — ", res$cond_class, " vs rest (", res$cond_var, ")"),
                  style = "margin:12px 0 4px;"),
        shiny::plotOutput(ns("roc_curves_plot"), height = paste0(roc_curves_h(), "px")))
    })

    output$roc_curves_plot <- shiny::renderPlot({
      res <- rocauc_result(); shiny::req(!is.null(res))
      raw <- res$res
      roc_p   <- raw$roc_plot
      shiny::req(!is.null(roc_p), inherits(roc_p, c("gg", "ggplot")))
      roc_df <- roc_p$data
      shiny::req(!is.null(roc_df), nrow(roc_df) > 0)

      n_genes_shown <- length(unique(roc_df$Gene))
      ncol_shown    <- res$roc_ncol
      font          <- input$roc_fontsize

      # Use the stable colour map built at compute time so colours match the AUC barplot.
      col_vals <- res$stable_colors
      if (is.null(col_vals)) {
        genes_uniq <- unique(roc_df$Gene)
        col_vals   <- stats::setNames(rep_len(.pp_palette, length(genes_uniq)), genes_uniq)
      }

      p <- ggplot2::ggplot(roc_df, ggplot2::aes(
             x = .data$FPR, y = .data$TPR,
             color = .data$Gene, group = .data$Gene)) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::scale_color_manual(values = col_vals) +
        ggplot2::facet_wrap(~ .data$Gene, ncol = ncol_shown) +
        ggplot2::geom_abline(linetype = "dashed", color = "grey70") +
        ggplot2::labs(x = "False Positive Rate (1 - Specificity)",
                      y = "True Positive Rate (Sensitivity)",
                      title = paste0("ROC Curves — ", res$cond_class, " vs rest")) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "none",
                       strip.text  = ggplot2::element_text(size = font, face = "bold"),
                       axis.text   = ggplot2::element_text(size = font - 1),
                       axis.title  = ggplot2::element_text(size = font),
                       plot.title  = ggplot2::element_text(hjust = 0.5))

      n_rows_shown <- ceiling(n_genes_shown / max(1L, ncol_shown))
      roc_curves_h(max(350L, as.integer(n_rows_shown * 280 + 60)))
      print(p)
    }, height = function() roc_curves_h(), res = 150)

    output$rocauc_plot_ui <- shiny::renderUI({
      res <- rocauc_result()
      if (is.null(res)) return(NULL)
      shiny::tagList(
        shiny::h5(paste0("AUC per Gene — ", res$cond_class,
                         " vs rest (", res$cond_var, ")"), style = "margin:20px 0 4px;"),
        plotly::plotlyOutput(ns("rocauc_plot"), height = paste0(rocauc_h(), "px")))
    })

    output$rocauc_plot <- plotly::renderPlotly({
      res  <- rocauc_result(); shiny::req(!is.null(res))
      font <- input$roc_fontsize
      raw  <- res$res
      auc_df <- raw$auc_values
      shiny::req(!is.null(auc_df), is.data.frame(auc_df), nrow(auc_df) > 0)

      genes_ord   <- auc_df$Gene[order(auc_df$AUC)]
      auc_df$Gene <- factor(auc_df$Gene, levels = genes_ord)
      # Reuse the stable colour map from compute time — same colours as ROC curves.
      color_vals <- res$stable_colors
      if (is.null(color_vals))
        color_vals <- stats::setNames(rep_len(.pp_palette, length(genes_ord)), genes_ord)
      auc_df$hover <- paste0(as.character(auc_df$Gene),
                              "<br>AUC: ", sprintf("%.3f", auc_df$AUC))

      p <- ggplot2::ggplot(auc_df,
             ggplot2::aes(x = .data$AUC, y = .data$Gene,
                          fill = .data$Gene, text = .data$hover)) +
        ggplot2::geom_col() +
        ggplot2::scale_fill_manual(values = color_vals) +
        ggplot2::geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey50") +
        ggplot2::coord_cartesian(xlim = c(0, 1)) +
        ggplot2::labs(x = "AUC", y = "Gene",
                      title = paste0("AUC per Gene — ", res$cond_class, " vs rest")) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "none",
                       axis.text  = ggplot2::element_text(size = font),
                       axis.title = ggplot2::element_text(size = font + 1),
                       plot.title = ggplot2::element_text(hjust = 0.5))

      plotly::ggplotly(p, tooltip = "text") |> plotly::layout(height = rocauc_h())
    })

  })
}
