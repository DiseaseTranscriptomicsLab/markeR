# =============================================================================
# Gene Sets Tab – Shiny Module
# =============================================================================

# ---- Helpers -----------------------------------------------------------------

if (!exists("%||%", mode = "function"))
  `%||%` <- function(a, b) if (!is.null(a)) a else b

# Compute (x0,x1,y0,y1) segment coordinates for drawing a dendrogram.
# Leaves are placed at integers 1..n in the visual order given by hc$order.
.hclust_to_segments <- function(hc) {
  n <- length(hc$order)   # always set; hc$labels can be NULL for unlabelled dists
  if (n < 2L)
    return(data.frame(x0 = numeric(), x1 = numeric(), y0 = numeric(), y1 = numeric()))
  # Map original (1-based) label index → visual position 1..n
  lpos <- integer(n); lpos[hc$order] <- seq_len(n)
  node_x <- numeric(n - 1L)
  node_h <- hc$height
  x0 <- x1 <- y0 <- y1 <- numeric(3L * (n - 1L))
  for (i in seq_len(n - 1L)) {
    l <- hc$merge[i, 1L]; r <- hc$merge[i, 2L]
    xl <- if (l < 0L) lpos[[-l]] else node_x[[l]]
    xr <- if (r < 0L) lpos[[-r]] else node_x[[r]]
    yl <- if (l < 0L) 0 else node_h[[l]]
    yr <- if (r < 0L) 0 else node_h[[r]]
    node_x[[i]] <- (xl + xr) / 2L
    hi <- node_h[[i]]
    j  <- 3L * (i - 1L)
    x0[[j+1L]]<-xl; x1[[j+1L]]<-xl; y0[[j+1L]]<-yl; y1[[j+1L]]<-hi  # left arm
    x0[[j+2L]]<-xl; x1[[j+2L]]<-xr; y0[[j+2L]]<-hi; y1[[j+2L]]<-hi  # crossbar
    x0[[j+3L]]<-xr; x1[[j+3L]]<-xr; y0[[j+3L]]<-yr; y1[[j+3L]]<-hi  # right arm
  }
  data.frame(x0=x0, x1=x1, y0=y0, y1=y1)
}

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
                                 show_row_names = TRUE, show_col_names = TRUE,
                                 cluster_rows = TRUE, cluster_columns = TRUE) {
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
      cluster_rows      = cluster_rows,
      cluster_columns   = cluster_columns,
      show_row_names    = show_row_names,
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

# Helper: nav_panel title with a small ℹ tooltip after the name.
.tab_title <- function(name, info) {
  shiny::tagList(
    name,
    bslib::tooltip(
      shiny::icon("circle-info",
                  style = "color:#bbb; cursor:help; font-size:0.78em; margin-left:5px; vertical-align:middle;"),
      shiny::tags$div(style = "max-width:320px; font-size:0.83em; line-height:1.5; text-align:left;", info),
      placement = "bottom"
    )
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

  # Custom JS handler for tab visibility — bypasses bslib nav_hide/nav_show
  # timing issues by directly manipulating the DOM.
  # Registered at script-load time (no navset-init dependency).
  js_ready_script <- shiny::tags$script(shiny::HTML(sprintf(
    # markeR.setTabVisible: show or hide a bslib nav-tab button + handle active switching
    "
    Shiny.addCustomMessageHandler('markeR.setTabVisible', function(msg) {
      // msg: { navsetId, tabValue, visible }
      // Find the navset container by id (may be the <ul> or a wrapper)
      var root = document.getElementById(msg.navsetId);
      if (!root) {
        // Fall back: search the whole document for a [data-navs-id] or similar
        root = document;
      }

      // Find every nav-link with this data-value (handles <a> and <button>)
      var escapedVal = msg.tabValue.replace(/'/g, \"\\\\'\");
      var btns = root.querySelectorAll('[data-value=\"' + msg.tabValue + '\"]');
      // Filter to only nav-link elements (skip tab-pane divs)
      var navBtns = Array.from(btns).filter(function(el) {
        return el.classList.contains('nav-link');
      });

      navBtns.forEach(function(btn) {
        var li = btn.closest('.nav-item') || btn.parentElement;
        if (msg.visible) {
          li.style.removeProperty('display');
          btn.style.removeProperty('display');
        } else {
          // If hiding the active tab, switch to first remaining visible tab first
          if (btn.classList.contains('active')) {
            var container = btn.closest('.nav, ul') || root;
            var others = Array.from(
              container.querySelectorAll('.nav-link:not([data-value=\"' + msg.tabValue + '\"])')
            ).filter(function(el) {
              var p = el.closest('.nav-item') || el.parentElement;
              return p.style.display !== 'none';
            });
            if (others.length > 0) {
              if (typeof bootstrap !== 'undefined' && bootstrap.Tab) {
                new bootstrap.Tab(others[0]).show();
              } else {
                others[0].click();
              }
            }
          }
          li.style.display = 'none';
          btn.style.display = 'none';
        }
      });
    });

    // Signal server once connected (JS is fully ready)
    $(document).on('shiny:connected', function() {
      Shiny.setInputValue('%s', true, {priority: 'event'});
    });
    ",
    ns("js_initialized")
  )))

  bslib::layout_sidebar(

    # ---- Sidebar: shared gene/sample selection ------------------------------
    sidebar = bslib::sidebar(
      width = 310,
      shiny::div(
        style = "padding-bottom:6px;",
        shiny::h4("Gene Sets", style = "font-weight:700; color:#801F4F; margin-bottom:4px;"),
        shiny::p("Explore gene set composition, overlap with MSigDB, cross-set similarity, and per-gene expression profiles.",
                 style = "color:#6c757d; font-size:0.85em;"),
        shiny::p(style = "font-size:0.78em; color:#888; margin:0;",
          "All analyses use only the genes from the loaded gene sets, not all genes in the dataset."),
        shiny::hr()
      ),
      shiny::tags$small(shiny::tags$b("Active gene set")),
      shiny::helpText(style = "margin-bottom:4px; font-size:0.8em; color:#888;",
        shiny::icon("circle-info"), " Pairwise Overlap and Similarity use all loaded gene sets."),
      shiny::uiOutput(ns("gs_picker_ui")),
      shiny::uiOutput(ns("gs_overlap_warn_ui")),
      shiny::uiOutput(ns("gs_composition_ui"))
    ),

    # ---- Main panel ---------------------------------------------------------
    bslib::navset_card_tab(id = ns("gs_navset"),

      # ---- 1. Pairwise Overlap -----------------------------------------------
      bslib::nav_panel(
        .tab_title("Pairwise Overlap",
          "Computes pairwise gene overlap between all loaded gene sets using the Jaccard index or odds ratio. Useful for understanding redundancy and identifying highly similar or complementary signatures."),
        .tab_settings(
          # ── Metric (always visible) ───────────────────────────────────────────
          shiny::div(
            style = "margin-bottom:8px;",
            shiny::tags$label(
              style = "font-size:0.85em; font-weight:600;",
              "Metric ",
              bslib::tooltip(
                shiny::icon("circle-question", style = "color:#6c757d; cursor:help;"),
                shiny::tags$span(
                  shiny::tags$b("Jaccard index:"), " |A ∩ B| / |A ∪ B| — fraction of shared genes.",
                  " JI = 0: no overlap; JI = 1: identical sets.",
                  shiny::tags$br(),
                  shiny::tags$b("Odds ratio:"), " enrichment of overlap relative to chance.",
                  " OR = 1: chance overlap; OR > 1: enriched. Self-comparisons = Inf."
                ),
                placement = "right"
              )
            ),
            shiny::radioButtons(ns("pw_metric"), label = NULL,
              choices  = c("Jaccard index" = "jaccard", "Odds ratio" = "odds_ratio"),
              selected = "jaccard", inline = TRUE)
          ),
          # ── Display options (collapsible) ─────────────────────────────────────
          shiny::tags$details(
            style = "margin-bottom:10px;",
            shiny::tags$summary(
              style = "cursor:pointer; font-size:0.85em; color:#555; user-select:none;",
              "Display options"
            ),
            shiny::div(
              style = "padding:8px 0 4px;",
              shiny::numericInput(ns("pw_fontsize"), "Font size (pt):",
                                  value = 10, min = 6, max = 20, step = 1)
            )
          ),
          shiny::actionButton(ns("run_pairwise"), "Run Pairwise Overlap",
                              class = "btn-primary btn-sm", width = "100%")
        ),
        shiny::uiOutput(ns("pairwise_status_ui")),
        shiny::uiOutput(ns("pairwise_plot_ui")),
        shiny::uiOutput(ns("pairwise_table_hdr_ui")),
        DT::DTOutput(ns("pairwise_table"))
      ),

      # ---- 2. Similarity vs MSigDB ------------------------------------------
      bslib::nav_panel(
        .tab_title("Similarity vs MSigDB",
          "Compares your loaded gene sets against MSigDB collections or custom reference sets to find known pathways with high overlap. Even a strong discriminatory gene set may reflect a well-characterised biological pathway rather than a novel finding."),
        .tab_settings(
          # ── Metric (always visible, single compact row) ───────────────────────
          shiny::div(
            style = "margin-bottom:10px;",
            shiny::tags$label("Metric", style = "font-size:0.85em; font-weight:600;"),
            shiny::radioButtons(ns("gs_metric"), label = NULL,
              choices  = c("Jaccard index" = "jaccard", "Odds ratio" = "odds_ratio"),
              selected = "jaccard", inline = TRUE)
          ),
          # ── Organism + collection (always visible, two columns, aligned top) ──
          shiny::div(
            style = "display:grid; grid-template-columns:1fr 1fr; gap:12px; align-items:start; margin-bottom:8px;",
            shiny::div(
              shiny::tags$label("Organism", style = "font-size:0.85em; font-weight:600;"),
              shinyWidgets::pickerInput(ns("gs_organism"), label = NULL,
                choices  = .msigdb_organisms,
                selected = "Homo sapiens",
                options  = shinyWidgets::pickerOptions(liveSearch = TRUE)),
              shiny::uiOutput(ns("gs_db_species_ui")),
              shiny::uiOutput(ns("gs_organism_info_ui"))
            ),
            shiny::div(
              shiny::uiOutput(ns("gs_collection_ui"))
            )
          ),

          # ── Advanced options (collapsible) ────────────────────────────────────
          shiny::tags$details(
            style = "margin-bottom:6px;",
            shiny::tags$summary(
              style = "cursor:pointer; font-size:0.85em; color:#555; user-select:none;",
              "Advanced options"
            ),
            shiny::div(
              style = "padding:10px 0 4px;",
              shiny::textInput(ns("gs_subcollection"), "Sub-collection:",
                placeholder = "e.g. CP:REACTOME"),
              shiny::textAreaInput(ns("msig_subset"),
                label = shiny::tags$span("Restrict to specific pathways ",
                  bslib::tooltip(shiny::icon("circle-question", style="color:#aaa;cursor:help;font-size:0.85em;"),
                    "Enter exact MSigDB pathway names, comma-separated.", placement = "right")),
                placeholder = "HALLMARK_APOPTOSIS, HALLMARK_HYPOXIA, ...", rows = 2, width = "100%"),
              shiny::conditionalPanel(paste0("input['", ns("gs_metric"), "'] == 'jaccard'"),
                shiny::sliderInput(ns("jaccard_threshold"), "Min. Jaccard index:",
                                   min = 0, max = 1, value = 0, step = 0.01),
                shiny::uiOutput(ns("jaccard_size_warn_ui"))),
              shiny::conditionalPanel(paste0("input['", ns("gs_metric"), "'] == 'odds_ratio'"),
                shiny::div(
                  style = "display:grid; grid-template-columns:1fr 1fr; gap:10px;",
                  shiny::numericInput(ns("or_threshold"), "Min. odds ratio:", value = 2, min = 0, step = 0.5),
                  shiny::numericInput(ns("or_pval"), "Max adj. p-value:", value = 0.05, min = 0, max = 1, step = 0.01)
                ))
            )
          ),

          # ── Custom gene sets (collapsible) ────────────────────────────────────
          shiny::tags$details(
            style = "margin-bottom:6px;",
            shiny::tags$summary(
              style = "cursor:pointer; font-size:0.85em; color:#555; user-select:none;",
              "Custom reference gene sets"
            ),
            shiny::div(
              style = "padding:10px 0 4px;",
              shiny::tags$p(style = "font-size:0.8em; color:#888; margin-bottom:6px;",
                "Compare your signatures against your own gene sets in addition to (or instead of) MSigDB.",
                " Format: ", shiny::tags$code(">SetName"), " on one line, then gene symbols comma-separated on the next."),
              shiny::textAreaInput(ns("user_gs_paste"),
                label = NULL,
                placeholder = ">MySignature\nGENE1,GENE2,GENE3\n>AnotherSet\nGENEA,GENEB",
                rows = 4, width = "100%"),
              shiny::fileInput(ns("user_gs_file"), label = NULL,
                accept = ".txt",
                buttonLabel = shiny::tags$span(shiny::icon("upload"), " Upload .txt"),
                placeholder = "or upload a .txt file with the same format", width = "100%"),
              shiny::uiOutput(ns("user_gs_parse_ui"))
            )
          ),

          # ── Plot options (collapsible) ────────────────────────────────────────
          shiny::tags$details(
            style = "margin-bottom:10px;",
            shiny::tags$summary(
              style = "cursor:pointer; font-size:0.85em; color:#555; user-select:none;",
              "Plot options"
            ),
            shiny::div(
              style = "padding:10px 0 4px; display:grid; grid-template-columns:1fr 1fr; gap:10px;",
              shiny::numericInput(ns("sim_fontsize"), "Font size (pt):",
                                  value = 9, min = 6, max = 20, step = 1),
              shiny::div(
                shiny::tags$label(style = "font-size:0.85em;", "View mode"),
                shiny::checkboxInput(ns("sim_interactive"), "Interactive (zoom/pan)", value = FALSE),
                shiny::tags$small(style = "color:#aaa; font-size:0.78em;",
                  "Interactive is slower but lets you zoom into pathway clusters.")
              )
            )
          ),

          # ── Run button ────────────────────────────────────────────────────────
          shiny::actionButton(ns("run_similarity"), "Run Similarity",
                              class = "btn-primary btn-sm")
        ),
        shiny::uiOutput(ns("similarity_status_ui")),
        shiny::uiOutput(ns("similarity_plot_ui")),
        shiny::uiOutput(ns("similarity_table_hdr_ui")),
        shiny::uiOutput(ns("similarity_table_container"))
      ),

      # ---- 3. Correlation ---------------------------------------------------
      bslib::nav_panel(
        .tab_title("Correlation",
          "Displays pairwise correlations between selected genes, revealing co-expression patterns within the gene set. Highly correlated genes may be co-regulated or redundant, while uncorrelated genes contribute independently to the signature score."),
        .tab_settings(
          shiny::uiOutput(ns("corr_overlap_warn_ui")),
          # ── Primary settings (always visible) ────────────────────────────────
          shiny::div(
            style = "display:grid; grid-template-columns:1fr 1fr; gap:12px; align-items:start; margin-bottom:8px;",
            shiny::div(
              shiny::uiOutput(ns("corr_genes_ui")),
              shiny::uiOutput(ns("corr_gene_count_warn_ui"))
            ),
            shiny::div(
              shiny::selectInput(ns("corr_method"), "Correlation method:",
                choices = c("Spearman" = "spearman", "Pearson" = "pearson", "Kendall" = "kendall"),
                selected = "spearman")
            )
          ),
          # ── Advanced options (collapsible) ────────────────────────────────────
          shiny::tags$details(
            style = "margin-bottom:10px;",
            shiny::tags$summary(
              style = "cursor:pointer; font-size:0.85em; color:#555; user-select:none;",
              "Advanced options"
            ),
            shiny::div(
              style = "padding:10px 0 4px;",
              shiny::uiOutput(ns("corr_sep_ui")),
              shiny::div(
                style = "display:grid; grid-template-columns:1fr 1fr; gap:10px; align-items:start; margin-top:6px;",
                shiny::div(
                  shiny::checkboxInput(ns("corr_cluster_rows"), "Cluster genes (rows)",    value = TRUE),
                  shiny::checkboxInput(ns("corr_cluster_cols"), "Cluster genes (columns)", value = TRUE),
                  shiny::checkboxInput(ns("corr_show_row_names"), "Show row labels",       value = TRUE),
                  shiny::checkboxInput(ns("show_col_names"),      "Show column labels",    value = TRUE)
                ),
                shiny::div(
                  shiny::numericInput(ns("corr_fontsize"), "Font size (pt):",
                                      value = 10, min = 6, max = 20, step = 1),
                  shiny::checkboxInput(ns("corr_interactive"), "Interactive (zoom/pan)", value = FALSE),
                  shiny::tags$small(style = "color:#aaa; font-size:0.78em;",
                    "Interactive is slower but lets you zoom into co-expression clusters and inspect gene pairs.")
                )
              )
            )
          ),
          shiny::actionButton(ns("run_corr"), "Plot Correlation Heatmap",
                              class = "btn-primary btn-sm")
        ),
        shiny::uiOutput(ns("corr_status_ui")),
        shiny::uiOutput(ns("corr_plot_ui"))
      ),

      # ---- 5. Effect Size ---------------------------------------------------
      bslib::nav_panel(
        .tab_title("Effect Size",
          "Computes Cohen's d for each gene: the standardised difference in expression between the selected class and all others (one-vs-rest). Useful for identifying which individual genes drive the overall gene set signal and which contribute little."),
        .tab_settings(
          shiny::uiOutput(ns("cohend_overlap_warn_ui")),
          # ── Gene selection (always visible) ───────────────────────────────────
          shiny::tags$label("Gene selection", style = "font-size:0.85em; font-weight:600;"),
          shiny::radioButtons(ns("cohend_mode"), label = NULL,
            choices  = c("Specific genes" = "manual", "Top / random N" = "n_based"),
            selected = "n_based", inline = TRUE),
          shiny::conditionalPanel(
            condition = sprintf("input['%s'] === 'manual'", ns("cohend_mode")),
            shiny::uiOutput(ns("cohend_genes_ui"))
          ),
          shiny::conditionalPanel(
            condition = sprintf("input['%s'] === 'n_based'", ns("cohend_mode")),
            shiny::div(
              style = "margin-bottom:8px;",
              shiny::numericInput(ns("cohend_max_genes"), "Number of genes (N):",
                value = 30, min = 1, max = 500, step = 5, width = "50%"),
              shiny::radioButtons(ns("cohend_gene_sel"), "Criterion:",
                choices  = c("Top N by Cohen's d" = "top_d", "First N (gene set order)" = "first", "Random N" = "random"),
                selected = "top_d", inline = TRUE)
            )
          ),
          shiny::tags$hr(style = "margin:6px 0 10px;"),
          # ── Group comparison (always visible) ─────────────────────────────────
          shiny::tags$label(
            "Group comparison ",
            bslib::tooltip(
              shiny::icon("circle-question", style = "color:#aaa; cursor:help; font-size:0.85em;"),
              shiny::tags$span(
                shiny::tags$b("Cohen's d"), " measures the standardised difference in gene expression",
                " between the selected class and all other samples (one-vs-rest), using pooled standard deviation.",
                " It is always computed as a binary comparison regardless of how many levels the variable has.",
                " d = 0: no separation; d ≈ 0.5: moderate; d ≥ 0.8: large."
              ), placement = "right"
            ),
            style = "font-size:0.85em; font-weight:600;"
          ),
          shiny::div(
            style = "display:grid; grid-template-columns:1fr 1fr; gap:12px; align-items:start; margin-bottom:8px;",
            shiny::uiOutput(ns("cohend_cond_var_ui")),
            shiny::uiOutput(ns("cohend_cond_class_ui"))
          ),
          # ── Display options (collapsible) ─────────────────────────────────────
          shiny::tags$details(
            style = "margin-bottom:10px;",
            shiny::tags$summary(
              style = "cursor:pointer; font-size:0.85em; color:#555; user-select:none;",
              "Display options"
            ),
            shiny::div(
              style = "padding:8px 0 4px;",
              shiny::numericInput(ns("cohend_fontsize"), "Font size (pt):",
                                  value = 10, min = 6, max = 20, step = 1)
            )
          ),
          shiny::actionButton(ns("run_cohend"), "Plot Effect Size",
                              class = "btn-primary btn-sm")
        ),
        shiny::uiOutput(ns("cohend_status_ui")),
        shiny::uiOutput(ns("cohend_plot_ui"))
      ),

      # ---- 6. Gene PCA ------------------------------------------------------
      bslib::nav_panel(
        .tab_title("Gene PCA",
          "Performs PCA on the expression of the gene set genes to test whether they explain enough variance to visually separate groups of interest. A complementary, unsupervised view alongside quantitative metrics like AUC: if the genes drive the biology, samples should cluster by the colour variable."),
        .tab_settings(
          # ── Primary settings (always visible) ────────────────────────────────
          shiny::div(
            style = "display:grid; grid-template-columns:1fr 1fr; gap:12px; align-items:start; margin-bottom:8px;",
            shiny::div(
              shiny::uiOutput(ns("pca_genes_ui")),
              shiny::uiOutput(ns("pca_overlap_warn_ui")),
              shiny::tags$small(style = "color:#888; font-size:0.78em; display:block; margin-top:2px;",
                "All gene set genes with matching expression data are used.",
                " Higher coverage means the PCA better captures the biology of interest.")
            ),
            shiny::uiOutput(ns("pca_color_var_ui"))
          ),
          # ── Display options (collapsible) ─────────────────────────────────────
          shiny::tags$details(
            style = "margin-bottom:10px;",
            shiny::tags$summary(
              style = "cursor:pointer; font-size:0.85em; color:#555; user-select:none;",
              "Display options"
            ),
            shiny::div(
              style = "padding:8px 0 4px;",
              shiny::div(
                style = "display:grid; grid-template-columns:1fr 1fr; gap:8px;",
                shiny::selectInput(ns("pca_x_pc"), "X axis:", choices = paste0("PC", 1:10), selected = "PC1"),
                shiny::selectInput(ns("pca_y_pc"), "Y axis:", choices = paste0("PC", 1:10), selected = "PC2")
              ),
              shiny::numericInput(ns("pca_fontsize"), "Font size (pt):",
                                  value = 10, min = 6, max = 20, step = 1)
            )
          ),
          shiny::actionButton(ns("run_pca"), "Plot Gene PCA",
                              class = "btn-primary btn-sm")
        ),
        shiny::uiOutput(ns("pca_status_ui")),
        shiny::uiOutput(ns("pca_plots_ui"))
      ),

      # ---- 7. Discriminatory Power ------------------------------------------
      bslib::nav_panel(
        .tab_title("Discriminatory Power",
          shiny::tagList(
            "Evaluates each gene's ability to discriminate between two groups via ROC curves and AUC values.",
            shiny::tags$br(),
            "AUC = 0.5: no better than chance. AUC = 1.0: perfect separation.",
            " A raw AUC below 0.5 means the gene separates in the opposite direction and is still informative.",
            shiny::tags$br(),
            "Colours are consistent between the ROC curves and the AUC barplot.",
            shiny::tags$br(),
            shiny::tags$em("Tip: to explore per-sample expression for any gene, use the Expression tab.")
          )),
        .tab_settings(
          # ── Gene selection (always visible) ──────────────────────────────────
          shiny::uiOutput(ns("roc_overlap_warn_ui")),
          shiny::tags$small(style = "color:#888; font-size:0.78em; display:block; margin-bottom:6px;",
            "Genes shown are from the active gene set only."),
          shiny::tags$label("Gene selection", style = "font-size:0.85em; font-weight:600;"),
          shiny::radioButtons(ns("roc_mode"), label = NULL,
            choices  = c("Specific genes" = "manual",
                         "Top / random N" = "n_based"),
            selected = "n_based", inline = TRUE),
          shiny::conditionalPanel(
            condition = sprintf("input['%s'] === 'manual'", ns("roc_mode")),
            shiny::uiOutput(ns("roc_genes_ui"))
          ),
          shiny::conditionalPanel(
            condition = sprintf("input['%s'] === 'n_based'", ns("roc_mode")),
            shiny::numericInput(ns("roc_max_genes"), "Number of genes (N):",
              value = 20, min = 1, max = 500, step = 5, width = "50%"),
            shiny::radioButtons(ns("roc_gene_sel"), "Criterion:",
              choices  = c("Top N by AUC" = "top_auc", "Random N" = "random"),
              selected = "top_auc", inline = TRUE),
            shiny::tags$small(style = "color:#6c757d; display:block; margin-bottom:8px;",
              shiny::icon("circle-info"),
              " Top N by AUC ranks by ", shiny::tags$b("effective AUC"),
              " = max(AUC, 1 − AUC) — captures discriminatory genes in both directions.")
          ),
          # ── Separator ────────────────────────────────────────────────────────
          shiny::tags$hr(style = "margin:8px 0;"),
          # ── Group comparison (always visible) ─────────────────────────────────
          shiny::tags$p(style = "font-size:0.85em; font-weight:600; margin-bottom:4px;",
            "Group comparison"),
          shiny::div(
            style = "display:grid; grid-template-columns:1fr 1fr; gap:10px; align-items:start; margin-bottom:6px;",
            shiny::uiOutput(ns("roc_cond_var_ui")),
            shiny::uiOutput(ns("roc_cond_class_ui"))
          ),
          shiny::checkboxInput(ns("roc_invert_auc"),
            "Invert AUC < 0.5 genes (flip ROC direction)", value = TRUE),
          shiny::tags$small(style = "color:#6c757d; display:block; margin:-4px 0 8px;",
            "When checked: genes with raw AUC < 0.5 are flipped so the reported",
            " AUC = 1 − raw. An AUC = 0 becomes 1.0 (perfect in reverse direction)."),
          # ── Display options (collapsible) ─────────────────────────────────────
          shiny::tags$details(
            style = "margin-bottom:10px;",
            shiny::tags$summary(
              style = "cursor:pointer; font-size:0.85em; color:#555; user-select:none;",
              "Display options"
            ),
            shiny::div(
              style = "padding:8px 0 4px;",
              shiny::numericInput(ns("roc_fontsize"), "Font size (pt):",
                                  value = 10, min = 6, max = 20, step = 1)
            )
          ),
          shiny::actionButton(ns("run_rocauc"), "Plot Discriminatory Power",
                              class = "btn-primary btn-sm", width = "100%")
        ),
        shiny::uiOutput(ns("rocauc_status_ui")),
        shiny::uiOutput(ns("roc_curves_ui")),
        shiny::uiOutput(ns("rocauc_plot_ui"))
      ),

      # ---- 7. Expression (violins + heatmap) — last tab ---------------------
      bslib::nav_panel(
        .tab_title("Expression",
          "Displays expression levels of selected gene set genes across samples (violin plots and heatmap), annotated by a chosen metadata variable. Useful for inspecting individual gene behaviour and validating whether key genes show the expected pattern of change."),
        .tab_settings(
          shiny::uiOutput(ns("expr_overlap_warn_ui")),
          # ── Shared settings (always visible) ─────────────────────────────────
          shiny::div(
            style = "display:grid; grid-template-columns:1fr 1fr; gap:12px; align-items:start; margin-bottom:10px;",
            shiny::div(
              shiny::uiOutput(ns("expr_genes_ui")),
              shiny::uiOutput(ns("expr_gene_count_warn_ui"))
            ),
            shiny::div(
              shiny::uiOutput(ns("expr_group_var_ui")),
              shiny::uiOutput(ns("expr_level_warning_ui"))
            )
          ),
          # ── Violin plots (collapsible, open by default) ───────────────────────
          shiny::tags$details(
            open = NA,
            style = "margin-bottom:6px;",
            shiny::tags$summary(
              style = "cursor:pointer; font-size:0.85em; font-weight:600; color:#333; user-select:none; margin-bottom:6px;",
              "Violin plots"
            ),
            shiny::div(
              style = "padding:8px 0 4px; display:grid; grid-template-columns:1fr 1fr; gap:10px; align-items:end;",
              shiny::uiOutput(ns("expr_color_var_ui")),
              shiny::numericInput(ns("expr_fontsize"), "Font size (pt):",
                                  value = 10, min = 6, max = 20, step = 1)
            ),
            shiny::actionButton(ns("run_expression"), "Plot Violins",
                                class = "btn-primary btn-sm")
          ),
          shiny::tags$hr(style = "margin:4px 0 6px;"),
          # ── Expression heatmap (collapsible, open by default) ─────────────────
          shiny::tags$details(
            open = NA,
            style = "margin-bottom:4px;",
            shiny::tags$summary(
              style = "cursor:pointer; font-size:0.85em; font-weight:600; color:#333; user-select:none; margin-bottom:6px;",
              "Expression heatmap"
            ),
            shiny::div(
              style = "padding:8px 0 6px;",
              shiny::div(
                style = "display:grid; grid-template-columns:1fr 1fr; gap:10px; align-items:start; margin-bottom:6px;",
                shiny::div(
                  shiny::checkboxInput(ns("hm_cluster_rows"), "Cluster genes",   value = TRUE),
                  shiny::checkboxInput(ns("hm_cluster_cols"), "Cluster samples", value = TRUE)
                ),
                shiny::numericInput(ns("hm_fontsize"), "Font size (pt):",
                                    value = 10, min = 6, max = 20, step = 1)
              ),
              shiny::checkboxInput(ns("hm_interactive"), "Interactive view (zoom/pan)", value = FALSE),
              shiny::tags$small(style = "color:#aaa; font-size:0.78em; display:block; margin:-4px 0 8px;",
                "Interactive is slower but lets you zoom into gene clusters and hover over samples."),
              shiny::actionButton(ns("run_heatmap"), "Plot Heatmap",
                                  class = "btn-primary btn-sm")
            )
          )
        ),
        shiny::uiOutput(ns("violin_status_ui")),
        shiny::uiOutput(ns("expr_plot_ui")),
        shiny::uiOutput(ns("heatmap_status_ui")),
        shiny::uiOutput(ns("heatmap_plot_ui"))
      )
    ),
    js_ready_script
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
    user_gs_parsed    <- shiny::reactiveVal(NULL)   # user-defined ref gene sets for Similarity
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

    # ---- Overlap between active gene set and expression matrix ---------------
    # Uses RAW (non-normalised) gene names so that case differences between
    # organisms (mouse Title-case vs human ALL-CAPS) are flagged correctly.
    # toupper() would mask cross-species mismatches by making them look identical.
    gs_expr_overlap <- shiny::reactive({
      gs   <- get_gene_sets()
      sel  <- input$selected_gs
      expr <- get_expr()
      if (is.null(gs) || is.null(sel) || !sel %in% names(gs) || is.null(expr)) return(NULL)
      gs_genes   <- .gs_genes(gs[[sel]])              # raw — preserve case
      expr_genes <- rownames(expr)                     # raw — preserve case
      n_gs      <- length(gs_genes)
      n_overlap <- length(intersect(gs_genes, expr_genes))
      # Also try case-insensitive to distinguish "wrong case" from "truly absent"
      n_overlap_ci <- length(intersect(toupper(gs_genes), toupper(expr_genes)))
      list(n_overlap    = n_overlap,
           n_overlap_ci = n_overlap_ci,
           n_gs         = n_gs,
           pct    = if (n_gs > 0) round(100 * n_overlap    / n_gs) else 0L,
           pct_ci = if (n_gs > 0) round(100 * n_overlap_ci / n_gs) else 0L)
    })

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
      # Auto-run PCA immediately if data are ready
      gs <- get_gene_sets(); sel <- input$selected_gs
      if (!is.null(gs) && !is.null(sel) && sel %in% names(gs) &&
          !is.null(get_expr()) && !is.null(get_meta())) {
        genes_plot <- .gs_genes(gs[[sel]])
        cv <- shiny::isolate(input$pca_color_var)
        color_var <- if (!is.null(cv) && nzchar(trimws(cv %||% ""))) cv else NULL
        .run_pca_now(genes_plot, color_var)
      }
    }, ignoreInit = TRUE, ignoreNULL = TRUE)

    # ---- Tab visibility: show/hide based on available data ------------------
    # Uses a custom JS message handler (markeR.setTabVisible) instead of
    # bslib::nav_hide/nav_show to avoid bslib JS initialisation timing issues.
    .set_tab_visible <- function(tab_value, visible) {
      session$sendCustomMessage("markeR.setTabVisible", list(
        navsetId = ns("gs_navset"),
        tabValue = tab_value,
        visible  = visible
      ))
    }

    .apply_tab_visibility <- function() {
      gs           <- shiny::isolate(get_gene_sets())
      has_multi_gs <- !is.null(gs) && length(gs) > 1
      has_expr     <- !is.null(shiny::isolate(get_expr()))
      has_meta     <- !is.null(shiny::isolate(get_meta()))
      has_both     <- has_expr && has_meta

      .set_tab_visible("Pairwise Overlap", has_multi_gs)
      for (tab in c("Expression", "Correlation", "Effect Size", "Gene PCA", "Discriminatory Power")) {
        .set_tab_visible(tab, has_both)
      }
    }

    # Fire once after shiny:connected (JS handler is ready, DOM is rendered).
    shiny::observeEvent(input$js_initialized, {
      .apply_tab_visibility()
    }, once = TRUE, ignoreNULL = TRUE)

    # Re-apply whenever data availability changes.
    shiny::observe({
      .apply_tab_visibility()
    }) |> shiny::bindEvent(
      get_gene_sets(), get_expr(), get_meta(),
      ignoreInit = TRUE,
      ignoreNULL = FALSE
    )

    # ---- Reset violin/heatmap when grouping variable changes ----------------
    shiny::observeEvent(input$expr_group_var, {
      expression_result(NULL); heatmap_result(NULL)
    }, ignoreInit = TRUE, ignoreNULL = TRUE)
    shiny::observeEvent(input$pca_color_var, {
      pca_result(NULL)
      gs <- get_gene_sets(); sel <- shiny::isolate(input$selected_gs)
      if (!is.null(gs) && !is.null(sel) && sel %in% names(gs) &&
          !is.null(get_expr()) && !is.null(get_meta())) {
        genes_plot <- .gs_genes(gs[[sel]])
        cv <- input$pca_color_var
        color_var <- if (!is.null(cv) && nzchar(trimws(cv %||% ""))) cv else NULL
        .run_pca_now(genes_plot, color_var)
      }
    }, ignoreInit = TRUE, ignoreNULL = TRUE)
    shiny::observeEvent(input$expr_genes, {
      expression_result(NULL); heatmap_result(NULL)
    }, ignoreInit = TRUE, ignoreNULL = TRUE)
    shiny::observeEvent(input$corr_genes, {
      corr_result(NULL)
    }, ignoreInit = TRUE, ignoreNULL = TRUE)
    shiny::observeEvent(input$cohend_genes, {
      cohend_result(NULL)
    }, ignoreInit = TRUE, ignoreNULL = TRUE)
    shiny::observeEvent(input$pca_genes, {
      pca_result(NULL)
    }, ignoreInit = TRUE, ignoreNULL = TRUE)
    shiny::observeEvent(input$roc_genes, {
      rocauc_result(NULL)
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

    output$gs_composition_ui <- shiny::renderUI({
      gs  <- get_gene_sets()
      sel <- input$selected_gs
      if (is.null(gs) || length(gs) == 0 || is.null(sel) || !sel %in% names(gs)) return(NULL)
      entry <- gs[[sel]]
      if (is.data.frame(entry)) {
        df <- entry[, seq_len(min(2L, ncol(entry))), drop = FALSE]
        colnames(df)[1] <- "Gene"
        if (ncol(df) == 2) colnames(df)[2] <- "Direction"
      } else {
        df <- data.frame(Gene = as.character(entry), stringsAsFactors = FALSE)
      }
      shiny::tagList(
        shiny::tags$small(shiny::tags$b(
          sprintf("Composition (%d gene%s)", nrow(df), if (nrow(df) != 1) "s" else ""))),
        DT::datatable(df, rownames = FALSE, filter = "none",
          options = list(paging = FALSE, dom = "ft", scrollY = "260px",
                         scrollCollapse = TRUE, scrollX = FALSE,
                         columnDefs = list(list(className = "dt-left", targets = "_all"))))
      )
    })

    # Build the overlap alert div given the overlap list (used in both sidebar
    # and per-tab versions).  Returns NULL when overlap is fine (≥60% exact match).
    .overlap_alert <- function(ov, compact = FALSE) {
      if (is.null(ov)) return(NULL)
      style_base <- if (compact) "font-size:0.82em; padding:5px 9px; margin-bottom:6px;"
                    else         "font-size:0.82em; padding:6px 9px; margin-top:6px;"

      if (ov$n_overlap == 0 && ov$n_overlap_ci == 0) {
        # Truly no matching genes at all
        shiny::div(class = "alert alert-danger", style = style_base,
          shiny::icon("circle-xmark"), " ",
          shiny::strong("No gene overlap"), " with the expression matrix.",
          " Check that the gene ID type matches (e.g. symbol vs Ensembl)",
          " and that the correct organism is selected.")

      } else if (ov$n_overlap == 0 && ov$n_overlap_ci > 0) {
        # Case mismatch → organism mismatch (mouse Title-case vs human ALL-CAPS)
        shiny::div(class = "alert alert-danger", style = style_base,
          shiny::icon("circle-xmark"), " ",
          shiny::strong("Gene symbol case mismatch — likely organism mismatch."),
          sprintf(" 0 / %d gene(s) match exactly, but %d match after ignoring case.",
                  ov$n_gs, ov$n_overlap_ci),
          " Your gene set and expression matrix appear to use",
          " different species conventions (e.g. mouse Title-case vs human ALL-CAPS).",
          " Per-gene analyses will return empty results unless you fix the organism.")

      } else if (ov$pct < 25) {
        pct_label <- if (ov$pct == 0 && ov$n_overlap > 0) "<1%" else paste0(ov$pct, "%")
        shiny::div(class = "alert alert-warning", style = style_base,
          shiny::icon("triangle-exclamation"), " ",
          shiny::strong(sprintf("%d / %d genes (%s)", ov$n_overlap, ov$n_gs, pct_label)),
          " found in the expression matrix.",
          " Check organism and gene ID type — results may be unreliable.")

      } else if (ov$pct < 60) {
        shiny::div(class = "alert alert-info", style = style_base,
          shiny::icon("circle-info"), " ",
          shiny::strong(sprintf("%d / %d genes (%d%%)", ov$n_overlap, ov$n_gs, ov$pct)),
          " found in the expression matrix.")
      }
      # ≥60% exact overlap → no message
    }

    # Compact per-tab version
    .overlap_warn_inline <- function() {
      shiny::renderUI(.overlap_alert(gs_expr_overlap(), compact = TRUE))
    }
    output$expr_overlap_warn_ui   <- .overlap_warn_inline()
    output$corr_overlap_warn_ui   <- .overlap_warn_inline()
    output$cohend_overlap_warn_ui <- .overlap_warn_inline()
    output$pca_overlap_warn_ui    <- .overlap_warn_inline()
    output$roc_overlap_warn_ui    <- .overlap_warn_inline()

    # Sidebar version (slightly more padding, shown above composition table)
    output$gs_overlap_warn_ui <- shiny::renderUI(
      .overlap_alert(gs_expr_overlap(), compact = FALSE)
    )

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

    # In-tab gene-count warning (Expression tab) — shown when > 20 genes selected
    output$expr_gene_count_warn_ui <- shiny::renderUI({
      n <- length(input$expr_genes)
      if (n > 20)
        shiny::div(
          class = "alert alert-warning",
          style = "font-size:0.82em; padding:6px 10px; margin-bottom:8px;",
          shiny::icon("triangle-exclamation"), " ",
          shiny::strong(n), " genes currently selected.",
          " Violin plots and heatmaps with more than ~20 genes can be very slow",
          " or hard to read. Reduce the selection using the gene picker above.")
    })

    # ---- Per-tab: genes-to-plot pickers ------------------------------------
    .make_gene_picker <- function(input_id, label = "Genes to plot:") {
      shiny::renderUI({
        gs   <- get_gene_sets()
        sel  <- input$selected_gs
        expr <- get_expr()
        if (is.null(gs) || is.null(sel) || !sel %in% names(gs)) return(NULL)
        genes <- .gs_genes(gs[[sel]])

        # Split into "in data" / "not in data" if expression matrix is available
        if (!is.null(expr)) {
          genes_in  <- genes[genes %in% rownames(expr)]
          genes_out <- genes[!genes %in% rownames(expr)]
          n_in      <- length(genes_in)
          n_total   <- length(genes)
          count_note <- shiny::tags$small(
            style = "color:#6c757d; display:block; margin-top:2px;",
            sprintf("%d of %d gene set genes found in expression data.", n_in, n_total)
          )
          choices   <- if (length(genes_out) > 0)
            list("In expression data" = genes_in, "Not in expression data" = genes_out)
          else
            genes_in
          selected  <- genes_in   # pre-select only genes present in the data
        } else {
          choices   <- genes
          selected  <- genes
          count_note <- NULL
        }

        shiny::tagList(
          shinyWidgets::pickerInput(ns(input_id), label,
            choices  = choices,
            selected = selected,
            multiple = TRUE,
            options  = shinyWidgets::pickerOptions(
              actionsBox = TRUE, liveSearch = TRUE,
              selectedTextFormat = "count > 4",
              countSelectedText  = "{0} genes selected"
            )
          ),
          count_note
        )
      })
    }
    output$expr_genes_ui   <- .make_gene_picker("expr_genes")
    output$corr_genes_ui   <- .make_gene_picker("corr_genes")
    output$cohend_genes_ui <- .make_gene_picker("cohend_genes")
    output$pca_genes_ui    <- .make_gene_picker("pca_genes")
    output$roc_genes_ui    <- .make_gene_picker("roc_genes",  "Genes to analyse:")

    # ---- Per-tab: group variable (Expression only) -------------------------
    output$expr_group_var_ui <- shiny::renderUI({
      meta <- get_meta()
      if (is.null(meta)) return(NULL)
      cols <- setdiff(colnames(meta), colnames(meta)[1])
      shiny::selectInput(ns("expr_group_var"), "Grouping variable (x-axis / violin):",
                         choices = cols, selected = cols[1])
    })
    output$expr_level_warning_ui <- shiny::renderUI({
      meta <- get_meta()
      if (is.null(meta) || is.null(input$expr_group_var) ||
          !input$expr_group_var %in% colnames(meta)) return(NULL)
      n_lev <- length(unique(as.character(meta[[input$expr_group_var]])))
      if (n_lev > 10)
        shiny::div(class = "alert alert-warning",
                   style = "font-size:0.82em; padding:5px 9px; margin-bottom:4px;",
                   shiny::icon("triangle-exclamation"), " ", shiny::strong(n_lev),
                   " levels — plots may be hard to read.")
    })

    # ---- Per-tab: colour variable -------------------------------------------
    .make_color_picker <- function(input_id, label = "Colour by (optional):") {
      shiny::renderUI({
        meta <- get_meta()
        if (is.null(meta)) return(NULL)
        cols <- setdiff(colnames(meta), colnames(meta)[1])
        shiny::selectInput(ns(input_id), label,
                           choices = c("None" = "", cols), selected = "")
      })
    }
    output$expr_color_var_ui <- .make_color_picker("expr_color_var")
    output$pca_color_var_ui  <- .make_color_picker("pca_color_var", "Colour by:")

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
          NULL
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
      # Only show variables with ≥2 distinct non-NA values (a variable with 1 level
      # cannot define a contrast and causes "missing value where TRUE/FALSE needed").
      valid <- cols[sapply(cols, function(v) {
        length(unique(na.omit(as.character(meta[[v]])))) >= 2L
      })]
      if (length(valid) == 0)
        return(shiny::div(class = "alert alert-warning", style = "font-size:0.82em;",
          "No metadata variable with ≥2 unique values found."))
      shiny::selectInput(ns("cohend_cond_var"), "Grouping variable:",
                         choices = valid, selected = valid[1])
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
      valid <- cols[sapply(cols, function(v) {
        length(unique(na.omit(as.character(meta[[v]])))) >= 2L
      })]
      if (length(valid) == 0)
        return(shiny::div(class = "alert alert-warning", style = "font-size:0.82em;",
          "No metadata variable with ≥2 unique values found."))
      shiny::selectInput(ns("roc_cond_var"), "Grouping variable:",
                         choices = valid, selected = valid[1])
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

    # ---- Parse user-defined gene sets (paste or file upload) ----------------
    .parse_user_gs_text <- function(txt) {
      if (is.null(txt) || !nzchar(trimws(txt))) return(NULL)
      lines  <- strsplit(txt, "\n")[[1]]
      result <- list(); cur_name <- NULL; cur_genes <- NULL
      for (ln in lines) {
        ln <- trimws(ln)
        if (!nzchar(ln)) next
        if (startsWith(ln, ">")) {
          if (!is.null(cur_name) && length(cur_genes) > 0)
            result[[cur_name]] <- cur_genes
          cur_name  <- trimws(substring(ln, 2))
          cur_genes <- character(0)
        } else if (!is.null(cur_name)) {
          g <- trimws(unlist(strsplit(ln, "[,;\\s]+")))
          cur_genes <- c(cur_genes, g[nzchar(g)])
        }
      }
      if (!is.null(cur_name) && length(cur_genes) > 0)
        result[[cur_name]] <- cur_genes
      if (length(result) == 0) NULL else result
    }

    shiny::observeEvent(input$user_gs_paste, {
      parsed <- .parse_user_gs_text(input$user_gs_paste)
      user_gs_parsed(parsed)
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$user_gs_file, {
      f <- input$user_gs_file
      if (is.null(f)) return()
      txt <- paste(readLines(f$datapath, warn = FALSE), collapse = "\n")
      shiny::updateTextAreaInput(session, "user_gs_paste", value = txt)
      parsed <- .parse_user_gs_text(txt)
      user_gs_parsed(parsed)
    }, ignoreInit = TRUE)

    output$user_gs_parse_ui <- shiny::renderUI({
      parsed <- user_gs_parsed()
      if (is.null(parsed)) return(NULL)
      shiny::tags$small(style = "color:#28a745; display:block; margin:-4px 0 6px;",
        shiny::icon("circle-check"),
        sprintf(" %d gene set(s) loaded: %s", length(parsed),
                paste(names(parsed), collapse = ", ")))
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
          user_gs  <- shiny::isolate(user_gs_parsed())
          args <- list(signatures = gs_char, collection = collection,
                       subcollection = subcollection, msig_subset = msig_subset, metric = metric,
                       organism = organism, db_species = db_species)
          if (!is.null(user_gs) && length(user_gs) > 0)
            args$other_user_signatures <- user_gs
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
      shiny::req(expr_d, meta, gs, input$selected_gs, input$expr_group_var,
                 input$expr_genes, length(input$expr_genes) > 0)
      genes_plot <- shiny::isolate(input$expr_genes)
      group_var  <- shiny::isolate(input$expr_group_var)
      color_var  <- shiny::isolate(input$expr_color_var)
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
      shiny::req(expr_d, meta, gs, input$selected_gs, input$expr_group_var,
                 input$expr_genes, length(input$expr_genes) > 0)
      genes_plot <- shiny::isolate(input$expr_genes)
      group_var  <- shiny::isolate(input$expr_group_var)
      cl_rows    <- shiny::isolate(input$hm_cluster_rows)
      cl_cols    <- shiny::isolate(input$hm_cluster_cols)
      if (is.null(cl_rows)) cl_rows <- TRUE
      if (is.null(cl_cols)) cl_cols <- TRUE

      levs      <- unique(as.character(meta[[group_var]]))
      annot_col <- list(stats::setNames(.pp_palette[seq_along(levs)], levs))
      names(annot_col) <- group_var

      # Remove genes with zero variance (produce NaN after Z-scoring → hclust fails)
      genes_in_expr <- genes_plot[genes_plot %in% rownames(expr_d)]
      if (length(genes_in_expr) > 0) {
        gvars <- apply(expr_d[genes_in_expr, , drop = FALSE], 1, stats::var, na.rm = TRUE)
        genes_plot <- genes_in_expr[is.finite(gvars) & gvars > 0]
      }
      if (length(genes_plot) == 0) {
        shiny::showNotification("No genes with non-zero variance found in expression data.",
                                type = "error", duration = 8)
        return()
      }

      shiny::withProgress(message = "Generating expression heatmap...", value = 0.3, {
        result <- tryCatch({
          font_sz <- shiny::isolate(input$hm_fontsize)
          if (is.null(font_sz) || !is.numeric(font_sz)) font_sz <- 10
          hm <- ExpressionHeatmap(
            data              = expr_d,
            metadata          = meta,
            genes             = genes_plot,
            annotate.by       = group_var,
            annotation_colors = annot_col,
            colorlist         = list(low = "#4173B4", mid = "white", high = "#B44141"),
            cluster_rows      = cl_rows,
            cluster_columns   = cl_cols,
            title             = character(0),
            titlesize         = font_sz + 2,
            fontsize          = font_sz
          )
          shiny::incProgress(0.7)
          list(heatmap = hm, gs_name = shiny::isolate(input$selected_gs),
               n_genes = length(genes_plot), cl_rows = cl_rows, cl_cols = cl_cols,
               group_var     = group_var,
               annot_colors  = annot_col[[group_var]],
               sample_groups = stats::setNames(
                 as.character(meta[[group_var]]),
                 as.character(meta[[1L]])))
        }, error = function(e) {
          shiny::showNotification(paste("Expression heatmap failed:", conditionMessage(e)),
                                  type = "error", duration = 10); NULL
        })
      })
      if (!is.null(result)) {
        base_h  <- .dyn_h(result$n_genes, base = 300, per_item = 30)
        extra_h <- if (isTRUE(result$cl_cols)) 120L else 0L   # room for column dendrogram
        hm_h(base_h + extra_h)
      }
      heatmap_result(result)
    })

    shiny::observeEvent(input$run_corr, {
      expr_d <- get_expr(); meta <- get_meta()
      shiny::req(expr_d, meta, input$corr_genes, length(input$corr_genes) > 0)
      genes_plot   <- shiny::isolate(input$corr_genes)
      # Filter to genes present in the expression matrix with non-zero variance
      # (zero-variance genes produce NA correlations which crash hclust)
      genes_plot <- genes_plot[genes_plot %in% rownames(expr_d)]
      if (length(genes_plot) < 2) {
        shiny::showNotification("Need ≥2 genes with expression in the matrix for correlation.",
                                type = "error", duration = 8)
        return()
      }
      gvars <- apply(expr_d[genes_plot, , drop = FALSE], 1, stats::var, na.rm = TRUE)
      genes_plot <- genes_plot[is.finite(gvars) & gvars > 0]
      if (length(genes_plot) < 2) {
        shiny::showNotification("Need ≥2 genes with non-zero variance for correlation.",
                                type = "error", duration = 8)
        return()
      }
      sep          <- shiny::isolate(input$corr_sep)
      sep_by       <- if (!is.null(sep) && nzchar(sep)) sep else NULL
      method       <- shiny::isolate(input$corr_method)
      show_colnam  <- shiny::isolate(input$show_col_names)
      if (is.null(show_colnam)) show_colnam <- TRUE
      cl_rows <- shiny::isolate(input$corr_cluster_rows); if (is.null(cl_rows)) cl_rows <- TRUE
      cl_cols <- shiny::isolate(input$corr_cluster_cols); if (is.null(cl_cols)) cl_cols <- TRUE
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
               sep = sep_by, method = method, show_col_names = show_colnam,
               cluster_rows = cl_rows, cluster_columns = cl_cols)
        }, error = function(e) {
          shiny::showNotification(paste("Correlation failed:", conditionMessage(e)),
                                  type = "error", duration = 10); NULL
        })
      })
      if (!is.null(result)) {
        base_h  <- .dyn_h(result$n_genes, base = 300, per_item = 30)
        extra_h <- if (isTRUE(result$cluster_columns)) 120L else 0L
        corr_h(base_h + extra_h)
      }
      corr_result(result)
    })

    shiny::observeEvent(input$run_cohend, {
      expr_d <- get_expr(); meta <- get_meta()
      shiny::req(expr_d, meta, input$cohend_cond_var, input$cohend_cond_class)

      mode       <- shiny::isolate(input$cohend_mode) %||% "n_based"
      cond_var   <- shiny::isolate(input$cohend_cond_var)
      cond_class <- shiny::isolate(input$cohend_cond_class)

      # Gene universe depends on mode
      if (isTRUE(mode == "manual")) {
        shiny::req(input$cohend_genes, length(input$cohend_genes) > 0)
        genes_all <- shiny::isolate(input$cohend_genes)
        gene_sel  <- "manual"
        max_g     <- NULL
      } else {
        gs  <- get_gene_sets(); sel <- input$selected_gs
        shiny::req(gs, sel, sel %in% names(gs))
        genes_all <- .gs_genes(gs[[sel]])
        gene_sel  <- shiny::isolate(input$cohend_gene_sel) %||% "top_d"
        max_g     <- shiny::isolate(input$cohend_max_genes)
      }

      # Drop samples with NA in the condition variable
      valid_rows <- !is.na(meta[[cond_var]])
      if (sum(valid_rows) == 0) {
        shiny::showNotification("All samples have NA for the selected condition variable.",
                                type = "error", duration = 8)
        return()
      }
      meta_sub <- meta[valid_rows, , drop = FALSE]
      # Align expression columns to the filtered metadata samples
      shared_s <- intersect(colnames(expr_d), as.character(meta_sub[[1]]))
      meta_sub  <- meta_sub[as.character(meta_sub[[1]]) %in% shared_s, , drop = FALSE]
      expr_sub  <- expr_d[, shared_s, drop = FALSE]

      # Validate that both classes have ≥2 samples
      n_pos <- sum(as.character(meta_sub[[cond_var]]) %in% cond_class)
      n_neg <- nrow(meta_sub) - n_pos
      if (n_pos < 2 || n_neg < 2) {
        shiny::showNotification(
          sprintf("Need ≥2 samples in both groups. Got %d positive, %d negative.", n_pos, n_neg),
          type = "error", duration = 10)
        return()
      }

      # Keep only genes present in expression data and with non-zero variance
      genes_compute <- genes_all[genes_all %in% rownames(expr_sub)]
      if (length(genes_compute) == 0) {
        shiny::showNotification("None of the selected genes found in the expression matrix.",
                                type = "error", duration = 8)
        return()
      }
      gene_vars     <- apply(expr_sub[genes_compute, , drop = FALSE], 1, stats::var, na.rm = TRUE)
      genes_compute <- genes_compute[is.finite(gene_vars) & gene_vars > 0]
      if (length(genes_compute) == 0) {
        shiny::showNotification("All selected genes have zero variance across samples.",
                                type = "error", duration = 8)
        return()
      }

      # N-based sub-selection for non-top_d modes
      if (!isTRUE(gene_sel %in% c("top_d", "manual")) &&
          !is.null(max_g) && max_g > 0 && length(genes_compute) > max_g) {
        genes_compute <- if (isTRUE(gene_sel == "random"))
          sample(genes_compute, min(max_g, length(genes_compute)))
        else
          genes_compute[seq_len(min(max_g, length(genes_compute)))]
      }

      shiny::withProgress(message = "Computing effect size...", value = 0.3, {
        result <- tryCatch({
          res <- CohenD_IndividualGenes(
            data = expr_sub, metadata = meta_sub, genes = genes_compute,
            condition_var = cond_var, class = cond_class, group_var = NULL,
            params = list(colors = .pp_palette[1], limits = NULL))
          shiny::incProgress(0.6)

          # Post-filter for "top_d": keep only the N genes with largest |d|.
          d_df <- if (is.list(res) && is.data.frame(res$data)) res$data else NULL
          if (!is.null(d_df) && "CohensD" %in% colnames(d_df)) {
            d_valid <- d_df[!is.na(d_df$CohensD), , drop = FALSE]
            if (isTRUE(gene_sel == "top_d") && !is.null(max_g) && max_g > 0 &&
                nrow(d_valid) > max_g) {
              top_genes <- d_valid$Gene[order(abs(d_valid$CohensD), decreasing = TRUE)][seq_len(max_g)]
              res$data  <- d_valid[d_valid$Gene %in% top_genes, , drop = FALSE]
            } else {
              res$data <- d_valid
            }
          }

          shiny::incProgress(0.4)
          n_shown <- if (is.list(res) && is.data.frame(res$data)) nrow(res$data) else length(genes_compute)
          list(res = res, n_genes = n_shown, cond_var = cond_var, cond_class = cond_class,
               gene_sel = gene_sel)
        }, error = function(e) {
          shiny::showNotification(paste("Effect size failed:", conditionMessage(e)),
                                  type = "error", duration = 10); NULL
        })
      })
      if (!is.null(result)) cohend_h(.dyn_h(result$n_genes, base = 250, per_item = 28))
      cohend_result(result)
    })

    # Helper that runs PCA — shared by both auto-run and the button observer
    .run_pca_now <- function(genes_plot, color_var) {
      expr_d <- get_expr(); meta <- get_meta()
      if (is.null(expr_d) || is.null(meta) || length(genes_plot) == 0) return()
      shiny::withProgress(message = "Running PCA on genes...", value = 0.2, {
        result <- tryCatch({
          genes_use <- genes_plot[genes_plot %in% rownames(expr_d)]
          if (length(genes_use) < 2)
            stop("Need ≥2 genes present in the expression matrix for PCA.")
          pca_mat <- t(log2(as.matrix(expr_d[genes_use, , drop = FALSE]) + 1))
          pca_obj <- stats::prcomp(pca_mat, scale. = FALSE, center = TRUE)
          shiny::incProgress(0.5)
          var_exp <- round(100 * pca_obj$sdev^2 / sum(pca_obj$sdev^2), 2)
          n_pcs   <- min(10L, ncol(pca_obj$x))
          pca_df <- as.data.frame(pca_obj$x[, seq_len(n_pcs), drop = FALSE])
          pca_df$Sample <- rownames(pca_df)
          if (!is.null(meta)) {
            meta2 <- meta; colnames(meta2)[1] <- "Sample"
            pca_df <- merge(pca_df, meta2, by = "Sample", all.x = TRUE)
          }
          shiny::incProgress(0.3)
          list(pca_df = pca_df, var_exp = var_exp, n_pcs = n_pcs,
               color_var = color_var, n_genes = length(genes_use))
        }, error = function(e) {
          shiny::showNotification(paste("PCA failed:", conditionMessage(e)),
                                  type = "error", duration = 10); NULL
        })
      })
      pca_result(result)
      if (!is.null(result)) {
        pc_choices <- paste0("PC", seq_len(result$n_pcs))
        curr_x <- shiny::isolate(input$pca_x_pc) %||% "PC1"
        curr_y <- shiny::isolate(input$pca_y_pc) %||% "PC2"
        if (!curr_x %in% pc_choices) curr_x <- "PC1"
        if (!curr_y %in% pc_choices) curr_y <- if ("PC2" %in% pc_choices) "PC2" else "PC1"
        shiny::updateSelectInput(session, "pca_x_pc", choices = pc_choices, selected = curr_x)
        shiny::updateSelectInput(session, "pca_y_pc", choices = pc_choices, selected = curr_y)
      }
    }

    # Manual re-run button (e.g. after changing gene selection or colour variable)
    shiny::observeEvent(input$run_pca, {
      shiny::req(input$pca_genes, length(input$pca_genes) > 0)
      genes_plot <- shiny::isolate(input$pca_genes)
      cv         <- shiny::isolate(input$pca_color_var)
      color_var  <- if (!is.null(cv) && nzchar(trimws(cv %||% ""))) cv else NULL
      .run_pca_now(genes_plot, color_var)
    })

    shiny::observeEvent(input$run_rocauc, {
      expr_d <- get_expr(); meta <- get_meta()
      shiny::req(expr_d, meta, input$roc_cond_var, input$roc_cond_class)

      roc_mode <- shiny::isolate(input$roc_mode) %||% "n_based"
      cond_var   <- shiny::isolate(input$roc_cond_var)
      cond_class <- shiny::isolate(input$roc_cond_class)

      if (isTRUE(roc_mode == "manual")) {
        shiny::req(input$roc_genes, length(input$roc_genes) > 0)
        genes_plot <- shiny::isolate(input$roc_genes)
        gene_sel   <- "manual"
        max_g      <- NULL
      } else {
        gs  <- get_gene_sets(); sel <- input$selected_gs
        shiny::req(gs, sel, sel %in% names(gs))
        genes_plot <- .gs_genes(gs[[sel]])
        gene_sel   <- shiny::isolate(input$roc_gene_sel) %||% "top_auc"
        max_g      <- shiny::isolate(input$roc_max_genes)
      }

      # Drop NA samples for the condition variable
      valid_rows <- !is.na(meta[[cond_var]])
      if (sum(valid_rows) < 4) {
        shiny::showNotification("Too few non-NA samples for ROC analysis.", type = "error", duration = 8)
        return()
      }
      meta_roc <- meta[valid_rows, , drop = FALSE]
      # Align expression columns to filtered metadata
      shared_s  <- intersect(colnames(expr_d), as.character(meta_roc[[1]]))
      meta_roc  <- meta_roc[as.character(meta_roc[[1]]) %in% shared_s, , drop = FALSE]
      expr_roc  <- expr_d[, shared_s, drop = FALSE]

      # Keep only genes present in expression data with non-zero variance
      genes_plot <- genes_plot[genes_plot %in% rownames(expr_roc)]
      if (length(genes_plot) == 0) {
        shiny::showNotification("None of the selected genes found in the expression matrix.",
                                type = "error", duration = 8)
        return()
      }
      gene_vars  <- apply(expr_roc[genes_plot, , drop = FALSE], 1, stats::var, na.rm = TRUE)
      genes_plot <- genes_plot[is.finite(gene_vars) & gene_vars > 0]
      if (length(genes_plot) == 0) {
        shiny::showNotification("All selected genes have zero variance across samples.",
                                type = "error", duration = 8)
        return()
      }

      # For "top_auc": compute all genes first, then filter by AUC afterwards.
      # For "random": pre-sample before computing (expensive).
      # For "manual": use as-is.
      if (!isTRUE(gene_sel == "manual") &&
          gene_sel == "random" && !is.null(max_g) && max_g > 0 && length(genes_plot) > max_g) {
        genes_plot <- sample(genes_plot, min(max_g, length(genes_plot)))
      }
      n_genes    <- length(genes_plot)
      roc_ncol   <- max(1L, min(4L, as.integer(ceiling(sqrt(n_genes)))))
      # Build a stable gene→colour mapping in the original gene order so that
      # the same gene always gets the same colour in BOTH the ROC curves and AUC plots.
      stable_colors <- stats::setNames(rep_len(.pp_palette, length(genes_plot)), genes_plot)

      shiny::withProgress(message = "Computing discriminatory power...", value = 0.3, {
        tmp_pdf <- tempfile(fileext = ".pdf")
        grDevices::pdf(tmp_pdf, width = 8, height = 6)
        result <- tryCatch({
          res <- ROCandAUCplot(
            data = expr_roc, metadata = meta_roc, genes = genes_plot,
            condition_var = cond_var, class = cond_class,
            group_var = NULL, plot_type = "all",
            invert_auc = isTRUE(shiny::isolate(input$roc_invert_auc)),
            roc_params = list(ncol = roc_ncol))
          shiny::incProgress(0.7)

          # Post-filter for "top_auc": keep only genes with the highest effective AUC.
          # Always rank by pmax(AUC, 1 - AUC) so that genes with AUC near 0 are
          # treated as equally discriminatory as genes with AUC near 1, regardless
          # of whether the user has checked "Invert AUC". The invert checkbox only
          # controls how the AUC is *displayed*, not which genes are selected.
          if (gene_sel == "top_auc" && !is.null(max_g) && max_g > 0) {
            auc_df <- res$auc_values
            if (!is.null(auc_df) && is.data.frame(auc_df) && nrow(auc_df) > max_g) {
              effective_auc <- pmax(auc_df$AUC, 1 - auc_df$AUC)
              top_genes <- auc_df$Gene[order(effective_auc, decreasing = TRUE)][seq_len(max_g)]
              res$auc_values <- auc_df[auc_df$Gene %in% top_genes, , drop = FALSE]
              if (!is.null(res$roc_plot) && is.data.frame(res$roc_plot$data))
                res$roc_plot$data <- res$roc_plot$data[res$roc_plot$data$Gene %in% top_genes, ]
              n_genes <- length(top_genes)
              roc_ncol <- max(1L, min(4L, as.integer(ceiling(sqrt(n_genes)))))
              stable_colors <- stable_colors[top_genes]
            }
          }

          list(res = res, n_genes = n_genes, roc_ncol = roc_ncol,
               cond_var = cond_var, cond_class = cond_class,
               stable_colors = stable_colors)
        }, error = function(e) {
          shiny::showNotification(paste("Discriminatory power failed:", conditionMessage(e)),
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
      interactive <- isTRUE(input$sim_interactive)
      shiny::tagList(
        shiny::h5("Similarity Heatmap", style = "margin:12px 0 4px;"),
        if (interactive)
          plotly::plotlyOutput(ns("similarity_plot"), height = paste0(sim_h(), "px"))
        else
          shiny::plotOutput(ns("similarity_plot_static"), width = "100%", height = paste0(sim_h(), "px"))
      )
    })

    output$similarity_plot <- plotly::renderPlotly({
      shiny::req(isTRUE(input$sim_interactive))
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

    output$similarity_plot_static <- shiny::renderPlot({
      shiny::req(!isTRUE(input$sim_interactive))
      res  <- similarity_result(); shiny::req(!is.null(res))
      font <- input$sim_fontsize
      metric <- shiny::isolate(input$gs_metric)
      org_label <- if (!is.null(res$organism)) res$organism else "Homo sapiens"
      if (metric == "odds_ratio" && !is.null(res$data))
        .build_or_heatmap(res$data, font = font,
          title = paste0("Similarity vs MSigDB (Odds Ratio) — ", org_label))
      else if (!is.null(res$data))
        .build_ji_heatmap(res$data, font = font,
          title = paste0("Similarity vs MSigDB (Jaccard Index) — ", org_label))
      else shiny::req(FALSE)
    }, height = function() sim_h())

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
                   " Select a gene set and grouping variable above,",
                   " then click ", shiny::strong("Plot Violins"), "."))
      NULL
    })

    output$expr_plot_ui <- shiny::renderUI({
      res <- expression_result()
      if (is.null(res)) return(NULL)
      shiny::tagList(
        shiny::div(
          class = "d-flex align-items-center justify-content-between",
          style = "margin:12px 0 4px;",
          shiny::h5(paste("Per-gene Violin Plots –", res$gs_name), style = "margin:0;"),
          shiny::downloadButton(ns("download_violin"), "",
            icon = shiny::icon("download"), class = "btn-sm btn-outline-secondary")
        ),
        shiny::plotOutput(ns("violin_plot"), height = paste0(expr_h(), "px"))
      )
    })

    .render_violin <- function(res, font) {
      p <- res$violin$plot
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
    }

    output$violin_plot <- shiny::renderPlot({
      res <- expression_result(); shiny::req(!is.null(res), !is.null(res$violin))
      font <- if (!is.null(input$expr_fontsize)) input$expr_fontsize else 10L
      .render_violin(res, font)
    }, height = function() expr_h(), res = 150)

    output$download_violin <- shiny::downloadHandler(
      filename = function() paste0("violin_plots_", Sys.Date(), ".png"),
      content  = function(file) {
        res <- expression_result(); shiny::req(!is.null(res), !is.null(res$violin))
        font <- if (!is.null(input$expr_fontsize)) input$expr_fontsize else 10L
        h    <- max(500L, expr_h())
        n_g  <- length(input$expr_genes)
        w    <- max(800L, as.integer(n_g * 120 + 200))
        grDevices::png(file, width = w, height = h, res = 150)
        .render_violin(res, font)
        grDevices::dev.off()
      }
    )

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
      interactive <- isTRUE(input$hm_interactive)
      shiny::tagList(
        shiny::div(
          class = "d-flex align-items-center justify-content-between",
          style = "margin: 12px 0 6px;",
          shiny::h5(paste("Expression Heatmap —",
                          if (nchar(res$gs_name) > 55)
                            paste0(substr(res$gs_name, 1, 52), "…")
                          else res$gs_name),
                    style = "margin:0;", title = res$gs_name),
          if (!interactive)
            shiny::downloadButton(ns("download_heatmap"), "Download",
              icon = shiny::icon("download"), class = "btn-sm btn-outline-secondary")
        ),
        if (interactive)
          plotly::plotlyOutput(ns("heatmap_plotly"), height = paste0(hm_h(), "px"))
        else
          shiny::plotOutput(ns("heatmap_plot"), height = paste0(hm_h(), "px"))
      )
    })

    output$heatmap_plot <- shiny::renderPlot({
      res <- heatmap_result(); shiny::req(!is.null(res), !is.null(res$heatmap))
      hm <- res$heatmap
      ComplexHeatmap::draw(
        hm$plot,
        heatmap_legend_side    = hm$scale_side %||% "right",
        annotation_legend_side = hm$annot_side %||% "top"
      )
    }, height = function() hm_h(), res = 150)

    output$heatmap_plotly <- plotly::renderPlotly({
      res <- heatmap_result(); shiny::req(!is.null(res), !is.null(res$heatmap))
      mat <- res$heatmap$data
      shiny::req(is.matrix(mat) || is.data.frame(mat))
      mat <- as.matrix(mat)
      cl_r <- isTRUE(res$cl_rows); cl_c <- isTRUE(res$cl_cols)
      font_sz <- if (!is.null(input$hm_fontsize)) input$hm_fontsize else 10L

      # Remove non-finite rows before clustering
      finite_rows <- apply(mat, 1, function(r) all(is.finite(r)))
      if (!all(finite_rows)) mat <- mat[finite_rows, , drop = FALSE]
      shiny::req(nrow(mat) > 0)

      # Annotation info stored at compute time
      has_annot <- !is.null(res$group_var) && !is.null(res$annot_colors) &&
                   !is.null(res$sample_groups) && length(res$annot_colors) > 0

      # Cluster
      hc_r <- hc_c <- NULL
      if (cl_r && nrow(mat) > 1) { hc_r <- hclust(dist(mat));    mat <- mat[hc_r$order, , drop = FALSE] }
      if (cl_c && ncol(mat) > 1) { hc_c <- hclust(dist(t(mat))); mat <- mat[, hc_c$order, drop = FALSE] }
      nr <- nrow(mat); nc <- ncol(mat)
      # Gene-name labels are placed on the RIGHT side of the heatmap panel (side="right" below),
      # so the row dendrogram panel only needs to accommodate the tree lines, not the labels.
      rdend_w <- if (!is.null(hc_r)) 0.12 else 0
      cdend_h <- if (!is.null(hc_c)) 0.14 else 0
      heat_h  <- 1 - cdend_h   # annotation is embedded inside the heatmap, not a separate panel

      # Prepare annotation colorscale + codes (computed after column re-ordering)
      ann_codes <- NULL; ann_cscale <- NULL; ann_text <- NULL
      levs_ann  <- NULL; colors_ann <- NULL
      if (has_annot) {
        ann_vec    <- res$sample_groups[colnames(mat)]
        ann_vec[is.na(ann_vec)] <- "Unknown"
        levs_ann   <- names(res$annot_colors)
        colors_ann <- unname(res$annot_colors)
        n_l  <- length(levs_ann)
        ann_codes  <- match(ann_vec, levs_ann); ann_codes[is.na(ann_codes)] <- 1L
        # Solid discrete colorscale: two stops per level
        ann_cscale <- do.call(c, lapply(seq_along(levs_ann), function(i)
          list(list((i - 1) / n_l, colors_ann[i]), list(i / n_l, colors_ann[i]))))
        ann_text   <- matrix(paste0(res$group_var, ": ", ann_vec), nrow = 1L)
      }

      # y-axis: gene rows at y = 1..nr; annotation row at y = 0 (visually above genes)
      y_min      <- if (has_annot) -0.5 else 0.5
      y_tickvals <- if (has_annot) c(0L, seq_len(nr)) else seq_len(nr)
      y_ticktext <- if (has_annot) c("", rownames(mat)) else rownames(mat)

      # Hover text for main heatmap
      text_mat <- matrix(
        paste0("Gene: ", rep(rownames(mat), times = nc),
               "<br>Sample: ", rep(colnames(mat), each = nr),
               "<br>Z-score: ", sprintf("%.2f", as.vector(mat))),
        nrow = nr, ncol = nc)

      # Main heatmap — genes at y = 1..nr
      # height= set here (not in layout) to avoid the deprecation warning while
      # still allowing the figure to scale with the number of genes.
      p_heat <- plotly::plot_ly(
        height = hm_h(),
        z = mat, x = seq_len(nc), y = seq_len(nr),
        type = "heatmap",
        colorscale = list(c(0,"#4173B4"), c(0.5,"white"), c(1,"#B44141")),
        zmin = -2, zmax = 2,
        colorbar = list(title = "Z-score", titleside = "right"),
        text = text_mat, hovertemplate = "%{text}<extra></extra>",
        showscale = TRUE
      ) |> plotly::layout(
        xaxis = list(tickvals=seq_len(nc), ticktext=colnames(mat),
                     tickfont=list(size=font_sz), title=""),
        # side="right" moves gene labels to the right of the heatmap body,
        # completely away from the row dendrogram panel on the left.
        yaxis = list(tickvals=y_tickvals, ticktext=y_ticktext,
                     tickfont=list(size=font_sz), title="", side="right",
                     autorange=FALSE, range=c(nr + 0.5, y_min)))

      # Overlay annotation as a second heatmap trace at y = 0
      if (has_annot) {
        p_heat <- p_heat |> plotly::add_trace(
          z    = matrix(ann_codes, nrow = 1L), x = seq_len(nc), y = 0,
          type = "heatmap", colorscale = ann_cscale,
          zmin = 0.5, zmax = length(levs_ann) + 0.5,
          showscale = FALSE, text = ann_text,
          hovertemplate = "%{text}<extra></extra>",
          inherit = FALSE)
        # Invisible legend squares for each condition group
        for (i in seq_along(levs_ann))
          p_heat <- p_heat |> plotly::add_trace(
            x=NA_real_, y=NA_real_, type="scatter", mode="markers",
            marker=list(color=colors_ann[i], size=10, symbol="square"),
            name=levs_ann[i], showlegend=TRUE, hoverinfo="none",
            inherit=FALSE)
      }

      if (is.null(hc_r) && is.null(hc_c))
        return(p_heat |> plotly::layout(margin=list(l=20, r=160, b=80),
                                        showlegend=has_annot,
                                        legend=list(orientation="v", x=1.18, y=0.5)))

      # Dendrogram helpers
      .seg_trace <- function(segs, swap_xy=FALSE) {
        xs <- as.vector(rbind(segs$x0, segs$x1, NA_real_))
        ys <- as.vector(rbind(segs$y0, segs$y1, NA_real_))
        if (swap_xy) { tmp<-xs; xs<-ys; ys<-tmp }
        plotly::plot_ly(x=xs, y=ys, type="scatter", mode="lines",
                        color=I("grey40"), line=list(width=1),
                        showlegend=FALSE, hoverinfo="skip")
      }
      empty_p <- plotly::plot_ly(x=0, y=0, type="scatter", mode="markers",
        marker=list(color="rgba(0,0,0,0)", size=1), showlegend=FALSE, hoverinfo="none") |>
        plotly::layout(
          xaxis=list(showgrid=FALSE, zeroline=FALSE, showticklabels=FALSE, range=c(-1,1)),
          yaxis=list(showgrid=FALSE, zeroline=FALSE, showticklabels=FALSE, range=c(-1,1)))

      p_rdend <- NULL; p_cdend <- NULL
      if (!is.null(hc_r))
        # y-range must match heatmap exactly so gene leaves align correctly
        # (extends to y_min to cover the annotation row at y = 0 when present)
        p_rdend <- .seg_trace(.hclust_to_segments(hc_r), swap_xy=TRUE) |> plotly::layout(
          xaxis=list(autorange="reversed", showgrid=FALSE, zeroline=FALSE, showticklabels=FALSE),
          yaxis=list(range=c(nr + 0.5, y_min), showgrid=FALSE, zeroline=FALSE, showticklabels=FALSE))
      if (!is.null(hc_c))
        p_cdend <- .seg_trace(.hclust_to_segments(hc_c)) |> plotly::layout(
          xaxis=list(range=c(0.5, nc + 0.5), showgrid=FALSE, zeroline=FALSE, showticklabels=FALSE),
          yaxis=list(showgrid=FALSE, zeroline=FALSE, showticklabels=FALSE))

      # Subplot grid: rows = [cdend?][heat], cols = [rdend?][heat]
      # Annotation is now inside p_heat — max 2 rows × 2 cols
      has_rdend <- !is.null(p_rdend); has_cdend <- !is.null(p_cdend)
      row_heights <- c(if (has_cdend) cdend_h, heat_h)
      col_widths  <- c(if (has_rdend) rdend_w, 1 - rdend_w)
      n_rows <- length(row_heights)

      plot_list <- list()
      if (has_cdend) plot_list <- c(plot_list, if (has_rdend) list(empty_p, p_cdend) else list(p_cdend))
      plot_list <- c(plot_list, if (has_rdend) list(p_rdend, p_heat) else list(p_heat))

      sub_args <- list(nrows=n_rows, margin=0, shareX=FALSE, shareY=FALSE)
      if (length(col_widths) > 1) sub_args$widths  <- col_widths
      if (n_rows > 1)             sub_args$heights <- row_heights
      fig <- do.call(plotly::subplot, c(plot_list, sub_args)) |>
        plotly::layout(margin=list(r=160),
                       showlegend=has_annot,
                       legend=list(orientation="v", x=1.18, y=0.5))

      # Link dendrogram axes to the heatmap axes so zoom/pan stays in sync.
      # Axis numbering follows subplot plot_list order (1-based, left→right, top→bottom).
      # Layout:
      #   has_cdend & has_rdend → 2×2: [1=empty, 2=cdend | 3=rdend, 4=heat]
      #   !has_cdend & has_rdend → 1×2: [1=rdend, 2=heat]
      #   has_cdend & !has_rdend → 2×1: [1=cdend, 2=heat]  (nrows=2, 1 col)
      # Use plotly::layout() for deep-merge so existing axis properties are preserved.
      if (has_rdend && has_cdend) {
        fig <- fig |>
          plotly::layout(yaxis3 = list(matches = "y4"),
                         xaxis2 = list(matches = "x4"))
      } else if (has_rdend) {
        fig <- fig |> plotly::layout(yaxis = list(matches = "y2"))
      } else if (has_cdend) {
        fig <- fig |> plotly::layout(xaxis = list(matches = "x2"))
      }
      fig
    })

    output$download_heatmap <- shiny::downloadHandler(
      filename = function() paste0("expression_heatmap_", Sys.Date(), ".png"),
      content  = function(file) {
        res <- heatmap_result(); shiny::req(!is.null(res), !is.null(res$heatmap))
        hm <- res$heatmap
        grDevices::png(file, width = 1400, height = hm_h(), res = 150)
        ComplexHeatmap::draw(
          hm$plot,
          heatmap_legend_side    = hm$scale_side %||% "right",
          annotation_legend_side = hm$annot_side %||% "top"
        )
        grDevices::dev.off()
      }
    )

    # ---- Correlation --------------------------------------------------------

    output$corr_status_ui <- shiny::renderUI({
      if (is.null(corr_result()))
        return(shiny::div(class = "alert alert-info", style = "margin:10px 0;",
                   shiny::icon("circle-info"),
                   " Configure settings above and click ",
                   shiny::strong("Plot Correlation Heatmap"), "."))
      NULL
    })

    output$corr_gene_count_warn_ui <- shiny::renderUI({
      genes_sel <- input$corr_genes
      n <- if (!is.null(genes_sel)) length(genes_sel) else 0L
      if (n > 50)
        shiny::div(class = "alert alert-warning", style = "font-size:0.82em; padding:5px 9px; margin-bottom:4px;",
          shiny::icon("triangle-exclamation"), " ",
          sprintf("%d genes selected. Correlation heatmaps with >50 genes can be hard to interpret — consider selecting a subset.", n))
    })

    output$corr_plot_ui <- shiny::renderUI({
      res <- corr_result()
      if (is.null(res)) return(NULL)
      interactive <- isTRUE(input$corr_interactive)
      shiny::tagList(
        shiny::div(
          class = "d-flex align-items-center justify-content-between",
          style = "margin: 12px 0 6px;",
          shiny::span(shiny::strong("Pairwise Gene Correlation")),
          if (!interactive)
            shiny::downloadButton(ns("download_corr"), "Download",
              icon = shiny::icon("download"), class = "btn-sm btn-outline-secondary")
        ),
        if (interactive)
          plotly::plotlyOutput(ns("corr_plotly"), height = paste0(corr_h(), "px"))
        else
          shiny::plotOutput(ns("corr_plot"), height = paste0(corr_h(), "px"))
      )
    })

    # renderPlot reads input$corr_fontsize and input$show_col_names directly,
    # so changes to those inputs cause an instant re-draw WITHOUT re-computing
    # the correlation matrices (expensive step already done in run_corr).
    output$corr_plot <- shiny::renderPlot({
      res <- corr_result(); shiny::req(!is.null(res), !is.null(res$corrmat))
      font     <- if (!is.null(input$corr_fontsize))     input$corr_fontsize     else 10L
      show_row <- if (!is.null(input$corr_show_row_names)) input$corr_show_row_names else TRUE
      show_col <- if (!is.null(input$show_col_names))    input$show_col_names    else TRUE
      tryCatch(
        .build_corr_heatmap(res$corrmat, method = res$method,
                             font = font, show_row_names = show_row, show_col_names = show_col,
                             cluster_rows    = res$cluster_rows    %||% TRUE,
                             cluster_columns = res$cluster_columns %||% TRUE),
        error = function(e) shiny::showNotification(
          paste("Correlation plot failed:", e$message), type = "warning", duration = 8))
    }, height = function() corr_h(), res = 150)

    output$corr_plotly <- plotly::renderPlotly({
      res <- corr_result(); shiny::req(!is.null(res), !is.null(res$corrmat))
      cl_r <- res$cluster_rows    %||% TRUE
      cl_c <- res$cluster_columns %||% TRUE
      font_sz  <- if (!is.null(input$corr_fontsize))       input$corr_fontsize      else 10L
      show_row <- if (!is.null(input$corr_show_row_names)) input$corr_show_row_names else TRUE
      show_col <- if (!is.null(input$show_col_names))      input$show_col_names      else TRUE
      method_lbl <- switch(res$method %||% "spearman",
        spearman = "Spearman Corr.", pearson = "Pearson Corr.", kendall = "Kendall Corr.", "Correlation")

      # Helper: build one plotly heatmap for a single correlation matrix (no dendrogram)
      .make_one_corr <- function(mat_i, cond_title = NULL, show_y = TRUE, show_cbar = TRUE) {
        if (cl_r && nrow(mat_i) > 1) {
          hc_i <- hclust(dist(mat_i)); mat_i <- mat_i[hc_i$order, hc_i$order, drop = FALSE]
        }
        nr_i <- nrow(mat_i); nc_i <- ncol(mat_i)
        txt_i <- matrix(
          paste0(rep(rownames(mat_i), times=nc_i), " x ", rep(colnames(mat_i), each=nr_i),
                 "<br>", method_lbl, ": ", sprintf("%.3f", as.vector(mat_i))),
          nrow = nr_i, ncol = nc_i)
        p <- plotly::plot_ly(
          height = corr_h(),
          z = mat_i, x = seq_len(nc_i), y = seq_len(nr_i), type = "heatmap",
          colorscale = list(c(0,"#4173B4"), c(0.5,"white"), c(1,"#B44141")),
          zmin = -1, zmax = 1,
          colorbar = list(title = method_lbl), showscale = show_cbar,
          text = txt_i, hovertemplate = "%{text}<extra></extra>"
        ) |> plotly::layout(
          xaxis = list(tickvals=seq_len(nc_i), ticktext=colnames(mat_i),
                       tickfont=list(size=font_sz), title="", showticklabels=show_col),
          yaxis = list(tickvals=seq_len(nr_i), ticktext=rownames(mat_i),
                       tickfont=list(size=font_sz), title="",
                       autorange="reversed", showticklabels=show_y))
        if (!is.null(cond_title))
          p <- p |> plotly::layout(
            annotations = list(list(text=cond_title, x=0.5, y=1.04,
              xref="paper", yref="paper", showarrow=FALSE,
              font=list(size=font_sz+2, color="#333"))))
        p
      }

      corrmat_data <- res$corrmat

      # ---- Multiple conditions (separate.by was set) ----
      if (is.list(corrmat_data) && !is.matrix(corrmat_data) && !is.data.frame(corrmat_data)) {
        conds <- names(corrmat_data)
        n_c   <- length(conds)
        plots <- lapply(seq_along(conds), function(i)
          .make_one_corr(as.matrix(corrmat_data[[i]]),
                         cond_title = conds[i],
                         show_y     = (i == 1),
                         show_cbar  = (i == n_c)))
        if (n_c == 1L) return(plots[[1]])
        do.call(plotly::subplot, c(plots, list(
          nrows=1L, margin=0.04, shareX=FALSE, shareY=TRUE
        ))) |> plotly::layout(showlegend=FALSE)

      # ---- Single condition ----
      } else {
        mat <- if (is.matrix(corrmat_data) || is.data.frame(corrmat_data))
                 as.matrix(corrmat_data)
               else return(NULL)

        # For a symmetric matrix, row and col clustering are identical
        hc_r_obj <- hc_c_obj <- NULL
        if (cl_r && nrow(mat) > 1) {
          hc_r_obj <- hclust(dist(mat)); mat <- mat[hc_r_obj$order, hc_r_obj$order, drop=FALSE]
        } else if (cl_c && ncol(mat) > 1) {
          hc_c_obj <- hclust(dist(t(mat))); mat <- mat[, hc_c_obj$order, drop=FALSE]
        }
        nr <- nrow(mat); nc <- ncol(mat)
        hc_for_dend <- hc_r_obj %||% hc_c_obj
        # Gene labels on the right (side="right") so the row dendrogram panel
        # only needs to hold tree lines, not labels bleeding into it.
        rdend_w <- if (!is.null(hc_for_dend)) 0.12 else 0
        cdend_h <- if (!is.null(hc_for_dend)) 0.14 else 0

        text_mat <- matrix(
          paste0(rep(rownames(mat), times=nc), " x ", rep(colnames(mat), each=nr),
                 "<br>", method_lbl, ": ", sprintf("%.3f", as.vector(mat))),
          nrow = nr, ncol = nc)

        p_heat <- plotly::plot_ly(
          height = corr_h(),
          z = mat, x = seq_len(nc), y = seq_len(nr), type = "heatmap",
          colorscale = list(c(0,"#4173B4"), c(0.5,"white"), c(1,"#B44141")),
          zmin = -1, zmax = 1,
          colorbar = list(title = method_lbl),
          text = text_mat, hovertemplate = "%{text}<extra></extra>", showscale = TRUE
        ) |> plotly::layout(
          xaxis = list(tickvals=seq_len(nc), ticktext=colnames(mat),
                       tickfont=list(size=font_sz), title="", showticklabels=show_col),
          yaxis = list(tickvals=seq_len(nr), ticktext=rownames(mat),
                       tickfont=list(size=font_sz), title="", side="right",
                       autorange="reversed", showticklabels=show_row))

        if (is.null(hc_for_dend))
          return(p_heat |> plotly::layout(margin=list(l=20, r=160, b=120)))

        .seg_trace <- function(segs, swap_xy=FALSE) {
          xs <- as.vector(rbind(segs$x0, segs$x1, NA_real_))
          ys <- as.vector(rbind(segs$y0, segs$y1, NA_real_))
          if (swap_xy) { tmp<-xs; xs<-ys; ys<-tmp }
          plotly::plot_ly(x=xs, y=ys, type="scatter", mode="lines",
                          color=I("grey40"), line=list(width=1),
                          showlegend=FALSE, hoverinfo="skip")
        }
        empty_p <- plotly::plot_ly(x=0,y=0,type="scatter",mode="markers",
          marker=list(color="rgba(0,0,0,0)",size=1),showlegend=FALSE,hoverinfo="none") |>
          plotly::layout(
            xaxis=list(showgrid=FALSE,zeroline=FALSE,showticklabels=FALSE,range=c(-1,1)),
            yaxis=list(showgrid=FALSE,zeroline=FALSE,showticklabels=FALSE,range=c(-1,1)))

        segs <- .hclust_to_segments(hc_for_dend)
        p_rdend <- .seg_trace(segs, swap_xy=TRUE) |> plotly::layout(
          xaxis=list(autorange="reversed",showgrid=FALSE,zeroline=FALSE,showticklabels=FALSE),
          yaxis=list(range=c(nr+0.5,0.5),showgrid=FALSE,zeroline=FALSE,showticklabels=FALSE))
        p_cdend <- .seg_trace(segs) |> plotly::layout(
          xaxis=list(range=c(0.5,nc+0.5),showgrid=FALSE,zeroline=FALSE,showticklabels=FALSE),
          yaxis=list(showgrid=FALSE,zeroline=FALSE,showticklabels=FALSE))

        fig <- plotly::subplot(empty_p, p_cdend, p_rdend, p_heat,
                        nrows=2L, margin=0,
                        widths=c(rdend_w, 1-rdend_w),
                        heights=c(cdend_h, 1-cdend_h),
                        shareX=FALSE, shareY=FALSE) |>
          plotly::layout(margin=list(r=160), showlegend=FALSE)

        # Link dendrogram axes to heatmap axes for synchronised zoom/pan.
        # 2×2 layout: 1=empty, 2=cdend, 3=rdend, 4=heat
        # Use plotly::layout() deep-merge so existing axis properties are preserved.
        fig |>
          plotly::layout(yaxis3 = list(matches = "y4"),
                         xaxis2 = list(matches = "x4"))
      }
    })

    output$download_corr <- shiny::downloadHandler(
      filename = function() paste0("correlation_heatmap_", Sys.Date(), ".png"),
      content  = function(file) {
        res <- corr_result(); shiny::req(!is.null(res), !is.null(res$corrmat))
        font     <- if (!is.null(input$corr_fontsize))  input$corr_fontsize  else 10L
        show_col <- if (!is.null(input$show_col_names)) input$show_col_names else TRUE
        n_genes  <- res$n_genes %||% 10L
        # Generous sizing so dendrogram, labels, and legend all fit
        h <- max(800L, as.integer(n_genes * 45 + 200L))
        w <- max(h + 300L, 1000L)
        show_row <- if (!is.null(input$corr_show_row_names)) input$corr_show_row_names else TRUE
        grDevices::png(file, width = w, height = h, res = 150)
        .build_corr_heatmap(res$corrmat, method = res$method,
                             font = font, show_row_names = show_row, show_col_names = show_col,
                             cluster_rows    = res$cluster_rows    %||% TRUE,
                             cluster_columns = res$cluster_columns %||% TRUE)
        grDevices::dev.off()
      }
    )

    # ---- Cohen's d ----------------------------------------------------------

    output$cohend_status_ui <- shiny::renderUI({
      if (is.null(cohend_result()))
        return(shiny::div(class = "alert alert-info", style = "margin:10px 0;",
                   shiny::icon("circle-info"),
                   " Set condition variable and positive class, then click ",
                   shiny::strong("Plot Effect Size"), "."))
      NULL
    })

    output$cohend_plot_ui <- shiny::renderUI({
      res <- cohend_result()
      if (is.null(res)) return(NULL)
      shiny::tagList(
        shiny::h5(paste0("Effect Size (Cohen's d) per Gene — ", res$cond_class,
                         " vs rest (", res$cond_var, ")"), style = "margin:12px 0 4px;"),
        plotly::plotlyOutput(ns("cohend_plot"), height = paste0(cohend_h(), "px")))
    })

    output$cohend_plot <- plotly::renderPlotly({
      res  <- cohend_result(); shiny::req(!is.null(res))
      font <- if (!is.null(input$cohend_fontsize)) input$cohend_fontsize else 10L
      df   <- if (is.list(res$res) && is.data.frame(res$res$data)) res$res$data else NULL
      shiny::req(!is.null(df), "CohensD" %in% colnames(df), nrow(df) > 0)
      # Sort ascending so the largest effect is at the top of the horizontal bar chart
      df <- df[order(df$CohensD, decreasing = FALSE), , drop = FALSE]
      df$Gene <- factor(df$Gene, levels = df$Gene)
      col_vals <- stats::setNames(rep_len(.pp_palette, nrow(df)), as.character(df$Gene))
      plotly::plot_ly(
        data        = df,
        x           = ~CohensD,
        y           = ~Gene,
        type        = "bar",
        orientation = "h",
        marker      = list(color = col_vals[as.character(df$Gene)]),
        text        = ~paste0("Gene: ", Gene, "<br>Cohen's d: ", sprintf("%.3f", CohensD)),
        hoverinfo   = "text"
      ) |>
        plotly::layout(
          xaxis     = list(title = list(text = "Effect size (Cohen's d)", font = list(size = font + 1)),
                           tickfont = list(size = font)),
          yaxis     = list(title = "", tickfont = list(size = font), automargin = TRUE),
          showlegend = FALSE,
          margin    = list(l = 10),
          height    = cohend_h()
        )
    })

    # ---- Gene PCA -----------------------------------------------------------

    output$pca_status_ui <- shiny::renderUI({
      if (is.null(pca_result()))
        return(shiny::div(class = "alert alert-info", style = "margin:10px 0;",
                   shiny::icon("circle-info"),
                   " Click ", shiny::strong("Plot Gene PCA"), ".",
                   " Set 'Colour by' in the settings above to colour points by a metadata variable."))
      NULL
    })

    # PCA: combined side-by-side layout (scree left, scatter right)
    output$pca_plots_ui <- shiny::renderUI({
      res <- pca_result()
      if (is.null(res)) return(NULL)
      n_pcs <- res$n_pcs
      # Update PC selector choices to match available PCs
      pc_choices <- paste0("PC", seq_len(n_pcs))
      shiny::tagList(
        shiny::div(
          style = "display:grid; grid-template-columns:1fr 1.6fr; gap:16px; margin-top:12px; align-items:start;",
          # Left: scree plot
          shiny::div(
            shiny::h5("Variance Explained", style = "margin:0 0 4px;"),
            plotly::plotlyOutput(ns("pca_scree"), height = paste0(pca_h(), "px"))
          ),
          # Right: scatter
          shiny::div(
            shiny::h5(
              shiny::uiOutput(ns("pca_scatter_title_ui"), inline = TRUE),
              style = "margin:0 0 4px;"
            ),
            plotly::plotlyOutput(ns("pca_plot"), height = paste0(pca_h(), "px"))
          )
        )
      )
    })

    output$pca_scatter_title_ui <- shiny::renderUI({
      x_pc <- input$pca_x_pc %||% "PC1"
      y_pc <- input$pca_y_pc %||% "PC2"
      shiny::span(paste("Gene PCA —", x_pc, "vs", y_pc))
    })

    output$pca_scree <- plotly::renderPlotly({
      res <- pca_result(); shiny::req(!is.null(res))
      x_pc    <- input$pca_x_pc %||% "PC1"
      y_pc    <- input$pca_y_pc %||% "PC2"
      var_exp <- res$var_exp
      n_show  <- min(length(var_exp), 10L)
      active_pcs <- c(x_pc, y_pc)
      bar_colors <- ifelse(paste0("PC", seq_len(n_show)) %in% active_pcs,
                           "#EBB43E", "#d1d5db")
      df_s <- data.frame(
        PC  = factor(paste0("PC", seq_len(n_show)), levels = paste0("PC", seq_len(n_show))),
        Var = var_exp[seq_len(n_show)]
      )
      plotly::plot_ly(df_s) |>
        plotly::add_bars(x = ~PC, y = ~Var, name = "Variance (%)",
                         marker = list(color = bar_colors),
                         hovertemplate = "%{x}: %{y:.1f}%<extra></extra>") |>
        plotly::layout(
          xaxis      = list(title = "PC", tickfont = list(size = 10)),
          yaxis      = list(title = "Variance (%)"),
          showlegend = FALSE,
          height     = pca_h(),
          margin     = list(t = 40, b = 60),
          title      = list(text = "<span style='font-size:10px;color:grey'>Highlighted = active PCs</span>",
                            x = 0.5, xanchor = "center")
        )
    })

    output$pca_plot <- plotly::renderPlotly({
      res       <- pca_result(); shiny::req(!is.null(res))
      font      <- if (!is.null(input$pca_fontsize)) input$pca_fontsize else 10L
      pca_df    <- res$pca_df
      var_exp   <- res$var_exp
      color_var <- res$color_var
      x_pc      <- input$pca_x_pc %||% "PC1"
      y_pc      <- input$pca_y_pc %||% "PC2"
      x_idx     <- as.integer(sub("PC", "", x_pc))
      y_idx     <- as.integer(sub("PC", "", y_pc))
      shiny::req(!is.null(pca_df), x_pc %in% colnames(pca_df), y_pc %in% colnames(pca_df))

      x_lab <- sprintf("%s (%.1f%%)", x_pc, var_exp[x_idx])
      y_lab <- sprintf("%s (%.1f%%)", y_pc, var_exp[y_idx])

      sample_lbl <- if ("Sample" %in% colnames(pca_df)) as.character(pca_df$Sample) else
                      as.character(seq_len(nrow(pca_df)))

      pca_df$.x <- pca_df[[x_pc]]
      pca_df$.y <- pca_df[[y_pc]]

      if (!is.null(color_var) && color_var %in% colnames(pca_df)) {
        levs       <- unique(as.character(pca_df[[color_var]]))
        color_vals <- stats::setNames(rep_len(.pp_palette, length(levs)), levs)
        pca_df[[color_var]] <- factor(pca_df[[color_var]], levels = levs)
        hover_text <- paste0("Sample: ", sample_lbl,
                             "<br>", color_var, ": ", as.character(pca_df[[color_var]]))
        p <- ggplot2::ggplot(pca_df,
               ggplot2::aes(x = .data$.x, y = .data$.y,
                            fill = .data[[color_var]], text = hover_text)) +
          ggplot2::geom_point(size = 4, alpha = 0.8, shape = 21, color = "black") +
          ggplot2::scale_fill_manual(values = color_vals, name = color_var)
      } else {
        p <- ggplot2::ggplot(pca_df,
               ggplot2::aes(x = .data$.x, y = .data$.y,
                            text = paste0("Sample: ", sample_lbl))) +
          ggplot2::geom_point(size = 4, alpha = 0.8, shape = 21,
                              color = "black", fill = .pp_palette[1])
      }

      p <- p +
        ggplot2::labs(x = x_lab, y = y_lab) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dotted", color = "grey60") +
        ggplot2::geom_hline(yintercept = 0, linetype = "dotted", color = "grey60") +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text     = ggplot2::element_text(size = font),
                       axis.title    = ggplot2::element_text(size = font + 1),
                       legend.text   = ggplot2::element_text(size = font),
                       legend.position = "bottom")

      plotly::ggplotly(p, tooltip = "text") |>
        plotly::layout(height = pca_h(), legend = list(orientation = "h"))
    })

    # ---- ROC / AUC ----------------------------------------------------------

    output$rocauc_status_ui <- shiny::renderUI({
      if (is.null(rocauc_result()))
        return(shiny::div(class = "alert alert-info", style = "margin:10px 0;",
                   shiny::icon("circle-info"),
                   " Set condition variable and positive class, then click ",
                   shiny::strong("Plot Discriminatory Power"), ".",
                   " ROC curves and an AUC barplot will both be shown."))
      NULL
    })

    output$roc_curves_ui <- shiny::renderUI({
      res <- rocauc_result()
      if (is.null(res)) return(NULL)
      shiny::tagList(
        shiny::div(
          class = "d-flex align-items-center justify-content-between",
          style = "margin: 12px 0 6px;",
          shiny::h5(paste0("ROC Curves — ", res$cond_class, " vs rest (", res$cond_var, ")"),
                    style = "margin:0;"),
          shiny::downloadButton(ns("download_roc"), "Download",
            icon = shiny::icon("download"), class = "btn-sm btn-outline-secondary")
        ),
        shiny::plotOutput(ns("roc_curves_plot"), height = paste0(roc_curves_h(), "px"))
      )
    })

    output$roc_curves_plot <- shiny::renderPlot({
      res <- rocauc_result(); shiny::req(!is.null(res))
      raw <- res$res
      roc_p   <- raw$roc_plot
      shiny::req(!is.null(roc_p), inherits(roc_p, c("gg", "ggplot")))
      roc_df <- roc_p$data
      shiny::req(!is.null(roc_df), nrow(roc_df) > 0)

      n_genes_shown <- length(unique(roc_df$Gene))
      ncol_shown    <- max(1L, min(4L, as.integer(ceiling(sqrt(n_genes_shown)))))
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

    output$download_roc <- shiny::downloadHandler(
      filename = function() paste0("roc_curves_", Sys.Date(), ".png"),
      content  = function(file) {
        res <- rocauc_result(); shiny::req(!is.null(res))
        raw   <- res$res
        roc_p <- raw$roc_plot
        shiny::req(!is.null(roc_p), inherits(roc_p, c("gg", "ggplot")))
        h <- max(600L, roc_curves_h())
        grDevices::png(file, width = 1400, height = h, res = 150)
        print(roc_p)
        grDevices::dev.off()
      }
    )

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
        ggplot2::geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey50",
                             linewidth = 0.7) +
        { if (!isTRUE(input$roc_invert_auc))
            ggplot2::geom_vline(xintercept = 0, linetype = "solid",
                                color = "grey70", linewidth = 0.4) } +
        ggplot2::coord_cartesian(
          xlim = c(if (isTRUE(input$roc_invert_auc)) 0.5 else -0.03, 1)) +
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
