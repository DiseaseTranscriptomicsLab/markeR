# =============================================================================
# markeR_preprocessing.R  —  Sequential Preprocessing Tab (v6)
# =============================================================================
#' @importFrom shiny moduleServer NS reactiveValues reactiveVal reactive observe
#'   observeEvent req isolate renderPlot renderUI uiOutput plotOutput tagList
#'   tags div p h5 numericInput radioButtons selectInput checkboxInput
#'   checkboxGroupInput actionButton showNotification withProgress incProgress
#'   downloadButton downloadHandler fluidRow column HTML hr helpText validate
#'   need updateNumericInput showModal modalDialog modalButton
#' @importFrom bslib card card_header layout_sidebar sidebar
#' @importFrom ggplot2 ggplot aes geom_line geom_density geom_boxplot geom_vline
#'   geom_hline geom_point annotate theme_bw theme element_blank element_line
#'   element_text xlab ylab ggtitle labs margin
#' @importFrom edgeR DGEList calcNormFactors cpm
#' @importFrom limma voom lmFit
#' @importFrom stats prcomp model.matrix reformulate setNames median
#' @importFrom utils write.csv head tail

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !is.na(a[1])) a else b

# =============================================================================
# INTERNAL HELPERS
# =============================================================================

.pp_complexity_byexprgenes <- function(data, samples_to_plot,
                                        reads_in_percentage = FALSE,
                                        progress_fn = NULL) {
  n_s     <- length(samples_to_plot)
  df_list <- lapply(seq_along(samples_to_plot), function(i) {
    samp <- samples_to_plot[i]
    if (!is.null(progress_fn))
      progress_fn(i / n_s, paste0("Processing sample ", i, "/", n_s, " — ", samp))
    vals <- sort(data[, samp], decreasing = TRUE)
    data.frame(Sample = samp, CumulativeReads = cumsum(vals),
               GeneIndex = seq_along(vals), stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, df_list)
  if (!reads_in_percentage) {
    plt <- ggplot2::ggplot(df,
        ggplot2::aes(x = GeneIndex, y = log10(CumulativeReads + 1),
                     group = Sample, colour = Sample)) +
      ggplot2::geom_line(alpha = 0.75, linewidth = 0.7) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none",
                     panel.grid.minor = ggplot2::element_blank(),
                     panel.border     = ggplot2::element_blank(),
                     axis.line        = ggplot2::element_line(colour = "black")) +
      ggplot2::xlab("Gene index (desc.)") + ggplot2::ylab("Cumulative reads (log₁₀)")
  } else {
    totals <- tapply(df$CumulativeReads, df$Sample, max)
    df$Pct <- 100 * df$CumulativeReads / totals[df$Sample]
    plt <- ggplot2::ggplot(df,
        ggplot2::aes(x = GeneIndex, y = Pct, group = Sample, colour = Sample)) +
      ggplot2::geom_line(alpha = 0.75, linewidth = 0.7) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none",
                     panel.grid.minor = ggplot2::element_blank(),
                     panel.border     = ggplot2::element_blank(),
                     axis.line        = ggplot2::element_line(colour = "black")) +
      ggplot2::xlab("Gene index (desc.)") + ggplot2::ylab("Cumulative reads (%)")
  }
  list(plot = plt, data_to_plot = df)
}

.pp_complexity_samplingreads <- function(data, NumberGenes, NumberIterations,
                                          stepReads, samples_to_plot,
                                          reads_in_percentage = FALSE,
                                          progress_fn = NULL) {
  NumberGenes <- min(as.integer(NumberGenes), 1000L, nrow(data))
  gene_idx    <- sample(seq_len(nrow(data)), NumberGenes)
  n_s <- length(samples_to_plot)
  out <- data.frame()
  for (si in seq_along(samples_to_plot)) {
    samp <- samples_to_plot[si]
    if (!is.null(progress_fn))
      progress_fn(si / n_s, paste0("Sampling reads for sample ", si, "/", n_s, " — ", samp))
    sub     <- data.frame(gene = rownames(data)[gene_idx],
                           count = data[gene_idx, samp], stringsAsFactors = FALSE)
    vec     <- unlist(mapply(rep, sub$gene, ceiling(sub$count), SIMPLIFY = FALSE))
    libsize <- length(vec)
    if (libsize < 1L) next
    vec   <- sample(vec)
    steps <- unique(c(seq(1L, libsize, by = max(1L, as.integer(stepReads))), libsize))
    n_unique <- vapply(steps, function(n)
      mean(replicate(NumberIterations,
                     length(unique(sample(vec, n, replace = FALSE))))), numeric(1))
    out <- rbind(out, data.frame(Sample = samp, NumberReads = steps,
      PercentageReads = steps / libsize, NumberUniqueGenesDetected = n_unique,
      stringsAsFactors = FALSE))
  }
  if (nrow(out) == 0) return(NULL)
  if (reads_in_percentage) {
    plt <- ggplot2::ggplot(out,
      ggplot2::aes(x = PercentageReads, y = NumberUniqueGenesDetected,
                   colour = Sample, group = Sample)) +
      ggplot2::geom_line() + ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none") +
      ggplot2::xlab("Fraction of reads") + ggplot2::ylab("Unique genes detected")
  } else {
    plt <- ggplot2::ggplot(out,
      ggplot2::aes(x = NumberReads, y = NumberUniqueGenesDetected,
                   colour = Sample, group = Sample)) +
      ggplot2::geom_line() + ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none") +
      ggplot2::xlab("Number of reads") + ggplot2::ylab("Unique genes detected")
  }
  list(plot = plt, data_to_plot = out)
}

# ── PCA: compute once, render many times ──────────────────────────────────────

.pp_run_pca <- function(log_counts) {
  pca <- stats::prcomp(t(as.matrix(log_counts)), scale = FALSE, center = TRUE)
  list(
    df = as.data.frame(pca$x[, 1:min(2L, ncol(pca$x)), drop = FALSE]),
    ev = pca$sdev^2
  )
}

.pp_pca_from_df <- function(pca_df, ev, meta, color_by = NULL, title = "") {
  df     <- pca_df
  pc1pct <- round(100 * ev[1] / sum(ev), 1)
  pc2pct <- if (length(ev) >= 2L) round(100 * ev[2] / sum(ev), 1) else 0
  if (!is.null(meta) && !is.null(color_by) && color_by %in% colnames(meta)) {
    id_col   <- if ("SampleID" %in% colnames(meta)) "SampleID" else colnames(meta)[1]
    meta_ids <- as.character(meta[[id_col]])
    idx      <- match(rownames(df), meta_ids)
    df$ColorVar <- as.character(meta[[color_by]])[idx]
  } else {
    df$ColorVar <- rownames(df); color_by <- "Sample"
  }
  ggplot2::ggplot(df, ggplot2::aes(x = PC1, y = PC2, colour = ColorVar)) +
    ggplot2::geom_point(size = 3.5, alpha = 0.8) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", colour = "grey55") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted", colour = "grey55") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom", legend.box = "vertical",
                   legend.margin = ggplot2::margin(),
                   panel.border = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(colour = "black")) +
    ggplot2::xlab(paste0("PC1: ", pc1pct, "%")) +
    ggplot2::ylab(paste0("PC2: ", pc2pct, "%")) +
    ggplot2::ggtitle(title) + ggplot2::labs(colour = color_by)
}

# ── Gene ID type detection ─────────────────────────────────────────────────────

.pp_detect_gene_id_type <- function(gene_ids) {
  s <- utils::head(gene_ids[nchar(gene_ids) > 0], 30)
  if (length(s) == 0L) return("symbol")
  if (mean(grepl("^ENSG[0-9]",              s)) > 0.5) return("ensembl_human")
  if (mean(grepl("^ENSMUSG[0-9]",           s)) > 0.5) return("ensembl_mouse")
  if (mean(grepl("^ENS[A-Z]+G[0-9]",        s)) > 0.5) return("ensembl_other")
  if (mean(grepl("^(NM_|NR_|XM_|XR_)[0-9]",s)) > 0.5) return("refseq")
  if (mean(grepl("^[0-9]+$",                s)) > 0.5) return("entrez")
  return("symbol")
}

# ── Colour warning ─────────────────────────────────────────────────────────────

.pp_color_warn <- function(meta, color_by, n_thresh = 10L) {
  if (is.null(meta) || is.null(color_by) || !color_by %in% colnames(meta)) return(NULL)
  vals   <- meta[[color_by]]
  if (is.numeric(vals)) return(NULL)
  n_lvls <- length(unique(vals[!is.na(vals)]))
  if (n_lvls > n_thresh) {
    return(shiny::div(class = "alert alert-warning", role = "alert",
      style = "font-size:0.82em; padding:6px 10px; margin:6px 0 0 0;",
      shiny::HTML(paste0("⚠️ <b>", color_by, "</b> has <b>", n_lvls,
                          "</b> unique values — colours may be hard to distinguish. ",
                          "Choose a variable with ≤", n_thresh, " levels."))))
  }
  NULL
}

# ── Collapsible done card ─────────────────────────────────────────────────────

.pp_done_card <- function(step_num, step_name, summary_text, back_btn_id, plot_ui = NULL) {
  shiny::tags$details(
    style = "border:1px solid #dee2e6;border-radius:6px;margin-bottom:12px;background:#fff;",
    shiny::tags$summary(
      style = paste0("display:flex;align-items:center;justify-content:space-between;",
                     "padding:10px 14px;cursor:pointer;list-style:none;",
                     "background:#f8f9fa;border-radius:6px;user-select:none;"),
      shiny::tags$div(style = "display:flex;align-items:center;gap:8px;",
        shiny::tags$span(class = "pp-num pp-num-done",
          style = "width:20px;height:20px;font-size:0.72em;flex-shrink:0;", "✓"),
        shiny::tags$strong(paste0("Step ", step_num, " — ", step_name)),
        shiny::tags$span(style = "font-size:0.82em;color:#28a745;margin-left:4px;",
                          summary_text)
      ),
      shiny::tags$span(style = "font-size:0.78em;color:#6c757d;", "▼ expand to review")
    ),
    shiny::tags$div(style = "padding:10px 14px;border-top:1px solid #dee2e6;",
      if (!is.null(plot_ui)) shiny::div(style = "margin-bottom:10px;", plot_ui),
      shiny::actionButton(back_btn_id, paste0("↩ Edit Step ", step_num),
                           class = "btn-sm btn-outline-secondary")
    )
  )
}

# ── Remove suspicious samples section (Steps 1 & 3 only) ──────────────────────

.pp_remove_section_ui <- function(picker_id, btn_id) {
  shiny::tags$details(
    style = "margin-top:12px;border:1px dashed #d1d5db;border-radius:6px;padding:8px 12px;",
    shiny::tags$summary(
      style = "cursor:pointer;font-size:0.84em;color:#6b7280;user-select:none;",
      "⚠️ Remove suspicious samples (optional)"
    ),
    shiny::div(style = "margin-top:8px;",
      shiny::uiOutput(picker_id),
      shiny::actionButton(btn_id, "Remove selected & update data",
                           class = "btn-sm btn-outline-danger", style = "margin-top:6px;")
    )
  )
}

# ── Native listbox picker ──────────────────────────────────────────────────────

.pp_native_picker <- function(input_id, choices) {
  shiny::selectInput(inputId = input_id, label = "Select samples to remove:",
                      choices = choices, selected = NULL, multiple = TRUE,
                      selectize = FALSE, size = min(6L, max(3L, length(choices))),
                      width = "100%")
}

# ── Expand/download button ─────────────────────────────────────────────────────

.pp_expand_btn <- function(btn_id, label = "⤢ Expand / Download") {
  shiny::actionButton(btn_id, label,
    class = "btn-outline-secondary btn-sm",
    style = "font-size:0.78em;padding:2px 8px;margin-top:4px;")
}

# ── Plot modal dialog ──────────────────────────────────────────────────────────

.pp_plot_modal <- function(plot_output_id, dl_btn_id) {
  shiny::modalDialog(
    size = "l", easyClose = TRUE,
    shiny::plotOutput(plot_output_id, height = "550px"),
    footer = shiny::tagList(
      shiny::downloadButton(dl_btn_id, "⬇ Download PNG"),
      shiny::modalButton("Close")
    )
  )
}

# ── Slow-render warning ────────────────────────────────────────────────────────

.pp_slow_warn <- function(text = NULL) {
  msg <- text %||% paste0(
    "⏳ Rendering may take a moment for large datasets. ",
    "Progress is shown in the bar above — results will appear when ready.")
  shiny::div(class = "pp-slow-warn", msg)
}

# =============================================================================
# UI
# =============================================================================
#' @export
preprocessingUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::tags$style(shiny::HTML("
      .pp-num { display:inline-flex;align-items:center;justify-content:center;
        width:24px;height:24px;border-radius:50%;font-weight:700;font-size:0.8em;
        color:#fff;margin-right:6px;flex-shrink:0; }
      .pp-num-active { background:#EBB43E; }
      .pp-num-done   { background:#28a745; }
      .pp-num-locked { background:#9ca3af; }
      .pp-locked { filter:blur(3px);opacity:0.38;pointer-events:none!important;user-select:none; }
      .pp-step-card  { margin-bottom:12px; }
      .pp-hint { font-size:0.82em;color:#6c757d;font-style:italic;margin:4px 0 8px 0; }
      .pp-run-btn    { margin:10px 0 4px 0; }
      .pp-proceed-btn { margin-top:14px; }
      .pp-info-note { font-size:0.8em;color:#EBB43E;margin:2px 0 6px 0; }
      .pp-slow-warn { font-size:0.81em;color:#6b7280;background:#f9fafb;border-radius:4px;
        padding:5px 10px;margin-bottom:8px;border-left:3px solid #EBB43E; }
      details summary::-webkit-details-marker { display:none; }
      details > summary { outline:none; }
      .pp-main-scroll { overflow-y:auto;max-height:calc(100vh - 80px);padding-right:4px; }
      .pp-summary-card { background:#fff;border:2px solid #EBB43E;border-radius:8px;
        padding:16px 18px;margin-top:8px; }
      .pp-summary-row { display:flex;align-items:flex-start;gap:8px;padding:3px 0;
        font-size:0.84em;border-bottom:1px solid #f3f4f6; }
      .pp-summary-row:last-child { border-bottom:none; }
      .pp-summary-label { font-weight:600;color:#374151;min-width:130px;flex-shrink:0; }
      .pp-dl-btn { width:100%;margin-top:4px;font-size:0.78em; }
      .pp-geneid-banner { background:#eef2ff;border:1px solid #c7d2fe;border-radius:8px;
        padding:14px 16px;margin-bottom:12px; }
    ")),

    bslib::layout_sidebar(fill = FALSE,
      sidebar = bslib::sidebar(width = 220, open = TRUE,
        shiny::tags$h6("Progress", style = "font-weight:700;margin-bottom:8px;"),
        shiny::uiOutput(ns("sidebar_progress")),
        shiny::hr(style = "margin:10px 0;"),
        shiny::uiOutput(ns("sidebar_data_info")),
        shiny::hr(style = "margin:10px 0;"),
        shiny::actionButton(ns("reset_preprocessing"), "↺ Reset all",
          class = "btn-outline-danger btn-sm w-100", style = "margin-bottom:4px;",
          title = "Reset all preprocessing steps and restore original data"),
        shiny::actionButton(ns("skip_preprocessing"), "Skip (use data as-is)",
          class = "btn-outline-secondary btn-sm w-100",
          title = "Bypass preprocessing — imported data flows to all analyses"),
        shiny::uiOutput(ns("sidebar_downloads"))
      ),
      shiny::div(class = "pp-main-scroll", style = "padding:10px;",
        shiny::uiOutput(ns("gene_id_banner")),
        shiny::uiOutput(ns("step1_card")),
        shiny::uiOutput(ns("step2_card")),
        shiny::uiOutput(ns("step3_card")),
        shiny::uiOutput(ns("step4_card")),
        shiny::uiOutput(ns("pp_summary_card"))
      )
    )
  )
}

# =============================================================================
# SERVER
# =============================================================================
#' @export
preprocessingServer <- function(id, get_expr, get_meta, log_fn = NULL, get_log = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    ns   <- session$ns
    glog <- function(msg) if (!is.null(log_fn)) log_fn(msg)

    # ── Modal plot reactiveVal (shared across all plots) ───────────────────────
    modal_plot <- shiny::reactiveVal(NULL)

    # ── State ──────────────────────────────────────────────────────────────────
    rv <- shiny::reactiveValues(
      step              = 1L, active_step = 1L,
      data_s0           = NULL, data_s1 = NULL, data_s2 = NULL,
      norm_linear       = NULL, data_s3 = NULL,
      meta_s            = NULL, final_data = NULL, finalized = 0L,
      orig_n_genes      = NULL, orig_n_samples = NULL,
      s1_summary        = "", s2_summary = "", s3_summary = "",
      s1_run            = FALSE, s2_run = FALSE, s3_run = FALSE,
      s4_run            = FALSE, s4_bc_run = FALSE,
      # plot storage
      s1_plot_a         = NULL, s1_plot_b = NULL,
      s3_plot_before    = NULL, s3_plot_after = NULL,
      # PCA: df + eigenvalues stored separately so color change doesn't recompute
      s4_pca_df_before  = NULL, s4_pca_ev_before = NULL,
      s4_pca_df_after   = NULL, s4_pca_ev_after  = NULL,
      s4_bc_result      = NULL,
      s2_thr            = 1.0,
      # gene ID detection
      gene_id_type      = "symbol",
      gene_id_dismissed = FALSE,
      # downloads
      downloads_available = FALSE,
      pp_summary        = NULL
    )

    shiny::observe({
      e <- get_expr(); m <- get_meta(); shiny::req(e, m)
      if (shiny::isolate(rv$active_step) == 1L &&
          is.null(shiny::isolate(rv$data_s0))) {
        d0 <- as.data.frame(e)
        rv$data_s0 <- d0; rv$data_s1 <- d0; rv$meta_s <- m
        rv$orig_n_genes <- nrow(d0); rv$orig_n_samples <- ncol(d0)
        rv$gene_id_type <- .pp_detect_gene_id_type(rownames(d0))
        rv$gene_id_dismissed <- FALSE
      }
    })

    # =========================================================================
    # GENE ID DETECTION BANNER
    # =========================================================================

    output$gene_id_banner <- shiny::renderUI({
      shiny::req(rv$data_s0)
      id_type <- rv$gene_id_type
      if (id_type == "symbol" || rv$gene_id_dismissed) return(NULL)

      type_label <- switch(id_type,
        ensembl_human  = "Ensembl (human) — ENSG...",
        ensembl_mouse  = "Ensembl (mouse) — ENSMUSG...",
        ensembl_other  = "Ensembl (non-human/mouse)",
        refseq         = "RefSeq — NM_, NR_, ...",
        entrez         = "Entrez Gene IDs",
        "unknown format"
      )
      default_species <- if (id_type == "ensembl_mouse") "mouse" else "human"

      shiny::div(class = "pp-geneid-banner",
        shiny::tags$h6(style = "margin-bottom:6px;font-weight:700;color:#3730a3;",
          "\U0001f9ec Gene ID format detected: ", type_label),
        shiny::tags$p(style = "font-size:0.86em;margin-bottom:8px;",
          "Your gene identifiers appear to be in ",
          shiny::tags$b(type_label), " format. ",
          "Would you like to convert them to standard gene symbols (e.g. BRCA1, TP53)?",
          shiny::tags$br(),
          shiny::tags$span(style = "color:#6b7280;font-size:0.9em;",
            "Version numbers (e.g. .2 in ENSG00000141736.2) will be stripped automatically. ",
            "Requires org.Hs.eg.db / org.Mm.eg.db from Bioconductor.")),
        if (id_type %in% c("ensembl_other", "refseq", "entrez"))
          shiny::selectInput(ns("gene_id_species"), "Species:",
            choices  = c("Human (Homo sapiens)" = "human",
                          "Mouse (Mus musculus)" = "mouse"),
            selected = default_species, width = "280px"),
        shiny::div(style = "display:flex;gap:8px;margin-top:8px;",
          shiny::actionButton(ns("convert_gene_ids"), "✓ Convert to gene symbols",
                               class = "btn-primary btn-sm"),
          shiny::actionButton(ns("skip_gene_id"), "✗ Keep original IDs",
                               class = "btn-outline-secondary btn-sm")
        )
      )
    })

    shiny::observeEvent(input$skip_gene_id, {
      rv$gene_id_dismissed <- TRUE
      glog("Gene ID conversion skipped — keeping original identifiers.")
    })

    shiny::observeEvent(input$convert_gene_ids, {
      shiny::req(rv$data_s0)
      id_type <- shiny::isolate(rv$gene_id_type)
      species <- shiny::isolate(input$gene_id_species %||%
        if (id_type == "ensembl_mouse") "mouse" else "human")

      shiny::withProgress(message = "\U0001f9ec Converting gene IDs to symbols…", value = 0, {
        shiny::incProgress(0.1, "Loading annotation database…")

        pkg <- if (species == "mouse") "org.Mm.eg.db" else "org.Hs.eg.db"
        if (!requireNamespace(pkg, quietly = TRUE)) {
          shiny::showNotification(
            paste0("Package not installed. Run: BiocManager::install('", pkg, "')"),
            type = "error", duration = 12)
          return(NULL)
        }
        if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
          shiny::showNotification(
            "AnnotationDbi not installed. Run: BiocManager::install('AnnotationDbi')",
            type = "error", duration = 12)
          return(NULL)
        }

        gene_ids  <- rownames(rv$data_s0)
        ids_clean <- sub("\\.[0-9]+$", "", gene_ids)

        keytype <- switch(id_type,
          ensembl_human  = "ENSEMBL",
          ensembl_mouse  = "ENSEMBL",
          ensembl_other  = "ENSEMBL",
          refseq         = "REFSEQ",
          entrez         = "ENTREZID",
          "SYMBOL"
        )

        shiny::incProgress(0.3,
          paste0("Mapping ", length(ids_clean), " IDs to symbols using ", pkg, "…"))

        org_db <- getExportedValue(pkg, pkg)
        symbols <- tryCatch(
          suppressMessages(
            AnnotationDbi::mapIds(org_db, keys = ids_clean, column = "SYMBOL",
                                   keytype = keytype, multiVals = "first")
          ),
          error = function(e) {
            shiny::showNotification(paste("Conversion error:", e$message),
                                     type = "error", duration = 12)
            NULL
          }
        )
        if (is.null(symbols)) return(NULL)

        new_ids <- ifelse(is.na(symbols) | symbols == "",
                          ids_clean, as.character(symbols))

        # Deduplicate
        dups <- duplicated(new_ids)
        if (any(dups)) {
          dup_vals <- unique(new_ids[dups])
          for (dv in dup_vals) {
            idx <- which(new_ids == dv)
            new_ids[idx] <- paste0(dv, "_dup", seq_along(idx))
          }
        }

        shiny::incProgress(0.5, "Updating gene identifiers…")
        n_mapped <- sum(new_ids != ids_clean)

        rownames(rv$data_s0) <- new_ids
        rownames(rv$data_s1) <- new_ids
        rv$gene_id_dismissed <- TRUE

        glog(paste0("Gene ID Conversion: ", n_mapped, "/", length(gene_ids),
                     " IDs mapped from ", id_type, " to gene symbols (", species, ")."))
        shiny::showNotification(
          shiny::HTML(paste0(
            "✅ ", format(n_mapped, big.mark = ","), " IDs converted to gene symbols. ",
            "Unmapped IDs retain their (version-stripped) original identifier.")),
          type = "default", duration = 8)
        shiny::incProgress(0.1, "Done.")
      })
    })

    # =========================================================================
    # MODAL PLOT — shared output + download
    # =========================================================================

    output$modal_plot_render <- shiny::renderPlot({
      shiny::req(modal_plot())
      modal_plot()
    })

    output$dl_plot_png <- shiny::downloadHandler(
      filename = function() paste0("markeR_plot_", Sys.Date(), ".png"),
      content  = function(file) {
        shiny::req(modal_plot())
        ggplot2::ggsave(file, modal_plot(), width = 12, height = 8, dpi = 150, bg = "white")
      }
    )

    # ── Expand handlers ────────────────────────────────────────────────────────
    shiny::observeEvent(input$expand_s1a, {
      shiny::req(rv$s1_plot_a)
      modal_plot(rv$s1_plot_a$plot)
      shiny::showModal(.pp_plot_modal(ns("modal_plot_render"), ns("dl_plot_png")))
    })
    shiny::observeEvent(input$expand_s1b, {
      shiny::req(rv$s1_plot_b)
      modal_plot(rv$s1_plot_b$plot)
      shiny::showModal(.pp_plot_modal(ns("modal_plot_render"), ns("dl_plot_png")))
    })
    shiny::observeEvent(input$expand_s3_before, {
      shiny::req(rv$s3_plot_before)
      modal_plot(rv$s3_plot_before)
      shiny::showModal(.pp_plot_modal(ns("modal_plot_render"), ns("dl_plot_png")))
    })
    shiny::observeEvent(input$expand_s3_after, {
      shiny::req(rv$s3_plot_after)
      modal_plot(rv$s3_plot_after)
      shiny::showModal(.pp_plot_modal(ns("modal_plot_render"), ns("dl_plot_png")))
    })
    shiny::observeEvent(input$expand_pca_before, {
      shiny::req(rv$s4_pca_df_before, rv$s4_pca_ev_before)
      col <- shiny::isolate(input$s4_pca_color)
      plt <- .pp_pca_from_df(rv$s4_pca_df_before, rv$s4_pca_ev_before,
                              rv$meta_s, col, "PCA — before batch correction")
      modal_plot(plt)
      shiny::showModal(.pp_plot_modal(ns("modal_plot_render"), ns("dl_plot_png")))
    })
    shiny::observeEvent(input$expand_pca_after, {
      shiny::req(rv$s4_pca_df_after, rv$s4_pca_ev_after)
      col <- shiny::isolate(input$s4_pca_after_color %||% input$s4_pca_color)
      plt <- .pp_pca_from_df(rv$s4_pca_df_after, rv$s4_pca_ev_after,
                              rv$meta_s, col, "PCA — after batch correction")
      modal_plot(plt)
      shiny::showModal(.pp_plot_modal(ns("modal_plot_render"), ns("dl_plot_png")))
    })

    # =========================================================================
    # SIDEBAR
    # =========================================================================

    step_info <- list(
      list(n = 1L, label = "Sample QC",
           tip = paste0("Explore library size and read distribution across samples. ",
                        "Helps identify poorly sequenced samples or outliers. ",
                        "This step is OPTIONAL — you can skip straight to filtering.")),
      list(n = 2L, label = "Gene Filtering",
           tip = paste0("Remove lowly or non-expressed genes. ",
                        "Reduces noise and improves statistical power in downstream analyses. ",
                        "The threshold is the minimum mean log₂(count+1) a gene must have.")),
      list(n = 3L, label = "Normalisation",
           tip = paste0("Correct for differences in library size and composition between samples. ",
                        "TMM (edgeR) normalises all samples; the boxplot subset is for display only.")),
      list(n = 4L, label = "PCA & Batch",
           tip = paste0("Visualise major sources of variation with PCA. ",
                        "Optionally apply batch correction (voom/lmFit) to remove technical variation."))
    )

    output$sidebar_progress <- shiny::renderUI({
      step <- rv$step
      shiny::tagList(lapply(step_info, function(x) {
        if (x$n < step)       { cls <- "pp-num pp-num-done";   txt <- "✓" }
        else if (x$n == step) { cls <- "pp-num pp-num-active"; txt <- as.character(x$n) }
        else                  { cls <- "pp-num pp-num-locked";  txt <- as.character(x$n) }
        shiny::tags$div(
          style = "display:flex;align-items:center;justify-content:space-between;padding:4px 0;font-size:0.82em;",
          shiny::tags$div(style = "display:flex;align-items:center;",
            shiny::tags$span(class = cls,
              style = "width:18px;height:18px;font-size:0.7em;flex-shrink:0;", txt),
            shiny::tags$span(x$label,
              style = if (x$n == step) "font-weight:600;" else "")
          ),
          shiny::tags$abbr(
            title = x$tip,
            style = "cursor:help;font-size:0.82em;color:#9ca3af;text-decoration:none;border-bottom:1px dotted #9ca3af;",
            "ⓘ"
          )
        )
      }))
    })

    output$sidebar_data_info <- shiny::renderUI({
      e_curr <- if (!is.null(rv$data_s2)) rv$data_s2
                else if (!is.null(rv$data_s1)) rv$data_s1
                else get_expr()
      if (is.null(e_curr)) return(shiny::helpText("No data loaded.", style = "font-size:0.82em;"))
      cg <- nrow(e_curr); cs <- ncol(e_curr)
      og <- rv$orig_n_genes; os <- rv$orig_n_samples
      filtered <- !is.null(og) && (cg != og || cs != os)
      shiny::tags$div(style = "font-size:0.82em;color:#555;",
        shiny::tags$b("Current data:"), shiny::tags$br(),
        paste0(format(cg, big.mark = ","), " genes"), shiny::tags$br(),
        paste0(format(cs, big.mark = ","), " samples"),
        if (filtered) shiny::tags$div(
          style = "color:#9ca3af;font-size:0.82em;margin-top:3px;",
          paste0("(original: ", format(og, big.mark = ","), "×", format(os, big.mark = ","), ")")
        )
      )
    })

    output$sidebar_downloads <- shiny::renderUI({
      if (!rv$downloads_available) return(NULL)
      shiny::tagList(
        shiny::hr(style = "margin:10px 0;"),
        shiny::tags$div(style = "font-size:0.8em;font-weight:600;color:#555;margin-bottom:2px;",
                         "\U0001f4be Download results"),
        shiny::tags$p(style = "font-size:0.75em;color:#9ca3af;margin:0 0 4px 0;",
          "Optional — data is already available to all analyses."),
        shiny::downloadButton(ns("dl_expr"), "Expression data (.csv)",
                               class = "btn-outline-primary btn-sm pp-dl-btn"),
        shiny::downloadButton(ns("dl_meta"), "Metadata (.csv)",
                               class = "btn-outline-secondary btn-sm pp-dl-btn"),
        if (!is.null(get_log))
          shiny::downloadButton(ns("dl_log"), "Processing log (.txt)",
                                 class = "btn-outline-secondary btn-sm pp-dl-btn")
      )
    })

    output$dl_expr <- shiny::downloadHandler(
      filename = function() paste0("markeR_expr_preprocessed_", Sys.Date(), ".csv"),
      content  = function(file) {
        shiny::req(rv$final_data)
        utils::write.csv(rv$final_data, file)
      }
    )
    output$dl_meta <- shiny::downloadHandler(
      filename = function() paste0("markeR_metadata_preprocessed_", Sys.Date(), ".csv"),
      content  = function(file) {
        shiny::req(rv$meta_s)
        utils::write.csv(rv$meta_s, file, row.names = FALSE)
      }
    )
    output$dl_log <- shiny::downloadHandler(
      filename = function() paste0("markeR_log_", Sys.Date(), ".txt"),
      content  = function(file) {
        log <- if (!is.null(get_log)) get_log() else character(0)
        writeLines(log, file)
      }
    )

    shiny::observeEvent(input$skip_preprocessing, {
      d <- rv$data_s1
      if (is.null(d)) d <- tryCatch(as.data.frame(get_expr()), error = function(e) NULL)
      if (!is.null(d)) rv$final_data <- d
      rv$downloads_available <- TRUE
      glog(paste0("Preprocessing skipped — using imported data as-is (",
                   if (!is.null(d)) paste0(nrow(d), " genes × ", ncol(d), " samples")
                   else "data dimensions unknown", ")."))
      shiny::showNotification("Preprocessing skipped — using imported data as-is.",
                               type = "message", duration = 5)
    })

    shiny::observeEvent(input$reset_preprocessing, {
      rv$step <- 1L; rv$active_step <- 1L
      rv$data_s1 <- rv$data_s0; rv$data_s2 <- NULL
      rv$norm_linear <- NULL; rv$data_s3 <- NULL; rv$final_data <- NULL
      rv$meta_s <- get_meta()
      rv$s1_run <- FALSE; rv$s2_run <- FALSE; rv$s3_run <- FALSE
      rv$s4_run <- FALSE; rv$s4_bc_run <- FALSE; rv$s2_thr <- 1.0
      rv$s1_plot_a <- NULL; rv$s1_plot_b <- NULL
      rv$s3_plot_before <- NULL; rv$s3_plot_after <- NULL
      rv$s4_pca_df_before <- NULL; rv$s4_pca_ev_before <- NULL
      rv$s4_pca_df_after  <- NULL; rv$s4_pca_ev_after  <- NULL
      rv$s4_bc_result <- NULL
      rv$pp_summary <- NULL; rv$finalized <- 0L
      rv$downloads_available <- FALSE
      rv$gene_id_dismissed <- FALSE
      rv$gene_id_type <- .pp_detect_gene_id_type(rownames(rv$data_s0 %||% data.frame()))
      shiny::showNotification("Preprocessing reset to original data.", type = "message", duration = 4)
    })

    # =========================================================================
    # STEP 1 — Sample Complexity QC (optional, visualisation only)
    # =========================================================================

    .s1_plot_samples_calc <- shiny::reactive({
      e <- rv$data_s0; shiny::req(e)
      all_s  <- colnames(e)
      method <- input$s1_subset_method %||% "all"
      switch(method,
        all      = all_s,
        specific = { sel <- input$s1_specific_samples; if (length(sel) > 0L) intersect(sel, all_s) else all_s },
        random   = sample(all_s, min(as.integer(input$s1_n_samp %||% 10L), length(all_s))),
        topbottom = {
          n <- as.integer(input$s1_n_samp %||% 10L); tots <- log2(colSums(e) + 1); ord <- order(tots)
          nb <- floor(n / 2L); nt <- n - nb
          unique(c(all_s[utils::head(ord, nb)], all_s[utils::tail(ord, nt)]))
        },
        all_s
      )
    })

    output$step1_card <- shiny::renderUI({
      e <- rv$data_s0
      if (is.null(e)) {
        return(bslib::card(class = "pp-step-card",
          bslib::card_header(shiny::tags$div(style = "display:flex;align-items:center;",
            shiny::tags$span(class = "pp-num pp-num-active", "1"),
            shiny::tags$strong("Step 1 — Sample Complexity QC (Optional)"))),
          shiny::div(style = "padding:14px;",
            shiny::helpText("Load expression data first (Import tab)."))))
      }
      if (rv$step > 1L && rv$active_step != 1L) {
        return(.pp_done_card(1, "Sample Complexity QC", rv$s1_summary, ns("back_to_s1"),
          plot_ui = shiny::plotOutput(ns("plot_s1a"), height = "180px")))
      }
      controls <- shiny::fluidRow(
        shiny::column(6,
          shiny::tags$b("Plot type"),
          shiny::radioButtons(ns("s1_plot_type"), label = NULL,
            choices = c("Cumulative reads by gene index" = "byexpr",
                        "Unique genes vs reads sampled (slow)" = "sampling"),
            selected = "byexpr"),
          shiny::checkboxInput(ns("s1_pct"), "Show reads as % of total", FALSE),
          shiny::conditionalPanel(
            condition = sprintf("input['%s'] == 'sampling'", ns("s1_plot_type")),
            shiny::fluidRow(
              shiny::column(6, shiny::numericInput(ns("s1_n_genes"), "N genes (≤1000)", 500, min = 10, max = 1000, step = 50)),
              shiny::column(6, shiny::numericInput(ns("s1_n_iter"), "Iterations", 5, min = 1, max = 20))
            ),
            shiny::numericInput(ns("s1_step_reads"), "Read step size", 500, min = 10, max = 10000, step = 100),
            shiny::helpText("⚠️ Sampling is slow — use ≤1000 genes & ≤5 samples.",
                             style = "color:#856404;font-size:0.8em;")
          )
        ),
        shiny::column(6,
          shiny::tags$b("Samples to visualise"),
          shiny::div(class = "pp-info-note",
            "⚠️ This selection is for the plot only — it does not filter your data."),
          shiny::radioButtons(ns("s1_subset_method"), label = NULL,
            choices = c("All samples" = "all", "Specific samples" = "specific",
                        "Random N" = "random", "Top/Bottom N by lib-size" = "topbottom"),
            selected = "all"),
          shiny::conditionalPanel(
            condition = sprintf("input['%s'] == 'specific'", ns("s1_subset_method")),
            shiny::uiOutput(ns("s1_specific_picker"))),
          shiny::conditionalPanel(
            condition = sprintf("input['%s'] != 'all' && input['%s'] != 'specific'",
                                ns("s1_subset_method"), ns("s1_subset_method")),
            shiny::numericInput(ns("s1_n_samp"), "N samples", 10, min = 2, step = 1))
        )
      )
      bslib::card(class = "pp-step-card",
        bslib::card_header(shiny::tags$div(style = "display:flex;align-items:center;",
          shiny::tags$span(class = "pp-num pp-num-active", "1"),
          shiny::tags$strong("Step 1 — Sample Complexity QC"),
          shiny::tags$span(style = "margin-left:8px;font-size:0.78em;color:#6c757d;font-style:italic;",
                            "(optional — you can skip this step)"))),
        shiny::div(style = "padding:14px;",
          controls,
          shiny::hr(style = "margin:8px 0;"),
          shiny::uiOutput(ns("s1_run_section"))))
    })

    output$s1_run_section <- shiny::renderUI({
      run <- rv$s1_run
      n_s <- if (!is.null(rv$data_s0)) ncol(rv$data_s0) else 0L
      if (!run)
        return(shiny::div(class = "pp-run-btn",
          .pp_slow_warn(paste0(
            "⏳ Generating this plot processes each sample individually. ",
            "For datasets with many samples (>50), this may take a minute.")),
          shiny::actionButton(ns("s1_run_btn"), "\U0001f4ca Show Complexity Plot",
                               class = "btn-primary btn-sm", style = "margin-right:6px;"),
          shiny::actionButton(ns("skip_s1"), "Skip Step 1 →",
                               class = "btn-outline-secondary btn-sm")
        ))
      shiny::tagList(
        shiny::div(style = "margin-bottom:6px;",
          .pp_slow_warn(paste0(
            "⏳ Re-running will reprocess ", n_s, " sample(s). ",
            "Results will appear once computation completes.")),
          shiny::actionButton(ns("s1_run_btn"), "\U0001f504 Re-run Plot",
                               class = "btn-sm btn-outline-primary", style = "margin-right:6px;"),
          shiny::actionButton(ns("skip_s1"), "Skip & proceed →",
                               class = "btn-sm btn-outline-secondary")
        ),
        shiny::fluidRow(
          shiny::column(6,
            if (!is.null(rv$s1_plot_a))
              shiny::tagList(
                shiny::plotOutput(ns("plot_s1a"), height = "260px"),
                .pp_expand_btn(ns("expand_s1a"))
              )
            else
              shiny::div(style = "height:260px;display:flex;align-items:center;justify-content:center;color:#9ca3af;font-size:0.85em;",
                "Click 'Show Complexity Plot' above.")
          ),
          shiny::column(6,
            shiny::conditionalPanel(
              condition = sprintf("input['%s'] == 'sampling'", ns("s1_plot_type")),
              shiny::tagList(
                shiny::actionButton(ns("run_sampling"), "▶ Run sampling plot",
                                     class = "btn-sm btn-outline-secondary", style = "margin-bottom:6px;"),
                shiny::plotOutput(ns("plot_s1b"), height = "260px"),
                if (!is.null(rv$s1_plot_b)) .pp_expand_btn(ns("expand_s1b"))
              )
            ),
            shiny::conditionalPanel(
              condition = sprintf("input['%s'] != 'sampling'", ns("s1_plot_type")),
              shiny::div(
                style = "display:flex;align-items:center;justify-content:center;height:260px;color:#9ca3af;font-size:0.88em;text-align:center;",
                shiny::HTML("Select <em>Unique genes vs reads sampled</em><br>above for a second plot."))
            )
          )
        ),
        .pp_remove_section_ui(ns("s1_remove_picker"), ns("s1_remove_btn")),
        shiny::div(class = "pp-proceed-btn",
          shiny::actionButton(ns("proceed_s1"), "Proceed to Gene Filtering →",
                               class = "btn-warning btn-sm"))
      )
    })

    output$s1_specific_picker <- shiny::renderUI({
      e <- rv$data_s0; shiny::req(e)
      shiny::selectInput(ns("s1_specific_samples"), label = "Samples to show:",
                          choices = colnames(e), selected = colnames(e),
                          multiple = TRUE, selectize = TRUE, width = "100%")
    })

    output$s1_remove_picker <- shiny::renderUI({
      e <- rv$data_s1; shiny::req(e, rv$s1_run)
      .pp_native_picker(ns("s1_remove_samples"), colnames(e))
    })

    output$plot_s1a <- shiny::renderPlot({
      shiny::req(rv$s1_run, !is.null(rv$s1_plot_a))
      rv$s1_plot_a$plot
    })
    output$plot_s1b <- shiny::renderPlot({
      shiny::req(rv$s1_run)
      res <- rv$s1_plot_b
      shiny::validate(shiny::need(!is.null(res), "Click '▶ Run sampling plot' to generate."))
      res$plot
    })

    shiny::observeEvent(input$s1_run_btn, {
      shiny::req(rv$data_s0)
      samps <- shiny::isolate(.s1_plot_samples_calc())
      pct   <- shiny::isolate(isTRUE(input$s1_pct))
      n_s   <- length(samps)
      shiny::withProgress(
        message = paste0("\U0001f4ca Computing complexity plot — ", n_s, " sample(s)…"),
        value = 0, {
          res <- tryCatch(
            .pp_complexity_byexprgenes(rv$data_s0, samps, pct,
              progress_fn = function(frac, msg) shiny::incProgress(frac * 0.9, msg)),
            error = function(e) { shiny::showNotification(e$message, type = "error"); NULL }
          )
          shiny::incProgress(0.1, "Rendering plot…")
          rv$s1_plot_a <- res
          rv$s1_run    <- TRUE
        }
      )
    })

    shiny::observeEvent(input$run_sampling, {
      shiny::req(rv$data_s0)
      samps <- shiny::isolate(.s1_plot_samples_calc())
      ng    <- min(as.integer(shiny::isolate(input$s1_n_genes %||% 500L)), 1000L)
      ni    <- as.integer(shiny::isolate(input$s1_n_iter     %||% 5L))
      sr    <- as.integer(shiny::isolate(input$s1_step_reads %||% 500L))
      pct   <- shiny::isolate(isTRUE(input$s1_pct))
      shiny::withProgress(
        message = paste0("⏳ Sampling reads — ", length(samps), " sample(s), up to ", ng, " genes each…"),
        value = 0, {
          res <- tryCatch(
            .pp_complexity_samplingreads(rv$data_s0, ng, ni, sr, samps, pct,
              progress_fn = function(frac, msg) shiny::incProgress(frac * 0.9, msg)),
            error = function(e) {
              shiny::showNotification(e$message, type = "error", duration = 10); NULL
            }
          )
          shiny::incProgress(0.1, "Rendering…")
          rv$s1_plot_b <- res
        }
      )
    })

    shiny::observeEvent(input$s1_remove_btn, {
      to_rm <- input$s1_remove_samples; shiny::req(length(to_rm) >= 1L, rv$data_s1)
      keep <- setdiff(colnames(rv$data_s1), to_rm)
      if (length(keep) < 2L) {
        shiny::showNotification("Cannot remove: fewer than 2 samples would remain.", type = "error")
        return()
      }
      rv$data_s1 <- rv$data_s1[, keep, drop = FALSE]; .align_meta(keep)
      glog(paste0("Sample QC: removed ", length(to_rm), " sample(s): [",
                   paste(to_rm, collapse = ", "), "]. ", length(keep), " samples remain."))
      shiny::showNotification(
        paste0(length(to_rm), " sample(s) removed. Re-run to refresh the plot."),
        type = "warning", duration = 5)
    })

    shiny::observeEvent(input$skip_s1, { .go_to_s2() })
    shiny::observeEvent(input$proceed_s1, { .go_to_s2() })

    .go_to_s2 <- function() {
      shiny::req(rv$data_s1)
      n_s     <- ncol(rv$data_s1)
      removed <- (rv$orig_n_samples %||% n_s) - n_s
      rv$s1_summary <- paste0(n_s, " samples",
        if (removed > 0L) paste0(" (", removed, " removed)") else "")
      glog(paste0("Sample Complexity QC: ", n_s, " sample(s) passed to next step",
                   if (removed > 0L) paste0(" (", removed, " removed during QC)") else "",
                   "."))
      rv$step <- 2L; rv$active_step <- 2L; rv$s2_run <- TRUE; rv$s2_thr <- 1.0
    }

    shiny::observeEvent(input$back_to_s1, {
      rv$active_step <- 1L; rv$step <- 1L; rv$s1_run <- FALSE
      rv$data_s1 <- rv$data_s0; rv$meta_s <- get_meta()
      rv$data_s2 <- NULL; rv$norm_linear <- NULL; rv$data_s3 <- NULL
    })

    # =========================================================================
    # STEP 2 — Gene Filtering
    # =========================================================================

    shiny::observeEvent(input$s2_plot_click, {
      shiny::req(input$s2_plot_click)
      val <- round(as.numeric(input$s2_plot_click$x), 3)
      rv$s2_thr <- val; shiny::updateNumericInput(session, "s2_threshold_input", value = val)
    })
    shiny::observeEvent(input$s2_threshold_input, {
      val <- suppressWarnings(as.numeric(input$s2_threshold_input))
      if (!is.na(val)) rv$s2_thr <- val
    }, ignoreInit = TRUE)

    output$step2_card <- shiny::renderUI({
      if (rv$step < 2L) {
        return(shiny::div(class = "pp-locked",
          bslib::card(class = "pp-step-card",
            bslib::card_header(shiny::tags$div(style = "display:flex;align-items:center;",
              shiny::tags$span(class = "pp-num pp-num-locked", "2"),
              shiny::tags$strong("Step 2 — Filter Lowly Expressed Genes"))),
            shiny::div(style = "padding:14px;color:#9ca3af;",
              "Complete Step 1 (or skip it) to unlock."))))
      }
      if (rv$step > 2L && rv$active_step != 2L) {
        return(.pp_done_card(2, "Filter Lowly Expressed Genes", rv$s2_summary, ns("back_to_s2"),
          plot_ui = shiny::plotOutput(ns("plot_s2"), height = "180px")))
      }
      bslib::card(class = "pp-step-card",
        bslib::card_header(shiny::tags$div(style = "display:flex;align-items:center;",
          shiny::tags$span(class = "pp-num pp-num-active", "2"),
          shiny::tags$strong("Step 2 — Filter Lowly Expressed Genes"))),
        shiny::div(style = "padding:14px;",
          .pp_slow_warn(paste0(
            "⏳ Plots update as you move the threshold — ",
            "for very large gene sets, allow a moment for re-rendering.")),
          shiny::tags$p(class = "pp-hint",
            "\U0001f4a1 Click on the left plot to move the threshold, or type a value below."),
          shiny::fluidRow(
            shiny::column(4,
              shiny::numericInput(ns("s2_threshold_input"), "Threshold:",
                                   value = shiny::isolate(rv$s2_thr),
                                   min = 0, step = 0.1, width = "150px"),
              shiny::uiOutput(ns("s2_gene_count")),
              shiny::div(style = "margin-top:10px;",
                .pp_expand_btn(ns("expand_s2_all"), "⤢ Expand left"),
                .pp_expand_btn(ns("expand_s2_after"), "⤢ Expand right")
              )
            ),
            shiny::column(8,
              shiny::fluidRow(
                shiny::column(6,
                  shiny::tags$b(style = "font-size:0.85em;", "All genes"),
                  shiny::plotOutput(ns("plot_s2"), click = ns("s2_plot_click"), height = "250px")
                ),
                shiny::column(6,
                  shiny::tags$b(style = "font-size:0.85em;", "Genes passing filter"),
                  shiny::plotOutput(ns("plot_s2_after"), height = "250px")
                )
              )
            )
          ),
          shiny::div(class = "pp-proceed-btn",
            shiny::actionButton(ns("proceed_s2"),
              "Apply filter & proceed to Normalisation →",
              class = "btn-warning btn-sm"))
        )
      )
    })

    output$s2_gene_count <- shiny::renderUI({
      shiny::req(rv$step >= 2L, rv$data_s1)
      d <- rv$data_s1; thr <- rv$s2_thr; n_total <- nrow(d)
      n_keep <- sum(log2(rowMeans(as.matrix(d)) + 1) >= thr)
      shiny::tags$div(style = "font-size:0.83em;color:#555;margin-top:8px;line-height:1.7;",
        shiny::HTML(paste0(
          "<b>", format(n_keep, big.mark = ","), "</b> genes kept<br>",
          "<b>", format(n_total - n_keep, big.mark = ","), "</b> removed (",
          round(100 * (n_total - n_keep) / n_total, 1), "%)")))
    })

    output$plot_s2 <- shiny::renderPlot({
      shiny::req(rv$s2_run, rv$step >= 2L, rv$data_s1)
      d <- rv$data_s1; thr <- rv$s2_thr
      df <- data.frame(log_mean = log2(rowMeans(as.matrix(d)) + 1))
      ggplot2::ggplot(df, ggplot2::aes(x = log_mean)) +
        ggplot2::geom_density(fill = "#EBB43E", alpha = 0.45, colour = "#b8860b", linewidth = 0.8) +
        ggplot2::geom_vline(xintercept = as.numeric(thr), linetype = "dashed",
                             colour = "firebrick", linewidth = 1.2) +
        ggplot2::annotate("text", x = as.numeric(thr), y = Inf, vjust = 1.8, hjust = -0.1,
                           label = paste0("thr = ", round(thr, 3)), colour = "firebrick", size = 3.5) +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                       panel.border = ggplot2::element_blank(),
                       axis.line = ggplot2::element_line(colour = "black")) +
        ggplot2::xlab("Mean log₂(count + 1)") + ggplot2::ylab("Density") +
        ggplot2::ggtitle("All genes — click to set threshold")
    })

    output$plot_s2_after <- shiny::renderPlot({
      shiny::req(rv$s2_run, rv$step >= 2L, rv$data_s1)
      d <- rv$data_s1; thr <- rv$s2_thr
      log_means <- log2(rowMeans(as.matrix(d)) + 1)
      keep_vals <- log_means[log_means >= thr]
      if (length(keep_vals) == 0) {
        return(ggplot2::ggplot() +
          ggplot2::annotate("text", x = 0.5, y = 0.5,
                             label = "No genes pass this threshold", size = 4.5, colour = "#9ca3af") +
          ggplot2::theme_void())
      }
      df <- data.frame(log_mean = keep_vals)
      ggplot2::ggplot(df, ggplot2::aes(x = log_mean)) +
        ggplot2::geom_density(fill = "#28a745", alpha = 0.4, colour = "#155724", linewidth = 0.8) +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                       panel.border = ggplot2::element_blank(),
                       axis.line = ggplot2::element_line(colour = "black")) +
        ggplot2::xlab("Mean log₂(count + 1)") + ggplot2::ylab("Density") +
        ggplot2::ggtitle(paste0("Retained genes (thr ≥ ", round(thr, 2), ")"))
    })

    # Step 2 expand handlers (compute inline from current data)
    shiny::observeEvent(input$expand_s2_all, {
      shiny::req(rv$data_s1)
      d <- rv$data_s1; thr <- shiny::isolate(rv$s2_thr)
      df <- data.frame(log_mean = log2(rowMeans(as.matrix(d)) + 1))
      plt <- ggplot2::ggplot(df, ggplot2::aes(x = log_mean)) +
        ggplot2::geom_density(fill = "#EBB43E", alpha = 0.45, colour = "#b8860b", linewidth = 0.8) +
        ggplot2::geom_vline(xintercept = thr, linetype = "dashed", colour = "firebrick", linewidth = 1.2) +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), panel.border = ggplot2::element_blank(),
                       axis.line = ggplot2::element_line(colour = "black")) +
        ggplot2::xlab("Mean log₂(count + 1)") + ggplot2::ylab("Density") +
        ggplot2::ggtitle("All genes — expression density")
      modal_plot(plt)
      shiny::showModal(.pp_plot_modal(ns("modal_plot_render"), ns("dl_plot_png")))
    })
    shiny::observeEvent(input$expand_s2_after, {
      shiny::req(rv$data_s1)
      d <- rv$data_s1; thr <- shiny::isolate(rv$s2_thr)
      lm <- log2(rowMeans(as.matrix(d)) + 1); kv <- lm[lm >= thr]
      if (length(kv) == 0) return(NULL)
      df <- data.frame(log_mean = kv)
      plt <- ggplot2::ggplot(df, ggplot2::aes(x = log_mean)) +
        ggplot2::geom_density(fill = "#28a745", alpha = 0.4, colour = "#155724", linewidth = 0.8) +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), panel.border = ggplot2::element_blank(),
                       axis.line = ggplot2::element_line(colour = "black")) +
        ggplot2::xlab("Mean log₂(count + 1)") + ggplot2::ylab("Density") +
        ggplot2::ggtitle(paste0("Genes passing filter (thr ≥ ", round(thr, 2), ")"))
      modal_plot(plt)
      shiny::showModal(.pp_plot_modal(ns("modal_plot_render"), ns("dl_plot_png")))
    })

    shiny::observeEvent(input$proceed_s2, {
      shiny::req(rv$step >= 2L, rv$data_s1)
      d <- rv$data_s1; thr <- as.numeric(rv$s2_thr)
      keep <- rownames(d)[log2(rowMeans(as.matrix(d)) + 1) >= thr]
      if (length(keep) < 10L) {
        shiny::showNotification("Too few genes retained — lower the threshold.", type = "error")
        return()
      }
      rv$data_s2 <- d[keep, , drop = FALSE]
      n_rm  <- nrow(d) - length(keep)
      rv$s2_summary <- paste0(format(length(keep), big.mark = ","), " genes kept, ",
                                format(n_rm, big.mark = ","), " removed (thr=", thr, ")")
      glog(paste0("Gene Filtering: threshold=", thr, " — ",
                   format(length(keep), big.mark = ","), " genes retained, ",
                   format(n_rm, big.mark = ","), " removed."))
      rv$step <- 3L; rv$active_step <- 3L; rv$s3_run <- FALSE
    })

    shiny::observeEvent(input$back_to_s2, {
      rv$active_step <- 2L; rv$step <- 2L; rv$s2_run <- TRUE
      rv$data_s2 <- NULL; rv$norm_linear <- NULL; rv$data_s3 <- NULL
    })

    # =========================================================================
    # STEP 3 — Normalisation (before + after boxplots)
    # =========================================================================

    output$step3_card <- shiny::renderUI({
      if (rv$step < 3L) {
        return(shiny::div(class = "pp-locked",
          bslib::card(class = "pp-step-card",
            bslib::card_header(shiny::tags$div(style = "display:flex;align-items:center;",
              shiny::tags$span(class = "pp-num pp-num-locked", "3"),
              shiny::tags$strong("Step 3 — Normalisation"))),
            shiny::div(style = "padding:14px;color:#9ca3af;", "Complete Steps 1–2 to unlock."))))
      }
      if (rv$step > 3L && rv$active_step != 3L) {
        return(.pp_done_card(3, "Normalisation", rv$s3_summary, ns("back_to_s3"),
          plot_ui = shiny::plotOutput(ns("plot_s3_after"), height = "220px")))
      }
      meta_cols <- setdiff(colnames(rv$meta_s), "SampleID")
      controls <- shiny::fluidRow(
        shiny::column(4,
          shiny::tags$b("Method"),
          shiny::radioButtons(ns("s3_method"), label = NULL,
            choices = c("TMM (edgeR)" = "TMM", "VST (DESeq2)" = "DESeq2"), selected = "TMM"),
          shiny::conditionalPanel(sprintf("input['%s'] == 'DESeq2'", ns("s3_method")),
            shiny::helpText("Requires: BiocManager::install('DESeq2')",
                             style = "font-size:0.8em;color:#856404;"))
        ),
        shiny::column(4,
          shiny::tags$b("Samples in boxplot"),
          shiny::tags$p(style = "font-size:0.78em;color:#6b7280;margin:0 0 4px 0;",
            "TMM runs on ALL samples — this selection is for display only."),
          shiny::radioButtons(ns("s3_sample_display"), label = NULL,
            choiceValues = c("all", "random", "extreme"),
            choiceNames  = list(
              "All",
              "Random N",
              shiny::tags$span("Extreme N",
                shiny::tags$abbr(
                  style = "cursor:help;font-size:0.82em;color:#9ca3af;text-decoration:none;border-bottom:1px dotted #9ca3af;margin-left:3px;",
                  title = paste0("Shows the N/2 samples with the highest AND N/2 with the lowest ",
                                  "median log-count — useful to spot outliers without plotting all samples."),
                  "ⓘ"))
            ),
            selected = "all"),
          shiny::conditionalPanel(sprintf("input['%s'] != 'all'", ns("s3_sample_display")),
            shiny::numericInput(ns("s3_n_box"), "N", 20, min = 2, step = 1))
        ),
        shiny::column(4,
          shiny::tags$b("Colour samples by"),
          if (length(meta_cols) > 0L)
            shiny::selectInput(ns("s3_color_by"), label = NULL,
                                choices = meta_cols, selected = meta_cols[1])
          else
            shiny::helpText("No metadata columns."),
          shiny::uiOutput(ns("s3_color_warn"))
        )
      )
      bslib::card(class = "pp-step-card",
        bslib::card_header(shiny::tags$div(style = "display:flex;align-items:center;",
          shiny::tags$span(class = "pp-num pp-num-active", "3"),
          shiny::tags$strong("Step 3 — Normalisation & Quality Distribution"))),
        shiny::div(style = "padding:14px;",
          controls,
          shiny::hr(style = "margin:8px 0;"),
          shiny::uiOutput(ns("s3_run_section"))))
    })

    output$s3_color_warn <- shiny::renderUI({ .pp_color_warn(rv$meta_s, input$s3_color_by) })

    output$s3_run_section <- shiny::renderUI({
      n_g <- if (!is.null(rv$data_s2)) nrow(rv$data_s2) else "?"
      n_s <- if (!is.null(rv$data_s2)) ncol(rv$data_s2) else "?"
      if (!rv$s3_run)
        return(shiny::div(class = "pp-run-btn",
          .pp_slow_warn(paste0(
            "⏳ Normalisation runs on all ", n_s, " samples × ", n_g,
            " genes, then builds two boxplots. This may take a minute for large datasets.")),
          shiny::actionButton(ns("s3_run_btn"),
            "\U0001f4ca Run Normalisation & Show Boxplots", class = "btn-primary btn-sm")
        ))
      shiny::tagList(
        shiny::div(style = "margin-bottom:6px;",
          .pp_slow_warn(paste0(
            "⏳ Re-running will re-normalise all ", n_s, " samples and rebuild both boxplots.")),
          shiny::actionButton(ns("s3_run_btn"), "\U0001f504 Re-run",
                               class = "btn-sm btn-outline-primary")),
        shiny::fluidRow(
          shiny::column(6,
            shiny::tags$b(style = "font-size:0.85em;", "Before normalisation"),
            shiny::plotOutput(ns("plot_s3_before"), height = "300px"),
            .pp_expand_btn(ns("expand_s3_before"))
          ),
          shiny::column(6,
            shiny::tags$b(style = "font-size:0.85em;", "After normalisation"),
            shiny::plotOutput(ns("plot_s3_after"), height = "300px"),
            .pp_expand_btn(ns("expand_s3_after"))
          )
        ),
        .pp_remove_section_ui(ns("s3_remove_picker"), ns("s3_remove_btn")),
        shiny::div(class = "pp-proceed-btn",
          shiny::actionButton(ns("proceed_s3"),
            "Confirm normalisation & proceed to PCA →",
            class = "btn-warning btn-sm"))
      )
    })

    shiny::observeEvent(input$s3_run_btn, {
      shiny::req(rv$step >= 3L, rv$data_s2)
      d      <- as.matrix(rv$data_s2)
      method <- shiny::isolate(input$s3_method %||% "TMM")

      shiny::withProgress(
        message = paste0("⚙️ Normalising with ", method, " — ",
                          nrow(d), " genes × ", ncol(d), " samples…"),
        value = 0, {
          shiny::incProgress(0.1, paste0("Running ", method, " on all ", ncol(d), " samples…"))

          res <- if (method == "TMM") {
            dge <- edgeR::DGEList(counts = d); dge <- edgeR::calcNormFactors(dge)
            nl  <- edgeR::cpm(dge, log = FALSE)
            list(norm_linear = nl, log_norm = log2(nl + 1), method = "TMM")
          } else {
            if (!requireNamespace("DESeq2", quietly = TRUE)) {
              shiny::showNotification(
                "DESeq2 not installed: BiocManager::install('DESeq2')",
                type = "error", duration = 10)
              return(NULL)
            }
            int_d <- round(d); storage.mode(int_d) <- "integer"
            col_df <- if (!is.null(rv$meta_s)) rv$meta_s else data.frame(row.names = colnames(int_d))
            dds <- tryCatch(
              DESeq2::DESeqDataSetFromMatrix(int_d, colData = col_df, design = ~1),
              error = function(e) { shiny::showNotification(e$message, type = "error"); NULL })
            if (is.null(dds)) return(NULL)
            vst_obj <- tryCatch(
              DESeq2::vst(dds, blind = TRUE),
              error = function(e) { shiny::showNotification(e$message, type = "error"); NULL })
            if (is.null(vst_obj)) return(NULL)
            lg <- DESeq2::assay(vst_obj)
            list(norm_linear = 2^lg - 1, log_norm = lg, method = "DESeq2/VST")
          }

          if (is.null(res)) return(NULL)

          shiny::incProgress(0.3, "Building boxplots…")

          # Helper: build boxplot from log-count matrix
          .make_boxplot <- function(log_d, all_s, display, n_box, color_by, meta, title_str) {
            samps_plot <- switch(display,
              random  = sample(all_s, min(n_box, length(all_s))),
              extreme = {
                meds <- apply(log_d, 2, stats::median)
                ord  <- order(meds); nb <- floor(n_box / 2L); nt <- n_box - nb
                unique(c(all_s[utils::head(ord, nb)], all_s[utils::tail(ord, nt)]))
              },
              all_s
            )
            sub_d <- log_d[, samps_plot, drop = FALSE]
            df    <- do.call(rbind, lapply(samps_plot, function(s)
              data.frame(Sample = s, LogCount = sub_d[, s], stringsAsFactors = FALSE)))
            if (!is.null(color_by) && !is.null(meta) && color_by %in% colnames(meta)) {
              id_col <- if ("SampleID" %in% colnames(meta)) "SampleID" else colnames(meta)[1]
              cmap   <- stats::setNames(as.character(meta[[color_by]]),
                                        as.character(meta[[id_col]]))
              df$Group <- cmap[df$Sample]
            } else { df$Group <- df$Sample; color_by <- "Sample" }
            med_lvls  <- names(sort(tapply(df$LogCount, df$Sample, stats::median, na.rm = TRUE)))
            df$Sample <- factor(df$Sample, levels = med_lvls)
            ggplot2::ggplot(df, ggplot2::aes(x = Sample, y = LogCount, fill = Group)) +
              ggplot2::geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3, lwd = 0.3) +
              ggplot2::theme_bw() +
              ggplot2::theme(
                axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 6),
                panel.grid.minor = ggplot2::element_blank(),
                panel.border = ggplot2::element_blank(),
                axis.line = ggplot2::element_line(colour = "black"),
                legend.position = "bottom") +
              ggplot2::xlab("") +
              ggplot2::ylab("log₂(count + 1)") +
              ggplot2::ggtitle(paste0(title_str, " — ", length(samps_plot), " samples")) +
              ggplot2::labs(fill = color_by)
          }

          all_s    <- colnames(res$log_norm)
          display  <- shiny::isolate(input$s3_sample_display %||% "all")
          n_box    <- as.integer(shiny::isolate(input$s3_n_box %||% 20L))
          color_by <- shiny::isolate(input$s3_color_by)
          meta     <- rv$meta_s

          # BEFORE: raw log-counts from filtered (un-normalised) data
          log_before <- log2(d + 1)
          rv$s3_plot_before <- .make_boxplot(log_before, all_s, display, n_box,
                                              color_by, meta, "Raw counts")

          # AFTER: TMM/VST normalised
          rv$s3_plot_after  <- .make_boxplot(res$log_norm, all_s, display, n_box,
                                              color_by, meta, res$method)

          rv$norm_linear <- res$norm_linear
          rv$data_s3     <- as.data.frame(res$log_norm)

          shiny::incProgress(0.6, "Done.")
          rv$s3_run <- TRUE
        }
      )
    })

    output$plot_s3_before <- shiny::renderPlot({
      shiny::req(rv$s3_run, rv$step >= 3L)
      plt <- rv$s3_plot_before
      shiny::validate(shiny::need(!is.null(plt), "Click 'Run Normalisation' to generate."))
      plt
    })
    output$plot_s3_after <- shiny::renderPlot({
      shiny::req(rv$s3_run, rv$step >= 3L)
      plt <- rv$s3_plot_after
      shiny::validate(shiny::need(!is.null(plt), "Click 'Run Normalisation' to generate."))
      plt
    })

    output$s3_remove_picker <- shiny::renderUI({
      e <- rv$data_s2; shiny::req(e, rv$s3_run)
      .pp_native_picker(ns("s3_remove_samples"), colnames(e))
    })

    shiny::observeEvent(input$s3_remove_btn, {
      to_rm <- input$s3_remove_samples; shiny::req(length(to_rm) >= 1L, rv$data_s2)
      keep <- setdiff(colnames(rv$data_s2), to_rm)
      if (length(keep) < 2L) {
        shiny::showNotification("Cannot remove: fewer than 2 samples would remain.", type = "error")
        return()
      }
      rv$data_s2 <- rv$data_s2[, keep, drop = FALSE]
      if (!is.null(rv$data_s1)) rv$data_s1 <- rv$data_s1[, keep, drop = FALSE]
      .align_meta(keep)
      glog(paste0("Normalisation step: removed ", length(to_rm), " sample(s): [",
                   paste(to_rm, collapse = ", "), "]. Re-run to refresh boxplots."))
      shiny::showNotification(
        paste0(length(to_rm), " sample(s) removed. Re-run to refresh."),
        type = "warning", duration = 5)
    })

    shiny::observeEvent(input$proceed_s3, {
      shiny::req(rv$step >= 3L, rv$norm_linear, rv$data_s3)
      method <- shiny::isolate(input$s3_method %||% "TMM")
      rv$s3_summary <- paste0(method, " — ",
                                format(nrow(rv$data_s3), big.mark = ","),
                                " genes × ", ncol(rv$data_s3), " samples")
      glog(paste0("Normalisation: ", method, " applied to all samples — ",
                   format(nrow(rv$data_s3), big.mark = ","), " genes × ",
                   ncol(rv$data_s3), " samples."))
      rv$step <- 4L; rv$active_step <- 4L; rv$s4_run <- FALSE; rv$s4_bc_run <- FALSE
    })

    shiny::observeEvent(input$back_to_s3, {
      rv$active_step <- 3L; rv$step <- 3L; rv$s3_run <- FALSE
      rv$norm_linear <- NULL; rv$data_s3 <- NULL
    })

    # =========================================================================
    # STEP 4 — PCA & Optional Batch Correction
    # PCA computation is split from rendering: changing colour does NOT re-run PCA
    # =========================================================================

    output$step4_card <- shiny::renderUI({
      if (rv$step < 4L) {
        return(shiny::div(class = "pp-locked",
          bslib::card(class = "pp-step-card",
            bslib::card_header(shiny::tags$div(style = "display:flex;align-items:center;",
              shiny::tags$span(class = "pp-num pp-num-locked", "4"),
              shiny::tags$strong("Step 4 — PCA & Batch Correction (Optional)"))),
            shiny::div(style = "padding:14px;color:#9ca3af;",
              "Complete Steps 1–3 to unlock."))))
      }
      if (rv$step > 4L && rv$active_step != 4L) {
        return(.pp_done_card(4, "PCA & Batch Correction",
          if (!is.null(rv$pp_summary)) rv$pp_summary$bc_label else "done",
          ns("back_to_s4"),
          plot_ui = shiny::fluidRow(
            shiny::column(6, shiny::plotOutput(ns("pca_before"), height = "200px")),
            shiny::column(6, if (rv$s4_bc_run) shiny::plotOutput(ns("pca_after"), height = "200px"))
          )
        ))
      }
      meta_cols <- setdiff(colnames(rv$meta_s), "SampleID")
      controls <- shiny::fluidRow(
        shiny::column(4,
          shiny::tags$b("PCA colour variable"),
          if (length(meta_cols) > 0L)
            shiny::selectInput(ns("s4_pca_color"), label = NULL,
                                choices = meta_cols, selected = meta_cols[1])
          else shiny::helpText("No metadata variables."),
          shiny::tags$p(style = "font-size:0.78em;color:#6b7280;margin:2px 0 0 0;",
            "Changing colour re-draws the plot without re-running PCA."),
          shiny::uiOutput(ns("s4_color_warn"))
        ),
        shiny::column(8,
          shiny::checkboxInput(ns("s4_do_bc"),
            "Apply batch correction (voom / lmFit)", value = FALSE),
          shiny::uiOutput(ns("s4_bc_section"))
        )
      )
      bslib::card(class = "pp-step-card",
        bslib::card_header(shiny::tags$div(style = "display:flex;align-items:center;",
          shiny::tags$span(class = "pp-num pp-num-active", "4"),
          shiny::tags$strong("Step 4 — PCA & Batch Correction (Optional)"))),
        shiny::div(style = "padding:14px;",
          controls,
          shiny::hr(style = "margin:8px 0;"),
          shiny::uiOutput(ns("s4_run_section"))))
    })

    output$s4_color_warn <- shiny::renderUI({ .pp_color_warn(rv$meta_s, input$s4_pca_color) })

    output$s4_bc_section <- shiny::renderUI({
      if (!isTRUE(input$s4_do_bc)) return(NULL)
      meta_cols <- setdiff(colnames(rv$meta_s), "SampleID")
      shiny::tags$details(
        style = "border:1px solid #dee2e6;border-radius:6px;padding:8px 12px;margin-top:4px;",
        open = TRUE,
        shiny::tags$summary(style = "cursor:pointer;font-weight:600;font-size:0.88em;user-select:none;",
                             "Batch correction options"),
        shiny::div(style = "margin-top:8px;",
          shiny::fluidRow(
            shiny::column(6,
              shiny::tags$b(style = "font-size:0.88em;", "Batch variable(s) to remove"),
              shiny::helpText("Technical sources of variation (e.g. DatasetID, batch).",
                               style = "font-size:0.78em;margin-top:0;"),
              shiny::div(
                style = "max-height:130px;overflow-y:auto;border:1px solid #dee2e6;border-radius:4px;padding:4px 8px;",
                shiny::checkboxGroupInput(ns("s4_batch_cols"), label = NULL,
                                           choices = meta_cols, selected = NULL))
            ),
            shiny::column(6,
              shiny::tags$b(style = "font-size:0.88em;", "Biological variable(s) to retain"),
              shiny::helpText("Effects NOT to remove (e.g. Condition, CellType).",
                               style = "font-size:0.78em;margin-top:0;"),
              shiny::div(
                style = "max-height:130px;overflow-y:auto;border:1px solid #dee2e6;border-radius:4px;padding:4px 8px;",
                shiny::checkboxGroupInput(ns("s4_effect_cols"), label = NULL,
                                           choices = meta_cols, selected = NULL))
            )
          ),
          shiny::fluidRow(
            shiny::column(6,
              shiny::tags$b(style = "font-size:0.88em;", "Colour after-BC PCA by"),
              shiny::selectInput(ns("s4_pca_after_color"), label = NULL,
                                  choices = meta_cols,
                                  selected = if (length(meta_cols) > 0) meta_cols[1] else NULL)
            ),
            shiny::column(6, shiny::div(style = "padding-top:24px;",
              shiny::actionButton(ns("bc_run_btn"), "\U0001f4ca Run Batch Correction",
                                   class = "btn-sm btn-primary")))
          )
        )
      )
    })

    output$s4_run_section <- shiny::renderUI({
      n_g <- if (!is.null(rv$data_s3)) nrow(rv$data_s3) else "?"
      n_s <- if (!is.null(rv$data_s3)) ncol(rv$data_s3) else "?"
      if (!rv$s4_run)
        return(shiny::div(class = "pp-run-btn",
          .pp_slow_warn(paste0(
            "⏳ PCA decomposes the full expression matrix (", n_g, " genes × ", n_s,
            " samples). For large datasets this may take a minute.")),
          shiny::actionButton(ns("s4_run_btn"), "\U0001f4ca Show PCA",
                               class = "btn-primary btn-sm")
        ))

      bc_done <- rv$s4_bc_run && isTRUE(input$s4_do_bc) && !is.null(rv$s4_bc_result)

      finalize_ui <- if (bc_done) {
        shiny::tagList(
          shiny::tags$p(style = "font-size:0.84em;color:#555;margin-top:10px;",
            "Choose which data to use for all downstream analyses:"),
          shiny::div(style = "display:flex;gap:8px;flex-wrap:wrap;margin-top:4px;",
            shiny::actionButton(ns("finalize_bc"),
              "✓ Use batch-corrected data", class = "btn-warning btn-sm"),
            shiny::actionButton(ns("finalize_norm"),
              "✓ Use normalised data (without batch correction)",
              class = "btn-warning btn-sm")
          )
        )
      } else {
        shiny::div(style = "display:flex;align-items:center;gap:12px;margin-top:12px;",
          shiny::actionButton(ns("finalize_btn"), "✓ Finalise & apply preprocessing",
                               class = "btn-success btn-sm"),
          shiny::tags$span(style = "font-size:0.82em;color:#6c757d;",
                            "Updates expression data used by all downstream analyses."))
      }

      shiny::tagList(
        shiny::div(style = "margin-bottom:6px;",
          .pp_slow_warn(paste0(
            "⏳ Re-running PCA recomputes the decomposition. ",
            "Changing the colour variable below only redraws the plot — no recomputation.")),
          shiny::actionButton(ns("s4_run_btn"), "\U0001f504 Re-run PCA",
                               class = "btn-sm btn-outline-primary")
        ),
        shiny::fluidRow(
          shiny::column(6,
            shiny::plotOutput(ns("pca_before"), height = "330px"),
            .pp_expand_btn(ns("expand_pca_before"))
          ),
          shiny::column(6,
            shiny::uiOutput(ns("pca_after_panel")))
        ),
        shiny::hr(style = "margin:10px 0;"),
        finalize_ui
      )
    })

    # PCA COMPUTATION — runs on button click, stores df + eigenvalues
    shiny::observeEvent(input$s4_run_btn, {
      shiny::req(rv$step >= 4L, rv$data_s3)
      n_g <- nrow(rv$data_s3); n_s <- ncol(rv$data_s3)
      shiny::withProgress(
        message = paste0("\U0001f4ca Computing PCA — ",
                          format(n_g, big.mark = ","), " genes × ", n_s, " samples…"),
        value = 0, {
          shiny::incProgress(0.4, "Decomposing variance…")
          result <- tryCatch(
            .pp_run_pca(rv$data_s3),
            error = function(e) {
              shiny::showNotification(paste("PCA error:", e$message), type = "error"); NULL
            }
          )
          if (!is.null(result)) {
            rv$s4_pca_df_before <- result$df
            rv$s4_pca_ev_before <- result$ev
          }
          shiny::incProgress(0.6, "Done.")
          rv$s4_run <- TRUE
        }
      )
    })

    # PCA RENDERING — reacts to colour input (no recomputation)
    output$pca_before <- shiny::renderPlot({
      shiny::req(rv$s4_run, rv$s4_pca_df_before, rv$s4_pca_ev_before)
      col   <- input$s4_pca_color  # reactive — triggers re-render only
      title <- if (isTRUE(input$s4_do_bc) && rv$s4_bc_run)
                 "PCA — before batch correction" else "PCA"
      .pp_pca_from_df(rv$s4_pca_df_before, rv$s4_pca_ev_before, rv$meta_s, col, title)
    })

    output$pca_after_panel <- shiny::renderUI({
      shiny::req(rv$s4_run)
      if (!isTRUE(input$s4_do_bc))
        return(shiny::div(
          style = "display:flex;align-items:center;justify-content:center;height:330px;color:#9ca3af;font-size:0.88em;text-align:center;",
          shiny::HTML("Enable batch correction above<br>and click 'Run Batch Correction'.")))
      if (!rv$s4_bc_run)
        return(shiny::div(
          style = "display:flex;align-items:center;justify-content:center;height:330px;color:#9ca3af;font-size:0.88em;text-align:center;",
          shiny::HTML("Select batch and biological variables,<br>then click 'Run Batch Correction'.")))
      shiny::tagList(
        shiny::plotOutput(ns("pca_after"), height = "330px"),
        .pp_expand_btn(ns("expand_pca_after"))
      )
    })

    # BATCH CORRECTION
    shiny::observeEvent(input$bc_run_btn, {
      shiny::req(rv$step >= 4L, rv$norm_linear, rv$meta_s)
      batch_cols  <- shiny::isolate(input$s4_batch_cols)
      effect_cols <- shiny::isolate(input$s4_effect_cols)
      if (length(batch_cols) < 1L) {
        shiny::showNotification("Select at least one batch variable.", type = "warning"); return()
      }
      if (length(effect_cols) < 1L) {
        shiny::showNotification("Select at least one biological variable.", type = "warning"); return()
      }
      norm_counts <- as.matrix(rv$norm_linear); meta <- rv$meta_s
      n_s <- ncol(norm_counts)
      shiny::withProgress(
        message = paste0("⚙️ Running batch correction — ",
                          n_s, " samples, batch: [", paste(batch_cols, collapse = ", "), "]…"),
        value = 0, {
          shiny::incProgress(0.15, "Building design matrices…")
          corrdf <- tryCatch({
            bc_form  <- stats::reformulate(paste(batch_cols,  collapse = " + "), intercept = FALSE)
            eff_form <- stats::reformulate(paste(effect_cols, collapse = " + "), intercept = FALSE)
            mmbatch  <- stats::model.matrix(bc_form,  data = meta)
            mmkeep   <- stats::model.matrix(eff_form, data = meta)
            mm       <- cbind(mmkeep, mmbatch)
            shiny::incProgress(0.30, "Running voom + lmFit…")
            D0  <- edgeR::DGEList(norm_counts); D0 <- edgeR::calcNormFactors(D0)
            y   <- limma::voom(D0, mm, plot = FALSE)
            fit <- limma::lmFit(y, mm)
            n_keep <- ncol(mmkeep)
            beta   <- fit$coefficients[, -(1:n_keep), drop = FALSE]; beta[is.na(beta)] <- 0
            corrcounts <- as.matrix(y$E) - beta %*% t(mmbatch)
            offset <- apply(corrcounts, 1, min) - apply(log2(norm_counts + 1), 1, min)
            as.data.frame(2^(corrcounts - offset))
          }, error = function(e) {
            shiny::showNotification(paste("Batch correction failed:", e$message),
                                     type = "error", duration = 12)
            NULL
          })

          if (!is.null(corrdf)) {
            shiny::incProgress(0.25, "Computing post-correction PCA…")
            col_after <- shiny::isolate(input$s4_pca_after_color %||% input$s4_pca_color)
            log_bc    <- log2(as.matrix(corrdf) + 1)
            pca_after <- tryCatch(
              .pp_run_pca(log_bc),
              error = function(e) {
                shiny::showNotification(paste("Post-BC PCA error:", e$message), type = "error"); NULL
              }
            )
            rv$s4_bc_result    <- corrdf
            rv$s4_pca_df_after <- if (!is.null(pca_after)) pca_after$df else NULL
            rv$s4_pca_ev_after <- if (!is.null(pca_after)) pca_after$ev else NULL
            rv$s4_bc_run       <- TRUE
            glog(paste0("Batch Correction: removed [", paste(batch_cols, collapse = ", "),
                         "], retained [", paste(effect_cols, collapse = ", "), "]. ",
                         n_s, " samples."))
          }
          shiny::incProgress(0.30, "Done.")
        }
      )
    })

    # AFTER-BC PCA RENDERING — reacts to colour input (no recomputation)
    output$pca_after <- shiny::renderPlot({
      shiny::req(rv$s4_run, rv$s4_bc_run,
                 rv$s4_pca_df_after, rv$s4_pca_ev_after)
      col <- input$s4_pca_after_color %||% input$s4_pca_color
      .pp_pca_from_df(rv$s4_pca_df_after, rv$s4_pca_ev_after,
                       rv$meta_s, col, "PCA — after batch correction")
    })

    # ── Finalise ──────────────────────────────────────────────────────────────
    .do_finalize <- function(use_bc) {
      shiny::req(rv$step >= 4L)
      if (use_bc) {
        bc <- rv$s4_bc_result
        if (is.null(bc)) {
          shiny::showNotification("Batch correction not completed.", type = "error"); return()
        }
        rv$final_data <- bc
      } else {
        rv$final_data <- as.data.frame(rv$data_s3)
      }
      n_g  <- nrow(rv$final_data); n_s <- ncol(rv$final_data)
      orig_g <- rv$orig_n_genes %||% n_g; orig_s <- rv$orig_n_samples %||% n_s
      genes_removed   <- orig_g - n_g
      genes_kept      <- n_g
      samples_removed <- orig_s - n_s
      norm_method     <- shiny::isolate(input$s3_method %||% "TMM")
      bc_vars         <- if (use_bc) shiny::isolate(input$s4_batch_cols) else NULL

      rv$pp_summary <- list(
        n_genes         = genes_kept,
        genes_kept      = genes_kept,
        genes_removed   = genes_removed,
        n_samples       = n_s,
        samples_removed = samples_removed,
        norm_method     = norm_method,
        bc_applied      = use_bc,
        bc_vars         = bc_vars,
        bc_effect       = if (use_bc) shiny::isolate(input$s4_effect_cols) else NULL,
        bc_label        = if (use_bc)
                            paste0("batch-corrected ([", paste(bc_vars, collapse = ", "), "])")
                          else "normalised (no batch correction)"
      )
      glog(paste0("Preprocessing finalised",
                   if (use_bc)
                     paste0(" WITH batch correction [", paste(bc_vars, collapse = ", "), "]")
                   else " (no batch correction)",
                   " — ", format(genes_kept, big.mark = ","), " genes kept",
                   if (genes_removed > 0L) paste0(", ", format(genes_removed, big.mark = ","), " removed") else "",
                   " — ", n_s, " samples."))
      rv$downloads_available <- TRUE
      rv$finalized   <- rv$finalized + 1L
      rv$step        <- 5L; rv$active_step <- 5L
      shiny::showNotification(
        shiny::HTML(paste0(
          "✅ Preprocessing complete — ",
          format(genes_kept, big.mark = ","), " genes × ", n_s, " samples. ",
          "Expression data updated for all analyses.")),
        type = "default", duration = 8)
    }

    shiny::observeEvent(input$finalize_btn,  { .do_finalize(FALSE) })
    shiny::observeEvent(input$finalize_bc,   { .do_finalize(TRUE)  })
    shiny::observeEvent(input$finalize_norm, { .do_finalize(FALSE) })

    shiny::observeEvent(input$back_to_s4, {
      rv$active_step <- 4L; rv$step <- 4L
      rv$s4_run <- FALSE; rv$s4_bc_run <- FALSE
      rv$final_data <- NULL; rv$finalized <- 0L
      rv$pp_summary <- NULL; rv$downloads_available <- FALSE
    })

    # =========================================================================
    # SUMMARY CARD
    # =========================================================================

    output$pp_summary_card <- shiny::renderUI({
      shiny::req(rv$finalized > 0L, !is.null(rv$pp_summary))
      s <- rv$pp_summary
      .row <- function(label, value)
        shiny::tags$div(class = "pp-summary-row",
          shiny::tags$span(class = "pp-summary-label", label),
          shiny::tags$span(value))

      shiny::div(class = "pp-summary-card",
        shiny::tags$h6(style = "color:#EBB43E;font-weight:700;margin-bottom:12px;",
                        "✅ Preprocessing Summary"),
        .row("Final dataset:",
             paste0(format(s$genes_kept, big.mark = ","), " genes × ", s$n_samples, " samples")),
        .row("Genes kept:",   format(s$genes_kept, big.mark = ",")),
        if (s$genes_removed > 0L)
          .row("Genes removed:",
               paste0(format(s$genes_removed, big.mark = ","), " (below expression threshold)")),
        if (s$samples_removed > 0L)
          .row("Samples removed:", as.character(s$samples_removed)),
        .row("Normalisation:", s$norm_method),
        .row("Batch correction:",
             if (s$bc_applied)
               paste0("✓ Applied — batch: [", paste(s$bc_vars, collapse = ", "),
                      "] | retained: [", paste(s$bc_effect, collapse = ", "), "]")
             else "Not applied"),
        .row("Data used downstream:",
             if (s$bc_applied) "Batch-corrected expression data" else "Normalised expression data"),
        shiny::div(style = "margin-top:12px;font-size:0.82em;color:#6c757d;",
          "Download the preprocessed data and metadata from the sidebar (optional). ",
          "The processing log is available in the bar at the bottom of the page.")
      )
    })

    # =========================================================================
    # HELPERS
    # =========================================================================

    .align_meta <- function(keep) {
      meta <- rv$meta_s; if (is.null(meta)) return()
      id_col <- if ("SampleID" %in% colnames(meta)) "SampleID" else colnames(meta)[1]
      idx    <- match(keep, as.character(meta[[id_col]])); valid <- !is.na(idx)
      if (sum(valid) >= 1L) rv$meta_s <- meta[idx[valid], , drop = FALSE]
    }

    list(final_data = shiny::reactive(rv$final_data),
         finalized  = shiny::reactive(rv$finalized))
  })
}
