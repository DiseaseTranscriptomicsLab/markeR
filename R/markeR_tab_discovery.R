# =============================================================================
# Discovery Tab – Shiny Module
#
# Explores associations between all (selected) metadata variables and the gene
# set signal, using either scoring or GSEA methods.
#
# Primary functions (NOT modified here - only their returned objects are laid
# out and sized):
#   VariableAssociation()        – umbrella wrapper (Score or GSEA path)
#   Score_VariableAssociation()  – score path. Returns, per call:
#     $Overall             data.frame(Variable, Cohen_f, P_Value)
#     $Contrasts           data.frame(Variable, Contrast, Group1, Group2,
#                                      CohenD, PValue, padj) or NULL if every
#                                      tested variable is numeric
#     $plot_overall        one ggplot: Cohen's f lollipop, one row per variable
#     $plot_contrasts      one ggplot: Cohen's d lollipop, facetted by Variable
#     $plot_distributions  named list of ggplots, one per variable (density or
#                                      scatter depending on variable type)
#   GSEA_VariableAssociation()   – GSEA path. Returns:
#     $data                data.frame(pathway, NES, padj, stat_used, Contrast)
#                                      flattened across every tested variable
#     $plot                one ggplot: NES lollipop, one row per contrast
#
# This tab lays those pieces out as separate, properly-sized cards instead of
# printing the pre-baked combined figure, and adds a "Method Summary" tab that
# runs all four methods (log2-median, ranking, ssGSEA, GSEA) for the same
# variables/gene set and shows them side by side as a heatmap (effect size +
# significance), similar in spirit to a "gene set x method" overview figure.
#
# Data contract (reactive accessors from main server):
#   get_expr()      – data.frame: genes × samples (normalised counts)
#   get_meta()      – data.frame: samples × variables (first col = SampleID)
#   get_gene_sets() – named list: character vectors or 2-col data.frames
# =============================================================================

# ── UI -----------------------------------------------------------------------

#' @importFrom bslib layout_sidebar sidebar navset_card_tab nav_panel card card_header
#' @importFrom shiny NS radioButtons sliderInput selectInput numericInput actionButton
#'   uiOutput plotOutput hr h4 p checkboxGroupInput
#' @importFrom shinyWidgets pickerInput pickerOptions
#' @importFrom DT DTOutput
if (!exists("%||%", mode = "function"))
  `%||%` <- function(a, b) if (!is.null(a)) a else b

# ---- Dynamic height helpers (mirrors the *_h() reactiveVal pattern used in
# the Benchmarking tab: height grows with how much content a plot has to show,
# instead of a single fixed pixel value that gets cramped or wastes space). --
.disc_h_overall <- function(n_vars)
  as.integer(min(1400L, max(260L, 220L + n_vars * 42L)))

# Distribution plots are always laid out at most 2 per row (any more and the
# per-variable density/scatter panels get too narrow to read).
.DISC_DIST_NCOL <- 2L

.disc_h_dist <- function(n_vars, ncol_dist = .DISC_DIST_NCOL)
  as.integer(min(3200L, max(320L, 120L + ceiling(n_vars / ncol_dist) * 300L)))

.disc_h_contrasts <- function(n_rows, n_groups)
  as.integer(min(5000L, max(400L, 160L + n_rows * 42L + n_groups * 55L)))

.disc_h_gsea <- function(n_rows)
  as.integer(min(4000L, max(320L, 150L + n_rows * 32L)))

.disc_h_summary <- function(n_vars)
  as.integer(min(2200L, max(420L, 260L + n_vars * 48L)))

discoveryUI <- function(id) {
  ns <- shiny::NS(id)

  bslib::layout_sidebar(
    # fillable = FALSE at this top level (not just on the nav_panels below):
    # otherwise the main content area still participates in bslib's flexbox
    # "fill" chain and constrains card height to the viewport, which is what
    # was forcing scrollbars inside individual plot cards instead of letting
    # the page itself scroll.
    fillable = FALSE,

    sidebar = bslib::sidebar(
      width = 330,

      shiny::div(
        style = "padding-bottom:10px;",
        shiny::h4("Discovery Mode", style = "font-weight:700; color:#1E497B; margin-bottom:4px;"),
        shiny::p(
          "Systematically test all metadata variables for association with",
          " the gene set signal.",
          style = "color:#6c757d; font-size:0.85em;"
        ),
        shiny::div(
          class = "alert alert-warning", style = "font-size:0.82em; padding:8px 10px; margin-bottom:0;",
          shiny::icon("circle-exclamation"),
          " Discovery Mode is for ", shiny::tags$b("hypothesis generation"),
          ": it quickly screens many variables so you can spot promising",
          " candidates. If something here looks interesting, take it to the",
          shiny::tags$b(" Benchmarking"), " tab for a closer, more rigorous",
          " look. And however it holds up there, it should still be",
          " confirmed in an independent dataset before you rely on it."
        ),
        shiny::hr()
      ),

      # ---- Gene set selector (single) ---------------------------------------
      shiny::uiOutput(ns("disc_gs_picker_ui")),

      # ---- Metadata variables to test ---------------------------------------
      shiny::uiOutput(ns("disc_vars_ui")),
      shiny::uiOutput(ns("disc_vars_warn_ui")),

      shiny::hr(),

      # ---- Mode-specific controls: only the ones relevant to the active
      # main tab are shown, so it's unambiguous which button/setting feeds
      # which output. Shared settings (contrast mode, significance, display)
      # stay visible for both, since both modes use them. -------------------
      shiny::conditionalPanel(
        condition = sprintf("input['%s'] == 'single'", ns("disc_active_tab")),
        shiny::tags$label(
          "Association method ", .method_tooltip,
          style = "font-size:0.85em; font-weight:600;"
        ),
        shiny::radioButtons(
          ns("disc_method"), label = NULL,
          choices  = c(
            "Log-median score" = "logmedian",
            "Ranking score"    = "ranking",
            "ssGSEA score"     = "ssGSEA",
            "GSEA enrichment"  = "GSEA"
          ),
          selected = "logmedian"
        )
      ),
      shiny::conditionalPanel(
        condition = sprintf("input['%s'] == 'summary'", ns("disc_active_tab")),
        shiny::div(
          class = "alert alert-secondary", style = "font-size:0.8em; padding:8px 10px;",
          shiny::icon("circle-info"),
          " Runs ", shiny::strong("all four methods"),
          " (log2-median, ranking, ssGSEA, GSEA) automatically for the",
          " variables and gene set selected above - no method choice needed here."
        )
      ),

      .bm_collapsible_settings(
        .bm_contrast_mode_radio(ns, "disc_mode", selected = "simple"),
        shiny::tags$hr(style = "margin:8px 0;"),
        shiny::sliderInput(
          ns("disc_sig"), "Significance threshold (adjusted p-value):",
          min = 0.001, max = 0.2, value = 0.05, step = 0.001
        ),
        shiny::selectInput(
          ns("disc_padj"), "p-adjust method:",
          choices  = c("BH", "bonferroni", "holm", "BY", "fdr", "none"),
          selected = "BH"
        ),
        shiny::tags$hr(style = "margin:8px 0;"),
        shiny::tags$details(
          shiny::tags$summary(
            style = paste("cursor:pointer; font-size:0.85em; font-weight:600; color:#555;",
                          "user-select:none; list-style:none;"),
            "Display options"
          ),
          shiny::div(
            style = "display:grid; grid-template-columns:1fr 1fr; gap:8px; margin-top:8px;",
            shiny::numericInput(ns("disc_pointsize"), "Point size:",
                                value = 4, min = 1, max = 10, step = 1),
            shiny::numericInput(ns("disc_wraplabels"), "Contrast wrap (chars):",
                                value = 32, min = 8, max = 60, step = 2),
            shiny::numericInput(ns("disc_titlesize"), "Title size (pt):",
                                value = 8, min = 6, max = 26, step = 1),
            shiny::numericInput(ns("disc_labsize"), "Label size (pt):",
                                value = 8, min = 6, max = 22, step = 1)
          )
        )
      ),

      shiny::hr(),

      shiny::conditionalPanel(
        condition = sprintf("input['%s'] == 'single'", ns("disc_active_tab")),
        shiny::actionButton(
          ns("run_discovery"), "Run Selected Method",
          class = "btn-primary", width = "100%"
        )
      ),
      shiny::conditionalPanel(
        condition = sprintf("input['%s'] == 'summary'", ns("disc_active_tab")),
        shiny::actionButton(
          ns("run_summary"), "Build Method Summary",
          class = "btn-primary", width = "100%"
        ),
        shiny::tags$p(
          style = "color:#888; font-size:0.78em; margin:6px 0 0;",
          "Slower than a single method - GSEA is run once per variable."
        )
      )
    ),

    # ── Main panel ----------------------------------------------------------
    # `fillable = FALSE` on every tab: without it, bslib fits card content to
    # the viewport height and scrolls *inside* the card. Since these plots are
    # sized to their content (which can be tall), we want the normal page
    # scrollbar instead, matching the fix already used for the Benchmarking
    # tab's FPR Simulation panel.
    bslib::navset_card_tab(
      id = ns("disc_active_tab"),

      # ---- Tab 1: Single Method Explorer -------------------------------------
      bslib::nav_panel(
        title = "Single Method Explorer", value = "single", fillable = FALSE,

        shiny::uiOutput(ns("disc_status_ui")),
        shiny::uiOutput(ns("disc_overall_card_ui")),
        shiny::uiOutput(ns("disc_dist_card_ui")),
        shiny::uiOutput(ns("disc_contrasts_card_ui")),
        shiny::uiOutput(ns("disc_gsea_card_ui")),

        shiny::tags$details(
          style = "margin-top:16px;",
          shiny::tags$summary(
            style = paste("cursor:pointer; font-size:0.9em; font-weight:600; color:#555;",
                          "user-select:none; padding:6px 2px;"),
            shiny::icon("table", style = "color:#2D6A4F; font-size:0.9em;"),
            " Results table & insights"
          ),
          shiny::div(
            style = "padding-top:10px;",
            bslib::card(
              full_screen = TRUE,
              bslib::card_header("Variable Association Results"),
              shiny::uiOutput(ns("disc_table_status_ui")),
              DT::DTOutput(ns("disc_table"))
            )
          )
        )
      ),

      # ---- Tab 2: Method Summary ---------------------------------------------
      bslib::nav_panel(
        title = "Method Summary (all methods)", value = "summary", fillable = FALSE,

        shiny::uiOutput(ns("disc_summary_status_ui")),
        shiny::uiOutput(ns("disc_summary_card_ui")),

        shiny::tags$details(
          style = "margin-top:16px;",
          shiny::tags$summary(
            style = paste("cursor:pointer; font-size:0.9em; font-weight:600; color:#555;",
                          "user-select:none; padding:6px 2px;"),
            shiny::icon("table", style = "color:#2D6A4F; font-size:0.9em;"),
            " Underlying values"
          ),
          shiny::div(
            style = "padding-top:10px;",
            bslib::card(
              full_screen = TRUE,
              bslib::card_header("Score methods - Cohen's f"),
              DT::DTOutput(ns("disc_summary_table_score"))
            ),
            bslib::card(
              full_screen = TRUE,
              bslib::card_header("GSEA - significant contrasts per variable"),
              DT::DTOutput(ns("disc_summary_table_gsea"))
            ),
            bslib::card(
              full_screen = TRUE,
              bslib::card_header("GSEA - contrast-level detail"),
              shiny::div(
                class = "text-muted", style = "font-size:0.78em; padding:0 6px 6px;",
                shiny::icon("circle-info"),
                " Every individual contrast that was considered for each variable,",
                " with its gene set, NES, and adjusted p-value."
              ),
              DT::DTOutput(ns("disc_summary_table_gsea_detail"))
            )
          )
        )
      )
    )
  )
}

# ── Server -------------------------------------------------------------------

#' @importFrom shiny moduleServer reactive reactiveVal observeEvent req
#'   renderPlot renderUI withProgress showNotification isolate
#' @importFrom DT renderDT datatable
discoveryServer <- function(id, get_expr, get_meta, get_gene_sets) {

  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ── Cached results --------------------------------------------------------
    discovery_result <- shiny::reactiveVal(NULL)
    summary_result    <- shiny::reactiveVal(NULL)

    ov_h   <- shiny::reactiveVal(300L)
    dist_h <- shiny::reactiveVal(400L)
    cont_h <- shiny::reactiveVal(350L)
    gsea_h <- shiny::reactiveVal(400L)
    sum_h  <- shiny::reactiveVal(300L)

    # ── Dynamic UI: gene set selector ------------------------------------------
    output$disc_gs_picker_ui <- shiny::renderUI({
      gs <- get_gene_sets()
      shiny::req(!is.null(gs), length(gs) > 0)
      shinyWidgets::pickerInput(
        ns("disc_gs"), "Gene set to test:",
        choices  = names(gs), selected = names(gs)[1],
        options  = shinyWidgets::pickerOptions(liveSearch = TRUE)
      )
    })

    # ── Dynamic UI: variable multi-selector -------------------------------------
    # Only variables with more than one observed level can be tested at all
    # (a constant column has no variance for lm()/cor.test() to work with),
    # so those are excluded from the choices entirely rather than just left
    # unselected. Variables with many levels (> 6) are still offered and
    # selected, but flagged - contrasts/lollipops get unwieldy and slow with
    # a lot of levels, especially in "All combinations" contrast mode.
    output$disc_vars_ui <- shiny::renderUI({
      meta <- get_meta()
      shiny::req(!is.null(meta))
      cols <- setdiff(colnames(meta), colnames(meta)[1])

      level_counts <- vapply(cols, function(cl) {
        length(unique(stats::na.omit(meta[[cl]])))
      }, integer(1))

      usable   <- cols[level_counts > 1]
      excluded <- setdiff(cols, usable)

      shiny::tagList(
        shinyWidgets::pickerInput(
          ns("disc_vars"), "Metadata variables to test:",
          # Nothing pre-selected - the user chooses which variables to test.
          choices  = usable, selected = character(0), multiple = TRUE,
          options  = shinyWidgets::pickerOptions(
            actionsBox = TRUE, liveSearch = TRUE,
            selectedTextFormat = "count > 3"
          )
        ),
        if (length(excluded) > 0) {
          shiny::tags$p(
            style = "color:#888; font-size:0.78em; margin:2px 0 4px;",
            shiny::icon("circle-info"),
            sprintf(" Excluded (only 1 level, no variance to test): %s.",
                    paste(excluded, collapse = ", "))
          )
        }
      )
    })

    # Separate output (depends on the live selection, not just the metadata)
    # so the "too many levels" warning only calls out variables the user has
    # actually selected to test, not every eligible variable.
    output$disc_vars_warn_ui <- shiny::renderUI({
      meta <- get_meta()
      shiny::req(!is.null(meta))
      sel <- input$disc_vars
      shiny::req(!is.null(sel), length(sel) > 0)

      level_counts <- vapply(sel, function(cl) {
        length(unique(stats::na.omit(meta[[cl]])))
      }, integer(1))
      highlvl <- sel[level_counts > 6]

      if (length(highlvl) > 0) {
        shiny::tags$p(
          style = "color:#b45309; font-size:0.78em; margin:2px 0 0;",
          shiny::icon("triangle-exclamation"),
          sprintf(" %s %s more than 6 levels - contrasts/plots may get large",
                  paste(highlvl, collapse = ", "),
                  if (length(highlvl) > 1) "have" else "has"),
          " and slow, especially in \"All combinations\" mode."
        )
      }
    })

    # ── Shared settings reader ---------------------------------------------
    .disc_settings <- function() {
      list(
        gs        = shiny::isolate(input$disc_gs),
        vars      = shiny::isolate(input$disc_vars),
        mode      = shiny::isolate(input$disc_mode)      %||% "simple",
        sig_thr   = shiny::isolate(input$disc_sig)        %||% 0.05,
        padj_m    = shiny::isolate(input$disc_padj)       %||% "BH",
        pointsize = shiny::isolate(input$disc_pointsize)  %||% 4,
        wraplab   = shiny::isolate(input$disc_wraplabels) %||% 32,
        titlesize = shiny::isolate(input$disc_titlesize)  %||% 8,
        labsize   = shiny::isolate(input$disc_labsize)    %||% 8
      )
    }

    # Blanket font-size override applied on top of whatever the (unmodified)
    # core functions returned, so the Display-option sliders always take
    # effect regardless of which sub-plot is being shown.
    .disc_font_theme <- function(labsize, titlesize) {
      ggplot2::theme(
        axis.text     = ggplot2::element_text(size = labsize),
        axis.title    = ggplot2::element_text(size = labsize + 1),
        strip.text    = ggplot2::element_text(size = labsize + 1, face = "bold"),
        legend.text   = ggplot2::element_text(size = labsize),
        legend.title  = ggplot2::element_text(size = labsize),
        # Shrink the legend keys (the coloured rectangles/points themselves)
        # along with the text - with many levels, the default (fixed, larger)
        # key size is what makes the legend overlap/get cropped even after
        # the text itself is made smaller.
        legend.key.size    = ggplot2::unit(max(0.28, labsize / 20), "cm"),
        legend.spacing.y   = ggplot2::unit(1, "pt"),
        plot.title    = ggplot2::element_text(size = titlesize),
        plot.subtitle = ggplot2::element_text(size = titlesize - 2)
      )
    }

    # ── Run: single-method Discovery -------------------------------------------
    shiny::observeEvent(input$run_discovery, {
      expr <- get_expr(); meta <- get_meta(); gs <- get_gene_sets()
      s    <- .disc_settings()
      shiny::req(expr, meta, gs, s$gs, s$vars)

      method          <- shiny::isolate(input$disc_method) %||% "logmedian"
      gene_set_single <- gs[s$gs]

      shiny::withProgress(message = "Running variable association analysis...", value = 0, {
        result <- tryCatch({
          shiny::incProgress(0.2, detail = paste("Method:", method))

          res <- VariableAssociation(
            method          = method, data = expr, metadata = meta,
            cols            = s$vars, gene_set = gene_set_single,
            mode            = s$mode, sig_threshold = s$sig_thr,
            p.adjust.method = s$padj_m, widthlabels = s$wraplab,
            labsize         = s$labsize, titlesize = s$titlesize,
            pointSize       = s$pointsize, printplt = FALSE
          )

          shiny::incProgress(0.8, detail = "Finalising...")
          list(result = res, method = method, vars = s$vars, mode = s$mode,
              gs_name = s$gs, sig_thr = s$sig_thr,
              labsize = s$labsize, titlesize = s$titlesize)
        }, error = function(e) {
          shiny::showNotification(
            paste("Discovery analysis failed:", conditionMessage(e)),
            type = "error", duration = 12)
          NULL
        })
      })

      discovery_result(result)
      if (!is.null(result)) {
        r <- result$result
        if (method == "GSEA") {
          gsea_h(.disc_h_gsea(nrow(r$data)))
        } else {
          ov_h(.disc_h_overall(nrow(r$Overall)))
          dist_h(.disc_h_dist(length(r$plot_distributions)))
          if (!is.null(r$Contrasts) && nrow(r$Contrasts) > 0) {
            cont_h(.disc_h_contrasts(nrow(r$Contrasts),
                                     length(unique(r$Contrasts$Variable))))
          }
        }
      }
    })

    # ── Renders: Status banner --------------------------------------------------
    output$disc_status_ui <- shiny::renderUI({
      if (is.null(discovery_result())) {
        shiny::div(
          class = "alert alert-info", style = "margin:10px 0;",
          shiny::icon("circle-info"),
          " Select a gene set and the variables to test, then click ",
          shiny::strong("Run Selected Method"), "."
        )
      }
    })

    # ── Card: Effect size overview (score methods only) ------------------------
    output$disc_overall_card_ui <- shiny::renderUI({
      res_wrap <- discovery_result()
      shiny::req(!is.null(res_wrap), res_wrap$method != "GSEA")
      bslib::card(
        full_screen = TRUE,
        .bm_card_header("Effect Size Overview (Cohen's f per variable)", "dl_disc_overall", ns),
        shiny::plotOutput(ns("disc_overall_plot"), height = paste0(ov_h(), "px"))
      )
    })

    output$disc_overall_plot <- shiny::renderPlot({
      res_wrap <- discovery_result(); shiny::req(!is.null(res_wrap), res_wrap$method != "GSEA")
      tryCatch(
        print(res_wrap$result$plot_overall +
                .disc_font_theme(res_wrap$labsize, res_wrap$titlesize)),
        error = function(e)
          shiny::showNotification(paste("Plot error:", e$message), type = "warning", duration = 8))
    }, height = function() ov_h(), res = 150)

    output$dl_disc_overall <- shiny::downloadHandler(
      filename = function() paste0("discovery_effect_overview_", Sys.Date(), ".png"),
      content  = function(file) {
        res_wrap <- discovery_result(); shiny::req(!is.null(res_wrap))
        ggplot2::ggsave(file, plot = res_wrap$result$plot_overall +
                          .disc_font_theme(res_wrap$labsize, res_wrap$titlesize),
                        width = 9, height = max(4, ov_h() / 96), dpi = 150, units = "in")
      }
    )

    # ── Card: Score distributions (score methods only) --------------------------
    output$disc_dist_card_ui <- shiny::renderUI({
      res_wrap <- discovery_result()
      shiny::req(!is.null(res_wrap), res_wrap$method != "GSEA")
      bslib::card(
        full_screen = TRUE,
        .bm_card_header("Score Distributions per Variable", "dl_disc_dist", ns),
        shiny::plotOutput(ns("disc_dist_plot"), height = paste0(dist_h(), "px"))
      )
    })

    .disc_dist_combined <- function(res_wrap) {
      plist <- res_wrap$result$plot_distributions
      plist <- lapply(plist, function(p) p + .disc_font_theme(res_wrap$labsize, res_wrap$titlesize))
      n <- length(plist)
      ncol_d <- max(1L, min(.DISC_DIST_NCOL, n))
      nrow_d <- ceiling(n / ncol_d)
      ggpubr::ggarrange(plotlist = plist, ncol = ncol_d, nrow = nrow_d, align = "hv")
    }

    output$disc_dist_plot <- shiny::renderPlot({
      res_wrap <- discovery_result(); shiny::req(!is.null(res_wrap), res_wrap$method != "GSEA")
      tryCatch(print(.disc_dist_combined(res_wrap)),
               error = function(e)
                 shiny::showNotification(paste("Plot error:", e$message), type = "warning", duration = 8))
    }, height = function() dist_h(), res = 150)

    output$dl_disc_dist <- shiny::downloadHandler(
      filename = function() paste0("discovery_score_distributions_", Sys.Date(), ".png"),
      content  = function(file) {
        res_wrap <- discovery_result(); shiny::req(!is.null(res_wrap))
        ggplot2::ggsave(file, plot = .disc_dist_combined(res_wrap),
                        width = 12, height = max(4, dist_h() / 96), dpi = 150, units = "in")
      }
    )

    # ── Card: Pairwise contrasts (score methods, categorical variables only) ---
    output$disc_contrasts_card_ui <- shiny::renderUI({
      res_wrap <- discovery_result()
      shiny::req(!is.null(res_wrap), res_wrap$method != "GSEA",
                !is.null(res_wrap$result$Contrasts), nrow(res_wrap$result$Contrasts) > 0)
      bslib::card(
        full_screen = TRUE,
        .bm_card_header("Pairwise Contrasts (Cohen's d)", "dl_disc_contrasts", ns),
        shiny::plotOutput(ns("disc_contrasts_plot"), height = paste0(cont_h(), "px"))
      )
    })

    output$disc_contrasts_plot <- shiny::renderPlot({
      res_wrap <- discovery_result()
      shiny::req(!is.null(res_wrap), res_wrap$method != "GSEA", !is.null(res_wrap$result$plot_contrasts))
      tryCatch(
        print(res_wrap$result$plot_contrasts +
                .disc_font_theme(res_wrap$labsize, res_wrap$titlesize)),
        error = function(e)
          shiny::showNotification(paste("Plot error:", e$message), type = "warning", duration = 8))
    }, height = function() cont_h(), res = 150)

    output$dl_disc_contrasts <- shiny::downloadHandler(
      filename = function() paste0("discovery_contrasts_", Sys.Date(), ".png"),
      content  = function(file) {
        res_wrap <- discovery_result(); shiny::req(!is.null(res_wrap))
        ggplot2::ggsave(file, plot = res_wrap$result$plot_contrasts +
                          .disc_font_theme(res_wrap$labsize, res_wrap$titlesize),
                        width = 9, height = max(4, cont_h() / 96), dpi = 150, units = "in")
      }
    )

    # ── Card: GSEA lollipop (GSEA method only) ----------------------------------
    output$disc_gsea_card_ui <- shiny::renderUI({
      res_wrap <- discovery_result()
      shiny::req(!is.null(res_wrap), res_wrap$method == "GSEA")
      bslib::card(
        full_screen = TRUE,
        .bm_card_header("GSEA Enrichment Across Contrasts", "dl_disc_gsea", ns),
        shiny::div(
          class = "text-muted", style = "font-size:0.78em; padding:4px 6px 8px;",
          shiny::icon("circle-info"),
          " One row per contrast, across every selected variable.",
          " Dashed lollipops/grey points mark gene sets ranked by the",
          " B-statistic (no direction) with negative NES (generally uninformative)."
        ),
        shiny::plotOutput(ns("disc_gsea_plot"), height = paste0(gsea_h(), "px"))
      )
    })

    output$disc_gsea_plot <- shiny::renderPlot({
      res_wrap <- discovery_result(); shiny::req(!is.null(res_wrap), res_wrap$method == "GSEA")
      tryCatch(
        print(res_wrap$result$plot + .disc_font_theme(res_wrap$labsize, res_wrap$titlesize)),
        error = function(e)
          shiny::showNotification(paste("Plot error:", e$message), type = "warning", duration = 8))
    }, height = function() gsea_h(), res = 150)

    output$dl_disc_gsea <- shiny::downloadHandler(
      filename = function() paste0("discovery_gsea_lollipop_", Sys.Date(), ".png"),
      content  = function(file) {
        res_wrap <- discovery_result(); shiny::req(!is.null(res_wrap))
        ggplot2::ggsave(file, plot = res_wrap$result$plot +
                          .disc_font_theme(res_wrap$labsize, res_wrap$titlesize),
                        width = 9, height = max(4, gsea_h() / 96), dpi = 150, units = "in")
      }
    )

    # ── Run: Method Summary (all four methods, same variables/gene set) --------
    shiny::observeEvent(input$run_summary, {
      expr <- get_expr(); meta <- get_meta(); gs <- get_gene_sets()
      s    <- .disc_settings()
      shiny::req(expr, meta, gs, s$gs, s$vars)

      gene_set_single <- gs[s$gs]
      score_methods   <- c("logmedian" = "Log2-median", "ranking" = "Ranking",
                           "ssGSEA" = "ssGSEA")

      shiny::withProgress(message = "Building method summary...", value = 0, {
        result <- tryCatch({

          score_rows <- lapply(names(score_methods), function(m) {
            shiny::incProgress(0.15, detail = paste("Scoring:", score_methods[[m]]))
            res_m <- Score_VariableAssociation(
              data = expr, metadata = meta, cols = s$vars, method = m,
              gene_set = gene_set_single, mode = s$mode,
              sig_threshold = s$sig_thr, p.adjust.method = s$padj_m,
              printplt = FALSE
            )
            ov <- res_m$Overall
            data.frame(Variable = ov$Variable, Method = unname(score_methods[m]),
                      Effect = ov$Cohen_f, PValue = ov$P_Value,
                      stringsAsFactors = FALSE)
          })
          score_df <- do.call(rbind, score_rows)

          # For GSEA, a variable can involve several contrasts (e.g. a 3-level
          # categorical variable tested pairwise has 3). Rather than picking
          # one arbitrary "strongest" contrast, summarise as: how many of
          # that variable's contrasts were significant, out of how many were
          # tested - comparable across variables and doesn't hide behind a
          # single potentially-unrepresentative NES value.
          gsea_pairs <- lapply(s$vars, function(v) {
            shiny::incProgress(0.55 / length(s$vars), detail = paste("GSEA:", v))
            tryCatch({
              res_g <- GSEA_VariableAssociation(
                data = expr, metadata = meta, cols = v, gene_set = gene_set_single,
                mode = s$mode, sig_threshold = s$sig_thr,
                p.adjust.method = s$padj_m, printplt = FALSE
              )
              # GSEA_VariableAssociation()$data is (or is derived from) a
              # data.table (fgsea's own output type). A data.table's `[`
              # treats a single unnamed argument as a row-filter *expression*
              # evaluated in its own scope, not a column-selecting vector like
              # a data.frame's `[` does - so d[, detail_cols, drop = FALSE]
              # below would try to find a literal column called "detail_cols"
              # instead of using the variable's value. Coerce to a plain
              # data.frame first to get ordinary `[` semantics.
              d <- as.data.frame(res_g$data)
              n_total <- nrow(d)
              n_sig   <- sum(d$padj < s$sig_thr, na.rm = TRUE)
              summary_row <- data.frame(Variable = v, n_sig = n_sig, n_total = n_total,
                                        frac_sig = n_sig / n_total, stringsAsFactors = FALSE)
              # Keep the underlying per-contrast values too (contrast tested,
              # its NES, and its adjusted p-value), not just the aggregate
              # fraction-significant summary above.
              detail_cols <- intersect(c("Contrast", "pathway", "NES", "padj", "stat_used"),
                                       colnames(d))
              detail_row <- cbind(Variable = v, d[, detail_cols, drop = FALSE])
              list(summary = summary_row, detail = detail_row)
            }, error = function(e) {
              shiny::showNotification(
                paste0("GSEA summary skipped for '", v, "': ", conditionMessage(e)),
                type = "warning", duration = 8)
              NULL
            })
          })
          gsea_pairs <- Filter(Negate(is.null), gsea_pairs)
          gsea_df        <- do.call(rbind, lapply(gsea_pairs, `[[`, "summary"))
          gsea_detail_df <- do.call(rbind, lapply(gsea_pairs, `[[`, "detail"))

          list(score_df = score_df, gsea_df = gsea_df, gsea_detail_df = gsea_detail_df,
              vars = s$vars, gs_name = s$gs, sig_thr = s$sig_thr,
              labsize = s$labsize, titlesize = s$titlesize)

        }, error = function(e) {
          shiny::showNotification(paste("Method summary failed:", conditionMessage(e)),
                                  type = "error", duration = 12)
          NULL
        })
      })

      summary_result(result)
      if (!is.null(result)) sum_h(.disc_h_summary(length(result$vars)))
    })

    output$disc_summary_status_ui <- shiny::renderUI({
      if (is.null(summary_result())) {
        shiny::div(
          class = "alert alert-info", style = "margin:10px 0;",
          shiny::icon("circle-info"),
          " Click ", shiny::strong("Build Method Summary"),
          " in the sidebar to compare all four methods for the selected variables."
        )
      }
    })

    # labsize/titlesize are taken as explicit arguments (read live from
    # input$ by the callers below) rather than from `res`, so changing the
    # Display-option font sliders re-styles this plot instantly without
    # re-running the (expensive) per-method scoring/GSEA computation.
    .build_disc_summary_plot <- function(res, labsize = 8, titlesize = 8) {
      vars_order <- rev(res$vars)

      score_df <- res$score_df
      score_df$Variable <- factor(score_df$Variable, levels = vars_order)
      score_df$Method   <- factor(score_df$Method,
                                  levels = c("Log2-median", "Ranking", "ssGSEA"))

      p_score <- ggplot2::ggplot(score_df,
          ggplot2::aes(x = .data$Method, y = .data$Variable, fill = .data$Effect)) +
        ggplot2::geom_tile(color = "white", linewidth = 0.6) +
        ggplot2::geom_text(
          ggplot2::aes(label = ifelse(.data$PValue < res$sig_thr,
                                      sprintf("%.2f*", .data$Effect),
                                      sprintf("%.2f", .data$Effect))),
          size = labsize / 3, color = "black") +
        ggplot2::scale_fill_distiller(palette = "BuGn", direction = 1,
                                      name = "Cohen's f") +
        ggplot2::labs(x = NULL, y = NULL) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1, size = labsize),
          axis.text.y  = ggplot2::element_text(size = labsize),
          legend.position = "bottom",
          legend.text  = ggplot2::element_text(size = labsize),
          legend.title = ggplot2::element_text(size = labsize),
          panel.grid = ggplot2::element_blank())

      p_gsea <- NULL
      if (!is.null(res$gsea_df) && nrow(res$gsea_df) > 0) {
        gsea_df <- res$gsea_df
        gsea_df$Variable <- factor(gsea_df$Variable, levels = vars_order)
        gsea_df$Method   <- "GSEA"

        # A variable can involve several contrasts (e.g. a 3-level categorical
        # variable tested pairwise has 3); fill/label show what fraction of
        # those contrasts were significant, rather than any single one.
        p_gsea <- ggplot2::ggplot(gsea_df,
            ggplot2::aes(x = .data$Method, y = .data$Variable, fill = .data$frac_sig)) +
          ggplot2::geom_tile(color = "white", linewidth = 0.6) +
          ggplot2::geom_text(
            ggplot2::aes(label = sprintf("%d/%d", .data$n_sig, .data$n_total)),
            size = labsize / 3, color = "black") +
          ggplot2::scale_fill_distiller(palette = "RdPu", direction = 1,
                                        limits = c(0, 1),
                                        name = "Fraction of\ncontrasts significant") +
          ggplot2::labs(x = NULL, y = NULL) +
          ggplot2::theme_bw() +
          ggplot2::theme(
            axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1, size = labsize),
            axis.text.y  = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            legend.position = "bottom",
            legend.text  = ggplot2::element_text(size = labsize),
            legend.title = ggplot2::element_text(size = labsize),
            panel.grid = ggplot2::element_blank())
      }

      if (!is.null(p_gsea)) {
        combined <- ggpubr::ggarrange(p_score, p_gsea, ncol = 2, widths = c(0.72, 0.28), align = "h")
      } else {
        combined <- p_score
      }
      ggpubr::annotate_figure(combined,
        top = grid::textGrob(
          paste0("Method summary - ", res$gs_name, "  (* = adjusted p-value < ", res$sig_thr, ")"),
          gp = grid::gpar(cex = 1.1, fontsize = titlesize)))
    }

    output$disc_summary_card_ui <- shiny::renderUI({
      res <- summary_result(); shiny::req(!is.null(res))
      bslib::card(
        full_screen = TRUE,
        .bm_card_header("All-Methods Summary", "dl_disc_summary", ns),
        shiny::div(
          class = "text-muted", style = "font-size:0.78em; padding:4px 6px 8px;",
          shiny::icon("circle-info"),
          " Left block: Cohen's f from each score-based method (* = adjusted",
          " p-value below threshold). Right column: fraction of that",
          " variable's GSEA contrasts that were significant - e.g. a 3-level",
          " categorical variable tested pairwise has 3 contrasts; \"2/3\"",
          " means 2 of them passed the significance threshold."
        ),
        shiny::plotOutput(ns("disc_summary_plot"), height = paste0(sum_h(), "px"))
      )
    })

    output$disc_summary_plot <- shiny::renderPlot({
      res <- summary_result(); shiny::req(!is.null(res))
      lab_s   <- input$disc_labsize   %||% 8
      title_s <- input$disc_titlesize %||% 8
      tryCatch(print(.build_disc_summary_plot(res, lab_s, title_s)),
               error = function(e)
                 shiny::showNotification(paste("Summary plot error:", e$message),
                                         type = "warning", duration = 8))
    }, height = function() sum_h(), res = 150)

    output$dl_disc_summary <- shiny::downloadHandler(
      filename = function() paste0("discovery_method_summary_", Sys.Date(), ".png"),
      content  = function(file) {
        res <- summary_result(); shiny::req(!is.null(res))
        lab_s   <- input$disc_labsize   %||% 8
        title_s <- input$disc_titlesize %||% 8
        ggplot2::ggsave(file, plot = .build_disc_summary_plot(res, lab_s, title_s),
                        width = 10, height = max(4, sum_h() / 96), dpi = 150, units = "in")
      }
    )

    output$disc_summary_table_score <- DT::renderDT({
      res <- summary_result(); shiny::req(!is.null(res), !is.null(res$score_df))
      df <- res$score_df
      df$Effect <- round(df$Effect, 4)
      df$PValue <- round(df$PValue, 4)
      DT::datatable(df, rownames = FALSE, filter = "top",
                    options = list(pageLength = 10, scrollX = TRUE))
    })

    output$disc_summary_table_gsea <- DT::renderDT({
      res <- summary_result(); shiny::req(!is.null(res), !is.null(res$gsea_df))
      df <- res$gsea_df
      df$frac_sig <- round(df$frac_sig, 3)
      DT::datatable(df, rownames = FALSE, filter = "top",
                    options = list(pageLength = 10, scrollX = TRUE))
    })

    output$disc_summary_table_gsea_detail <- DT::renderDT({
      res <- summary_result(); shiny::req(!is.null(res), !is.null(res$gsea_detail_df))
      df <- res$gsea_detail_df
      if ("NES" %in% colnames(df))  df$NES  <- round(df$NES, 4)
      if ("padj" %in% colnames(df)) df$padj <- round(df$padj, 4)
      DT::datatable(df, rownames = FALSE, filter = "top", extensions = "Buttons",
                    options = list(pageLength = 15, dom = "Bfrtip",
                                  buttons = c("csv", "excel"), scrollX = TRUE))
    })

    # ── Renders: Results table ---------------------------------------------------
    output$disc_table_status_ui <- shiny::renderUI({
      res_wrap <- discovery_result()
      if (is.null(res_wrap)) {
        return(shiny::div(class = "alert alert-info",
                          "Run Discovery to populate this table."))
      }
      method_label <- switch(res_wrap$method %||% "",
        "logmedian" = "Log-median score", "ranking" = "Ranking score",
        "ssGSEA"    = "ssGSEA score",     "GSEA"    = "GSEA enrichment",
        res_wrap$method %||% "unknown")
      vars_str <- paste(res_wrap$vars %||% "-", collapse = ", ")
      gs_str   <- res_wrap$gs_name %||% "-"
      sig_str  <- if (!is.null(res_wrap$sig_thr)) as.character(res_wrap$sig_thr) else "0.05"
      shiny::div(
        style = "margin-bottom:6px; font-size:0.85em;",
        shiny::icon("circle-info", style = "color:#0d6efd; margin-right:4px;"),
        shiny::tags$span(style = "color:#555; margin-right:6px;", "Results filtered by:"),
        shiny::tags$span(class = "badge bg-primary me-1",   style = "font-weight:500;",
          paste0("Gene set: ", gs_str)),
        shiny::tags$span(class = "badge bg-secondary me-1", style = "font-weight:500;",
          paste0("Method: ", method_label)),
        shiny::tags$span(class = "badge bg-secondary me-1", style = "font-weight:500;",
          paste0("Variables: ", vars_str)),
        shiny::tags$span(class = "badge bg-warning text-dark me-1", style = "font-weight:500;",
          paste0("Sig. threshold (adj. p-value): ", sig_str))
      )
    })

    output$disc_table <- DT::renderDT({
      res_wrap <- discovery_result()
      shiny::req(!is.null(res_wrap))

      res <- res_wrap$result

      df <- NULL
      if (res_wrap$method == "GSEA") {
        df <- res$data
      } else if (!is.null(res$Contrasts) && nrow(res$Contrasts) > 0) {
        df <- res$Contrasts
      } else {
        df <- res$Overall
      }

      shiny::req(!is.null(df))

      num_cols <- sapply(df, is.numeric)
      df[num_cols] <- lapply(df[num_cols], round, digits = 4)

      sig_thr <- res_wrap$sig_thr
      pval_col_match <- intersect(c("P_Value", "PValue", "padj"), colnames(df))

      dt <- DT::datatable(
        df, rownames = FALSE, filter = "top", extensions = "Buttons",
        options = list(pageLength = 20, dom = "Bfrtip",
                      buttons = c("csv", "excel"), scrollX = TRUE)
      )
      if (length(pval_col_match) > 0) {
        dt <- dt |> DT::formatStyle(
          columns    = pval_col_match[1],
          target     = "row",
          backgroundColor = DT::styleInterval(sig_thr, c("#fff3cd", "white"))
        )
      }
      dt
    })

    # Note: the "Statistical Insights" card (gene set/method/variable summary,
    # significant-hits list, and the hypothesis-generation disclaimer) was
    # removed as a separate section - the disclaimer now lives as a prominent
    # alert at the top of the sidebar instead, since it's important enough to
    # see before running anything.

  })
}
