# =============================================================================
# Discovery Tab – Shiny Module
#
# Explores associations between all (selected) metadata variables and the gene
# set signal, using either scoring or GSEA methods.
#
# Primary functions:
#   VariableAssociation()  – umbrella wrapper (Score or GSEA path)
#   Score_VariableAssociation() – score path → lollipop per variable
#   GSEA_VariableAssociation()  – GSEA path → lollipop per variable
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

.dyn_height_disc <- function(n, base = 300, per_item = 40, max_h = 3000) {
  as.integer(min(max_h, base + n * per_item))
}

discoveryUI <- function(id) {
  ns <- shiny::NS(id)

  bslib::layout_sidebar(

    sidebar = bslib::sidebar(
      width = 310,

      shiny::div(
        style = "padding-bottom:10px;",
        shiny::h4("Discovery Mode", style = "font-weight:700; color:#1E497B; margin-bottom:4px;"),
        shiny::p(
          "Systematically test all metadata variables for association with",
          " the gene set signal.",
          " Use for hypothesis generation once a robust gene set has been established.",
          style = "color:#6c757d; font-size:0.85em;"
        ),
        shiny::hr()
      ),

      # ---- Gene set selector (single) ---------------------------------------
      shiny::uiOutput(ns("disc_gs_picker_ui")),

      # ---- Method -----------------------------------------------------------
      shiny::radioButtons(
        ns("disc_method"),
        "Association method:",
        choices  = c(
          "Log-median score" = "logmedian",
          "Ranking score"    = "ranking",
          "ssGSEA score"     = "ssGSEA",
          "GSEA enrichment"  = "GSEA"
        ),
        selected = "logmedian"
      ),

      # ---- Metadata variables to test ---------------------------------------
      shiny::uiOutput(ns("disc_vars_ui")),

      # ---- Contrast mode ----------------------------------------------------
      shiny::radioButtons(
        ns("disc_mode"),
        "Contrast mode (categorical variables):",
        choices  = c(
          "Simple pairwise"   = "simple",
          "One vs rest"       = "medium",
          "All combinations"  = "extensive"
        ),
        selected = "simple"
      ),

      # ---- Significance settings -------------------------------------------
      shiny::sliderInput(
        ns("disc_sig"),
        "Significance threshold (padj):",
        min   = 0.001, max = 0.2,
        value = 0.05,  step = 0.001
      ),

      shiny::selectInput(
        ns("disc_padj"),
        "p-adjust method:",
        choices  = c("BH", "bonferroni", "holm", "BY", "fdr", "none"),
        selected = "BH"
      ),

      shiny::numericInput(ns("disc_fontsize"), "Axis font size (pt):",
                          value = 10, min = 6, max = 20, step = 1),

      shiny::hr(),

      shiny::actionButton(
        ns("run_discovery"),
        "Run Discovery",
        class = "btn-primary",
        width  = "100%"
      )
    ),

    # ── Main panel ----------------------------------------------------------
    bslib::navset_card_tab(

      # ---- Tab 1: Association Plot ------------------------------------------
      bslib::nav_panel(
        "Variable Associations",

        shiny::uiOutput(ns("disc_status_ui")),
        shiny::uiOutput(ns("disc_main_plot_ui"))
      ),

      # ---- Tab 2: Results Table ---------------------------------------------
      bslib::nav_panel(
        "Results Table",

        bslib::card(
          full_screen = TRUE,
          bslib::card_header("Variable Association Results"),
          shiny::uiOutput(ns("disc_table_status_ui")),
          DT::DTOutput(ns("disc_table"))
        )
      ),

      # ---- Tab 3: Insights --------------------------------------------------
      bslib::nav_panel(
        "Insights",
        bslib::card(
          bslib::card_header("Statistical Insights"),
          shiny::uiOutput(ns("disc_insights_ui"))
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

    # ── Cached result --------------------------------------------------------
    discovery_result <- shiny::reactiveVal(NULL)
    disc_h           <- shiny::reactiveVal(500L)

    # ── Dynamic UI: gene set selector ----------------------------------------
    output$disc_gs_picker_ui <- shiny::renderUI({
      gs <- get_gene_sets()
      shiny::req(!is.null(gs), length(gs) > 0)
      shinyWidgets::pickerInput(
        ns("disc_gs"),
        "Gene set to test:",
        choices  = names(gs),
        selected = names(gs)[1],
        options  = shinyWidgets::pickerOptions(liveSearch = TRUE)
      )
    })

    # ── Dynamic UI: variable multi-selector ----------------------------------
    output$disc_vars_ui <- shiny::renderUI({
      meta <- get_meta()
      shiny::req(!is.null(meta))
      # Exclude first column (SampleID)
      cols <- setdiff(colnames(meta), colnames(meta)[1])
      shinyWidgets::pickerInput(
        ns("disc_vars"),
        "Metadata variables to test:",
        choices  = cols,
        selected = cols,          # all selected by default
        multiple = TRUE,
        options  = shinyWidgets::pickerOptions(
          actionsBox     = TRUE,
          liveSearch     = TRUE,
          selectedTextFormat = "count > 3"
        )
      )
    })

    # ── Run Discovery --------------------------------------------------------
    shiny::observeEvent(input$run_discovery, {
      expr <- get_expr()
      meta <- get_meta()
      gs   <- get_gene_sets()
      shiny::req(expr, meta, gs, input$disc_gs, input$disc_vars)

      selected_gs   <- shiny::isolate(input$disc_gs)
      selected_vars <- shiny::isolate(input$disc_vars)
      method        <- shiny::isolate(input$disc_method)
      mode          <- shiny::isolate(input$disc_mode)
      sig_thr       <- shiny::isolate(input$disc_sig)
      padj_m        <- shiny::isolate(input$disc_padj)

      # Subset to the single selected gene set
      gene_set_single <- gs[selected_gs]

      shiny::withProgress(
        message = "Running variable association analysis…",
        value   = 0,
        {
          result <- tryCatch({
            shiny::incProgress(0.2, detail = paste("Method:", method))

            res <- VariableAssociation(
              method          = method,
              data            = expr,
              metadata        = meta,
              cols            = selected_vars,
              gene_set        = gene_set_single,
              mode            = mode,
              sig_threshold   = sig_thr,
              p.adjust.method = padj_m,
              printplt        = FALSE
            )

            shiny::incProgress(0.8, detail = "Finalising…")
            list(
              result     = res,
              method     = method,
              vars       = selected_vars,
              gs_name    = selected_gs,
              sig_thr    = sig_thr
            )

          }, error = function(e) {
            shiny::showNotification(
              paste("Discovery analysis failed:", conditionMessage(e)),
              type = "error", duration = 12
            )
            NULL
          })
        }
      )

      discovery_result(result)

      # Update plot height once result is stored
      if (!is.null(result)) {
        n_cont <- length(result$vars) *
          if (result$method == "GSEA") 3L else 1L
        disc_h(.dyn_height_disc(n_cont, base = 300, per_item = 45))
      }
    })

    # ── Renders: Status banner -----------------------------------------------
    output$disc_status_ui <- shiny::renderUI({
      if (is.null(discovery_result())) {
        shiny::div(
          class = "alert alert-info",
          style = "margin:10px 0;",
          shiny::icon("circle-info"),
          " Select a gene set and the variables to test, then click ",
          shiny::strong("Run Discovery"), "."
        )
      }
    })

    # ── Dynamic plot container (height driven by disc_h reactiveVal)
    output$disc_main_plot_ui <- shiny::renderUI({
      shiny::req(!is.null(discovery_result()))
      shiny::tagList(
        shiny::h5("Association Plot", style = "margin:12px 0 4px;"),
        shiny::plotOutput(ns("disc_main_plot"), height = paste0(disc_h(), "px"))
      )
    })

    # ── Renders: Main association plot ---------------------------------------
    output$disc_main_plot <- shiny::renderPlot({
      res_wrap <- discovery_result()
      shiny::req(!is.null(res_wrap))

      res <- res_wrap$result


      font <- input$disc_fontsize
      tryCatch({
        p <- NULL
        if (is.list(res) && !inherits(res, "gg")) {
          p <- res$plot %||% res$plt %||%
            (if (length(res) > 0 && inherits(res[[1]], c("gg","ggplot"))) res[[1]])
        } else if (inherits(res, c("gg", "ggplot"))) {
          p <- res
        }
        if (!is.null(p) && inherits(p, c("gg", "ggplot"))) {
          print(p + ggplot2::theme(
            axis.text  = ggplot2::element_text(size = font),
            axis.title = ggplot2::element_text(size = font + 1),
            strip.text = ggplot2::element_text(size = font + 1, face = "bold"),
            legend.text = ggplot2::element_text(size = font)
          ))
        } else if (!is.null(p)) {
          print(p)
        } else {
          shiny::showNotification("Plot structure not recognised.",
                                  type = "warning", duration = 6)
        }
      }, error = function(e) {
        shiny::showNotification(paste("Plot rendering failed:", e$message),
                                type = "warning", duration = 8)
      })
    }, height = function() disc_h(), res = 150)

    # ── Renders: Results table -----------------------------------------------
    output$disc_table_status_ui <- shiny::renderUI({
      res_wrap <- discovery_result()
      if (is.null(res_wrap)) {
        return(shiny::div(class = "alert alert-info",
                          "Run Discovery to populate this table."))
      }
      # Show filter info for the last completed run
      method_label <- switch(res_wrap$method %||% "",
        "t.test"   = "t-test",
        "wilcox"   = "Wilcoxon",
        "lm"       = "Linear Model",
        "kruskal"  = "Kruskal-Wallis",
        res_wrap$method %||% "unknown")
      vars_str <- paste(res_wrap$vars %||% "—", collapse = ", ")
      gs_str   <- res_wrap$gs_name %||% "—"
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
          paste0("Sig. threshold (padj): ", sig_str))
      )
    })

    output$disc_table <- DT::renderDT({
      res_wrap <- discovery_result()
      shiny::req(!is.null(res_wrap))

      res <- res_wrap$result

      # Extract tabular data from the result list
      df <- NULL
      if (is.list(res) && !is.null(res$data)) {
        df <- res$data
      } else if (is.list(res) && !is.null(res$results)) {
        df <- res$results
      } else if (is.data.frame(res)) {
        df <- res
      } else if (is.list(res)) {
        # Try to bind any data frames found at the top level
        dfs <- Filter(is.data.frame, res)
        if (length(dfs) > 0) df <- do.call(rbind, dfs)
      }

      shiny::req(!is.null(df))

      num_cols <- sapply(df, is.numeric)
      df[num_cols] <- lapply(df[num_cols], round, digits = 4)

      # Highlight significant rows
      sig_thr <- res_wrap$sig_thr
      DT::datatable(
        df,
        rownames   = FALSE,
        filter     = "top",
        extensions = "Buttons",
        options    = list(
          pageLength = 20,
          dom        = "Bfrtip",
          buttons    = c("csv", "excel"),
          scrollX    = TRUE
        )
      ) |> DT::formatStyle(
        columns    = if ("P_Value" %in% colnames(df)) "P_Value"
                     else if ("pval" %in% colnames(df)) "pval"
                     else character(0),
        target     = "row",
        backgroundColor = DT::styleInterval(
          sig_thr,
          c("#fff3cd", "white")    # highlight rows with p < threshold
        )
      )
    })

    # ── Insights -------------------------------------------------------------
    output$disc_insights_ui <- shiny::renderUI({

      res_wrap <- discovery_result()

      if (is.null(res_wrap)) {
        return(shiny::p("Run Discovery to generate insights.",
                        style = "color:#6c757d;"))
      }

      res     <- res_wrap$result
      method  <- res_wrap$method
      vars    <- res_wrap$vars
      gs_name <- res_wrap$gs_name
      sig_thr <- res_wrap$sig_thr

      # Build insight text
      items <- list(
        shiny::h5("Discovery Analysis", style = "color:#1E497B;"),
        shiny::p(
          "Gene set: ", shiny::strong(gs_name), "  |  ",
          "Method: ", shiny::strong(method), "  |  ",
          shiny::strong(length(vars)), " variables tested."
        )
      )

      # Try to summarise significant hits
      sig_summary <- tryCatch({
        df <- NULL
        if (is.list(res) && !is.null(res$data)) df <- res$data
        else if (is.data.frame(res)) df <- res

        if (!is.null(df)) {
          pval_col <- intersect(c("P_Value", "pval", "padj"), colnames(df))[1]
          if (!is.na(pval_col)) {
            n_sig <- sum(df[[pval_col]] < sig_thr, na.rm = TRUE)
            top <- df[!is.na(df[[pval_col]]) & df[[pval_col]] < sig_thr, ]
            top <- top[order(top[[pval_col]]), ][seq_len(min(5, nrow(top))), ]

            list(
              shiny::p(
                shiny::strong(n_sig), " variable(s) showed significant
                association (padj < ", sig_thr, ")."
              ),
              if (nrow(top) > 0) {
                var_col <- intersect(c("Variable", "contrast", "Contrast"),
                                     colnames(top))[1]
                shiny::tags$ul(lapply(seq_len(nrow(top)), function(i) {
                  shiny::tags$li(
                    if (!is.na(var_col))
                      sprintf("%s  (p = %.3e)", top[[var_col]][i], top[[pval_col]][i])
                    else
                      sprintf("Row %d  (p = %.3e)", i, top[[pval_col]][i])
                  )
                }))
              }
            )
          }
        }
      }, error = function(e) NULL)

      if (!is.null(sig_summary)) items <- c(items, sig_summary)

      items <- c(items, list(
        shiny::hr(),
        shiny::p(
          style = "color:#6c757d; font-size:0.9em;",
          shiny::icon("circle-exclamation"),
          " Discovery Mode is intended for hypothesis generation. All findings
           should be validated with independent datasets. Treat associations
           as exploratory unless confirmed by orthogonal evidence."
        )
      ))

      do.call(shiny::tagList, items)
    })

  })
}
