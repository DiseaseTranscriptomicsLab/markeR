# =============================================================================
# Benchmarking Tab - Shiny Module
#
# Two primary modules + one diagnostic:
#   Score Analysis  - PlotScores / ROC_Scores / AUC_Scores
#   Enrichment      - calculateDE -> runGSEA -> plotGSEAenrichment /
#                     plotNESlollipop / plotCombinedGSEA
#   FPR Simulation  - FPR_Simulation (score-based Cohen's d vs null distribution)
#
# Data contract (reactive accessors from main server):
#   get_expr()      - data.frame: genes x samples (normalised counts)
#   get_meta()      - data.frame: samples x variables (first col = SampleID)
#   get_gene_sets() - named list: char vectors or 2-col data.frames
# =============================================================================

# ---- Helpers ----------------------------------------------------------------

.bm_cat_cols <- function(meta) {
  if (is.null(meta)) return(character(0))
  cols <- setdiff(colnames(meta), colnames(meta)[1])
  cols[sapply(meta[cols], function(x) is.character(x) || is.factor(x))]
}

# Contrast-mode tooltip (reused across tabs)
.contrast_mode_tooltip <- bslib::tooltip(
  shiny::icon("circle-question", style = "color:#aaa; cursor:help; font-size:0.85em;"),
  shiny::tags$div(
    shiny::tags$b("Simple pairwise:"),
    " all group-vs-group comparisons (A vs B, A vs C, B vs C ...).",
    shiny::tags$br(),
    shiny::tags$b("One vs rest:"),
    " each group compared to the union of all others (A vs B+C+D ...).",
    shiny::tags$br(),
    shiny::tags$b("All combinations:"),
    " every algebraic combination, including multi-group contrasts like (A+B) vs (C+D).",
    shiny::tags$br(),
    shiny::tags$em(
      "Tip: start with Simple pairwise. Use broader modes only if you need",
      " to capture more complex group structures."
    )
  ),
  placement = "right"
)

# Scoring method tooltip
.method_tooltip <- bslib::tooltip(
  shiny::icon("circle-question", style = "color:#aaa; cursor:help; font-size:0.85em;"),
  shiny::tags$div(
    shiny::tags$b("Log-median:"),
    " log2-transforms expression and median-centres each gene across samples.",
    " Fast and interpretable; score = mean of centred values.",
    shiny::tags$br(),
    shiny::tags$b("Ranking:"),
    " non-parametric; scores based on within-sample rank of gene set members.",
    shiny::tags$br(),
    shiny::tags$b("ssGSEA:"),
    " single-sample enrichment score (Barbie et al. 2009). More sensitive to gene set size."
  ),
  placement = "right"
)

# Contrast mode radio (reused in score + fpr)
.bm_contrast_mode_radio <- function(ns, input_id, selected = "simple") {
  shiny::div(
    shiny::tags$label(
      "Contrast mode ", .contrast_mode_tooltip,
      style = "font-size:0.85em; font-weight:600;"
    ),
    shiny::radioButtons(ns(input_id), label = NULL,
      choices  = c("Simple pairwise"  = "simple",
                   "One vs rest"       = "medium",
                   "All combinations"  = "extensive"),
      selected = selected, inline = TRUE)
  )
}

# Collapsible settings wrapper (starts open; user can collapse after configuring)
.bm_collapsible_settings <- function(...) {
  shiny::tags$details(
    open = "open",
    shiny::tags$summary(
      style = paste("cursor:pointer; font-size:0.82em; font-weight:600; color:#555;",
                    "user-select:none; padding:4px 2px 8px; list-style:none;",
                    "display:flex; align-items:center; gap:5px;"),
      shiny::icon("sliders", style = "color:#2D6A4F; font-size:0.9em;"),
      " Settings"
    ),
    shiny::div(
      style = paste("background:#f8f9fa; border:1px solid #dee2e6; border-radius:6px;",
                    "padding:12px 14px; margin-bottom:14px;"),
      ...
    )
  )
}

# Card header with download button
.bm_card_header <- function(label, dl_btn_id, ns) {
  bslib::card_header(
    class = "d-flex align-items-center justify-content-between",
    label,
    shiny::downloadButton(ns(dl_btn_id), "",
      icon  = shiny::icon("download"),
      class = "btn-sm btn-outline-secondary py-0 px-2",
      style = "font-size:0.75em;")
  )
}

# ---- UI ---------------------------------------------------------------------

#' @importFrom bslib layout_sidebar sidebar navset_card_tab navset_tab nav_panel tooltip
#' @importFrom shiny NS radioButtons selectInput numericInput checkboxInput
#'   actionButton uiOutput div tags conditionalPanel downloadButton
benchmarkingUI <- function(id) {
  ns <- shiny::NS(id)

  bslib::layout_sidebar(

    # =========================================================================
    # SIDEBAR
    # =========================================================================
    sidebar = bslib::sidebar(
      width = 290,

      shiny::div(
        style = "padding-bottom:6px;",
        # Matches the "markeR" wordmark colour used for this section on the About
        # page (black), distinct from Discovery Mode's blue - the two used to
        # share the same blue, making them hard to tell apart at a glance.
        shiny::h4("Benchmarking Mode", style = "font-weight:700; color:#020201; margin-bottom:4px;"),
        shiny::div(
          style = "color:#6c757d; font-size:0.85em; line-height:1.5;",
          shiny::tags$b("Score Analysis:"),
          " computes per-sample gene set scores (log-median, ranking, ssGSEA) and",
          " tests how well they separate groups of a chosen metadata variable, via",
          " score distributions, ROC/AUC classification, a combined Cohen's d/f",
          " summary across all three methods, and a null-distribution FPR check.",
          shiny::tags$br(), shiny::tags$br(),
          shiny::tags$b("Enrichment Analysis (GSEA):"),
          " uses limma linear models to rank genes per contrast, then fgsea tests",
          " whether each gene set is coordinately shifted towards the top or bottom",
          " of that ranking, reporting a Normalised Enrichment Score (NES) per",
          " gene set × contrast."
        ),
        shiny::hr()
      ),

      shiny::tags$label("Gene sets:", style = "font-size:0.85em; font-weight:600;"),
      shiny::radioButtons(ns("bm_gs_scope"), label = NULL,
        choices  = c("All loaded" = "all", "Select specific" = "selected"),
        selected = "all", inline = TRUE),
      shiny::conditionalPanel(
        condition = sprintf("input['%s'] === 'selected'", ns("bm_gs_scope")),
        shiny::uiOutput(ns("bm_gs_picker_ui"))
      )
    ),

    # =========================================================================
    # MAIN PANEL
    # =========================================================================
    bslib::navset_card_tab(

      # =======================================================================
      # TAB 1: SCORE ANALYSIS
      # =======================================================================
      bslib::nav_panel(
        .tab_title("Score Analysis",
          shiny::tagList(
            "Computes per-sample gene set scores (log-median, ranking, ssGSEA)",
            " and visualises them against a metadata variable.",
            shiny::tags$br(),
            shiny::tags$b("Score Distributions:"),
            " violin / scatter plots of scores per group for a chosen method. Fast.",
            shiny::tags$br(),
            shiny::tags$b("Classification (ROC/AUC):"),
            " ROC curves and AUC heatmap showing how well scores separate groups. Medium speed.",
            shiny::tags$br(),
            shiny::tags$b("All Methods Summary:"),
            " Cohen's d/f heatmap and volcano comparing all three methods at once.",
            shiny::tags$br(),
            shiny::tags$b("Score FPR Simulation:"),
            " compares observed Cohen's d to a null distribution from random gene sets.",
            shiny::tags$br(),
            shiny::tags$em(
              "Run each sub-analysis independently.",
              " Variable and contrast settings below apply to all sub-analyses.")
          )
        ),

        # ---- Shared settings (variable + contrast + contrast filter) ----------
        .bm_collapsible_settings(
          shiny::uiOutput(ns("score_var_ui")),
          shiny::uiOutput(ns("score_cond_ui")),
          shiny::uiOutput(ns("score_contrasts_ui"))
        ),

        # ---- Sub-tabs -------------------------------------------------------
        bslib::navset_tab(

          # -- A: Score Distributions -----------------------------------------
          bslib::nav_panel(
            "Score Distributions",
            .bm_collapsible_settings(
              shiny::tags$label(
                "Scoring method ", .method_tooltip,
                style = "font-size:0.85em; font-weight:600;"
              ),
              shiny::radioButtons(ns("score_method"), label = NULL,
                choices  = c("Log-median" = "logmedian",
                             "Ranking"    = "ranking",
                             "ssGSEA"     = "ssGSEA"),
                selected = "logmedian", inline = TRUE),
              shiny::uiOutput(ns("score_color_var_ui")),
              shiny::tags$details(
                style = "margin-bottom:8px;",
                shiny::tags$summary(
                  style = "cursor:pointer; font-size:0.85em; color:#555; user-select:none;",
                  shiny::icon("chevron-right", style = "font-size:0.7em;"),
                  " Advanced options"),
                shiny::div(
                  style = "padding:8px 0 4px;",
                  shiny::checkboxInput(ns("connect_groups"),
                    "Connect paired samples across groups", value = FALSE),
                  shiny::tags$small(
                    style = "color:#888; display:block; margin:-4px 0 6px; font-size:0.78em;",
                    "Draws lines between the same sample across groups.",
                    " Set ", shiny::tags$b("Colour by"), " to the sample / patient identifier",
                    " so each line connects the correct pair."),
                  shiny::tags$hr(style = "margin:6px 0;"),
                  shiny::uiOutput(ns("score_cohend_ui")),
                  shiny::tags$hr(style = "margin:6px 0;"),
                  shiny::checkboxInput(ns("pvalcalc"),
                    "Show t-test p-value in subtitle", value = FALSE),
                  shiny::tags$small(
                    style = "color:#888; display:block; margin:-4px 0 6px; font-size:0.78em;",
                    "Runs a t-test between Group A and Group B (set in the contrast settings above",
                    " under 'Focus Cohen's d') and adds the p-value to the plot subtitle.",
                    " Groups A and B must be specified.")
                )
              ),
              shiny::tags$details(
                style = "margin-bottom:10px;",
                shiny::tags$summary(
                  style = "cursor:pointer; font-size:0.85em; color:#555; user-select:none;",
                  shiny::icon("chevron-right", style = "font-size:0.7em;"),
                  " Display options"),
                shiny::div(
                  style = "padding:8px 0 4px;",
                  shiny::div(
                    style = "display:grid; grid-template-columns:1fr 1fr; gap:8px;",
                    shiny::numericInput(ns("score_pointsize"), "Point size:",
                                        value = 4, min = 1, max = 10, step = 1),
                    shiny::div(),   # spacer
                    bslib::tooltip(
                      shiny::numericInput(ns("score_width_title"), "Title wrap (chars):",
                                          value = 32, min = 10, max = 60, step = 2),
                      "Maximum characters per line in plot titles. Long gene set names are broken",
                      " at underscores or capital letters to fit this width.",
                      placement = "right"
                    ),
                    bslib::tooltip(
                      shiny::numericInput(ns("score_width_legend"), "Legend wrap (chars):",
                                          value = 32, min = 10, max = 60, step = 2),
                      "Maximum characters per line in legend entries. Affects how gene set names",
                      " are displayed in the plot legend.",
                      placement = "right"
                    ),
                    shiny::numericInput(ns("score_titlesize"), "Title size (pt):",
                                        value = 8, min = 6, max = 30, step = 1),
                    shiny::numericInput(ns("score_labsize"), "Label size (pt):",
                                        value = 8, min = 6, max = 26, step = 1)
                  )
                )
              ),
              shiny::actionButton(ns("run_scores"), "Plot Score Distributions",
                                  class = "btn-primary btn-sm", width = "100%")
            ),
            shiny::uiOutput(ns("score_dist_status_ui")),
            shiny::uiOutput(ns("score_dist_ui"))
          ),

          # -- B: Classification (ROC/AUC) ------------------------------------
          bslib::nav_panel(
            "Classification (ROC / AUC)",
            .bm_collapsible_settings(
              shiny::div(
                class = "alert alert-info",
                style = "font-size:0.82em; padding:6px 10px; margin-bottom:10px;",
                shiny::icon("circle-info"),
                " Requires a categorical variable.",
                " Uses the same variable and contrast mode set above."
              ),
              shiny::tags$label(
                "Scoring method ", .method_tooltip,
                style = "font-size:0.85em; font-weight:600;"
              ),
              shiny::radioButtons(ns("roc_method"), label = NULL,
                choices  = c("Log-median" = "logmedian",
                             "Ranking"    = "ranking",
                             "ssGSEA"     = "ssGSEA"),
                selected = "logmedian", inline = TRUE),
              shiny::tags$small(
                style = "color:#888; font-size:0.78em; display:block; margin-bottom:8px;",
                "To compare all three methods in one view, use the",
                shiny::tags$b(" All Methods Summary"), " sub-tab."
              ),
              shiny::tags$details(
                style = "margin-bottom:10px;",
                shiny::tags$summary(
                  style = "cursor:pointer; font-size:0.85em; color:#555; user-select:none;",
                  shiny::icon("chevron-right", style = "font-size:0.7em;"),
                  " Display options"),
                shiny::div(
                  style = "padding:8px 0 4px;",
                  shiny::div(
                    style = "display:grid; grid-template-columns:1fr 1fr; gap:8px;",
                    bslib::tooltip(
                      shiny::numericInput(ns("roc_width_title"), "Title wrap (chars):",
                                          value = 32, min = 10, max = 60, step = 2),
                      "Maximum characters per line for gene set and contrast titles on each ROC/AUC panel.",
                      " Lower values wrap more aggressively; increase if names are being cut off.",
                      placement = "right"
                    ),
                    shiny::numericInput(ns("roc_titlesize"), "Title size (pt):",
                                        value = 8, min = 6, max = 20, step = 1),
                    shiny::numericInput(ns("roc_labsize"), "Label size (pt):",
                                        value = 8, min = 6, max = 20, step = 1)
                  ),
                  bslib::tooltip(
                    shiny::checkboxInput(ns("roc_restore_direction"),
                      "Ignore direction (fold AUC to >= 0.5)",
                      value = FALSE),
                    "Off by default: shows the ROC curves and AUC values exactly as computed,",
                    " which preserves directionality - a curve can dip below the diagonal",
                    " (AUC < 0.5), meaning the score runs opposite to the group comparison.",
                    " When checked, direction is discarded: any AUC below 0.5 is displayed",
                    " as 1 - AUC instead, so it only reflects separation strength, not direction",
                    " (the curve itself is unchanged either way).",
                    placement = "right"
                  )
                )
              ),
              shiny::actionButton(ns("run_roc"), "Compute ROC Curves and AUC",
                                  class = "btn-primary btn-sm", width = "100%")
            ),
            shiny::uiOutput(ns("score_roc_status_ui")),
            shiny::uiOutput(ns("score_roc_ui")),
            shiny::uiOutput(ns("score_auc_ui"))
          ),

          # -- C: All Methods Summary -----------------------------------------
          bslib::nav_panel(
            "All Methods Summary",
            .bm_collapsible_settings(
              shiny::div(
                style = "font-size:0.82em; color:#555; margin-bottom:10px;",
                "Runs all three scoring methods and summarises discriminability as a",
                shiny::tags$b(" Cohen's d/f heatmap"),
                " (gene sets × contrasts) and a",
                shiny::tags$b(" volcano plot"),
                " (|effect size| vs significance).",
                " Use to identify which gene sets and contrasts are most consistently",
                " discriminatory regardless of scoring method."
              ),
              shiny::div(
                style = "display:grid; grid-template-columns:1fr 1fr; gap:8px;",
                shiny::numericInput(ns("sig_threshold"), "Significance threshold:",
                                    value = 0.05, min = 0.001, max = 0.5, step = 0.005),
                shiny::numericInput(ns("cohen_threshold"), "Effect size threshold:",
                                    value = 0.5, min = 0, max = 5, step = 0.1)
              ),
              shiny::tags$details(
                style = "margin-bottom:10px;",
                shiny::tags$summary(
                  style = "cursor:pointer; font-size:0.85em; color:#555; user-select:none;",
                  shiny::icon("chevron-right", style = "font-size:0.7em;"),
                  " Display options"),
                shiny::div(
                  style = "padding:8px 0 4px;",
                  shiny::div(
                    style = "font-size:0.78em; color:#888; margin-bottom:6px;",
                    "All the font-size settings below apply instantly on change,",
                    " no need to re-run."
                  ),
                  shiny::div(
                    style = "display:grid; grid-template-columns:1fr 1fr; gap:8px;",
                    bslib::tooltip(
                      shiny::numericInput(ns("all_width_title"), "Title wrap (chars):",
                                          value = 32, min = 10, max = 60, step = 2),
                      "Maximum characters per line for gene set names in heatmap rows and",
                      " volcano axis labels.",
                      placement = "right"
                    ),
                    bslib::tooltip(
                      shiny::numericInput(ns("all_width_legend"), "Legend wrap (chars):",
                                          value = 32, min = 10, max = 60, step = 2),
                      "Maximum characters per line for gene set names in the volcano legend.",
                      placement = "right"
                    ),
                    bslib::tooltip(
                      shiny::numericInput(ns("all_heatmap_textsize"), "Heatmap cell text size (mm):",
                                          value = 2, min = 1, max = 8, step = 0.5),
                      "Font size of the Cohen's d/f and p-value labels drawn inside each",
                      " heatmap cell. In millimetres, not points, because that's the unit",
                      " ggplot2's geom_text() uses (roughly pt ÷ 2.85).",
                      placement = "right"
                    ),
                    bslib::tooltip(
                      shiny::numericInput(ns("all_heatmap_titlesize"), "Heatmap title size (pt):",
                                          value = 8, min = 4, max = 20, step = 1),
                      "Font size of each gene set's panel title in the Cohen's d/f heatmap.",
                      placement = "right"
                    ),
                    bslib::tooltip(
                      shiny::numericInput(ns("all_heatmap_labsize"), "Heatmap axis label size (pt):",
                                          value = 8, min = 4, max = 20, step = 1),
                      "Font size of the contrast labels on the heatmap's axes.",
                      " Defaults to the same size as the heatmap title above.",
                      placement = "right"
                    ),
                    bslib::tooltip(
                      shiny::numericInput(ns("all_volcano_titlesize"), "Volcano title size (pt):",
                                          value = 8, min = 4, max = 20, step = 1),
                      "Font size of each contrast's facet title in the Effect Size Volcano plot.",
                      placement = "right"
                    ),
                    bslib::tooltip(
                      shiny::numericInput(ns("all_volcano_labsize"), "Volcano label size (pt):",
                                          value = 8, min = 4, max = 20, step = 1),
                      "Font size of the axis/legend text in the Effect Size Volcano plot.",
                      placement = "right"
                    )
                  )
                )
              ),
              shiny::actionButton(ns("run_all_methods"), "Run All Methods Summary",
                                  class = "btn-primary btn-sm", width = "100%")
            ),
            shiny::uiOutput(ns("score_all_status_ui")),
            shiny::uiOutput(ns("score_all_heatmap_ui")),
            shiny::uiOutput(ns("score_all_volcano_ui"))
          ),

          # -- D: Score FPR Simulation ----------------------------------------
          bslib::nav_panel(
            fillable = FALSE,
            .tab_title("Score FPR Simulation",
              shiny::tagList(
                "Evaluates whether a gene set's score-based discriminability is specific",
                " or could arise from any random set of genes of the same size.",
                shiny::tags$br(),
                shiny::tags$b("How it works:"),
                " gene set scores are computed (log-median, ranking, ssGSEA)",
                " and the observed Cohen's d or f is compared",
                " to a null distribution from random gene sets of the same size,",
                " drawn from the expression matrix.",
                shiny::tags$br(),
                "Select a contrast from the dropdown to view the plot for that comparison.",
                shiny::tags$br(),
                shiny::tags$em(
                  "Note: this tests score-level specificity, not gene-level discriminability.",
                  " Use in conjunction with the Classification and Effect Size tabs in Gene Sets.")
              )
            ),
            .bm_collapsible_settings(
              .bm_contrast_mode_radio(ns, "fpr_mode"),
              shiny::numericInput(ns("fpr_nsims"), "Random simulations per gene set:",
                                  value = 100, min = 10, max = 2000, step = 10),
              shiny::tags$small(
                style = "color:#888; display:block; margin:-4px 0 8px; font-size:0.78em;",
                shiny::icon("clock"),
                " Runtime scales with (gene sets) x (simulations) x 3 methods.",
                " This can take several minutes or longer for many gene sets or simulations.",
                " Start with 100 simulations; increase only for publication-quality results.",
                shiny::tags$br(),
                shiny::icon("lightbulb"),
                " Tip: to quickly test a single gene set of interest, set the gene set selector",
                " above to 'Select specific' and choose just that one gene set before running."
              ),
              shiny::tags$details(
                style = "margin-bottom:10px;",
                shiny::tags$summary(
                  style = "cursor:pointer; font-size:0.85em; color:#555; user-select:none;",
                  "Display options"),
                shiny::div(
                  style = "padding:8px 0 4px;",
                  shiny::numericInput(ns("fpr_fontsize"), "Font size (pt):",
                                      value = 10, min = 6, max = 20, step = 1)
                )
              ),
              shiny::actionButton(ns("run_fpr"), "Run Score FPR Simulation",
                                  class = "btn-primary btn-sm", width = "100%")
            ),
            shiny::uiOutput(ns("fpr_status_ui")),
            shiny::uiOutput(ns("fpr_plot_ui"))
          ),

          # -- E: Score Results Table -------------------------------------------
          bslib::nav_panel(
            "Score Results Table",
            .bm_collapsible_settings(
              shiny::tags$label("Scoring method(s):", style = "font-size:0.85em; font-weight:600;"),
              shiny::radioButtons(ns("scores_table_method"), label = NULL,
                choices  = c("Log-median" = "logmedian",
                             "Ranking"    = "ranking",
                             "ssGSEA"     = "ssGSEA",
                             "All methods" = "all"),
                selected = "all", inline = TRUE),
              shiny::actionButton(ns("run_scores_table"), "Compute Score Results Table",
                                  class = "btn-primary btn-sm", width = "100%")
            ),
            shiny::uiOutput(ns("scores_table_status_ui")),
            shiny::uiOutput(ns("scores_table_ui"))
          )
        )
      ),

      # =======================================================================
      # TAB 2: ENRICHMENT ANALYSIS (GSEA)
      # =======================================================================
      bslib::nav_panel(
        .tab_title("Enrichment Analysis",
          shiny::tagList(
            shiny::tags$b("Gene Set Enrichment Analysis (GSEA)"),
            " uses limma linear models to derive per-gene ranking statistics",
            " (t- or B-statistic) for each contrast,",
            " then fgsea tests whether each gene set is coordinately shifted",
            " towards the top or bottom of the ranked list.",
            shiny::tags$br(),
            "The output is a ", shiny::tags$b("Normalised Enrichment Score (NES)"),
            " per gene set × contrast: positive NES indicates enrichment among",
            " up-regulated genes; negative NES indicates enrichment among",
            " down-regulated genes.",
            shiny::tags$br(),
            shiny::tags$em(
              "Results appear in sub-tabs: Combined Volcano · NES Lollipop",
              " · Enrichment Plots · Results Table.")
          )
        ),

        .bm_collapsible_settings(

          # ---- Design specification ------------------------------------------
          shiny::div(
            style = "margin-bottom:6px;",
            shiny::tags$label(
              style = "font-size:0.82em; font-weight:700; color:#374151;",
              "Design specification "
            ),
            bslib::tooltip(
              shiny::icon("circle-question",
                style = "color:#aaa; cursor:help; font-size:0.82em;"),
              shiny::tags$div(
                shiny::tags$b("Automatic:"), " Pick a grouping variable. markeR builds",
                " the design matrix and all contrasts automatically. Best for simple",
                " case/control or multi-group comparisons.",
                shiny::tags$hr(style = "margin:4px 0;"),
                shiny::tags$b("Manual:"), " Paste an R formula (e.g. ",
                shiny::tags$code(
                  style = "background:rgba(255,255,255,0.2); color:#f0f0f0; padding:1px 4px; border-radius:3px;",
                  "~ 0 + Condition + Batch"),
                ") and define contrasts explicitly. Use when you need covariates,",
                " paired designs, or non-standard comparisons."
              ),
              placement = "right"
            )
          ),
          shiny::radioButtons(ns("gsea_design_mode"), label = NULL,
            choices  = c("Automatic" = "auto", "Manual (custom formula)" = "custom"),
            selected = "auto", inline = TRUE),
          shiny::uiOutput(ns("enrich_var_ui")),
          shiny::uiOutput(ns("enrich_custom_design_ui")),
          shiny::conditionalPanel(
            sprintf("input['%s'] === 'auto'", ns("gsea_design_mode")),
            .bm_contrast_mode_radio(ns, "enrich_mode")
          ),
          shiny::uiOutput(ns("enrich_contrasts_ui")),
          shiny::tags$hr(style = "margin:8px 0;"),

          # ---- Statistical options (collapsed) --------------------------------
          shiny::tags$details(
            style = "margin-bottom:10px;",
            shiny::tags$summary(
              style = "cursor:pointer; font-size:0.85em; color:#555; user-select:none;",
              shiny::icon("chevron-right", style = "font-size:0.7em;"),
              " Statistical options"),
            shiny::div(
              style = "padding:8px 0 4px;",
              shiny::tags$label(
                "Ranking statistic ",
                bslib::tooltip(
                  shiny::icon("circle-question",
                               style = "color:#aaa; cursor:help; font-size:0.85em;"),
                  shiny::tags$div(
                    shiny::tags$b("Auto (recommended):"),
                    " uses the t-statistic for directional gene sets (data frames with",
                    " +1/−1 weights) and the B-statistic for simple gene lists.",
                    shiny::tags$br(),
                    shiny::tags$b("t-statistic:"),
                    " ranks genes from most up-regulated to most down-regulated.",
                    " Positive NES = gene set enriched among up-regulated genes.",
                    shiny::tags$br(),
                    shiny::tags$b("B-statistic:"),
                    " ranks by strength of differential evidence regardless of direction."
                  ),
                  placement = "right"
                ),
                style = "font-size:0.82em; font-weight:600;"
              ),
              shiny::radioButtons(ns("gsea_stat"), label = NULL,
                choices  = c("Auto (recommended)" = "auto",
                             "t-statistic"         = "t",
                             "B-statistic"         = "B"),
                selected = "auto", inline = TRUE),
              shiny::div(
                style = "display:grid; grid-template-columns:1fr 1fr; gap:8px; margin-top:6px;",
                shiny::numericInput(ns("gsea_nperm"), "Permutations:",
                                    value = 1000, min = 100, max = 100000, step = 100),
                shiny::selectInput(ns("padj_method"), "p-adjust method:",
                  choices  = c("BH", "bonferroni", "holm", "BY", "fdr", "none"),
                  selected = "BH")
              ),
              shiny::checkboxInput(ns("contrast_correction"),
                "Apply FDR correction across all contrasts", value = FALSE),
              shiny::tags$small(
                style = "color:#6c757d; display:block; margin:-4px 0 4px; font-size:0.78em;",
                "When checked, p-values are corrected pooling all contrasts together",
                " (stricter). When unchecked, correction is per contrast."
              )
            )
          ),

          shiny::actionButton(ns("run_gsea"), "Run Enrichment Analysis (GSEA)",
                              class = "btn-primary btn-sm", width = "100%")
        ),
        shiny::uiOutput(ns("gsea_status_ui")),
        shiny::uiOutput(ns("gsea_results_ui"))
      ),

    )
  )
}


# ---- Server -----------------------------------------------------------------

#' @importFrom ggpubr ggarrange
benchmarkingServer <- function(id, get_expr, get_meta, get_gene_sets) {

  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Cached results
    score_dist_result <- shiny::reactiveVal(NULL)
    score_roc_result  <- shiny::reactiveVal(NULL)
    score_all_result  <- shiny::reactiveVal(NULL)
    gsea_result       <- shiny::reactiveVal(NULL)
    fpr_result        <- shiny::reactiveVal(NULL)
    scores_table_result <- shiny::reactiveVal(NULL)

    # Clear any previously computed plots/results when the underlying
    # expression data or metadata actually changes (e.g. new file loaded, or
    # preprocessing re-finalised) - stale results from the old dataset should
    # not linger on screen.
    shiny::observeEvent(list(get_expr(), get_meta()), {
      score_dist_result(NULL); score_roc_result(NULL); score_all_result(NULL)
      gsea_result(NULL); fpr_result(NULL); scores_table_result(NULL)
    }, ignoreInit = TRUE)

    # Dynamic plot heights
    score_h  <- shiny::reactiveVal(500L)
    roc_h    <- shiny::reactiveVal(600L)
    auc_h    <- shiny::reactiveVal(400L)
    all_h    <- shiny::reactiveVal(500L)
    all_vol_h <- shiny::reactiveVal(480L)
    nes_h      <- shiny::reactiveVal(500L)
    enrich_h   <- shiny::reactiveVal(500L)
    gsea_vol_h <- shiny::reactiveVal(500L)
    combined_h <- shiny::reactiveVal(650L)
    gsea_vol_w <- shiny::reactiveVal("100%")  # DEG volcano width (narrow when few contrasts)
    fpr_h    <- shiny::reactiveVal(500L)

    # ---- Gene-set scope (shared across all tabs) ----------------------------

    output$bm_gs_picker_ui <- shiny::renderUI({
      gs <- get_gene_sets()
      shiny::req(!is.null(gs), length(gs) > 0)
      shinyWidgets::pickerInput(
        ns("bm_gs_picker"), label = NULL,
        choices  = names(gs), selected = names(gs)[1],
        multiple = TRUE,
        options  = shinyWidgets::pickerOptions(
          actionsBox = TRUE, liveSearch = TRUE,
          selectedTextFormat = "count > 2",
          container = "body",   # render dropdown in body, not sidebar (avoids overflow clip)
          countSelectedText  = "{0} gene sets selected"
        )
      )
    })

    # GSEA (auto mode) builds contrast strings like "A-B" directly from a
    # metadata group's raw values, then hands them (plus that same metadata
    # column) to calculateDE(), which uses them as design-matrix/limma
    # contrast level names. limma::makeContrasts() requires those to be
    # syntactically valid R names (see help(make.names)) - a value like
    # "cellstate:growing" contains a colon (the formula-interaction operator)
    # and crashes it with "levels must be syntactically valid names". Fixed
    # here in the app layer (not in calculateDE() itself, which is a
    # general-purpose exported package function) by sanitising the group's
    # values before they're ever used to build contrasts or passed to
    # calculateDE. Used identically by both the contrast-checkbox UI and the
    # actual GSEA run so the two always agree on the same (sanitised) names.
    .sanitize_group_vals <- function(x) make.names(gsub(" ", "", as.character(x)))

    .active_gs <- function() {
      gs    <- get_gene_sets()
      scope <- input$bm_gs_scope
      if (is.null(scope) || scope == "all") return(gs)
      sel <- input$bm_gs_picker
      if (length(sel) == 0) return(gs)
      gs[intersect(sel, names(gs))]
    }

    # ---- Score Analysis: shared variable UI --------------------------------

    output$score_var_ui <- shiny::renderUI({
      meta <- get_meta(); shiny::req(!is.null(meta))
      cols <- setdiff(colnames(meta), colnames(meta)[1])
      if (length(cols) == 0)
        return(shiny::helpText("No metadata variables found."))
      shiny::selectInput(ns("score_var"), "Variable:",
                         choices  = c("None (density plot)" = "__none__", cols),
                         selected = cols[1])
    })

    score_var_is_cat <- shiny::reactive({
      meta <- get_meta(); var <- input$score_var
      if (is.null(meta) || is.null(var) || var == "__none__") return(FALSE)
      if (!var %in% colnames(meta)) return(FALSE)
      !is.numeric(meta[[var]])
    })

    output$score_cond_ui <- shiny::renderUI({
      meta <- get_meta(); var <- input$score_var
      if (is.null(meta) || is.null(var) || var == "__none__") return(NULL)
      shiny::req(var %in% colnames(meta))
      is_cat <- !is.numeric(meta[[var]])

      if (is_cat) {
        shiny::div(
          shiny::tags$hr(style = "margin:8px 0;"),
          .bm_contrast_mode_radio(ns, "score_mode")
        )
      } else {
        shiny::div(
          shiny::tags$hr(style = "margin:8px 0;"),
          shiny::tags$label("Correlation method:", style = "font-size:0.85em; font-weight:600;"),
          shiny::radioButtons(ns("cor_method"), label = NULL,
            choices  = c("Pearson"  = "pearson",
                         "Spearman" = "spearman",
                         "Kendall"  = "kendall"),
            selected = "pearson", inline = TRUE)
        )
      }
    })

    # ---- Contrast filter (Score Analysis shared) ---------------------------
    # Shows all available pairwise contrasts. The selected subset is used
    # as a display filter: FPR data is filtered post-hoc by contrast column;
    # the ROC contrast dropdown, the AUC heatmap, and the All Methods Summary
    # heatmap/volcano are all restricted to this selection.
    # Note: Score Distributions always shows every contrast because that
    # function builds contrasts internally and cannot be post-hoc filtered
    # without re-running.

    output$score_contrasts_ui <- shiny::renderUI({
      meta <- get_meta(); var <- input$score_var; mode <- input$score_mode %||% "simple"
      if (is.null(meta) || is.null(var) || var == "__none__") return(NULL)
      if (!var %in% colnames(meta)) return(NULL)
      if (is.numeric(meta[[var]])) return(NULL)
      uvals <- sort(unique(na.omit(as.character(meta[[var]]))))
      conts <- tryCatch(
        remove_division(generate_all_contrasts(uvals, mode = mode)),
        error = function(e) character(0))
      if (length(conts) == 0L) return(NULL)
      shiny::div(
        shiny::tags$hr(style = "margin:8px 0;"),
        shiny::tags$label("Contrasts to show:", style = "font-size:0.85em; font-weight:600;"),
        shiny::tags$small(
          style = "display:block; color:#6b7280; font-size:0.78em; margin-bottom:4px;",
          "Filters the FPR plots, the ROC contrast dropdown, the AUC heatmap,",
          " and the All Methods Summary heatmap/volcano.",
          " Score Distributions always shows every contrast."
        ),
        shinyWidgets::pickerInput(
          ns("score_selected_groups"), label = NULL,
          choices = conts, selected = conts, multiple = TRUE,
          options = shinyWidgets::pickerOptions(
            actionsBox = TRUE, liveSearch = TRUE,
            selectedTextFormat = "count > 3",
            countSelectedText = "{0} contrasts selected",
            container = "body"
          )
        )
      )
    })

    # ---- Focus Cohen's d UI (Score Distributions only) --------------------

    output$score_cohend_ui <- shiny::renderUI({
      meta <- get_meta(); var <- input$score_var
      if (is.null(meta) || is.null(var) || var == "__none__") return(NULL)
      if (is.numeric(meta[[var]])) return(NULL)   # not applicable for numeric variables
      levs <- sort(unique(as.character(na.omit(meta[[var]]))))
      shiny::tags$details(
        open = "open",
        style = "margin-bottom:6px;",
        shiny::tags$summary(
          style = "cursor:pointer; font-size:0.85em; color:#555; user-select:none;",
          shiny::icon("chevron-right", style = "font-size:0.7em;"),
          " Effect size (Cohen's d) options"
        ),
        shiny::div(
          style = "padding:6px 0 2px;",
          shiny::checkboxInput(ns("compute_cohen"),
            "Show effect size (Cohen's d/f) on plots", value = FALSE),
          shiny::tags$small(
            style = "color:#888; display:block; margin:-4px 0 8px; font-size:0.78em;",
            "Annotates each distribution plot with Cohen's d (categorical) or",
            " Cohen's f (numeric). Uncheck for faster rendering."),
          shiny::tags$hr(style = "margin:4px 0 8px;"),
          shiny::tags$small(
            style = "color:#888; display:block; margin-bottom:6px;",
            "Optionally focus the Cohen's d annotation on one specific pairwise comparison",
            " e.g. ", shiny::tags$em("Disease"), " vs. ", shiny::tags$em("Control"),
            ". Also required to compute the t-test p-value (below)."
          ),
          shiny::div(
            style = "display:grid; grid-template-columns:1fr 1fr; gap:8px;",
            shiny::selectInput(ns("cohend_group_a"), "Group A (positive):",
                               choices = levs, selected = levs[1]),
            shiny::selectInput(ns("cohend_group_b"), "Group B (reference):",
                               choices = c("All others" = "__rest__", levs),
                               selected = "__rest__")
          )
        )
      )
    })

    # ---- Colour variable renderer (Score Distributions) --------------------

    output$score_color_var_ui <- shiny::renderUI({
      meta <- get_meta(); shiny::req(!is.null(meta))
      cols <- setdiff(colnames(meta), colnames(meta)[1])
      shiny::selectInput(
        ns("score_color_var"), "Colour by:",
        choices  = c("None" = "", cols),
        selected = "")
    })

    # ---- Shared param extractor -------------------------------------------

    .score_shared_params <- function() {
      meta      <- get_meta()
      var       <- shiny::isolate(input$score_var)
      is_cat    <- shiny::isolate(score_var_is_cat())
      var_arg   <- if (!is.null(var) && var != "__none__") var else NULL
      color_var <- shiny::isolate(input$score_color_var)
      color_arg <- if (!is.null(color_var) && nzchar(color_var)) color_var else NULL
      mode      <- shiny::isolate(input$score_mode) %||% "simple"
      pvalcalc  <- shiny::isolate(input$pvalcalc) %||% FALSE
      connect_g <- shiny::isolate(input$connect_groups) %||% FALSE
      cor_m     <- shiny::isolate(input$cor_method) %||% "pearson"

      cond_cohend <- NULL
      if (is_cat && !is.null(var_arg) && !is.null(shiny::isolate(input$cohend_group_a))) {
        grp_a <- shiny::isolate(input$cohend_group_a)
        grp_b <- shiny::isolate(input$cohend_group_b)
        if (!is.null(grp_a) && nzchar(grp_a)) {
          all_levs <- sort(unique(as.character(na.omit(meta[[var_arg]]))))
          b_levs   <- if (is.null(grp_b) || grp_b == "__rest__") setdiff(all_levs, grp_a)
                      else grp_b
          cond_cohend <- list(A = grp_a, B = b_levs)
        }
      }

      list(var_arg = var_arg, is_cat = is_cat, mode = mode,
           color_arg = color_arg,
           pvalcalc = pvalcalc, connect_g = connect_g,
           cor_m = cor_m, cond_cohend = cond_cohend)
    }

    # (score_filter_meta removed - contrast filtering is now display-only via
    #  the FPR data filter and the ROC contrast dropdown. Observers no longer
    #  pre-filter metadata before passing to analysis functions.)

    # Helper: estimate n contrasts from variable + mode
    .est_n_cont <- function(meta, var, mode) {
      n <- length(unique(na.omit(meta[[var]])))
      switch(mode,
        simple    = max(1L, n * (n - 1L) %/% 2L),
        medium    = n,
        extensive = max(1L, n * (n - 1L) %/% 2L),
        1L)
    }

    # Helper: build the score distribution figure directly from
    # CalculateScores() output, replicating the *exact* per-signature
    # aesthetic of PlotScores_Categorical() / PlotScores_Numeric() - gene set
    # name as a plot title (not a facet strip), angled x-axis labels, and
    # Cohen's d/f (+ p-value) reported as a plot *subtitle*, exactly as those
    # exported functions do it - without modifying those functions. Each
    # signature becomes its own ggplot, arranged with ggpubr::ggarrange()
    # (same as the package functions) and rendered as a static plot.
    .build_score_dist_plots <- function(expr, meta, gs, method, var_arg, is_cat,
                                         color_arg, compute_cohen, cond_cohend,
                                         pvalcalc, connect_g, cor_m, pointsize,
                                         wid_title, ncol = 3L, color_values = NULL,
                                         titlesize = 16, labsize = 12) {
      scores_list <- CalculateScores(data = expr, metadata = meta,
                                     gene_sets = gs, method = method)
      sig_names <- names(scores_list)
      ncol_grid <- max(1L, min(ncol, length(sig_names)))
      nrow_grid <- ceiling(length(sig_names) / ncol_grid)

      gg_list <- list()

      # ---- Case 1: no grouping variable -> density plot per gene set --------
      if (is.null(var_arg) || is.null(meta)) {

        fill_col <- if (is.null(color_values)) "#ECBD78" else unname(color_values)[1]

        for (sig in sig_names) {
          df <- scores_list[[sig]]
          wrapped_title <- wrap_title(sig, wid_title)

          p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$score)) +
            ggplot2::geom_density(fill = fill_col, alpha = 0.5) +
            ggplot2::geom_rug(ggplot2::aes(x = .data$score), color = fill_col, sides = "b",
                              alpha = 0.8, linewidth = .5, length = grid::unit(0.035, "npc")) +
            ggplot2::theme_classic() +
            ggplot2::labs(title = wrapped_title, x = "", y = "") +
            ggplot2::theme(
              axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = labsize - .5),
              axis.text.y = ggplot2::element_text(size = labsize - .5),
              plot.title  = ggplot2::element_text(hjust = 0.5, size = titlesize - 1))

          gg_list[[sig]] <- p
        }

        xlab <- switch(method,
                       ssGSEA    = "ssGSEA Enrichment Score",
                       logmedian = "Normalised Signature Score",
                       ranking   = "Signature Genes' Ranking")
        ylab <- "Density"
        tooltip_vars <- c("x")

      } else if (is_cat) {

        # ---- Case 2: categorical variable -> violin + jitter ---------------
        for (sig in sig_names) {
          df <- scores_list[[sig]]
          df[[var_arg]] <- factor(df[[var_arg]],
                                  levels = sort(unique(as.character(df[[var_arg]]))))
          wrapped_title <- wrap_title(sig, wid_title)

          p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[var_arg]], y = .data$score))

          if (!is.null(color_arg)) {
            p <- p + ggplot2::geom_jitter(
              ggplot2::aes(color = .data[[color_arg]], text = .data$sample),
              width = 0.2, height = 0, size = pointsize, alpha = 0.6)
          } else {
            p <- p + ggplot2::geom_jitter(
              ggplot2::aes(text = .data$sample),
              width = 0.2, height = 0, size = pointsize, alpha = 0.6, color = "#5264B6")
          }

          p <- p +
            ggplot2::geom_violin(alpha = 0.4, scale = "width") +
            ggplot2::stat_summary(fun = stats::median, fun.min = stats::median,
                                  fun.max = stats::median, geom = "crossbar", width = 0.25,
                                  position = ggplot2::position_dodge(width = 0.13))

          if (isTRUE(connect_g) && !is.null(color_arg)) {
            p <- p + ggplot2::stat_summary(
              ggplot2::aes(group = .data[[color_arg]], color = .data[[color_arg]]),
              fun = stats::median, geom = "line", linewidth = 1, alpha = 0.75,
              show.legend = FALSE)
          }

          # ---- Cohen's d/f (+ p-value) subtitle, exactly as computed in
          # PlotScores_Categorical() -----------------------------------------
          subtitle <- NULL
          if (isTRUE(compute_cohen)) {
            subtitle <- tryCatch({
              if (!is.null(cond_cohend)) {

                x <- df[df[[var_arg]] %in% cond_cohend$A, "score", drop = TRUE]
                y <- df[df[[var_arg]] %in% cond_cohend$B, "score", drop = TRUE]
                cohen_d_results <- cohen_d(x, y)

                if (isTRUE(pvalcalc)) {
                  pv <- stats::t.test(x, y, var.equal = TRUE)$p.value
                  line1 <- wrap_title(paste0("Cohen's d = ",
                                             format(signif(cohen_d_results, digits = 3),
                                                    scientific = FALSE)), width = wid_title)
                  line2 <- wrap_title(paste0("p = ", format(signif(pv, digits = 3),
                                                            scientific = TRUE)), width = wid_title)
                  paste(line1, line2, sep = "\n")
                } else {
                  wrap_title(paste0("Cohen's d = ",
                                    format(signif(cohen_d_results, digits = 3),
                                           scientific = FALSE)), width = wid_title)
                }

              } else if (length(levels(df[[var_arg]])) == 2) {

                g1 <- levels(df[[var_arg]])[1]; g2 <- levels(df[[var_arg]])[2]
                x <- df[df[[var_arg]] == g1, "score", drop = TRUE]
                y <- df[df[[var_arg]] == g2, "score", drop = TRUE]
                cohen_d_results <- cohen_d(x, y)

                if (isTRUE(pvalcalc)) {
                  ttest_results <- rstatix::t_test(df, formula =
                                                     as.formula(paste("score ~", var_arg)))
                  pv <- ttest_results$p[1]
                  line1 <- wrap_title(paste0("Cohen's d = ",
                                             format(signif(cohen_d_results, digits = 3),
                                                    scientific = FALSE)), width = wid_title)
                  line2 <- wrap_title(paste0("p = ", format(signif(pv, digits = 3),
                                                            scientific = TRUE)), width = wid_title)
                  paste(line1, line2, sep = "\n")
                } else {
                  wrap_title(paste0("Cohen's d = ",
                                    format(signif(cohen_d_results, digits = 3),
                                           scientific = FALSE)), width = wid_title)
                }

              } else if (length(levels(df[[var_arg]])) > 2) {

                type  <- identify_variable_type(df, var_arg)[var_arg]
                model <- lm(score ~ get(var_arg), data = df)
                results_var <- compute_cohens_f_pval(model, type)

                if (isTRUE(pvalcalc)) {
                  line1 <- wrap_title(paste0("Cohen's f = ",
                                             format(signif(results_var["Cohen_f"], digits = 3),
                                                    scientific = FALSE)), width = wid_title)
                  line2 <- wrap_title(paste0("p = ",
                                             format(signif(results_var["P_Value"], digits = 3),
                                                    scientific = TRUE)), width = wid_title)
                  paste(line1, line2, sep = "\n")
                } else {
                  wrap_title(paste0("Cohen's f = ",
                                    format(signif(results_var["Cohen_f"], digits = 3),
                                           scientific = FALSE)), width = wid_title)
                }

              } else NULL
            }, error = function(e) NULL)
          }

          p <- p + ggplot2::theme_bw() +
            ggplot2::theme(
              axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = labsize),
              axis.text.y = ggplot2::element_text(size = labsize),
              plot.title    = ggplot2::element_text(hjust = 0.5, size = titlesize - 1),
              plot.subtitle = ggplot2::element_text(hjust = 0.5, size = titlesize - 1.5,
                                                    face = "italic")) +
            ggplot2::labs(title = wrapped_title, subtitle = subtitle, color = "", x = "", y = "")

          if (!is.null(color_arg) && !is.null(color_values)) {
            p <- p + ggplot2::scale_color_manual(values = color_values)
          } else if (!is.null(color_arg)) {
            p <- p + ggplot2::scale_color_brewer(palette = "Set3")
          }

          gg_list[[sig]] <- p
        }

        xlab <- var_arg
        ylab <- paste0("Gene Set's Score (", method, ")")
        tooltip_vars <- c("x", "y", "text", "colour")

      } else {

        # ---- Case 3: numeric variable -> scatter + linear fit --------------
        pt_col <- if (is.null(color_values)) "#5264B6" else unname(color_values)[1]

        for (sig in sig_names) {
          df <- scores_list[[sig]]
          wrapped_title <- wrap_title(sig, wid_title)

          p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[var_arg]], y = .data$score,
                                                text = .data$sample)) +
            ggplot2::geom_point(size = pointsize, alpha = 0.5, color = pt_col) +
            ggplot2::geom_smooth(method = "lm", col = "black", se = FALSE, linewidth = 1)

          subtitle <- tryCatch({
            if (isTRUE(compute_cohen)) {
              type  <- identify_variable_type(df, var_arg)[var_arg]
              model <- lm(score ~ get(var_arg), data = df)
              results_var <- compute_cohens_f_pval(model, type)

              if (isTRUE(pvalcalc)) {
                line1 <- wrap_title(paste0("Cohen's f = ",
                                           format(signif(results_var["Cohen_f"], digits = 3),
                                                  scientific = FALSE)), width = wid_title)
                line2 <- wrap_title(paste0("p = ",
                                           format(signif(results_var["P_Value"], digits = 3),
                                                  scientific = TRUE)), width = wid_title)
                paste(line1, line2, sep = "\n")
              } else {
                wrap_title(paste0("Cohen's f = ",
                                  format(signif(results_var["Cohen_f"], digits = 3),
                                         scientific = FALSE)), width = wid_title)
              }
            } else NULL
          }, error = function(e) NULL)

          p <- p + ggplot2::theme_bw() +
            ggplot2::theme(
              axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = labsize),
              axis.text.y = ggplot2::element_text(size = labsize),
              plot.title    = ggplot2::element_text(hjust = 0.5, size = titlesize - 1),
              plot.subtitle = ggplot2::element_text(hjust = 0.5, size = titlesize - 1.5,
                                                    face = "italic")) +
            ggplot2::labs(title = wrapped_title, subtitle = subtitle, x = "", y = "")

          gg_list[[sig]] <- p
        }

        xlab <- var_arg
        ylab <- paste0("Gene Set's Score (", method, ")")
        tooltip_vars <- c("x", "y", "text")
      }

      # ---- Combined static figure: identical layout to the original
      # PlotScores_Categorical()/PlotScores_Numeric() (ggarrange + shared
      # outer axis labels, common legend, per-panel title + subtitle). -----
      combined_static <- ggpubr::ggarrange(plotlist = gg_list, ncol = ncol_grid,
                                           nrow = nrow_grid, common.legend = TRUE,
                                           align = "h")
      combined_static <- ggpubr::annotate_figure(combined_static,
        left = grid::textGrob(ylab, rot = 90, vjust = 1,
                              gp = grid::gpar(cex = 1.3, fontsize = labsize)),
        bottom = grid::textGrob(xlab,
                                gp = grid::gpar(cex = 1.3, fontsize = labsize)))

      list(gg = combined_static)
    }

    # ---- Run: Score Distributions -----------------------------------------

    shiny::observeEvent(input$run_scores, {
      expr <- get_expr(); meta <- get_meta()
      gs   <- .active_gs(); shiny::req(expr, meta, gs, length(gs) > 0)

      p    <- .score_shared_params()
      meth <- shiny::isolate(input$score_method) %||% "logmedian"
      do_cohen <- shiny::isolate(input$compute_cohen) %||% TRUE

      # Build palette-matched ColorValues.
      # PlotScores applies ColorValues as scale_color_manual - the names must
      # match whichever column is actually mapped to color:
      #  • ColorVariable set → name by ColorVariable levels
      #  • Variable categorical → name by Variable (grouping) levels
      score_color_vals <- if (!is.null(p$color_arg)) {
        col_levs <- levels(as.factor(meta[[p$color_arg]]))
        stats::setNames(rep_len(.pp_palette, length(col_levs)), col_levs)
      } else if (!is.null(p$var_arg) && p$is_cat) {
        levs <- levels(as.factor(meta[[p$var_arg]]))
        stats::setNames(rep_len(.pp_palette, length(levs)), levs)
      } else NULL

      # Violins get hard to read when a grouping variable has many levels
      # crammed into a 3-column grid, so drop to 2 columns once there are
      # more than 3 levels (numeric/no-grouping cases keep 3, since those are
      # density/scatter plots without this problem).
      n_levels <- if (isTRUE(p$is_cat) && !is.null(p$var_arg))
        length(unique(stats::na.omit(meta[[p$var_arg]])))
      else NA_integer_
      ncol_use <- if (!is.na(n_levels) && n_levels > 3L) 2L else 3L

      shiny::withProgress(message = "Computing score distributions...", value = 0, {
        result <- tryCatch({
          shiny::incProgress(0.4, detail = paste("Scoring:", meth))
          plot_res <- .build_score_dist_plots(
            expr = expr, meta = meta, gs = gs,
            method = meth, var_arg = p$var_arg, is_cat = p$is_cat,
            color_arg = p$color_arg, compute_cohen = do_cohen,
            cond_cohend = p$cond_cohend, pvalcalc = p$pvalcalc,
            connect_g = p$connect_g, cor_m = p$cor_m,
            pointsize = shiny::isolate(input$score_pointsize) %||% 4L,
            wid_title = shiny::isolate(input$score_width_title) %||% 32L,
            titlesize = shiny::isolate(input$score_titlesize) %||% 8L,
            labsize   = shiny::isolate(input$score_labsize)   %||% 8L,
            ncol = ncol_use, color_values = score_color_vals
          )
          shiny::incProgress(0.6, detail = "Done.")
          list(plot_res = plot_res, n_gs = length(gs), ncol_use = ncol_use,
               score_var = shiny::isolate(input$score_var),
               is_cat = p$is_cat, method = meth)
        }, error = function(e) {
          shiny::showNotification(paste("Score plot failed:", conditionMessage(e)),
                                  type = "error", duration = 12); NULL
        })
      })
      if (!is.null(result))
        score_h(as.integer(min(5000L, 480L + ceiling(result$n_gs / result$ncol_use) * 380L)))
      score_dist_result(result)
    })

    # Helper: build the AUC heatmap grid (one panel per gene set) directly
    # from the raw ROCAUC_Scores_Calculate() output, with an optional
    # contrasts filter and a fixed-column grid. Mirrors AUC_Scores() without
    # modifying the exported package function (and avoids recomputing the
    # ROC/AUC data a second time, since the caller already has it).

    # pROC::roc() is called by ROCAUC_Scores_Calculate() (exported, untouched)
    # with its default direction = "auto", which picks whichever direction
    # gives AUC >= 0.5 - meaning the stored $ROC curve is ALREADY direction-
    # corrected, and its AUC essentially never drops below 0.5 on its own.
    # The manual `ifelse(auc < 0.5, 1 - auc, auc)` in that function is mostly
    # a no-op because of this. To show genuinely "raw" curves/AUC that can
    # dip below the diagonal, we rebuild the roc object here with a fixed
    # direction ("<": higher score assumed to indicate the case group) using
    # the same predictor/response pROC already stored - no core function
    # changes and no recomputation of scores needed.
    .raw_roc <- function(roc_obj) {
      tryCatch(
        pROC::roc(response = roc_obj$response, predictor = roc_obj$predictor,
                  levels = roc_obj$levels, direction = "<", quiet = TRUE),
        error = function(e) roc_obj
      )
    }

    .build_auc_heatmap <- function(auc_list, contrasts = NULL,
                                   ncol = NULL, nrow = NULL, limits = NULL,
                                   widthTitle = 22, titlesize = 12, labsize = 12,
                                   ColorValues = c("#F9F4AE", "#B44141"),
                                   title = NULL, restore_direction = FALSE) {
      heatmaps <- list()
      max_row_lbl <- 0L

      for (signature_name in names(auc_list[[1]])) {
        auc_matrix <- matrix(nrow = length(auc_list[[1]][[signature_name]]),
                             ncol = length(auc_list))
        rownames(auc_matrix) <- names(auc_list[[1]][[signature_name]])  # Contrasts
        colnames(auc_matrix) <- names(auc_list)                        # Methods

        for (method_name in names(auc_list)) {
          for (contrast_name in names(auc_list[[method_name]][[signature_name]])) {
            entry <- auc_list[[method_name]][[signature_name]][[contrast_name]]
            # By default show the genuinely raw (fixed-direction) AUC via
            # .raw_roc(); only when restore_direction is TRUE do we use the
            # pROC "auto"-direction value ROCAUC_Scores_Calculate() stored
            # in $AUC (always >= 0.5 by construction - see .raw_roc() above).
            auc_matrix[contrast_name, method_name] <- if (restore_direction)
              entry$AUC
            else
              as.numeric(pROC::auc(.raw_roc(entry$ROC)))
          }
        }

        # Restrict to the requested contrasts, if given (falls back to all if
        # none of the requested names match, e.g. a stale selection).
        if (!is.null(contrasts) && length(contrasts) > 0) {
          keep_rows <- intersect(contrasts, rownames(auc_matrix))
          if (length(keep_rows) > 0) auc_matrix <- auc_matrix[keep_rows, , drop = FALSE]
        }

        long_data <- as.data.frame(as.table(auc_matrix))
        colnames(long_data) <- c("Contrast", "Method", "AUC")
        long_data$label <- sprintf("%.2f", long_data$AUC)

        signature_title <- wrap_title(signature_name, widthTitle)
        # Raw (un-flipped) AUC can fall below 0.5, so the colour scale needs
        # the full 0-1 range in that case; the restored-direction AUC is
        # always >= 0.5, so a tighter 0.5-1 range shows more contrast.
        lims <- if (is.null(limits)) (if (restore_direction) c(0.5, 1) else c(0, 1)) else limits
        max_row_lbl <- max(max_row_lbl, nchar(as.character(long_data$Contrast)), na.rm = TRUE)

        # A diverging scale centred on white: an AUC of 0 and an AUC of 1 both
        # mean "strong separation" (just in opposite directions), so both
        # ends use the same colour and the same visual weight. When direction
        # is folded to >= 0.5 (restore_direction), only the white-to-red half
        # is relevant since AUC can no longer go below 0.5.
        hi_col   <- ColorValues[length(ColorValues)]
        mid_col  <- "white"
        fill_scale <- if (restore_direction)
          ggplot2::scale_fill_gradient(low = mid_col, high = hi_col, limits = lims, name = "AUC")
        else
          ggplot2::scale_fill_gradient2(low = hi_col, mid = mid_col, high = hi_col,
                                        midpoint = 0.5, limits = lims, name = "AUC")

        p <- ggplot2::ggplot(long_data, ggplot2::aes(x = .data$Method, y = .data$Contrast,
                                                     fill = .data$AUC)) +
          ggplot2::geom_tile() +
          ggplot2::geom_text(ggplot2::aes(label = .data$label), color = "black",
                             size = labsize / 3) +
          fill_scale +
          ggplot2::labs(title = signature_title, x = NULL, y = NULL, fill = "AUC") +
          ggplot2::theme_bw() +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = labsize),
                        axis.text.y = ggplot2::element_text(size = labsize),
                        plot.title = ggplot2::element_text(hjust = 0.5, size = titlesize))

        heatmaps[[signature_name]] <- p
      }

      num_signatures <- length(heatmaps)
      if (is.null(nrow) && is.null(ncol)) {
        ncol <- min(4L, num_signatures)
        nrow <- ceiling(num_signatures / ncol)
      } else if (is.null(nrow)) {
        nrow <- ceiling(num_signatures / ncol)
      } else if (is.null(ncol)) {
        ncol <- ceiling(num_signatures / nrow)
      }

      for (i in seq_along(heatmaps)) {
        col_idx <- (i - 1) %% ncol + 1
        heatmaps[[i]] <- heatmaps[[i]] +
          ggplot2::theme(
            axis.text.y  = if (col_idx == 1) ggplot2::element_text(size = labsize) else ggplot2::element_blank(),
            axis.ticks.y = if (col_idx == 1) ggplot2::element_line() else ggplot2::element_blank(),
            plot.margin  = ggplot2::margin(4, 0, 0, 0)
          )
      }

      # The first column alone carries the row (contrast) labels, so it needs
      # a wider share of the grid than the others - scaled to how long those
      # labels actually are rather than a flat, often-insufficient guess.
      first_col_w <- max(1.3, min(3.5, 1 + max_row_lbl * 0.05))
      widths <- c(first_col_w, rep(1, ncol - 1))
      plt <- ggpubr::ggarrange(
        plotlist = heatmaps, ncol = ncol, nrow = nrow,
        common.legend = TRUE, legend = "right", align = "h", widths = widths
      )
      if (!is.null(title))
        plt <- ggpubr::annotate_figure(plt, top = grid::textGrob(
          title, gp = grid::gpar(cex = 1.3, fontsize = titlesize + 2)))

      list(plt = plt, ncol = ncol, nrow = nrow)
    }

    # ---- Run: ROC / AUC ---------------------------------------------------

    shiny::observeEvent(input$run_roc, {
      expr <- get_expr(); meta <- get_meta()
      gs   <- .active_gs(); shiny::req(expr, meta, gs, length(gs) > 0)
      p    <- .score_shared_params()
      meth <- shiny::isolate(input$roc_method) %||% "logmedian"

      if (!p$is_cat) {
        shiny::showNotification(
          "Classification requires a categorical variable. Please select one above.",
          type = "warning", duration = 8)
        return(NULL)
      }

      # Use the same contrast logic as ROC_Scores for an accurate panel count
      n_cont <- tryCatch({
        levs_var <- unique(na.omit(meta[[p$var_arg]]))
        length(remove_division(generate_all_contrasts(levs_var, mode = p$mode)))
      }, error = function(e) .est_n_cont(meta, p$var_arg, p$mode))

      # AUC heatmap only shows the contrasts selected in "Contrasts to show" above.
      sel_conts <- shiny::isolate(input$score_selected_groups)

      shiny::withProgress(message = "Computing ROC curves and AUC...", value = 0, {
        result <- tryCatch({
          roc_wid_title <- shiny::isolate(input$roc_width_title) %||% 32L
          roc_titlesize <- shiny::isolate(input$roc_titlesize)   %||% 8L
          roc_labsize   <- shiny::isolate(input$roc_labsize)     %||% 8L

          shiny::incProgress(0.4, detail = "ROC curves...")
          # Use the raw calculate function so we can build per-contrast plots reactively
          roc_raw <- ROCAUC_Scores_Calculate(
            data = expr, metadata = meta, gene_sets = gs,
            method = meth, variable = p$var_arg, mode = p$mode)

          shiny::incProgress(0.6, detail = "Done.")
          auc_ncol <- min(4L, length(gs))

          # Extract contrast names from raw data (they're in [[method]][[sig]][[contrast]])
          roc_contrasts <- names(roc_raw[[1]][[1]])

          # The AUC heatmap itself is built live inside the render function
          # below (reading font-size and the "restore directionality" checkbox
          # live), so toggling those doesn't require re-running this (the
          # expensive part - ROC/AUC computation) again.
          list(roc_raw      = roc_raw,
               roc_contrasts = roc_contrasts,
               auc_ncol      = auc_ncol,
               n_gs          = length(gs),
               n_cont        = n_cont,
               roc_wid_title = roc_wid_title,
               roc_titlesize = roc_titlesize,
               roc_labsize   = roc_labsize)
        }, error = function(e) {
          shiny::showNotification(paste("ROC/AUC failed:", conditionMessage(e)),
                                  type = "error", duration = 12); NULL
        })
      })
      if (!is.null(result)) {
        # Same schema as the Score Distributions grid: fixed columns (up to 3),
        # height grows with the number of rows needed, not linearly with n_gs.
        n_cols_roc <- min(3L, result$n_gs)
        n_rows_roc <- ceiling(result$n_gs / n_cols_roc)
        roc_h(as.integer(min(6000L, max(400L, n_rows_roc * 420L))))

        n_rows_auc <- ceiling(result$n_gs / result$auc_ncol)
        auc_h(as.integer(min(3000L, max(300L, n_rows_auc * 300L + 100L))))
      }
      score_roc_result(result)
    })

    # Helper: subset a cohenlist (output of CohenD_allConditions /
    # CohenF_allConditions) down to a chosen set of contrasts (columns).
    # Falls back to the full cohenlist if no contrasts match/are given.
    .filter_cohenlist_contrasts <- function(cohenlist, contrasts) {
      if (is.null(contrasts) || length(contrasts) == 0) return(cohenlist)
      lapply(cohenlist, function(sig_res) {
        keep <- intersect(contrasts, colnames(sig_res[[1]]))
        if (length(keep) == 0) return(sig_res)
        lapply(sig_res, function(df) df[, keep, drop = FALSE])
      })
    }

    # Helper: rebuild the Cohen's d/f heatmap grid from cohenlist data, with an
    # adjustable cell text size. Mirrors Heatmap_Cohen() (which hard-codes the
    # cell label size) without modifying the exported package function.
    .build_cohen_heatmap <- function(cohenlist, ncol = NULL, nrow = NULL,
                                     limits = NULL, widthTitle = 22,
                                     titlesize = 12, textsize = 3, labsize = 12,
                                     ColorValues = NULL, title = NULL) {
      cohentype <- if ("CohenD" %in% names(cohenlist[[1]])) "d"
                   else if ("CohenF" %in% names(cohenlist[[1]])) "f"
                   else stop("Error: cohenlist format not valid.")

      if (!is.null(ColorValues)) {
        ColorValues <- if (!is.null(ColorValues[["heatmap"]])) ColorValues[["heatmap"]]
                       else ColorValues[[1]]
      }
      if (is.null(ColorValues)) ColorValues <- c("#F9F4AE", "#B44141")

      heatmaps <- list()
      for (signature in names(cohenlist)) {
        cohen_mat <- if (cohentype == "d") t(as.matrix(cohenlist[[signature]]$CohenD))
                     else t(as.matrix(cohenlist[[signature]]$CohenF))
        p_value_mat <- t(as.matrix(cohenlist[[signature]]$padj))

        long_data <- data.frame(
          Var1 = rep(rownames(cohen_mat), times = ncol(cohen_mat)),
          Var2 = rep(colnames(cohen_mat), each = nrow(cohen_mat)),
          Cohen = abs(as.vector(cohen_mat)),
          PValue = as.vector(p_value_mat),
          stringsAsFactors = FALSE
        )
        long_data$label <- paste0(sprintf("%.2f", long_data$Cohen), "\n(",
                                  format.pval(long_data$PValue, digits = 1), ")")

        signature_title <- wrap_title(signature, widthTitle)
        lims <- if (is.null(limits)) c(0, max(long_data$Cohen, na.rm = TRUE)) else limits

        p <- ggplot2::ggplot(long_data, ggplot2::aes(x = .data$Var2, y = .data$Var1,
                                                     fill = .data$Cohen)) +
          ggplot2::geom_tile() +
          ggplot2::geom_text(ggplot2::aes(label = .data$label), color = "black",
                             size = textsize) +
          ggplot2::scale_fill_gradientn(colors = ColorValues, limits = lims) +
          ggplot2::labs(title = signature_title, x = NULL, y = NULL,
                       fill = if (cohentype == "d") "|Cohen's d|" else "|Cohen's f|") +
          ggplot2::theme_bw() +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = labsize),
                        axis.text.y = ggplot2::element_text(size = labsize),
                        plot.title = ggplot2::element_text(hjust = 0.5, size = titlesize))

        heatmaps[[signature]] <- p
      }

      num_signatures <- length(heatmaps)
      if (is.null(nrow) && is.null(ncol)) {
        ncol <- min(3L, num_signatures)
        nrow <- ceiling(num_signatures / ncol)
      } else if (is.null(nrow)) {
        nrow <- ceiling(num_signatures / ncol)
      } else if (is.null(ncol)) {
        ncol <- ceiling(num_signatures / nrow)
      }

      for (i in seq_along(heatmaps)) {
        col_idx <- (i - 1) %% ncol + 1
        heatmaps[[i]] <- heatmaps[[i]] +
          ggplot2::theme(
            axis.text.y  = if (col_idx == 1) ggplot2::element_text(size = labsize) else ggplot2::element_blank(),
            axis.ticks.y = if (col_idx == 1) ggplot2::element_line() else ggplot2::element_blank(),
            plot.margin  = ggplot2::margin(4, 0, 0, 0)
          )
      }

      widths <- c(1.5, rep(1, ncol - 1))
      plt <- ggpubr::ggarrange(
        plotlist = heatmaps, ncol = ncol, nrow = nrow,
        common.legend = TRUE, legend = "right", align = "h", widths = widths
      )
      plt <- ggpubr::annotate_figure(plt, top = grid::textGrob(
        title, gp = grid::gpar(cex = 1.3, fontsize = titlesize + 2)))

      list(plt = plt, data = cohenlist, ncol = ncol, nrow = nrow)
    }

    # Helper: rebuild the Cohen's d/f volcano from cohenlist data, with the
    # facet grid capped at a given number of columns. Mirrors Volcano_Cohen()
    # (whose facet_wrap() ignores its own ncol/nrow arguments) without
    # modifying the exported package function.
    .build_cohen_volcano <- function(cohenlist, ncol = NULL, nrow = NULL,
                                     titlesize = 12, labsize = 12, ColorValues = NULL,
                                     title = NULL, widthlegend = 22,
                                     pointSize = 3, sig_threshold = 0.05,
                                     cohen_threshold = 0.5,
                                     colorPalette = "Set3") {
      cohentype <- if ("CohenD" %in% names(cohenlist[[1]])) "d"
                   else if ("CohenF" %in% names(cohenlist[[1]])) "f"
                   else stop("Error: cohenlist format not valid.")

      rows <- list()
      for (signature in names(cohenlist)) {
        sig_data <- cohenlist[[signature]]
        cohen_mat <- if (cohentype == "d") sig_data$CohenD else sig_data$CohenF
        padj_mat  <- sig_data$padj

        for (method in rownames(cohen_mat)) {
          for (contrast in colnames(cohen_mat)) {
            rows[[length(rows) + 1]] <- data.frame(
              signature = signature, contrast = contrast, method = method,
              cohen = cohen_mat[method, contrast],
              padj  = padj_mat[method, contrast],
              stringsAsFactors = FALSE
            )
          }
        }
      }
      final_df <- do.call(rbind, rows)
      final_df$signature <- vapply(final_df$signature,
                                   function(x) wrap_title(x, widthlegend),
                                   character(1))

      if (is.null(ColorValues)) {
        ColorValues <- grDevices::colorRampPalette(
          RColorBrewer::brewer.pal(12, colorPalette))(length(unique(final_df$signature)))
      } else {
        ColorValues <- if (!is.null(ColorValues[["volcano"]])) ColorValues[["volcano"]]
                       else ColorValues[[2]]
      }

      ggplot2::ggplot(final_df, ggplot2::aes(x = abs(.data$cohen),
                                             y = -log10(.data$padj),
                                             shape = .data$method)) +
        ggplot2::geom_point(colour = "black", size = pointSize) +
        ggplot2::geom_point(ggplot2::aes(colour = .data$signature),
                            size = pointSize - 1.5) +
        ggplot2::facet_wrap(. ~ .data$contrast, scales = "free", ncol = ncol, nrow = nrow) +
        ggplot2::geom_hline(yintercept = -log10(sig_threshold),
                            linetype = "dashed", color = "black", linewidth = 0.5) +
        ggplot2::geom_vline(xintercept = cohen_threshold,
                            linetype = "dashed", color = "black", linewidth = 0.5) +
        ggplot2::scale_color_manual(values = ColorValues) +
        ggplot2::scale_shape_manual(
          values = 15:(15 + length(unique(final_df$method)) - 1)) +
        ggplot2::labs(
          x = if (cohentype == "d") "|Cohen's d|" else "|Cohen's f|",
          y = "-log10(Adj. p-value)", color = "Signature", shape = "Method"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          legend.position  = "right",
          axis.text        = ggplot2::element_text(size = labsize),
          axis.title       = ggplot2::element_text(size = labsize + 1),
          legend.text      = ggplot2::element_text(size = labsize),
          legend.title     = ggplot2::element_text(size = labsize),
          strip.text       = ggplot2::element_text(size = titlesize, face = "bold"),
          plot.title       = ggplot2::element_text(hjust = 0.5, size = titlesize + 1, face = "bold"),
          strip.background = ggplot2::element_rect(fill = "white")
        ) +
        ggplot2::ggtitle(if (!is.null(title)) title else
          if (cohentype == "d") "Cohen's d Volcano Plot" else "Cohen's f Volcano Plot") +
        ggplot2::scale_x_continuous(limits = c(0, NA)) +
        ggplot2::scale_y_continuous(limits = c(0, NA))
    }

    # ---- Run: All Methods Summary ------------------------------------------

    # Only the (expensive) Cohen's d/f computation lives in this observer -
    # font/text sizes are read live inside the render functions below instead
    # of being "baked in" here, so tweaking them doesn't require re-running
    # the whole computation again.
    shiny::observeEvent(input$run_all_methods, {
      expr <- get_expr(); meta <- get_meta()
      gs   <- .active_gs(); shiny::req(expr, meta, gs, length(gs) > 0)
      p <- .score_shared_params()

      # Only the currently selected contrasts (from the "Contrasts to show"
      # picker above) are included in the heatmap and volcano.
      sel_conts <- shiny::isolate(input$score_selected_groups)

      shiny::withProgress(message = "Running all scoring methods...", value = 0, {
        result <- tryCatch({
          shiny::incProgress(0.6, detail = "Computing Cohen's d/f across all methods...")

          cohentype <- if (p$is_cat) "d" else "f"
          cohenlist_full <- if (cohentype == "f") {
            CohenF_allConditions(data = expr, metadata = meta, gene_sets = gs,
                                 variable = p$var_arg)
          } else {
            CohenD_allConditions(data = expr, metadata = meta, gene_sets = gs,
                                 variable = p$var_arg, mode = p$mode)
          }
          cohenlist <- .filter_cohenlist_contrasts(cohenlist_full, sel_conts)
          n_cont    <- ncol(cohenlist[[1]][[1]])

          # Grid layout (ncol/nrow) only depends on the number of gene sets,
          # not on any font/display setting, so it's safe to compute once here.
          num_sigs  <- length(cohenlist)
          heat_ncol <- min(3L, num_sigs)
          heat_nrow <- ceiling(num_sigs / heat_ncol)
          shiny::incProgress(0.4, detail = "Done.")

          list(cohenlist  = cohenlist,
               gs_names   = names(gs),
               n_gs       = length(gs),
               n_cont     = n_cont,
               heat_ncol  = heat_ncol, heat_nrow = heat_nrow,
               cohentype  = cohentype)
        }, error = function(e) {
          shiny::showNotification(paste("All methods failed:", conditionMessage(e)),
                                  type = "error", duration = 12); NULL
        })
      })
      if (!is.null(result)) {
        # Heatmap: height grows with the number of gene-set-panel rows.
        all_h(as.integer(min(4000L, max(350L, result$heat_nrow * 280L + 220L))))
        # Volcano: fixed 2-column facet grid (by contrast); height grows with rows.
        vol_ncol <- min(2L, result$n_cont)
        vol_nrow <- ceiling(result$n_cont / max(1L, vol_ncol))
        all_vol_h(as.integer(min(3000L, max(420L, vol_nrow * 380L + 120L))))
      }
      score_all_result(result)
    })

    # ---- Run: GSEA ---------------------------------------------------------

    shiny::observeEvent(input$run_gsea, {
      expr <- get_expr(); meta <- get_meta()
      gs   <- .active_gs(); shiny::req(expr, meta, gs, length(gs) > 0)

      design_mode <- shiny::isolate(input$gsea_design_mode) %||% "auto"
      nperm       <- shiny::isolate(input$gsea_nperm)
      padj_m      <- shiny::isolate(input$padj_method)
      cont_corr   <- as.logical(shiny::isolate(input$contrast_correction) %||% FALSE)
      stat_use    <- {
        s <- shiny::isolate(input$gsea_stat) %||% "auto"
        if (s == "auto") NULL else s
      }

      shiny::withProgress(message = "Running GSEA...", value = 0, {
        result <- tryCatch({
          shiny::incProgress(0.1, detail = "Building design/contrasts...")

          if (identical(design_mode, "custom")) {
            # ---- Custom formula path ----------------------------------------
            formula_str  <- trimws(shiny::isolate(input$gsea_formula) %||% "")
            custom_conts <- trimws(shiny::isolate(input$gsea_custom_contrasts) %||% "")
            shiny::req(nchar(formula_str) > 0)
            mmat <- tryCatch(
              model.matrix(stats::as.formula(formula_str), data = meta),
              error = function(e) stop("Formula error: ", e$message))
            # Strip the variable-name prefix from column names so contrasts like
            # "Senescent-Proliferative" work instead of needing the full
            # "ConditionSenescent-ConditionProliferative" form.
            # This mirrors what calculateDE does internally in auto mode.
            vars_in_formula <- tryCatch(
              all.vars(stats::as.formula(formula_str)),
              error = function(e) character(0))
            for (v in vars_in_formula)
              colnames(mmat) <- sub(paste0("^", v), "", colnames(mmat))
            colnames(mmat) <- gsub(" ", "", colnames(mmat))
            contrasts_vec <- if (nchar(custom_conts) > 0L)
              trimws(strsplit(custom_conts, "\n")[[1L]])
            else
              NULL   # calculateDE returns all columns
            enrich_var  <- formula_str
            enrich_mode <- "custom"
            shiny::incProgress(0.2, detail = "Fitting linear models (custom)...")
            DEGs <- calculateDE(data = expr, metadata = NULL,
                                modelmat = mmat, contrasts = contrasts_vec)
          } else {
            # ---- Auto path --------------------------------------------------
            enrich_var  <- shiny::isolate(input$enrich_var)
            enrich_mode <- shiny::isolate(input$enrich_mode) %||% "simple"
            shiny::req(enrich_var)
            # Sanitize the group's values in this local copy of `meta` BEFORE
            # building contrasts or calling calculateDE(), so the design
            # matrix calculateDE() builds internally (from this same `meta`)
            # and the contrast strings built from `uvals` below always refer
            # to the same, valid names - see .sanitize_group_vals() above.
            raw_vals  <- as.character(meta[[enrich_var]])
            safe_vals <- .sanitize_group_vals(raw_vals)
            if (!identical(raw_vals, safe_vals)) {
              meta[[enrich_var]] <- safe_vals
              shiny::showNotification(
                paste0("Some values of '", enrich_var, "' contained characters not valid ",
                       "in R names (e.g. ':') and were adjusted for this analysis (e.g. '.' ",
                       "instead of ':'). This only affects labels in the GSEA results, not your data."),
                type = "warning", duration = 10)
            }
            uvals <- unique(safe_vals)
            all_contrasts <- remove_division(
              generate_all_contrasts(uvals, mode = enrich_mode))
            # Filter to user-selected contrasts
            sel <- shiny::isolate(input$gsea_selected_contrasts)
            contrasts_vec <- if (!is.null(sel) && length(sel) > 0L)
              intersect(all_contrasts, sel)
            else
              all_contrasts
            if (length(contrasts_vec) == 0L)
              stop("No contrasts selected. Please select at least one contrast.")
            shiny::incProgress(0.2, detail = "Fitting linear models...")
            DEGs <- calculateDE(data = expr, metadata = meta,
                                variables = enrich_var, contrasts = contrasts_vec)
          }

          shiny::incProgress(0.5, detail = sprintf(
            "Running fgsea (%d permutations)...", nperm))
          gsea_res <- runGSEA(DEGList = DEGs, gene_sets = gs,
                              stat = stat_use, ContrastCorrection = cont_corr,
                              nPermSimple = nperm, p.adjust.method = padj_m)
          shiny::incProgress(0.2, detail = "Done.")
          n_cont <- length(DEGs); n_sets <- length(gs)
          list(gsea_res = gsea_res, DEGs = DEGs, gene_sets = gs,
               contrasts = names(DEGs), var = enrich_var,
               n_sets = n_sets, n_cont = n_cont,
               run_params = list(enrich_var = enrich_var, enrich_mode = enrich_mode,
                                 nperm = nperm, padj_m = padj_m,
                                 n_gene_sets = n_sets))
        }, error = function(e) {
          shiny::showNotification(paste("GSEA failed:", conditionMessage(e)),
                                  type = "error", duration = 12); NULL
        })
      })
      if (!is.null(result)) {
        # NES lollipop is rendered as one full-width row per contrast (ncol = 1,
        # see output$nes_lollipop below), so the total figure height scales
        # linearly with the number of contrasts (rows stacked) and, within each
        # row, with the number of gene sets on the y-axis (~22px per gene set
        # plus a fixed allowance for the title, the two stat_used facet strips,
        # and axis text).
        nes_h(as.integer(min(10000L, max(400L,
          result$n_cont * (result$n_sets * 34L + 180L)))))
        enrich_ncol_est <- min(2L, max(1L, result$n_sets * result$n_cont))
        enrich_nrow_est <- ceiling((result$n_sets * result$n_cont) / enrich_ncol_est)
        enrich_h(as.integer(min(14000L, max(460L, enrich_nrow_est * 400L + 60L))))
        gsea_vol_h(as.integer(min(10000L, max(460L, result$n_sets * 470L + 70L))))
        # The legend now sits at the bottom (wrapped to a few rows, capped
        # regardless of gene set count - see .build_combined_gsea_gg()), so
        # height mainly needs a generous, fairly fixed base for the two
        # side-by-side scatter panels themselves plus that bottom legend
        # block, rather than scaling with the number of gene sets.
        combined_h(as.integer(min(3000L, max(650L,
          650L + max(0L, result$n_cont - 6L) * 12L))))
      }
      gsea_result(result)
    })

    # Helper: build one small violin panel per (gene set x contrast)
    # combination directly from the $data returned by FPR_Simulation(), so we
    # can lay panels out in a flat grid with the contrast shown as a subtitle
    # (FPR_Simulation()'s own plot facets contrasts as a strip label inside
    # each gene-set panel, which is harder to read and cannot be restyled
    # without modifying the exported function).
    .build_fpr_plots <- function(fpr_data, ylab, pointSize = 2,
                                 titlesize = 12, labsize = 10, widthTitle = 32) {
      color_values <- c(Original = "#68B393", Simulated = "#666666")
      plot_list <- list()

      for (sig in names(fpr_data)) {
        df_sig <- fpr_data[[sig]]
        for (ct in unique(df_sig$contrast)) {
          df_ct <- df_sig[df_sig$contrast == ct, ]
          lvls  <- levels(df_ct$method)
          if (is.null(lvls)) lvls <- unique(as.character(df_ct$method))

          q_data <- data.frame()
          for (mt in lvls) {
            subset_df <- df_ct[df_ct$method == mt, ]
            if (nrow(subset_df) == 0) next
            q95  <- stats::quantile(subset_df$cohen, 0.95, na.rm = TRUE)
            xpos <- which(lvls == mt)
            q_data <- rbind(q_data, data.frame(
              method = mt, q_high = q95, xmin = xpos - 0.3, xmax = xpos + 0.3))
          }

          all_max <- stats::aggregate(cohen ~ method, data = df_ct, FUN = max)
          original_df <- df_ct[df_ct$type == "Original", ]
          all_max$FPR <- NA
          for (i in seq_len(nrow(all_max))) {
            match_idx <- which(original_df$method == all_max$method[i])
            if (length(match_idx) > 0) all_max$FPR[i] <- original_df$FPR[match_idx[1]]
          }
          all_max$label  <- sprintf("FPR=%.2f", all_max$FPR)
          all_max$y      <- all_max$cohen + 0.3
          all_max$method <- factor(all_max$method, levels = lvls)

          p <- ggplot2::ggplot() +
            ggplot2::geom_jitter(
              data = df_ct[df_ct$type == "Simulated", ],
              ggplot2::aes(y = .data$cohen, x = .data$method, color = .data$type),
              width = 0.3, height = 0, size = pointSize, alpha = 0.5) +
            ggplot2::geom_violin(
              data = df_ct, ggplot2::aes(y = .data$cohen, x = .data$method),
              fill = "#F0F0F0", color = "black", alpha = 0.5) +
            ggplot2::geom_jitter(
              data = df_ct[df_ct$type == "Original", ],
              ggplot2::aes(y = .data$cohen, x = .data$method, color = .data$type),
              width = 0.3, height = 0, size = pointSize, alpha = 1) +
            ggplot2::geom_text(
              data = all_max,
              ggplot2::aes(x = .data$method, y = .data$y, label = .data$label),
              size = labsize / 3, inherit.aes = FALSE) +
            ggplot2::geom_segment(
              data = q_data,
              ggplot2::aes(x = .data$xmin, xend = .data$xmax, y = .data$q_high, yend = .data$q_high),
              linetype = "dashed", color = "red", inherit.aes = FALSE) +
            ggplot2::labs(title = wrap_title(sig, widthTitle),
                         subtitle = wrap_title(as.character(ct), widthTitle),
                         y = ylab, x = "Method", color = "") +
            ggplot2::theme_classic() +
            ggplot2::theme(
              plot.title      = ggplot2::element_text(hjust = 0.5, size = titlesize),
              plot.subtitle   = ggplot2::element_text(hjust = 0.5, size = titlesize - 1.5),
              axis.text       = ggplot2::element_text(size = labsize),
              axis.title      = ggplot2::element_text(size = labsize),
              legend.text     = ggplot2::element_text(size = labsize),
              legend.title    = ggplot2::element_text(size = labsize)) +
            ggplot2::scale_color_manual(values = color_values)

          plot_list[[paste(sig, ct, sep = " | ")]] <- p
        }
      }
      plot_list
    }

    # ---- Run: FPR Simulation -----------------------------------------------

    shiny::observeEvent(input$run_fpr, {
      expr <- get_expr(); meta <- get_meta()
      gs   <- .active_gs(); shiny::req(expr, meta, gs, length(gs) > 0, input$score_var)

      fpr_var   <- shiny::isolate(input$score_var)
      fpr_nsims <- shiny::isolate(input$fpr_nsims)
      fpr_mode  <- shiny::isolate(input$fpr_mode)

      gs_names  <- names(gs); n_gs <- length(gs_names)

      font_s <- shiny::isolate(input$fpr_fontsize) %||% 10L

      shiny::withProgress(
        message = sprintf("Score FPR Simulation (%d gene sets, %d simulations each)...",
                          n_gs, fpr_nsims),
        value = 0.05, {
          result <- tryCatch({
            shiny::setProgress(value = 0.1, detail = "Running simulations...")
            # We only use $data from FPR_Simulation() (one data frame per gene
            # set, with a "contrast" column); the grid layout and per-panel
            # styling are built entirely in the app (see .build_fpr_plots(),
            # used by the renderer below) from that data.
            res_all <- FPR_Simulation(
              data                = expr,
              metadata            = meta,
              original_signatures = gs,
              Variable            = fpr_var,
              number_of_sims      = fpr_nsims,
              mode                = fpr_mode,
              titlesize           = font_s,
              widthTitle          = 32L
            )
            shiny::setProgress(value = 1, detail = "Done.")
            # Determine cohen type so the plot y-axis label is correct
            var_type   <- identify_variable_type(meta, fpr_var)
            cohentype  <- if (var_type == "Numeric" || fpr_mode == "none") "f" else "d"
            list(fpr_res = res_all, var = fpr_var, nsims = fpr_nsims,
                 n_gs = n_gs, cohentype = cohentype)
          }, error = function(e) {
            shiny::showNotification(paste("FPR failed:", conditionMessage(e)),
                                    type = "error", duration = 12); NULL
          })
        })
      if (!is.null(result)) {
        # Total panels = one per (gene set x contrast); same grid schema as
        # the Score Distributions panel (fixed columns, height grows by row).
        ncol_fpr    <- 3L
        n_contrasts <- length(unique(unlist(
          lapply(result$fpr_res$data, function(df) unique(df$contrast)))))
        n_panels    <- max(1L, result$n_gs * max(1L, n_contrasts))
        n_rows_fpr  <- ceiling(n_panels / ncol_fpr)
        fpr_h(as.integer(min(8000L, max(500L, 150L + n_rows_fpr * 380L))))
      }
      fpr_result(result)
    })

    # ---- Variable selectors ------------------------------------------------

    output$enrich_var_ui <- shiny::renderUI({
      meta <- get_meta(); shiny::req(!is.null(meta))
      design_mode <- input$gsea_design_mode %||% "auto"
      if (identical(design_mode, "custom")) return(NULL)
      cat_cols <- .bm_cat_cols(meta)
      if (length(cat_cols) == 0)
        return(shiny::div(class = "alert alert-warning",
                          style = "font-size:0.82em; padding:5px 9px;",
                          "No categorical metadata variables found."))
      shiny::selectInput(ns("enrich_var"), "Contrast variable:",
                         choices = cat_cols, selected = cat_cols[1])
    })

    # Custom design matrix UI
    output$enrich_custom_design_ui <- shiny::renderUI({
      design_mode <- input$gsea_design_mode %||% "auto"
      if (!identical(design_mode, "custom")) return(NULL)
      shiny::div(
        shiny::tags$small(
          style = "color:#6b7280; display:block; margin-bottom:6px; font-size:0.8em;",
          "Paste an R formula for ", shiny::tags$code("model.matrix()"),
          " (e.g. ", shiny::tags$code("~ 0 + Condition"), " or ",
          shiny::tags$code("~ 0 + Condition + Batch"), ").",
          " Column names must match your metadata. Leave Contrasts blank to return all columns."
        ),
        shiny::textInput(ns("gsea_formula"), "Model formula:",
                         value = "~ 0 + Condition", width = "100%"),
        shiny::uiOutput(ns("gsea_design_cols_ui")),
        shiny::textAreaInput(ns("gsea_custom_contrasts"),
          "Contrasts (one per line, e.g. Senescent - Proliferative):",
          value = "", rows = 3, width = "100%"),
        shiny::tags$small(
          style = "color:#6b7280; font-size:0.78em;",
          "Use the group names shown above (variable prefix is stripped automatically).",
          " Leave blank to return all design matrix columns as results."
        )
      )
    })

    # Show design matrix column names for custom mode (so user knows what to type)
    output$gsea_design_cols_ui <- shiny::renderUI({
      design_mode <- input$gsea_design_mode %||% "auto"
      if (!identical(design_mode, "custom")) return(NULL)
      meta <- get_meta(); formula_str <- trimws(input$gsea_formula %||% "")
      if (is.null(meta) || nchar(formula_str) == 0L) return(NULL)
      cols <- tryCatch({
        mmat <- model.matrix(stats::as.formula(formula_str), data = meta)
        vars_in_formula <- all.vars(stats::as.formula(formula_str))
        cn <- colnames(mmat)
        for (v in vars_in_formula) cn <- sub(paste0("^", v), "", cn)
        cn <- gsub(" ", "", cn)
        cn[cn != "(Intercept)"]
      }, error = function(e) character(0))
      if (length(cols) == 0L) return(NULL)
      shiny::div(
        style = "background:#f0f4ff; border:1px solid #c7d2fe; border-radius:4px; padding:5px 8px; margin:4px 0 6px; font-size:0.78em;",
        shiny::tags$b("Available group names: "),
        paste(cols, collapse = ", ")
      )
    })

    # Available contrasts checkboxes (auto mode only)
    output$enrich_contrasts_ui <- shiny::renderUI({
      design_mode <- input$gsea_design_mode %||% "auto"
      if (!identical(design_mode, "auto")) return(NULL)
      meta <- get_meta(); var <- input$enrich_var; mode <- input$enrich_mode
      shiny::req(!is.null(meta), !is.null(var), nchar(var) > 0)
      contrasts_all <- tryCatch({
        uvals <- unique(.sanitize_group_vals(meta[[var]]))
        remove_division(generate_all_contrasts(uvals, mode = mode %||% "simple"))
      }, error = function(e) character(0))
      if (length(contrasts_all) == 0L) return(NULL)
      shiny::div(
        shiny::tags$label(
          style = "font-size:0.82em; font-weight:600; color:#374151; display:block; margin-bottom:4px;",
          "Contrasts to test ",
          bslib::tooltip(
            shiny::icon("circle-question",
              style = "color:#aaa; cursor:help; font-size:0.85em;"),
            "All contrasts are pre-selected. Remove any you don't need.",
            " Fewer contrasts = faster computation and less multiple-testing correction.",
            placement = "right"
          )
        ),
        # Same widget as the "Contrasts to show" picker in Score Analysis, for
        # a consistent contrast-selection experience across tabs.
        shinyWidgets::pickerInput(
          ns("gsea_selected_contrasts"), label = NULL,
          choices = contrasts_all, selected = contrasts_all, multiple = TRUE,
          options = shinyWidgets::pickerOptions(
            actionsBox = TRUE, liveSearch = TRUE,
            selectedTextFormat = "count > 3",
            countSelectedText = "{0} contrasts selected",
            container = "body"
          )
        )
      )
    })

    output$fpr_var_ui <- shiny::renderUI({
      meta <- get_meta(); shiny::req(!is.null(meta))
      cat_cols <- .bm_cat_cols(meta)
      if (length(cat_cols) == 0)
        return(shiny::div(class = "alert alert-warning",
                          style = "font-size:0.82em; padding:5px 9px;",
                          "No categorical metadata variables found."))
      shiny::selectInput(ns("fpr_var"), "Grouping variable:",
                         choices = cat_cols, selected = cat_cols[1])
    })

    # ---- Renders: Score Distributions --------------------------------------

    output$score_dist_status_ui <- shiny::renderUI({
      if (is.null(score_dist_result()))
        shiny::div(class = "alert alert-info", style = "margin:10px 0;",
                   shiny::icon("circle-info"),
                   " Set the variable above, select a scoring method,",
                   " then click ", shiny::strong("Plot Score Distributions"), ".")
    })

    output$score_dist_ui <- shiny::renderUI({
      shiny::req(!is.null(score_dist_result()))
      res   <- score_dist_result()
      label <- if (res$is_cat) "Score Distributions per Group"
               else if (!is.null(res$score_var) && res$score_var == "__none__")
                 "Score Density"
               else "Score vs Continuous Variable"
      bslib::card(
        full_screen = TRUE,
        .bm_card_header(label, "dl_score_dist", ns),
        shiny::plotOutput(ns("score_dist_plot"), height = paste0(score_h(), "px"))
      )
    })

    output$score_dist_plot <- shiny::renderPlot({
      res <- score_dist_result(); shiny::req(!is.null(res), !is.null(res$plot_res$gg))
      tryCatch(print(res$plot_res$gg),
               error = function(e)
                 shiny::showNotification(paste("Score plot error:", e$message),
                                         type = "warning", duration = 8))
    }, height = function() score_h(), res = 150)

    output$dl_score_dist <- shiny::downloadHandler(
      filename = function() paste0("score_distributions_", Sys.Date(), ".png"),
      content  = function(file) {
        res <- score_dist_result(); shiny::req(!is.null(res), !is.null(res$plot_res$gg))
        h_in <- score_h() / 96
        ggplot2::ggsave(file, plot = res$plot_res$gg,
                        width = 12, height = max(4, h_in),
                        dpi = 150, units = "in")
      }
    )

    # ---- Renders: ROC / AUC ------------------------------------------------

    # ---- Helper: build ROC plot for one contrast from raw data ----------------
    .build_roc_contrast_plot <- function(roc_raw, sel_contrast,
                                         titlesize = 11L, labsize = 11L, wid_title = 32L,
                                         restore_direction = FALSE) {
      roc_colors <- c(logmedian = "#3E5587", ssGSEA = "#B65285", ranking = "#B68C52")
      method_names <- names(roc_raw)
      sig_names    <- names(roc_raw[[1]])

      per_sig <- lapply(sig_names, function(sig) {
        combined_df <- data.frame()
        auc_values  <- list()
        for (meth_name in method_names) {
          roc_data <- roc_raw[[meth_name]][[sig]][[sel_contrast]]
          if (is.null(roc_data)) next
          # By default use .raw_roc() (fixed direction, can dip below the
          # diagonal / AUC < 0.5). When restore_direction is TRUE, use the
          # pROC "auto"-direction curve ROCAUC_Scores_Calculate() computed
          # (always >= 0.5 by construction), matching its stored $AUC.
          curve <- if (restore_direction) roc_data$ROC else .raw_roc(roc_data$ROC)
          combined_df <- rbind(combined_df, data.frame(
            FPR    = rev(1 - curve$specificities),
            TPR    = rev(curve$sensitivities),
            Method = meth_name
          ))
          auc_values[[meth_name]] <- if (restore_direction)
            roc_data$AUC
          else
            as.numeric(pROC::auc(curve))
        }
        if (nrow(combined_df) == 0L) return(NULL)

        # Ensure color mapping works even for a single method
        cols_use <- roc_colors[names(auc_values)]
        cols_use[is.na(cols_use)] <- "#3E5587"

        auc_texts <- data.frame(
          Method = names(auc_values),
          AUC    = unlist(auc_values),
          x      = 1,
          y      = seq(0.05, 0.30, length.out = length(auc_values))
        )

        p <- ggplot2::ggplot(combined_df,
                             ggplot2::aes(x = .data$FPR, y = .data$TPR,
                                          color = .data$Method)) +
          ggplot2::geom_line(linewidth = 1) +
          ggplot2::scale_color_manual(values = cols_use) +
          ggplot2::geom_abline(linetype = "dashed", color = "gray") +
          ggplot2::geom_label(
            data = auc_texts,
            ggplot2::aes(x = .data$x, y = .data$y,
                         label = paste0("AUC ", .data$Method, " = ",
                                        round(.data$AUC, 2)),
                         color = .data$Method),
            size = labsize / 3, vjust = 0, hjust = 1, inherit.aes = FALSE, fill = "white") +
          ggplot2::labs(
            title    = wrap_title(sig, wid_title),
            subtitle = wrap_title(sel_contrast, wid_title),
            x = "False Positive Rate", y = "True Positive Rate") +
          ggplot2::theme_classic() +
          ggplot2::theme(
            legend.position = "none",
            axis.text     = ggplot2::element_text(size = labsize),
            axis.title    = ggplot2::element_text(size = labsize + 1),
            plot.title    = ggplot2::element_text(hjust = 0.5, size = titlesize),
            plot.subtitle = ggplot2::element_text(hjust = 0.5, size = titlesize - 1.5))
        p
      })
      per_sig <- Filter(Negate(is.null), per_sig)
      if (length(per_sig) == 0L) return(NULL)
      if (length(per_sig) == 1L) return(per_sig[[1]])
      ncol_p <- min(3L, length(per_sig))
      ggpubr::ggarrange(plotlist = per_sig, ncol = ncol_p,
                        nrow = ceiling(length(per_sig) / ncol_p), align = "h")
    }

    output$score_roc_status_ui <- shiny::renderUI({
      if (is.null(score_roc_result()))
        shiny::div(class = "alert alert-info", style = "margin:10px 0;",
                   shiny::icon("circle-info"),
                   " Set a categorical variable above, then click",
                   shiny::strong(" Compute ROC Curves and AUC"), ".")
    })

    output$score_roc_ui <- shiny::renderUI({
      res <- score_roc_result()
      shiny::req(!is.null(res), !is.null(res$roc_raw))
      all_conts <- res$roc_contrasts
      # Filter available contrasts by the "Contrasts to show" picker
      sel_conts <- input$score_selected_groups
      conts <- if (!is.null(sel_conts) && length(sel_conts) > 0L)
        intersect(all_conts, sel_conts) else all_conts
      if (length(conts) == 0L) conts <- all_conts
      bslib::card(
        full_screen = TRUE,
        .bm_card_header("ROC Curves", "dl_roc", ns),
        shiny::div(
          class = "text-muted", style = "font-size:0.78em; padding:4px 6px 8px;",
          shiny::icon("circle-info"),
          " Each panel shows the ROC curve for one gene set.",
          " Dashed diagonal = chance (AUC = 0.5). Curves above the diagonal indicate discriminatory power.",
          " By default the AUC shown matches the curve exactly, preserving directionality,",
          " even if it dips below 0.5. Tick \"Ignore direction\" in Display options above",
          " to instead show 1 - AUC whenever AUC < 0.5.",
          " The AUC heatmap below shows only the contrasts selected above (\"Contrasts to show\")."
        ),
        shiny::div(
          style = "padding:4px 6px 2px;",
          shiny::selectInput(ns("roc_sel_contrast"), "Contrast:",
                             choices = conts, selected = conts[1], width = "100%")
        ),
        shiny::plotOutput(ns("score_roc_plot"), height = paste0(roc_h(), "px"))
      )
    })

    output$score_roc_plot <- shiny::renderPlot({
      res <- score_roc_result()
      shiny::req(!is.null(res), !is.null(res$roc_raw))
      sel_ct <- input$roc_sel_contrast %||% res$roc_contrasts[1]
      shiny::req(nzchar(sel_ct %||% ""))
      tryCatch({
        p <- .build_roc_contrast_plot(
          res$roc_raw, sel_ct,
          titlesize = input$roc_titlesize %||% res$roc_titlesize,
          labsize   = input$roc_labsize   %||% res$roc_labsize,
          wid_title = res$roc_wid_title,
          restore_direction = isTRUE(input$roc_restore_direction))
        if (!is.null(p)) print(p)
      }, error = function(e)
        shiny::showNotification(paste("ROC plot error:", e$message),
                                type = "warning", duration = 8))
    }, height = function() roc_h(), res = 150)

    output$dl_roc <- shiny::downloadHandler(
      filename = function() paste0("roc_curves_", Sys.Date(), ".png"),
      content  = function(file) {
        res <- score_roc_result(); shiny::req(!is.null(res), !is.null(res$roc_raw))
        sel_ct <- shiny::isolate(input$roc_sel_contrast) %||% res$roc_contrasts[1]
        p    <- .build_roc_contrast_plot(res$roc_raw, sel_ct,
                                          titlesize   = input$roc_titlesize %||% res$roc_titlesize,
                                          labsize     = input$roc_labsize   %||% res$roc_labsize,
                                          wid_title   = res$roc_wid_title,
                                          restore_direction = isTRUE(input$roc_restore_direction))
        h_in <- roc_h() / 96
        ggplot2::ggsave(file, plot = p,
                        width = 14, height = max(5, h_in),
                        dpi = 150, units = "in")
      }
    )

    output$score_auc_ui <- shiny::renderUI({
      shiny::req(!is.null(score_roc_result()), !is.null(score_roc_result()$roc_raw))
      bslib::card(
        full_screen = TRUE,
        .bm_card_header("AUC Heatmap (contrasts × methods)", "dl_auc", ns),
        shiny::div(
          class = "text-muted", style = "font-size:0.78em; padding:4px 6px 8px;",
          shiny::icon("circle-info"),
          " Numbers are ", shiny::tags$b("AUC values"),
          " (Area Under the ROC Curve). AUC = 0.5: no better than chance.",
          " AUC = 1.0: perfect separation. AUC < 0.5: reversed direction (still informative -",
          " tick \"Ignore direction\" above to instead fold it to 1 - AUC, hiding which way it points)."
        ),
        shiny::plotOutput(ns("score_auc_plot"), height = paste0(auc_h(), "px"))
      )
    })

    # Built live from the raw ROC/AUC data, reading font-size and the "restore
    # directionality" checkbox directly from input$, so toggling either is
    # instant and never requires re-running the (expensive) ROC/AUC computation.
    output$score_auc_plot <- shiny::renderPlot({
      res <- score_roc_result(); shiny::req(!is.null(res), !is.null(res$roc_raw))
      tryCatch({
        auc_heat <- .build_auc_heatmap(
          res$roc_raw, contrasts = input$score_selected_groups,
          ncol = res$auc_ncol,
          titlesize = input$roc_titlesize %||% res$roc_titlesize,
          labsize   = input$roc_labsize   %||% res$roc_labsize,
          restore_direction = isTRUE(input$roc_restore_direction))
        print(auc_heat$plt)
      }, error = function(e)
        shiny::showNotification(paste("AUC plot error:", e$message),
                                type = "warning", duration = 8))
    }, height = function() auc_h(), res = 150)

    output$dl_auc <- shiny::downloadHandler(
      filename = function() paste0("auc_heatmap_", Sys.Date(), ".png"),
      content  = function(file) {
        res <- score_roc_result(); shiny::req(!is.null(res), !is.null(res$roc_raw))
        auc_heat <- .build_auc_heatmap(
          res$roc_raw, contrasts = input$score_selected_groups,
          ncol = res$auc_ncol,
          titlesize = input$roc_titlesize %||% res$roc_titlesize,
          labsize   = input$roc_labsize   %||% res$roc_labsize,
          restore_direction = isTRUE(input$roc_restore_direction))
        h_in <- auc_h() / 96
        ggplot2::ggsave(file, plot = auc_heat$plt,
                        width = 10, height = max(4, h_in),
                        dpi = 150, units = "in")
      }
    )

    # ---- Renders: All Methods Summary --------------------------------------

    output$score_all_status_ui <- shiny::renderUI({
      if (is.null(score_all_result()))
        shiny::div(class = "alert alert-info", style = "margin:10px 0;",
                   shiny::icon("circle-info"),
                   " Click ", shiny::strong("Run All Methods Summary"),
                   "; all three scoring methods will be computed.")
    })

    output$score_all_heatmap_ui <- shiny::renderUI({
      res <- score_all_result(); shiny::req(!is.null(res))
      p_test <- if (identical(res$cohentype, "f"))
        "a linear model F-test (score ~ variable)."
      else
        "a two-sample t-test (equal variances assumed)."
      bslib::card(
        full_screen = TRUE,
        .bm_card_header("Cohen's d/f Summary Heatmap", "dl_all_heatmap", ns),
        shiny::div(
          class = "text-muted", style = "font-size:0.78em; padding:4px 6px 8px;",
          shiny::icon("circle-info"),
          " Numbers are ",
          shiny::tags$b("Cohen's d"),
          " (categorical variable: standardised mean difference between groups) or ",
          shiny::tags$b("Cohen's f"),
          " (numeric variable). |d| > 0.2 small, |d| > 0.5 medium, |d| > 0.8 large.",
          " Only the contrasts selected above (\"Contrasts to show\") are included.",
          shiny::tags$br(),
          shiny::tags$b("P-value shown in brackets: "), " computed per gene set/method via ", p_test,
          " Adjusted (BH) across gene sets and contrasts, within each method."
        ),
        shiny::plotOutput(ns("score_all_heatmap"), height = paste0(all_h(), "px"))
      )
    })

    # ColorValues must be a named list:
    #   [[1]]/"heatmap" → 2-colour gradient
    #   [[2]]/"volcano" → named per-signature vector for Volcano_Cohen
    # Volcano_Cohen wraps signature names with wrap_title(width=widthlegend)
    # before building the data frame, so our color names must use the same
    # wrapping to avoid mismatches that produce grey points. Cheap to rebuild
    # on every render (no stats involved), so it can read the live legend-wrap
    # width input directly.
    .build_all_gs_colors <- function(gs_names, wid_legend) {
      gs_wrapped <- vapply(gs_names, wrap_title, character(1L), width = wid_legend)
      list(
        heatmap = c("#F9F4AE", "#B44141"),
        volcano = stats::setNames(rep_len(.pp_palette, length(gs_names)), gs_wrapped)
      )
    }

    output$score_all_heatmap <- shiny::renderPlot({
      res <- score_all_result(); shiny::req(!is.null(res), !is.null(res$cohenlist))
      gs_colors <- .build_all_gs_colors(res$gs_names, input$all_width_legend %||% 32L)
      tryCatch(
        print(.build_cohen_heatmap(
          cohenlist   = res$cohenlist,
          ncol        = res$heat_ncol, nrow = res$heat_nrow,
          widthTitle  = input$all_width_title  %||% 32L,
          textsize    = input$all_heatmap_textsize %||% 2,
          titlesize   = input$all_heatmap_titlesize %||% 8,
          labsize     = input$all_heatmap_labsize   %||% input$all_heatmap_titlesize %||% 8,
          ColorValues = gs_colors
        )$plt),
        error = function(e)
          shiny::showNotification(paste("Heatmap error:", e$message),
                                  type = "warning", duration = 8))
    }, height = function() all_h(), res = 150)

    output$dl_all_heatmap <- shiny::downloadHandler(
      filename = function() paste0("cohen_heatmap_", Sys.Date(), ".png"),
      content  = function(file) {
        res <- score_all_result(); shiny::req(!is.null(res), !is.null(res$cohenlist))
        gs_colors <- .build_all_gs_colors(res$gs_names, input$all_width_legend %||% 32L)
        heat_plt <- .build_cohen_heatmap(
          cohenlist   = res$cohenlist,
          ncol        = res$heat_ncol, nrow = res$heat_nrow,
          widthTitle  = input$all_width_title  %||% 32L,
          textsize    = input$all_heatmap_textsize %||% 2,
          titlesize   = input$all_heatmap_titlesize %||% 8,
          labsize     = input$all_heatmap_labsize   %||% input$all_heatmap_titlesize %||% 8,
          ColorValues = gs_colors
        )$plt
        h_in <- all_h() / 96
        ggplot2::ggsave(file, plot = heat_plt,
                        width = 12, height = max(4, h_in),
                        dpi = 150, units = "in")
      }
    )

    # Tableau-10 + extra colours, shared by the static combined GSEA volcano
    # below (kept as a named constant since it's reused there).
    .GSEA_COLOR_PALETTE <- c(
      "#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F",
      "#EDC948","#B07AA1","#FF9DA7","#9C755F","#BAB0AC",
      "#1F77B4","#FF7F0E","#2CA02C","#D62728","#9467BD",
      "#8C564B","#E377C2","#7F7F7F","#BCBD22","#17BECF"
    )

    output$score_all_volcano_ui <- shiny::renderUI({
      res <- score_all_result()
      shiny::req(!is.null(res), !is.null(res$cohenlist))
      p_test <- if (identical(res$cohentype, "f"))
        "a linear model F-test (score ~ variable)."
      else
        "a two-sample t-test (equal variances assumed)."
      bslib::card(
        full_screen = TRUE,
        bslib::card_header(
          class = "d-flex align-items-center justify-content-between",
          "Effect Size Volcano",
          shiny::downloadButton(ns("dl_all_volcano"), "",
            icon  = shiny::icon("download"),
            class = "btn-sm btn-outline-secondary py-0 px-2",
            style = "font-size:0.75em;")
        ),
        shiny::div(
          class = "text-muted", style = "font-size:0.78em; padding:4px 6px 8px;",
          shiny::icon("circle-info"),
          " X-axis: |Cohen's d| (effect size). Y-axis: −log₁₀(adjusted p-value), from ", p_test,
          " Shape = scoring method. Colour = gene set. One facet per selected contrast",
          " (up to 2 columns, additional rows added as needed).",
          " Dashed lines mark the significance and effect size thresholds above."
        ),
        shiny::plotOutput(ns("score_all_volcano"), height = paste0(all_vol_h(), "px"))
      )
    })

    .build_all_volcano_live <- function(res) {
      gs_colors <- .build_all_gs_colors(res$gs_names, input$all_width_legend %||% 32L)
      vol_ncol  <- min(2L, res$n_cont)
      .build_cohen_volcano(
        cohenlist       = res$cohenlist,
        ColorValues     = gs_colors,
        widthlegend     = input$all_width_legend %||% 32L,
        pointSize       = 4,
        titlesize       = input$all_volcano_titlesize %||% 8,
        labsize         = input$all_volcano_labsize    %||% 8,
        sig_threshold   = input$sig_threshold   %||% 0.05,
        cohen_threshold = input$cohen_threshold %||% 0.5,
        ncol            = vol_ncol
      )
    }

    output$score_all_volcano <- shiny::renderPlot({
      res <- score_all_result(); shiny::req(!is.null(res), !is.null(res$cohenlist))
      tryCatch(print(.build_all_volcano_live(res)),
               error = function(e)
                 shiny::showNotification(paste("Volcano error:", e$message),
                                         type = "warning", duration = 8))
    }, height = function() all_vol_h(), res = 150)

    output$dl_all_volcano <- shiny::downloadHandler(
      filename = function() paste0("effect_volcano_", Sys.Date(), ".png"),
      content  = function(file) {
        res <- score_all_result(); shiny::req(!is.null(res), !is.null(res$cohenlist))
        ggplot2::ggsave(file, plot = .build_all_volcano_live(res),
                        width = 14, height = 7, dpi = 150, units = "in")
      }
    )

    # ---- Renders: GSEA -----------------------------------------------------

    output$gsea_status_ui <- shiny::renderUI({
      if (is.null(gsea_result()))
        shiny::div(class = "alert alert-info", style = "margin:10px 0;",
                   shiny::icon("circle-info"),
                   " Select a contrast variable and click ",
                   shiny::strong("Run GSEA"), ".")
    })

    output$gsea_results_ui <- shiny::renderUI({
      res <- gsea_result(); shiny::req(!is.null(res))
      bslib::navset_tab(

        # 0. DEG Volcano plots
        bslib::nav_panel("DEG Volcanos",
          bslib::card(
            full_screen = TRUE,
            .bm_card_header("logFC vs −log₁₀(adjusted p-value)", "dl_gsea_deg_volcano", ns),
            shiny::div(
              class = "text-muted", style = "font-size:0.78em; padding:4px 6px 8px;",
              shiny::icon("circle-info"),
              " Volcano plot of differential expression for each contrast.",
              " Points are coloured by gene set membership when gene sets are provided.",
              " Thresholds mark significance and fold-change cutoffs."
            ),
            # Inline colour key - only relevant when gene set highlighting is on
            if (!is.null(res$gene_sets) && length(res$gene_sets) > 0L) shiny::div(
              style = "font-size:0.8em; padding:2px 6px 8px; display:flex; gap:14px; flex-wrap:wrap;",
              shiny::tags$span(shiny::tags$b(style="color:#B7B7B7;", "■"), " Background"),
              shiny::tags$span(shiny::tags$b(style="color:#038C65;", "■"),
                " Annotated UP in gene set",
                shiny::tags$i(style="font-size:0.9em; color:#666;", "(not observed DE)")),
              shiny::tags$span(shiny::tags$b(style="color:#8C0303;", "■"),
                " Annotated DOWN in gene set",
                shiny::tags$i(style="font-size:0.9em; color:#666;", "(not observed DE)")),
              shiny::tags$span(shiny::tags$b(style="color:#05254A;", "■"),
                " In gene set (no direction annotated)"),
              shiny::tags$span(shiny::tags$b(style="color:#6489B4;", "■"), " Significant (threshold)")
            ),
            .bm_collapsible_settings(
              shiny::div(
                style = "display:grid; grid-template-columns:1fr 1fr; gap:8px;",
                shiny::numericInput(ns("vol_threshold_x"), "logFC threshold:",
                                    value = 1, min = 0, max = 10, step = 0.5),
                bslib::tooltip(
                  shiny::numericInput(ns("vol_threshold_y"), "Adjusted p-value threshold (BH):",
                                      value = 0.05, min = 0, max = 1, step = 0.01),
                  "BH (Benjamini-Hochberg) adjusted p-value from the limma moderated t-test",
                  " (adj.P.Val), used both to colour significant genes and to pick which",
                  " genes are eligible for text labels below.",
                  placement = "right"
                ),
                bslib::tooltip(
                  shiny::numericInput(ns("vol_top_n"), "Label top N genes by |logFC| (0 = none):",
                                      value = 0, min = 0, max = 50, step = 1),
                  "With \"Gene set members\" highlighted: labels the N most extreme members",
                  " of each gene set (optionally narrowed by the thresholds above).",
                  " With \"Significant genes\" highlighted: labels the N most extreme genes",
                  " among those passing the thresholds above.",
                  placement = "right"
                ),
                shiny::numericInput(ns("vol_pointsize"), "Point size:",
                                    value = 2, min = 1, max = 8, step = 0.5),
                shiny::numericInput(ns("vol_labsize"), "Label size (pt):",
                                    value = 8, min = 5, max = 16, step = 1)
              ),
              shiny::radioButtons(ns("vol_highlight_mode"), "Highlight:",
                choices  = c("Gene set members" = "gene_sets",
                             "Significant genes (threshold)" = "threshold"),
                selected = "gene_sets", inline = TRUE)
            ),
            shiny::div(
              style = paste0("max-width:", gsea_vol_w(), ";"),
              shiny::plotOutput(ns("gsea_deg_volcano_plot"),
                                width = "100%",
                                height = paste0(gsea_vol_h(), "px"))
            )
          )
        ),

        # 1. Combined Volcano
        bslib::nav_panel("Combined Volcano",
          bslib::card(
            full_screen = TRUE,
            .bm_card_header("NES vs −log₁₀(adjusted p-value)", "dl_gsea_combined", ns),
            shiny::div(
              class = "text-muted", style = "font-size:0.78em; padding:4px 6px 8px;",
              shiny::icon("circle-info"),
              " Each point is one gene set x contrast. Colour = gene set; shape = contrast.",
              shiny::tags$br(),
              shiny::tags$b("Altered panel (B-statistic):"),
              " gene sets without specified direction, ranked by probability of differential expression.",
              " Positive NES = enriched among the most altered genes.",
              " Negative NES = enriched among the least altered genes (generally uninformative).",
              shiny::tags$br(),
              shiny::tags$b("Enriched/Depleted panel (t-statistic):"),
              " gene sets with specified direction, ranked by moderated t-statistic.",
              " Positive NES = enriched among up-regulated genes;",
              " negative NES = enriched among down-regulated genes.",
              shiny::tags$br(),
              " y-axis: −log₁₀(adjusted p-value); dashed line marks the significance threshold."
            ),
            shiny::tags$details(
              style = "margin-bottom:10px;",
              shiny::tags$summary(
                style = "cursor:pointer; font-size:0.85em; color:#555; user-select:none;",
                shiny::icon("chevron-right", style = "font-size:0.7em;"),
                " Display options"),
              shiny::div(
                style = "padding:8px 0 4px; display:grid; grid-template-columns:1fr 1fr; gap:8px;",
                shiny::numericInput(ns("combined_pointsize"), "Point size:",
                                    value = 3, min = 1, max = 16, step = 1),
                shiny::numericInput(ns("combined_sig_threshold"),
                                    "Significance threshold (adjusted p-value):",
                                    value = 0.05, min = 0.001, max = 0.5, step = 0.01),
                bslib::tooltip(
                  shiny::numericInput(ns("combined_width_legend"), "Legend wrap (chars):",
                                      value = 32, min = 8, max = 60, step = 2),
                  "Maximum characters per line in the gene set legend.",
                  placement = "right"
                ),
                bslib::tooltip(
                  shiny::numericInput(ns("combined_labsize"), "Label size (pt):",
                                      value = 8, min = 6, max = 20, step = 1),
                  "Font size of axis text, axis titles, legend text, and panel strip labels.",
                  placement = "right"
                )
              )
            ),
            shiny::plotOutput(ns("gsea_combined_plot"), height = paste0(combined_h(), "px"))
          )
        ),

        # 2. NES Lollipop
        bslib::nav_panel("NES Lollipop",
          bslib::card(
            full_screen = TRUE,
            .bm_card_header("NES Lollipop", "dl_nes_lollipop", ns),
            shiny::div(
              class = "text-muted", style = "font-size:0.78em; padding:4px 6px 8px;",
              shiny::icon("circle-info"),
              " Lollipop length = Normalised Enrichment Score (NES) magnitude.",
              " Two panels are shown:",
              shiny::tags$br(),
              shiny::tags$b("Altered (B-statistic):"),
              " gene sets without specified direction are ranked by the B-statistic",
              " (probability of differential expression regardless of direction).",
              " A positive NES indicates enrichment among the most altered genes.",
              " A negative NES indicates enrichment among the least altered genes (generally uninformative).",
              shiny::tags$br(),
              shiny::tags$b("Enriched/Depleted (t-statistic):"),
              " gene sets with specified direction are ranked by the moderated t-statistic",
              " (multiplied by -1 for putatively down-regulated genes).",
              " A positive NES indicates enrichment among up-regulated genes;",
              " a negative NES indicates enrichment among down-regulated genes.",
              shiny::tags$br(),
              shiny::tags$b("Point colour:"),
              " adjusted p-value (red = significant, white = not significant).",
              shiny::tags$b(" Dashed line:"), " significance threshold. One panel per contrast."
            ),
            shiny::tags$details(
              style = "margin-bottom:10px;",
              shiny::tags$summary(
                style = "cursor:pointer; font-size:0.85em; color:#555; user-select:none;",
                shiny::icon("chevron-right", style = "font-size:0.7em;"),
                " Display options"),
              shiny::div(
                style = "padding:8px 0 4px; display:grid; grid-template-columns:1fr 1fr; gap:8px;",
                shiny::numericInput(ns("nes_titlesize"), "Title size (pt):",
                                    value = 8, min = 6, max = 20, step = 1),
                bslib::tooltip(
                  shiny::numericInput(ns("nes_labsize"), "Label size (pt):",
                                      value = 8, min = 6, max = 20, step = 1),
                  "Font size of axis text, axis titles, legend text, and gene set (y-axis) labels.",
                  placement = "right"
                ),
                shiny::numericInput(ns("nes_sig_threshold"),
                                    "Significance threshold (adjusted p-value):",
                                    value = 0.05, min = 0.001, max = 0.5, step = 0.01),
                bslib::tooltip(
                  shiny::numericInput(ns("nes_width_legend"), "Gene set label wrap (chars):",
                                      value = 32, min = 8, max = 60, step = 2),
                  "Maximum characters per line for gene set names on the y-axis.",
                  placement = "right"
                )
              )
            ),
            shiny::plotOutput(ns("nes_lollipop"), height = paste0(nes_h(), "px"))
          )
        ),

        # 3. Enrichment Plots
        bslib::nav_panel("Enrichment Plots",
          bslib::card(
            full_screen = TRUE,
            .bm_card_header("Enrichment Plots", "dl_gsea_enrich", ns),
            shiny::div(
              class = "text-muted", style = "font-size:0.78em; padding:4px 6px 8px;",
              shiny::icon("circle-info"),
              " Running enrichment score along the ranked gene list for each gene set × contrast.",
              " The peak of the curve indicates where the gene set is most concentrated.",
              " Bars at the bottom mark the positions of gene set members in the ranked list.",
              " Reported NES and adjusted p-values are from fgsea.",
              shiny::tags$br(),
              shiny::tags$b("Two panels per contrast: "),
              shiny::tags$b("Enriched/Depleted"),
              " uses the t-statistic (genes ranked from most up- to most down-regulated;",
              " positive NES = enrichment among up-regulated genes). ",
              shiny::tags$b("Altered"),
              " uses the B-statistic (genes ranked by strength of differential evidence,",
              " regardless of direction; positive NES = gene set concentrated among",
              " the most strongly altered genes)."
            ),
            shiny::tags$details(
              style = "margin-bottom:10px;",
              shiny::tags$summary(
                style = "cursor:pointer; font-size:0.85em; color:#555; user-select:none;",
                shiny::icon("chevron-right", style = "font-size:0.7em;"),
                " Display options"),
              shiny::div(
                style = "padding:8px 0 4px; display:grid; grid-template-columns:1fr 1fr; gap:8px;",
                shiny::numericInput(ns("enrich_titlesize"), "Title size (pt):",
                                    value = 9, min = 6, max = 20, step = 1),
                bslib::tooltip(
                  shiny::numericInput(ns("enrich_labsize"), "Label size (pt):",
                                      value = 9, min = 6, max = 20, step = 1),
                  "Font size of axis text and axis titles on each enrichment panel.",
                  placement = "right"
                ),
                bslib::tooltip(
                  shiny::numericInput(ns("enrich_width_title"), "Title wrap (chars):",
                                      value = 38, min = 10, max = 60, step = 2),
                  "Maximum characters per line in each enrichment panel's title.",
                  placement = "right"
                )
              )
            ),
            shiny::plotOutput(ns("gsea_enrich_plot"), height = paste0(enrich_h(), "px"))
          )
        ),

        # 4. Results Table
        bslib::nav_panel("Results Table",
          shiny::uiOutput(ns("gsea_table_hdr_ui")),
          shiny::div(
            class = "text-muted", style = "font-size:0.78em; padding:4px 6px 8px;",
            shiny::icon("circle-info"),
            " NES > 0: gene set enriched among up-regulated genes for this contrast.",
            " NES < 0: enriched among down-regulated genes.",
            " leadingEdge: the subset of gene set members that drove the enrichment."
          ),
          DT::DTOutput(ns("gsea_table"))
        )
      )
    })

    # ---- GSEA: DEG Volcano --------------------------------------------------

    # Helper: extract character vector of all gene names from a gene set entry
    .gs_gene_names <- function(gs_entry) {
      if (is.data.frame(gs_entry)) as.character(gs_entry[[1]])
      else as.character(gs_entry)
    }

    # Helper: extract per-gene direction (1 = up, -1 = down, NA = unknown) for
    # a gene set entry. Only data-frame gene sets (Gene, Signal) carry a
    # direction; plain vectors have no direction information.
    .gs_gene_directions <- function(gs_entry) {
      if (is.data.frame(gs_entry)) as.numeric(gs_entry[[2]])
      else rep(NA_real_, length(gs_entry))
    }

    # Helper: build one DEG volcano ggplot per contrast.
    #
    # Highlight mode drives BOTH which points are coloured as "of interest"
    # AND which genes are eligible for text labels - the two were previously
    # decoupled (labels for "Gene set members" mode used to fall back to
    # plotVolcano()'s own top-N-by-logFC labelling across *all* genes,
    # ignoring gene set membership entirely; "Significant genes" mode never
    # labelled anything). Now:
    #   - "Gene set members": points/labels restricted to genes in each set;
    #     labels are the top N (by |logFC|) members of that set, optionally
    #     narrowed further by the logFC / padj thresholds.
    #   - "Significant genes (threshold)": points/labels restricted to genes
    #     passing the logFC / padj thresholds; labels are the top N (by
    #     |logFC|) among those.
    .build_deg_volcano_plots <- function(res, thr_x, thr_y, n_top, pts, lsz, hl_mode) {
      use_gs    <- (hl_mode == "gene_sets")
      genes_arg <- if (use_gs && !is.null(res$gene_sets)) res$gene_sets else NULL
      # thr_y is already a raw p-value (user enters 0.05, not 1.3)
      thr_y_raw <- if (thr_y > 0 && thr_y < 1) thr_y else NULL
      thr_x_arg <- if (thr_x > 0) thr_x else NULL
      contrast_names <- names(res$DEGs)

      if (!use_gs) {
        # --- Significant genes (threshold) mode ---
        lapply(contrast_names, function(cn) {
          df <- as.data.frame(res$DEGs[[cn]])
          padj_col <- if ("adj.P.Val" %in% colnames(df)) "adj.P.Val" else
                      if ("padj"      %in% colnames(df)) "padj"      else "adj.P.Val"
          lfc_col  <- if ("logFC"     %in% colnames(df)) "logFC"     else
                      if ("log2FoldChange" %in% colnames(df)) "log2FoldChange" else "logFC"
          df$gene_ <- rownames(df)
          df$lfc_  <- df[[lfc_col]]
          df$padj_ <- df[[padj_col]]
          df$logp_ <- -log10(pmax(df$padj_, 1e-300))
          sig <- rep(TRUE, nrow(df))
          if (!is.null(thr_y_raw)) sig <- sig & (df$padj_ <= thr_y_raw)
          if (!is.null(thr_x_arg)) sig <- sig & (abs(df$lfc_) >= thr_x_arg)
          if (is.null(thr_y_raw) && is.null(thr_x_arg)) sig <- rep(FALSE, nrow(df))
          df$.sig <- sig

          gg <- ggplot2::ggplot(df, ggplot2::aes(x = .data$lfc_, y = .data$logp_)) +
            ggplot2::geom_point(data = df[!df$.sig, , drop = FALSE],
                                color = "#B7B7B7", size = pts, alpha = 0.4) +
            ggplot2::geom_point(data = df[df$.sig,  , drop = FALSE],
                                color = "#6489B4", size = pts, alpha = 0.7) +
            ggplot2::labs(title = cn, x = lfc_col,
                          y = "-log10(Adj. p-value)") +
            ggplot2::theme_bw(base_size = lsz)
          if (!is.null(thr_y_raw))
            gg <- gg + ggplot2::geom_hline(yintercept = -log10(thr_y_raw),
                                           linetype = "dashed", color = "#555555")
          if (!is.null(thr_x_arg))
            gg <- gg + ggplot2::geom_vline(xintercept = c(-thr_x_arg, thr_x_arg),
                                           linetype = "dashed", color = "#555555")

          # Label the most extreme (by |logFC|) genes among the significant ones.
          if (n_top > 0L) {
            df_sig <- df[df$.sig, , drop = FALSE]
            if (nrow(df_sig) > 0L) {
              df_sig <- df_sig[order(abs(df_sig$lfc_), decreasing = TRUE), ]
              df_top <- head(df_sig, n_top)
              gg <- gg + ggrepel::geom_text_repel(
                data = df_top,
                ggplot2::aes(x = .data$lfc_, y = .data$logp_, label = .data$gene_),
                size = lsz / 3, max.overlaps = n_top * 2L, force = 10,
                color = "#243B53", fontface = "bold")
            }
          }
          gg
        })
      } else {
        # --- Gene set members mode: colour + labels restricted to set members ---
        lapply(contrast_names, function(cn) {
          df <- as.data.frame(res$DEGs[[cn]])
          lfc_col  <- if ("logFC"     %in% colnames(df)) "logFC"
                      else if ("log2FoldChange" %in% colnames(df)) "log2FoldChange"
                      else "logFC"
          padj_col <- if ("adj.P.Val" %in% colnames(df)) "adj.P.Val"
                      else if ("padj" %in% colnames(df)) "padj"
                      else "adj.P.Val"
          df$gene_ <- rownames(df)
          df$lfc_  <- df[[lfc_col]]
          df$padj_ <- df[[padj_col]]
          df$logp_ <- -log10(pmax(df[[padj_col]], 1e-300))

          gs_plots <- lapply(names(genes_arg), function(gsn) {
            gs_genes <- .gs_gene_names(genes_arg[[gsn]])
            gs_dirs  <- .gs_gene_directions(genes_arg[[gsn]])
            dir_map  <- stats::setNames(gs_dirs, gs_genes)[!duplicated(gs_genes)]
            in_gs   <- df$gene_ %in% gs_genes
            gene_dir <- unname(dir_map[df$gene_])
            upreg   <- in_gs & !is.na(gene_dir) & gene_dir == 1
            downreg <- in_gs & !is.na(gene_dir) & gene_dir == -1
            noinfo  <- in_gs & !upreg & !downreg

            gg <- ggplot2::ggplot(df, ggplot2::aes(x = .data$lfc_, y = .data$logp_)) +
              ggplot2::geom_point(data = df[!in_gs, , drop = FALSE],
                                  color = "#B7B7B7", size = pts, alpha = 0.4) +
              ggplot2::geom_point(data = df[noinfo, , drop = FALSE],
                                  color = "#05254A", size = pts, alpha = 0.8) +
              ggplot2::geom_point(data = df[upreg, , drop = FALSE],
                                  color = "#038C65", size = pts, alpha = 0.8) +
              ggplot2::geom_point(data = df[downreg, , drop = FALSE],
                                  color = "#8C0303", size = pts, alpha = 0.8) +
              ggplot2::labs(title    = gsn,
                            subtitle = cn,
                            x = lfc_col,
                            y = "-log10(Adj. p-value)") +
              ggplot2::theme_bw(base_size = lsz)

            # Label the most extreme members of THIS gene set only, optionally
            # narrowed further by the logFC / padj thresholds.
            if (n_top > 0L) {
              df_gs <- df[in_gs, , drop = FALSE]
              if (!is.null(thr_x_arg)) df_gs <- df_gs[abs(df_gs$lfc_) >= thr_x_arg, , drop = FALSE]
              if (!is.null(thr_y_raw)) df_gs <- df_gs[df_gs$padj_ <= thr_y_raw, , drop = FALSE]
              if (nrow(df_gs) > 0L) {
                df_gs  <- df_gs[order(abs(df_gs$lfc_), decreasing = TRUE), ]
                df_top <- head(df_gs, n_top)
                gg <- gg + ggrepel::geom_text_repel(
                  data = df_top,
                  ggplot2::aes(x = .data$lfc_, y = .data$logp_, label = .data$gene_),
                  size = lsz / 3, max.overlaps = n_top * 2L, force = 10,
                  color = "#05254A", fontface = "bold"
                )
              }
            }
            gg
          })
          if (length(gs_plots) == 1L) gs_plots[[1]]
          else ggpubr::ggarrange(plotlist = gs_plots,
                                 nrow = length(gs_plots), ncol = 1L, align = "h")
        })
      }
    }

    output$gsea_deg_volcano_plot <- shiny::renderPlot({
      res <- gsea_result(); shiny::req(!is.null(res), !is.null(res$DEGs))
      thr_x    <- input$vol_threshold_x    %||% 1
      thr_y    <- input$vol_threshold_y    %||% 0.05  # raw padj threshold
      n_top    <- as.integer(input$vol_top_n   %||% 0L)
      pts      <- input$vol_pointsize      %||% 2
      lsz      <- input$vol_labsize        %||% 8L
      hl_mode  <- input$vol_highlight_mode %||% "gene_sets"

      contrast_names <- names(res$DEGs)
      n_cont <- length(contrast_names)
      n_gs   <- if (hl_mode == "gene_sets" && !is.null(res$gene_sets))
        length(res$gene_sets) else 1L
      n_rows_vol <- n_gs; n_cols_vol <- n_cont
      gsea_vol_h(as.integer(min(10000L, max(460L, n_rows_vol * 470L + 70L))))
      gsea_vol_w(if (n_cols_vol >= 3L) "100%"
                 else paste0(min(100L, n_cols_vol * 50L), "%"))

      tryCatch({
        per_plots <- .build_deg_volcano_plots(res, thr_x, thr_y, n_top, pts, lsz, hl_mode)
        combined_vol <- if (length(per_plots) == 1L) {
          per_plots[[1]]
        } else {
          ggpubr::ggarrange(plotlist = per_plots,
                            ncol = length(per_plots), nrow = 1L, align = "h")
        }
        print(combined_vol)
      },
      error = function(e)
        shiny::showNotification(paste("Volcano error:", e$message),
                                type = "warning", duration = 8))
    }, height = function() gsea_vol_h(), res = 150)

    output$dl_gsea_deg_volcano <- shiny::downloadHandler(
      filename = function() paste0("deg_volcano_", Sys.Date(), ".png"),
      content  = function(file) {
        res     <- gsea_result(); shiny::req(!is.null(res), !is.null(res$DEGs))
        thr_x   <- input$vol_threshold_x    %||% 1
        thr_y   <- input$vol_threshold_y    %||% 0.05
        n_top   <- as.integer(input$vol_top_n   %||% 0L)
        pts     <- input$vol_pointsize      %||% 2
        lsz     <- input$vol_labsize        %||% 8L
        hl_mode <- input$vol_highlight_mode %||% "gene_sets"

        per_plots <- .build_deg_volcano_plots(res, thr_x, thr_y, n_top, pts, lsz, hl_mode)
        combined_vol <- if (length(per_plots) == 1L) per_plots[[1]] else
          ggpubr::ggarrange(plotlist = per_plots,
                            ncol = length(per_plots), nrow = 1L, align = "h")
        h_in <- gsea_vol_h() / 150
        ggplot2::ggsave(file, plot = combined_vol,
                        width = 10, height = max(4, h_in), dpi = 150, units = "in")
      }
    )

    # ---- GSEA: Enrichment Plots ---------------------------------------------

    # App-only wrapper: plotGSEAenrichment() (exported, untouched) only lets
    # us control `titlesize` (title + subtitle), not axis text/title size. We
    # call it with grid = FALSE to get the per (gene set x contrast) ggplot
    # list, extend each panel's theme with a configurable label size, then
    # arrange the grid ourselves (mirroring the exported function's own
    # nrow/ncol auto-layout logic when neither is supplied).
    .build_gsea_enrich_grid <- function(gsea_res, DEGs, gene_sets, titlesize,
                                        widthTitle, labsize, ncol = NULL, nrow = NULL) {
      plot_list <- plotGSEAenrichment(
        GSEA_results = gsea_res, DEGList = DEGs, gene_sets = gene_sets,
        titlesize = titlesize, widthTitle = widthTitle, grid = FALSE)
      plot_list <- lapply(plot_list, function(p)
        p + ggplot2::theme(
          axis.text  = ggplot2::element_text(size = labsize),
          axis.title = ggplot2::element_text(size = labsize + 1)
        ))
      n <- length(plot_list)
      if (n == 0L) return(NULL)
      if (n == 1L) return(plot_list[[1]])
      if (is.null(ncol) && is.null(nrow)) {
        ncol <- ceiling(sqrt(n)); nrow <- ceiling(n / ncol)
      } else if (is.null(ncol)) ncol <- ceiling(n / nrow)
      else if (is.null(nrow)) nrow <- ceiling(n / ncol)
      ggpubr::ggarrange(plotlist = plot_list, ncol = ncol, nrow = nrow, align = "hv")
    }

    output$gsea_enrich_plot <- shiny::renderPlot({
      res <- gsea_result(); shiny::req(!is.null(res))
      font_s  <- input$enrich_titlesize   %||% 9L
      lab_s   <- input$enrich_labsize     %||% 9L
      w_title <- input$enrich_width_title %||% 38L
      # Panels are one per (gene set x contrast) pair, flattened into a grid
      # capped at 2 columns; height scales with the number of rows that
      # produces. Each panel needs real vertical room (title, subtitle,
      # running-score curve, and the gene-set rug at the bottom), so the
      # per-row allowance is generous - too little and every panel becomes
      # illegibly squished.
      n_total  <- length(res$gsea_res) * length(res$gene_sets)
      ncol_enr <- min(2L, max(1L, n_total))
      nrow_enr <- ceiling(n_total / ncol_enr)
      enrich_h(as.integer(min(14000L, max(460L, nrow_enr * 400L + 60L))))
      tryCatch(
        print(.build_gsea_enrich_grid(
          res$gsea_res, res$DEGs, res$gene_sets, font_s, w_title, lab_s,
          nrow = nrow_enr, ncol = ncol_enr)),
        error = function(e)
          shiny::showNotification(paste("Enrichment plot error:", e$message),
                                  type = "warning", duration = 8))
    }, height = function() enrich_h(), res = 150)

    # ---- GSEA: NES Lollipop -------------------------------------------------

    # App-only wrapper: same idea as above - plotNESlollipop() (exported,
    # untouched) only exposes `titlesize`, not axis/legend text size.
    .build_nes_lollipop_grid <- function(gsea_res, sig_thr, titlesize, widthlabels,
                                         labsize, ncol = NULL, nrow = NULL) {
      plot_list <- plotNESlollipop(
        GSEA_results = gsea_res, sig_threshold = sig_thr,
        titlesize = titlesize, widthlabels = widthlabels, grid = FALSE)
      plot_list <- lapply(plot_list, function(p)
        p + ggplot2::theme(
          axis.text    = ggplot2::element_text(size = labsize),
          axis.title   = ggplot2::element_text(size = labsize + 1),
          legend.text  = ggplot2::element_text(size = labsize),
          legend.title = ggplot2::element_text(size = labsize),
          strip.text   = ggplot2::element_text(size = labsize)
        ))
      n <- length(plot_list)
      if (n == 0L) return(NULL)
      if (n == 1L) return(plot_list[[1]])
      if (is.null(ncol) && is.null(nrow)) {
        ncol <- ceiling(sqrt(n)); nrow <- ceiling(n / ncol)
      } else if (is.null(ncol)) ncol <- ceiling(n / nrow)
      else if (is.null(nrow)) nrow <- ceiling(n / ncol)
      ggpubr::ggarrange(plotlist = plot_list, ncol = ncol, nrow = nrow, align = "hv",
                        common.legend = TRUE, legend = "right")
    }

    output$nes_lollipop <- shiny::renderPlot({
      res <- gsea_result(); shiny::req(!is.null(res))
      sig_thr  <- input$nes_sig_threshold  %||% 0.05
      font_s   <- input$nes_titlesize      %||% 8L
      lab_s    <- input$nes_labsize        %||% 8L
      w_labels <- input$nes_width_legend   %||% 32L
      tryCatch(
        print(.build_nes_lollipop_grid(
          res$gsea_res, sig_thr, font_s, w_labels, lab_s, ncol = 1L)),
        error = function(e)
          shiny::showNotification(paste("Lollipop error:", e$message),
                                  type = "warning", duration = 8))
    }, height = function() nes_h(), res = 150)

    # ---- GSEA: Combined Volcano ----------------------------------------------

    # Build the combined GSEA volcano as a single static ggplot: one point per
    # gene set x contrast, coloured by gene set and shaped by contrast, with
    # both shown as proper legends.
    .build_combined_gsea_gg <- function(gsea_res, sig_thr, pt_size, w_legend, labsize = 8) {
      combined_data <- do.call(rbind, lapply(names(gsea_res), function(cn) {
        df <- gsea_res[[cn]]; df$contrast <- cn; df
      }))
      combined_data$logpadj <- -log10(combined_data$padj)
      combined_data$pathway_w <- vapply(combined_data$pathway,
        function(x) wrap_title(x, w_legend), character(1))

      n_paths <- length(unique(combined_data$pathway_w))
      pal <- grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(min(12L, max(3L, n_paths)), "Set3")
      )(n_paths)

      n_cont  <- length(unique(combined_data$contrast))
      shapes  <- 15:(15 - 1L + n_cont)

      ggplot2::ggplot(combined_data,
          ggplot2::aes(x = NES, y = logpadj,
                       colour = factor(pathway_w),
                       shape  = contrast)) +
        # A slightly larger black point drawn first, then the coloured point on
        # top, gives every point a thin black outline (the plotting shapes used
        # here, 15+, are solid/non-fillable so this "shadow" trick stands in
        # for a real stroke colour).
        ggplot2::geom_point(size = pt_size + 1, color = "black", show.legend = FALSE) +
        ggplot2::geom_point(size = pt_size) +
        ggplot2::geom_hline(yintercept = -log10(sig_thr),
                            linetype = "dashed", color = "black", linewidth = 0.5) +
        ggplot2::geom_vline(xintercept = 0,
                            linetype = "dashed", color = "black", linewidth = 0.5) +
        ggplot2::scale_colour_manual(values = pal, name = "Gene set") +
        ggplot2::scale_shape_manual(values = shapes, name = "Contrast") +
        # With a facet_grid of two side-by-side panels, a right-hand legend
        # column has very little horizontal room - a "Gene set" legend with
        # many entries stacked in a single narrow column can end up taller
        # than the plot device itself, clipping its own title. Moving the
        # legend to the bottom lets it wrap across the plot's full width
        # instead, and stacking the two legends (colour, then shape) with
        # explicit spacing between them keeps them from crowding each other.
        ggplot2::guides(
          colour = ggplot2::guide_legend(title.position = "top", nrow = 3, order = 1),
          shape  = ggplot2::guide_legend(title.position = "top", nrow = 1, order = 2)
        ) +
        ggplot2::labs(x = "Normalized Enrichment Score (NES)",
                      y = "-log10(adj. p-value)") +
        ggplot2::theme_bw() +
        ggplot2::facet_grid(. ~ stat_used,
          labeller = ggplot2::labeller(
            stat_used = c("t" = "Enriched/Depleted", "B" = "Altered")),
          scales = "free", switch = "y") +
        ggplot2::theme(strip.background = ggplot2::element_rect(fill = "white"),
                       legend.position  = "bottom",
                       legend.box       = "vertical",
                       legend.spacing.y = ggplot2::unit(8, "pt"),
                       legend.margin    = ggplot2::margin(t = 4, b = 2),
                       plot.margin   = ggplot2::margin(t = 8, r = 14, b = 8, l = 8),
                       axis.text     = ggplot2::element_text(size = labsize),
                       axis.title    = ggplot2::element_text(size = labsize + 1),
                       legend.text   = ggplot2::element_text(size = labsize),
                       legend.title  = ggplot2::element_text(size = labsize, hjust = 0),
                       strip.text    = ggplot2::element_text(size = labsize + 1, face = "bold"))
    }

    output$gsea_combined_plot <- shiny::renderPlot({
      res <- gsea_result(); shiny::req(!is.null(res), !is.null(res$gsea_res))
      sig_thr  <- input$combined_sig_threshold %||% 0.05
      pt_size  <- input$combined_pointsize     %||% 3L
      w_legend <- input$combined_width_legend  %||% 32L
      lab_s    <- input$combined_labsize       %||% 8L
      tryCatch(
        print(.build_combined_gsea_gg(res$gsea_res, sig_thr, pt_size, w_legend, lab_s)),
        error = function(e)
          shiny::showNotification(paste("Combined GSEA error:", e$message),
                                  type = "warning", duration = 8))
    }, height = function() combined_h(), res = 150)

    # ---- GSEA: Download handlers --------------------------------------------

    output$dl_gsea_enrich <- shiny::downloadHandler(
      filename = function() paste0("gsea_enrichment_", Sys.Date(), ".png"),
      content  = function(file) {
        res <- gsea_result(); shiny::req(!is.null(res))
        font_s  <- input$enrich_titlesize   %||% 9L
        lab_s   <- input$enrich_labsize     %||% 9L
        w_title <- input$enrich_width_title %||% 38L
        n_total  <- length(res$gsea_res) * length(res$gene_sets)
        ncol_enr <- min(2L, max(1L, n_total))
        nrow_enr <- ceiling(n_total / ncol_enr)
        h_in    <- enrich_h() / 96
        ggplot2::ggsave(file,
          plot = .build_gsea_enrich_grid(
            res$gsea_res, res$DEGs, res$gene_sets, font_s, w_title, lab_s,
            nrow = nrow_enr, ncol = ncol_enr),
          width = 14, height = max(5, h_in), dpi = 150, units = "in")
      }
    )

    output$dl_nes_lollipop <- shiny::downloadHandler(
      filename = function() paste0("nes_lollipop_", Sys.Date(), ".png"),
      content  = function(file) {
        res <- gsea_result(); shiny::req(!is.null(res))
        sig_thr  <- input$nes_sig_threshold %||% 0.05
        font_s   <- input$nes_titlesize     %||% 8L
        lab_s    <- input$nes_labsize       %||% 8L
        w_labels <- input$nes_width_legend  %||% 32L
        h_in <- nes_h() / 96
        ggplot2::ggsave(file,
          plot = .build_nes_lollipop_grid(res$gsea_res, sig_thr, font_s, w_labels,
                                          lab_s, ncol = 1L),
          width = 12, height = max(4, h_in), dpi = 150, units = "in")
      }
    )

    output$dl_gsea_combined <- shiny::downloadHandler(
      filename = function() paste0("gsea_combined_", Sys.Date(), ".png"),
      content  = function(file) {
        res <- gsea_result(); shiny::req(!is.null(res), !is.null(res$gsea_res))
        sig_thr  <- input$combined_sig_threshold %||% 0.05
        pt_size  <- input$combined_pointsize     %||% 3L
        w_legend <- input$combined_width_legend  %||% 32L
        lab_s    <- input$combined_labsize       %||% 8L
        h_in <- combined_h() / 96
        ggplot2::ggsave(file,
          plot = .build_combined_gsea_gg(res$gsea_res, sig_thr, pt_size, w_legend, lab_s),
          width = 12, height = max(5, h_in), dpi = 150, units = "in")
      }
    )

    output$gsea_table_hdr_ui <- shiny::renderUI({
      res <- gsea_result(); shiny::req(!is.null(res))
      p <- res$run_params; shiny::req(!is.null(p))
      n_rows <- sum(vapply(res$gsea_res, nrow, integer(1)))
      shiny::div(
        style = "margin-bottom:6px; font-size:0.85em; padding:4px 6px;",
        shiny::icon("circle-info", style = "color:#0d6efd; margin-right:4px;"),
        shiny::tags$span(style = "color:#555; margin-right:6px;",
          sprintf("%d result(s) across %d contrast(s):", n_rows, res$n_cont)),
        shiny::tags$span(class = "badge bg-primary me-1",  paste0("Variable: ", p$enrich_var)),
        shiny::tags$span(class = "badge bg-secondary me-1", paste0("Mode: ", p$enrich_mode)),
        shiny::tags$span(class = "badge bg-secondary me-1", paste0("p-adjust method: ", p$padj_m)),
        shiny::tags$span(class = "badge bg-secondary me-1", paste0("Gene sets: ", p$n_gene_sets))
      )
    })

    output$gsea_table <- DT::renderDT({
      res <- gsea_result(); shiny::req(!is.null(res))
      combined <- do.call(rbind, lapply(names(res$gsea_res), function(ct) {
        df <- res$gsea_res[[ct]]; df$contrast <- ct
        if ("leadingEdge" %in% colnames(df))
          df$leadingEdge <- vapply(df$leadingEdge,
                                   function(x) paste(x, collapse = ", "), character(1))
        df
      }))
      # fgsea() results are data.tables, and rbind() of data.tables returns a
      # data.table too. Coerce to a plain data.frame before column-selecting
      # below: a data.table's `[` treats a logical vector as a ROW filter
      # (matched against nrow), not a column selector like data.frame's `[`
      # does, which throws "i evaluates to a logical vector length <ncol> but
      # there are <nrow> rows" the moment ncol != nrow.
      combined <- as.data.frame(combined)
      combined[sapply(combined, is.numeric)] <-
        lapply(combined[sapply(combined, is.numeric)], round, 4)
      DT::datatable(combined, rownames = FALSE, filter = "top",
                    extensions = "Buttons",
                    options = list(pageLength = 15, dom = "Bfrtip",
                                   buttons = c("csv", "excel"), scrollX = TRUE))
    })

    # ---- Renders: FPR ------------------------------------------------------

    output$fpr_status_ui <- shiny::renderUI({
      if (is.null(fpr_result()))
        shiny::div(class = "alert alert-info", style = "margin:10px 0;",
                   shiny::icon("circle-info"),
                   " Set the variable above (in Score Analysis shared settings), select a simulation count, then click ",
                   shiny::strong("Run Score FPR Simulation"), ".")
    })

    # Helper: combine per (gene set x contrast) panels (built by
    # .build_fpr_plots(), defined above near the FPR run observer) into a
    # single grid, filtered to the currently selected contrasts.
    .build_fpr_grid <- function(data_list, fontsize = 10L, cohentype = "d",
                                ncol_g = 3L, sel_conts = NULL) {
      if (!is.null(sel_conts) && length(sel_conts) > 0L) {
        data_list <- lapply(data_list, function(df)
          df[df$contrast %in% sel_conts, , drop = FALSE])
        data_list <- Filter(function(df) nrow(df) > 0L, data_list)
      }
      if (length(data_list) == 0L) return(NULL)

      ylab <- if (cohentype == "d") "|Cohen's d|" else "|Cohen's f|"
      plot_list <- .build_fpr_plots(data_list, ylab = ylab,
                                    titlesize = fontsize, labsize = fontsize)
      n <- length(plot_list)
      if (n == 0L) return(NULL)
      if (n == 1L) return(plot_list[[1]])

      nrow_g <- ceiling(n / ncol_g)
      ggpubr::ggarrange(plotlist = plot_list, ncol = ncol_g, nrow = nrow_g,
                        common.legend = TRUE, legend = "top")
    }

    output$fpr_plot_ui <- shiny::renderUI({
      res <- fpr_result(); shiny::req(!is.null(res))
      shiny::tagList(
        shiny::tags$style(
          ".bm-fpr-card { overflow: visible !important; height: auto !important;",
          " max-height: none !important; }",
          ".bm-fpr-card .card-body { overflow: visible !important; }"
        ),
        bslib::card(
          class = "bm-fpr-card",
          full_screen = TRUE,
          .bm_card_header(
            sprintf("Score FPR: Cohen's %s vs null (%d simulations), %s",
                    res$cohentype %||% "d", res$nsims, res$var),
            "dl_fpr", ns),
          shiny::div(
            style = "overflow: visible;",
            shiny::plotOutput(ns("fpr_plot"), height = paste0(fpr_h(), "px"))
          )
        )
      )
    })

    output$fpr_plot <- shiny::renderPlot({
      res <- fpr_result(); shiny::req(!is.null(res), !is.null(res$fpr_res$data))
      fontsize  <- input$fpr_fontsize %||% 10L
      cohentype <- res$cohentype %||% "d"
      # Reactively filter cached FPR data by selected contrasts (no re-run needed).
      # sel_conts is a character vector of contrast strings (e.g. "A-B", "A-C").
      sel_conts <- input$score_selected_groups
      tryCatch({
        p <- .build_fpr_grid(res$fpr_res$data, fontsize, cohentype,
                             ncol_g = 3L, sel_conts = sel_conts)
        shiny::req(!is.null(p))
        print(p)
      }, error = function(e)
        shiny::showNotification(paste("FPR plot error:", e$message),
                                type = "warning", duration = 8))
    }, height = function() fpr_h(), res = 150)

    output$dl_fpr <- shiny::downloadHandler(
      filename = function() paste0("fpr_simulation_", Sys.Date(), ".png"),
      content  = function(file) {
        res <- fpr_result(); shiny::req(!is.null(res), !is.null(res$fpr_res$data))
        fontsize  <- input$fpr_fontsize %||% 10L
        cohentype <- res$cohentype %||% "d"
        sel_conts <- input$score_selected_groups
        p    <- .build_fpr_grid(res$fpr_res$data, fontsize, cohentype,
                                ncol_g = 3L, sel_conts = sel_conts)
        h_in <- fpr_h() / 150
        ggplot2::ggsave(file, plot = p,
                        width = 12, height = max(4, h_in),
                        dpi = 150, units = "in")
      }
    )

    # ---- Run + Render: Score Results Table ---------------------------------

    shiny::observeEvent(input$run_scores_table, {
      expr <- get_expr(); meta <- get_meta()
      gs   <- .active_gs(); shiny::req(expr, meta, gs, length(gs) > 0)
      meth <- shiny::isolate(input$scores_table_method) %||% "all"

      shiny::withProgress(message = "Computing score results table...", value = 0.3, {
        result <- tryCatch({
          scores_list <- CalculateScores(data = expr, metadata = meta,
                                         gene_sets = gs, method = meth)
          shiny::incProgress(0.6)

          # Normalise to a nested list of method -> signature -> data.frame,
          # regardless of whether one method or "all" were requested, then
          # flatten into a single long table (one row per sample x gene set
          # x method) with a Method / GeneSet column prepended.
          by_method <- if (identical(meth, "all")) scores_list
                       else stats::setNames(list(scores_list), meth)

          long_df <- data.frame()
          for (m in names(by_method)) {
            for (sig in names(by_method[[m]])) {
              df <- by_method[[m]][[sig]]
              if (!is.data.frame(df) || nrow(df) == 0) next
              df$Method  <- m
              df$GeneSet <- sig
              long_df <- rbind(long_df, df)
            }
          }
          shiny::req(nrow(long_df) > 0)
          front <- c("Method", "GeneSet", "sample", "score")
          front <- front[front %in% colnames(long_df)]
          long_df <- long_df[, c(front, setdiff(colnames(long_df), front)), drop = FALSE]
          shiny::incProgress(0.1)
          long_df
        }, error = function(e) {
          shiny::showNotification(paste("Score results table failed:", conditionMessage(e)),
                                  type = "error", duration = 12); NULL
        })
      })
      scores_table_result(result)
    })

    output$scores_table_status_ui <- shiny::renderUI({
      if (is.null(scores_table_result()))
        shiny::div(class = "alert alert-info", style = "margin:10px 0;",
                   shiny::icon("circle-info"),
                   " Choose method(s) above, then click",
                   shiny::strong(" Compute Score Results Table"), ".")
    })

    output$scores_table_ui <- shiny::renderUI({
      res <- scores_table_result(); shiny::req(!is.null(res))
      bslib::card(
        full_screen = TRUE,
        shiny::div(
          class = "d-flex align-items-center justify-content-between",
          style = "margin: 4px 0 6px;",
          shiny::h5("Score Results", style = "margin:0;"),
          shiny::downloadButton(ns("dl_scores_table"), "Download",
            icon = shiny::icon("download"), class = "btn-sm btn-outline-secondary")
        ),
        shiny::div(
          class = "text-muted", style = "font-size:0.78em; padding:0 6px 8px;",
          shiny::icon("circle-info"),
          " One row per sample x gene set x scoring method, with all sample metadata included.",
          " Use the column filters below the header to narrow down results."
        ),
        DT::DTOutput(ns("scores_table"))
      )
    })

    output$scores_table <- DT::renderDT({
      res <- scores_table_result(); shiny::req(!is.null(res), is.data.frame(res))
      df <- res
      if ("score" %in% colnames(df)) df$score <- round(df$score, 4)
      DT::datatable(
        df, rownames = FALSE, filter = "top", extensions = "Buttons",
        options = list(pageLength = 15, dom = "Bfrtip",
                      buttons = c("csv", "excel"), scrollX = TRUE)
      )
    })

    output$dl_scores_table <- shiny::downloadHandler(
      filename = function() paste0("score_results_table_", Sys.Date(), ".csv"),
      content  = function(file) {
        res <- scores_table_result(); shiny::req(!is.null(res), is.data.frame(res))
        utils::write.csv(res, file, row.names = FALSE)
      }
    )

  })
}
