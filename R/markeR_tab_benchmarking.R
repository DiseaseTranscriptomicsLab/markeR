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
        shiny::h4("Benchmarking Mode", style = "font-weight:700; color:#1E497B; margin-bottom:4px;"),
        shiny::p(
          "Evaluate gene set discriminatory power using score distributions,",
          " ROC/AUC curves, GSEA enrichment, and null-model FPR simulations.",
          style = "color:#6c757d; font-size:0.85em;"
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
      ),

      shiny::hr(),

      shiny::div(
        style = "font-size:0.78em; color:#888; line-height:1.5;",
        shiny::tags$b("Score Analysis:"),
        " score distributions, ROC/AUC classification, all-methods summary,",
        " and FPR simulation (sub-tabs).",
        shiny::tags$br(),
        shiny::tags$b("Enrichment (GSEA):"),
        " gene-set-level enrichment via limma + fgsea."
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
                    )
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
                                        value = 11, min = 6, max = 20, step = 1)
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

          # ---- Display options (collapsed) ------------------------------------
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
                shiny::numericInput(ns("gsea_sig_threshold"), "Significance threshold (padj):",
                                    value = 0.05, min = 0.001, max = 0.5, step = 0.01),
                shiny::numericInput(ns("gsea_point_size"), "Point size:",
                                    value = 7, min = 1, max = 16, step = 1),
                bslib::tooltip(
                  shiny::numericInput(ns("gsea_title_size"), "Title size (pt):",
                                      value = 10, min = 6, max = 20, step = 1),
                  "Font size for plot titles and axis labels.",
                  placement = "right"
                ),
                bslib::tooltip(
                  shiny::numericInput(ns("gsea_width_title"), "Title wrap (chars):",
                                      value = 32, min = 10, max = 60, step = 2),
                  "Maximum characters per line in enrichment plot titles.",
                  placement = "right"
                ),
                bslib::tooltip(
                  shiny::numericInput(ns("gsea_width_legend"), "Legend wrap (chars):",
                                      value = 32, min = 8, max = 60, step = 2),
                  "Maximum characters per line in legend and lollipop labels.",
                  placement = "right"
                )
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

    # Dynamic plot heights
    score_h  <- shiny::reactiveVal(500L)
    roc_h    <- shiny::reactiveVal(600L)
    auc_h    <- shiny::reactiveVal(400L)
    all_h    <- shiny::reactiveVal(500L)
    nes_h      <- shiny::reactiveVal(500L)
    enrich_h   <- shiny::reactiveVal(500L)
    gsea_vol_h <- shiny::reactiveVal(500L)
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
    # the ROC contrast dropdown is populated from this list.
    # Note: Score Distributions and All Methods always show every contrast
    # because those functions build contrasts internally and cannot be
    # post-hoc filtered without re-running.

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
          "Filters the FPR plots and the ROC contrast dropdown.",
          " Score Distributions and All Methods always compute every contrast."
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

    # (score_filter_meta removed — contrast filtering is now display-only via
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

    # ---- Run: Score Distributions -----------------------------------------

    shiny::observeEvent(input$run_scores, {
      expr <- get_expr(); meta <- get_meta()
      gs   <- .active_gs(); shiny::req(expr, meta, gs, length(gs) > 0)

      p    <- .score_shared_params()
      meth <- shiny::isolate(input$score_method) %||% "logmedian"
      do_cohen <- shiny::isolate(input$compute_cohen) %||% TRUE

      # Build palette-matched ColorValues.
      # PlotScores applies ColorValues as scale_color_manual — the names must
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

      shiny::withProgress(message = "Computing score distributions...", value = 0, {
        result <- tryCatch({
          shiny::incProgress(0.4, detail = paste("Scoring:", meth))
          plot_res <- PlotScores(
            data          = expr, metadata = meta, gene_sets = gs,
            method        = meth, Variable = p$var_arg,
            ColorVariable = p$color_arg,
            ColorValues   = score_color_vals,
            mode          = p$mode, compute_cohen = do_cohen,
            cond_cohend   = p$cond_cohend, pvalcalc = p$pvalcalc,
            ConnectGroups = p$connect_g, cor = p$cor_m,
            pointSize     = shiny::isolate(input$score_pointsize) %||% 4L,
            widthTitle    = shiny::isolate(input$score_width_title)  %||% 32L,
            widthlegend   = shiny::isolate(input$score_width_legend) %||% 32L,
            ncol          = 3L
          )
          shiny::incProgress(0.6, detail = "Done.")
          list(plot_res = plot_res, n_gs = length(gs),
               score_var = shiny::isolate(input$score_var),
               is_cat = p$is_cat, method = meth)
        }, error = function(e) {
          shiny::showNotification(paste("Score plot failed:", conditionMessage(e)),
                                  type = "error", duration = 12); NULL
        })
      })
      if (!is.null(result))
        score_h(as.integer(min(4000L, 400L + ceiling(result$n_gs / 3L) * 400L)))
      score_dist_result(result)
    })

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

      shiny::withProgress(message = "Computing ROC curves and AUC...", value = 0, {
        result <- tryCatch({
          roc_wid_title <- shiny::isolate(input$roc_width_title) %||% 32L
          roc_titlesize <- shiny::isolate(input$roc_titlesize)   %||% 11L

          shiny::incProgress(0.4, detail = "ROC curves...")
          # Use the raw calculate function so we can build per-contrast plots reactively
          roc_raw <- ROCAUC_Scores_Calculate(
            data = expr, metadata = meta, gene_sets = gs,
            method = meth, variable = p$var_arg, mode = p$mode)

          shiny::incProgress(0.4, detail = "AUC heatmap...")
          auc_res <- AUC_Scores(data = expr, metadata = meta, gene_sets = gs,
                                method = meth, variable = p$var_arg, mode = p$mode,
                                titlesize = roc_titlesize)
          shiny::incProgress(0.2, detail = "Done.")

          # Extract contrast names from raw data (they're in [[method]][[sig]][[contrast]])
          roc_contrasts <- names(roc_raw[[1]][[1]])

          list(roc_raw      = roc_raw,
               roc_contrasts = roc_contrasts,
               auc_res       = auc_res,
               n_gs          = length(gs),
               n_cont        = n_cont,
               roc_wid_title = roc_wid_title,
               roc_titlesize = roc_titlesize)
        }, error = function(e) {
          shiny::showNotification(paste("ROC/AUC failed:", conditionMessage(e)),
                                  type = "error", duration = 12); NULL
        })
      })
      if (!is.null(result)) {
        # Height: one contrast displayed at a time => n_gs panels, 3 cols
        n_cols_roc <- min(3L, result$n_gs)
        n_rows_roc <- ceiling(result$n_gs / n_cols_roc)
        roc_h(as.integer(min(6000L, max(400L, n_rows_roc * 450L))))
        auc_h(as.integer(min(3000L, max(300L, result$n_gs * 250L))))
      }
      score_roc_result(result)
    })

    # ---- Run: All Methods Summary ------------------------------------------

    shiny::observeEvent(input$run_all_methods, {
      expr <- get_expr(); meta <- get_meta()
      gs   <- .active_gs(); shiny::req(expr, meta, gs, length(gs) > 0)
      p <- .score_shared_params()

      # Read display params (used both in PlotScores and for color name matching)
      all_wid_title  <- shiny::isolate(input$all_width_title)  %||% 32L
      all_wid_legend <- shiny::isolate(input$all_width_legend) %||% 32L

      # ColorValues must be a named list:
      #   [[1]]/"heatmap" → 2-colour gradient for Heatmap_Cohen
      #   [[2]]/"volcano" → named per-signature vector for Volcano_Cohen
      # Volcano_Cohen wraps signature names with wrap_title(width=widthlegend)
      # before building the data frame, so our color names must use the same
      # wrapping to avoid mismatches that produce grey points.
      gs_names    <- names(gs)
      gs_wrapped  <- vapply(gs_names, wrap_title, character(1L), width = all_wid_legend)
      gs_colors   <- list(
        heatmap = c("#F9F4AE", "#B44141"),   # diverging gradient (default)
        volcano = stats::setNames(rep_len(.pp_palette, length(gs_names)), gs_wrapped)
      )

      shiny::withProgress(message = "Running all scoring methods...", value = 0, {
        result <- tryCatch({
          shiny::incProgress(0.6, detail = "Computing Cohen's d/f across all methods...")
          plot_res <- PlotScores(
            data            = expr, metadata = meta, gene_sets = gs,
            method          = "all", Variable = p$var_arg,
            ColorValues     = gs_colors,
            mode            = p$mode, compute_cohen = TRUE,
            cond_cohend     = p$cond_cohend,
            ConnectGroups   = p$connect_g, cor = p$cor_m,
            sig_threshold   = shiny::isolate(input$sig_threshold)   %||% 0.05,
            cohen_threshold = shiny::isolate(input$cohen_threshold)  %||% 0.5,
            widthTitle      = all_wid_title,
            widthlegend     = all_wid_legend
          )
          shiny::incProgress(0.4, detail = "Done.")
          list(plot_res = plot_res, n_gs = length(gs))
        }, error = function(e) {
          shiny::showNotification(paste("All methods failed:", conditionMessage(e)),
                                  type = "error", duration = 12); NULL
        })
      })
      if (!is.null(result))
        all_h(as.integer(min(4000L, max(350L, result$n_gs * 80L + 300L))))
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
            uvals <- gsub(" ", "", unique(as.character(meta[[enrich_var]])))
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
        nes_h(as.integer(min(5000L, max(350L, result$n_sets * result$n_cont * 45L))))
        enrich_h(as.integer(min(6000L, max(400L, result$n_sets * result$n_cont * 220L))))
        gsea_vol_h(as.integer(min(5000L, max(400L, result$n_cont * 350L))))
      }
      gsea_result(result)
    })

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
            # FPR_Simulation accepts ncol/nrow to arrange gene set panels in a grid.
            # Calling it once for all gene sets is the intended usage and avoids
            # ggarrange-of-annotate_figure issues that can collapse to 1 column.
            ncol_g <- 3L
            nrow_g <- ceiling(n_gs / ncol_g)
            shiny::setProgress(value = 0.1, detail = "Running simulations...")
            res_all <- FPR_Simulation(
              data                = expr,
              metadata            = meta,
              original_signatures = gs,
              Variable            = fpr_var,
              number_of_sims      = fpr_nsims,
              mode                = fpr_mode,
              titlesize           = font_s,
              widthTitle          = 32L,
              ncol                = ncol_g,
              nrow                = nrow_g
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
        ncol_fpr    <- 3L
        n_rows_fpr  <- ceiling(result$n_gs / ncol_fpr)
        # Each gene set plot stacks n_contrasts facets vertically, ~220px each
        n_contrasts <- length(unique(unlist(
          lapply(result$fpr_res$data, function(df) unique(df$contrast)))))
        h_per_gs <- n_contrasts * 220L + 100L
        fpr_h(as.integer(min(12000L, max(600L, n_rows_fpr * h_per_gs + 120L))))
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
        uvals <- gsub(" ", "", unique(as.character(meta[[var]])))
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
        shiny::selectizeInput(
          ns("gsea_selected_contrasts"),
          label   = NULL,
          choices  = contrasts_all,
          selected = contrasts_all,
          multiple = TRUE,
          width    = "100%",
          options  = list(
            plugins     = list("remove_button"),
            placeholder = "Select contrasts…",
            maxItems    = length(contrasts_all)
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
      res <- score_dist_result(); shiny::req(!is.null(res))
      tryCatch(print(res$plot_res),
               error = function(e)
                 shiny::showNotification(paste("Score plot error:", e$message),
                                         type = "warning", duration = 8))
    }, height = function() score_h(), res = 150)

    output$dl_score_dist <- shiny::downloadHandler(
      filename = function() paste0("score_distributions_", Sys.Date(), ".png"),
      content  = function(file) {
        res <- score_dist_result(); shiny::req(!is.null(res))
        h_in <- score_h() / 96
        ggplot2::ggsave(file, plot = res$plot_res,
                        width = 12, height = max(4, h_in),
                        dpi = 150, units = "in")
      }
    )

    # ---- Renders: ROC / AUC ------------------------------------------------

    # ---- Helper: build ROC plot for one contrast from raw data ----------------
    .build_roc_contrast_plot <- function(roc_raw, sel_contrast,
                                         titlesize = 11L, wid_title = 32L) {
      roc_colors <- c(logmedian = "#3E5587", ssGSEA = "#B65285", ranking = "#B68C52")
      method_names <- names(roc_raw)
      sig_names    <- names(roc_raw[[1]])

      per_sig <- lapply(sig_names, function(sig) {
        combined_df <- data.frame()
        auc_values  <- list()
        for (meth_name in method_names) {
          roc_data <- roc_raw[[meth_name]][[sig]][[sel_contrast]]
          if (is.null(roc_data)) next
          combined_df <- rbind(combined_df, data.frame(
            FPR    = rev(1 - roc_data$ROC$specificities),
            TPR    = rev(roc_data$ROC$sensitivities),
            Method = meth_name
          ))
          auc_values[[meth_name]] <- roc_data$AUC
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
            size = 3, vjust = 0, hjust = 1, inherit.aes = FALSE, fill = "white") +
          ggplot2::labs(
            title    = wrap_title(sig, wid_title),
            subtitle = wrap_title(sel_contrast, wid_title),
            x = "False Positive Rate", y = "True Positive Rate") +
          ggplot2::theme_classic() +
          ggplot2::theme(
            legend.position = "none",
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
          " The displayed AUC is 1 minus raw AUC when < 0.5, restoring directionality.",
          " The AUC heatmap below always shows all contrasts."
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
          titlesize = res$roc_titlesize, wid_title = res$roc_wid_title)
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
                                          titlesize   = res$roc_titlesize,
                                          wid_title   = res$roc_wid_title)
        h_in <- roc_h() / 96
        ggplot2::ggsave(file, plot = p,
                        width = 14, height = max(5, h_in),
                        dpi = 150, units = "in")
      }
    )

    output$score_auc_ui <- shiny::renderUI({
      shiny::req(!is.null(score_roc_result()), !is.null(score_roc_result()$auc_res))
      bslib::card(
        full_screen = TRUE,
        .bm_card_header("AUC Heatmap (contrasts × methods)", "dl_auc", ns),
        shiny::div(
          class = "text-muted", style = "font-size:0.78em; padding:4px 6px 8px;",
          shiny::icon("circle-info"),
          " Numbers are ", shiny::tags$b("AUC values"),
          " (Area Under the ROC Curve). AUC = 0.5: no better than chance.",
          " AUC = 1.0: perfect separation. AUC < 0.5: reversed direction (still informative)."
        ),
        shiny::plotOutput(ns("score_auc_plot"), height = paste0(auc_h(), "px"))
      )
    })

    output$score_auc_plot <- shiny::renderPlot({
      res <- score_roc_result(); shiny::req(!is.null(res), !is.null(res$auc_res))
      tryCatch(print(res$auc_res),
               error = function(e)
                 shiny::showNotification(paste("AUC plot error:", e$message),
                                         type = "warning", duration = 8))
    }, height = function() auc_h(), res = 150)

    output$dl_auc <- shiny::downloadHandler(
      filename = function() paste0("auc_heatmap_", Sys.Date(), ".png"),
      content  = function(file) {
        res <- score_roc_result(); shiny::req(!is.null(res), !is.null(res$auc_res))
        h_in <- auc_h() / 96
        ggplot2::ggsave(file, plot = res$auc_res,
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
      shiny::req(!is.null(score_all_result()))
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
          " (numeric variable). |d| > 0.2 small, |d| > 0.5 medium, |d| > 0.8 large."
        ),
        shiny::plotOutput(ns("score_all_heatmap"), height = paste0(all_h(), "px"))
      )
    })

    output$score_all_heatmap <- shiny::renderPlot({
      res <- score_all_result(); shiny::req(!is.null(res), !is.null(res$plot_res$heatmap))
      tryCatch(print(res$plot_res$heatmap),
               error = function(e)
                 shiny::showNotification(paste("Heatmap error:", e$message),
                                         type = "warning", duration = 8))
    }, height = function() all_h(), res = 150)

    output$dl_all_heatmap <- shiny::downloadHandler(
      filename = function() paste0("cohen_heatmap_", Sys.Date(), ".png"),
      content  = function(file) {
        res <- score_all_result(); shiny::req(!is.null(res), !is.null(res$plot_res$heatmap))
        h_in <- all_h() / 96
        ggplot2::ggsave(file, plot = res$plot_res$heatmap,
                        width = 12, height = max(4, h_in),
                        dpi = 150, units = "in")
      }
    )

    # Helper: clean ggplotly legend for Volcano_Cohen plots.
    # ggplotly names combined-aesthetic traces as "(A,B)"; this renames them,
    # deduplicates within each aesthetic, and groups into Method / Signature sections.
    # "Signature" header is only shown once (above the first sig entry).
    .clean_volcano_plotly <- function(plt) {
      methods_known    <- c("logmedian", "ranking", "ssGSEA")
      seen_methods     <- character(0)
      seen_sigs        <- character(0)
      sig_header_shown <- FALSE

      plt$x$data <- lapply(plt$x$data, function(trace) {
        nm <- if (is.null(trace$name)) "" else trace$name
        # ggplotly names combined traces "(A,B)"
        if (grepl("^\\(.*,.*\\)$", nm)) {
          inner   <- sub("^\\((.*)\\)$", "\\1", nm)
          comma_p <- regexpr(",", inner, fixed = TRUE)[1]
          part_a  <- trimws(substr(inner, 1L, comma_p - 1L))

          # Thin point outline for cleaner look
          if (!is.null(trace$marker)) {
            trace$marker$line <- list(width = 0.5,
              color = if (!is.null(trace$marker$line$color))
                        trace$marker$line$color
                      else "rgba(0,0,0,0.35)")
          }

          if (part_a %in% methods_known) {
            # Shape / method trace
            trace$name             <- part_a
            trace$legendgroup      <- "Method"
            trace$legendgrouptitle <- list(text = "<b>Method</b>")
            if (part_a %in% seen_methods) {
              trace$showlegend <- FALSE
            } else {
              seen_methods <<- c(seen_methods, part_a)
              trace$showlegend <- TRUE
            }
          } else {
            # Colour / signature trace — all in one group so header shows once
            trace$name        <- part_a
            trace$legendgroup <- "Signatures"
            # Show the "Signature" header only on the very first sig trace
            if (!sig_header_shown) {
              trace$legendgrouptitle <- list(text = "<b>Signature</b>")
              sig_header_shown <<- TRUE
            } else {
              trace$legendgrouptitle <- list(text = "")
            }
            if (part_a %in% seen_sigs) {
              trace$showlegend <- FALSE
            } else {
              seen_sigs <<- c(seen_sigs, part_a)
              trace$showlegend <- TRUE
            }
          }
        }
        trace
      })

      # Remove ggtitle (renders oddly in faceted ggplotly) and tidy legend
      plt |>
        plotly::layout(
          title  = list(text = ""),
          legend = list(
            tracegroupgap = 14,
            font          = list(size = 11),
            itemsizing    = "constant",
            orientation   = "v"
          )
        )
    }

    # ---- GSEA combined-volcano legend cleanup --------------------------------
    # plotCombinedGSEA maps colour=pathway, shape=contrast.
    # ggplotly creates traces named "(pathway,contrast)" — one per combination.
    # Sometimes pathway_val is a numeric index ("1","2"…) instead of the name;
    # we resolve those using known_pathways.
    #
    # Strategy:
    #   Pass 1: resolve pathway names; assign all real traces to "Gene Set" group;
    #           deduplicate by pathway; thin point outline.
    #   Pass 2: add one dummy trace per contrast (fixed plotly symbol by position)
    #           so the "Contrast" shape legend is always correct and complete.
    .GSEA_PLOTLY_SYMBOLS <- c("square", "circle", "diamond", "triangle-up",
                               "cross", "star", "hexagram", "pentagon")

    # Tableau-10 + extra colours — enough for up to 20 gene sets
    .GSEA_COLOR_PALETTE <- c(
      "#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F",
      "#EDC948","#B07AA1","#FF9DA7","#9C755F","#BAB0AC",
      "#1F77B4","#FF7F0E","#2CA02C","#D62728","#9467BD",
      "#8C564B","#E377C2","#7F7F7F","#BCBD22","#17BECF"
    )

    .clean_gsea_volcano_plotly <- function(plt, known_pathways, known_contrasts) {
      seen_pathways  <- character(0)
      pathway_colors <- character(0)   # named: pathway_val -> hex color
      pathways_sorted <- sort(known_pathways)

      # ── Pass 0: remove Layer 1 black background dots ───────────────────────────
      # plotCombinedGSEA uses TWO geom_point layers:
      #   Layer 1: geom_point(colour="black", size=PointSize)  -- large black "outline"
      #   Layer 2: geom_point(aes(colour=pathway), size=PointSize-2.5) -- coloured
      # ggplotly names Layer 2 as "(pathway, contrast)" and Layer 1 only by contrast.
      # Filtering keeps coloured data traces and non-marker traces (hlines, vlines).
      is_data_tr <- function(tr)
        grepl("^\\(.*,.*\\)$", if (is.null(tr$name)) "" else tr$name)
      is_line_tr <- function(tr)
        grepl("^lines", if (is.null(tr$mode)) "" else tr$mode)
      plt$x$data <- Filter(function(tr) is_data_tr(tr) || is_line_tr(tr),
                           plt$x$data)

      n_real <- length(plt$x$data)

      # ── Pass 1: resolve names, assign palette colors, set Gene set group ───────
      plt$x$data <- lapply(plt$x$data, function(trace) {
        nm <- if (is.null(trace$name)) "" else trace$name

        if (grepl("^\\(.*,.*\\)$", nm)) {
          inner  <- sub("^\\((.*)\\)$", "\\1", nm)
          part_a <- trimws(sub(",.*$", "", inner))
          part_b <- trimws(sub("^[^,]*,", "", inner))

          if (part_a %in% known_contrasts) {
            pathway_val <- part_b
          } else {
            pathway_val <- part_a
          }

          if (grepl("^\\d+$", trimws(pathway_val))) {
            idx <- as.integer(trimws(pathway_val))
            if (idx >= 1L && idx <= length(pathways_sorted))
              pathway_val <- pathways_sorted[idx]
          }

          if (!(pathway_val %in% names(pathway_colors))) {
            ci <- (length(pathway_colors) %% length(.GSEA_COLOR_PALETTE)) + 1L
            pathway_colors[[pathway_val]] <<- .GSEA_COLOR_PALETTE[ci]
          }
          pal_col <- pathway_colors[[pathway_val]]

          if (is.null(trace$marker)) trace$marker <- list()
          trace$marker$color <- pal_col

          trace$name             <- pathway_val
          trace$legendgroup      <- "Gene set"
          trace$legendgrouptitle <- list(text = "Gene set")
          trace$showlegend       <- !(pathway_val %in% seen_pathways)
          seen_pathways          <<- unique(c(seen_pathways, pathway_val))
        }
        trace
      })

      # ── Pass 2: dummy contrast traces in "Contrast" legend section ─────────────
      for (i in seq_along(known_contrasts)) {
        cv  <- known_contrasts[i]
        sym <- .GSEA_PLOTLY_SYMBOLS[((i - 1L) %% length(.GSEA_PLOTLY_SYMBOLS)) + 1L]
        plt <- plotly::add_trace(plt,
          x = 0, y = 0, type = "scatter", mode = "markers",
          visible = "legendonly",
          marker  = list(symbol = sym, size = 10,
                         color  = "rgba(70,70,70,0.9)",
                         line   = list(width = 0.5, color = "black")),
          name             = cv,
          legendgroup      = "Contrast",
          legendgrouptitle = list(text = "Contrast"),
          showlegend       = TRUE,
          hoverinfo        = "none",
          inherit          = FALSE
        )
      }

      # ── Pass 3: black outlines on coloured data traces ─────────────────────────
      if (n_real > 0L)
        plt <- plotly::style(plt,
          marker.line.width = 1.4,
          marker.line.color = "black",
          traces = seq_len(n_real)
        )

      plt |>
        plotly::layout(
          title  = list(text = ""),
          legend = list(tracegroupgap = 14, font = list(size = 11),
                        itemsizing = "constant")
        )
    }

    output$score_all_volcano_ui <- shiny::renderUI({
      shiny::req(!is.null(score_all_result()), !is.null(score_all_result()$plot_res$volcano))
      bslib::card(
        full_screen = TRUE,
        bslib::card_header(
          class = "d-flex align-items-center justify-content-between",
          "Effect Size Volcano (hover to identify gene sets)",
          shiny::downloadButton(ns("dl_all_volcano"), "",
            icon  = shiny::icon("download"),
            class = "btn-sm btn-outline-secondary py-0 px-2",
            style = "font-size:0.75em;")
        ),
        shiny::div(
          class = "text-muted", style = "font-size:0.78em; padding:4px 6px 8px;",
          shiny::icon("circle-info"),
          " X-axis: |Cohen's d| (effect size). Y-axis: −log₁₀(adjusted p-value).",
          " Shape = scoring method. Colour = gene set.",
          " Hover any point to see gene set, contrast, method and exact values.",
          " Dashed lines mark the significance / effect size thresholds above."
        ),
        plotly::plotlyOutput(ns("score_all_volcano"), height = "480px")
      )
    })

    output$score_all_volcano <- plotly::renderPlotly({
      res <- score_all_result(); shiny::req(!is.null(res), !is.null(res$plot_res$volcano))
      p_vol <- res$plot_res$volcano
      tryCatch({
        plt <- plotly::ggplotly(p_vol,
                                tooltip = c("x", "y", "colour", "shape"))
        .clean_volcano_plotly(plt)
      },
      error = function(e) {
        shiny::showNotification(paste("Volcano error:", e$message),
                                type = "warning", duration = 8)
        plotly::plotly_empty()
      })
    })

    output$dl_all_volcano <- shiny::downloadHandler(
      filename = function() paste0("effect_volcano_", Sys.Date(), ".png"),
      content  = function(file) {
        res <- score_all_result(); shiny::req(!is.null(res), !is.null(res$plot_res$volcano))
        ggplot2::ggsave(file, plot = res$plot_res$volcano,
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
            .bm_card_header("logFC vs −log₁₀(padj)", "dl_gsea_deg_volcano", ns),
            shiny::div(
              class = "text-muted", style = "font-size:0.78em; padding:4px 6px 8px;",
              shiny::icon("circle-info"),
              " Volcano plot of differential expression for each contrast.",
              " Points are coloured by gene set membership when gene sets are provided.",
              " Thresholds mark significance and fold-change cutoffs."
            ),
            # Inline colour key — only relevant when gene set highlighting is on
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
                shiny::numericInput(ns("vol_threshold_y"), "padj threshold:",
                                    value = 0.05, min = 0, max = 1, step = 0.01),
                shiny::numericInput(ns("vol_top_n"), "Annotate top/bottom N genes by logFC (0 = none):",
                                    value = 0, min = 0, max = 50, step = 1),
                shiny::checkboxInput(ns("vol_topn_gs_only"),
                                     "Restrict N labels to gene set members only",
                                     value = FALSE),
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

        # 1. Combined Volcano (interactive)
        bslib::nav_panel("Combined Volcano",
          bslib::card(
            full_screen = TRUE,
            .bm_card_header("NES vs −log₁₀(padj)", "dl_gsea_combined", ns),
            shiny::div(
              class = "text-muted", style = "font-size:0.78em; padding:4px 6px 8px;",
              shiny::icon("circle-info"),
              " Each point is one gene set x contrast. Colour = gene set; shape = contrast. Hover for details.",
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
              " y-axis: −log₁₀(padj); dashed line marks the significance threshold."
            ),
            plotly::plotlyOutput(ns("gsea_combined_plot"), height = "520px")
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
              " Reported NES and padj values are from fgsea.",
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

    output$gsea_deg_volcano_plot <- shiny::renderPlot({
      res <- gsea_result(); shiny::req(!is.null(res), !is.null(res$DEGs))
      thr_x    <- input$vol_threshold_x    %||% 1
      thr_y    <- input$vol_threshold_y    %||% 0.05  # raw padj threshold
      n_top    <- as.integer(input$vol_top_n   %||% 0L)
      pts      <- input$vol_pointsize      %||% 2
      lsz      <- input$vol_labsize        %||% 8L
      hl_mode  <- input$vol_highlight_mode %||% "gene_sets"
      gs_only  <- isTRUE(input$vol_topn_gs_only)

      # Map highlight mode to plotVolcano arguments
      # "gene_sets": color gene set members; no threshold-based coloring
      # "threshold": manual ggplot highlighting genes with padj <= thr_y and |logFC| >= thr_x
      use_gs    <- (hl_mode == "gene_sets")
      genes_arg <- if (use_gs && !is.null(res$gene_sets)) res$gene_sets else NULL
      # thr_y is already a raw p-value (user enters 0.05, not 1.3)
      thr_y_raw <- if (!use_gs && thr_y > 0 && thr_y < 1) thr_y else NULL
      thr_x_arg <- if (!use_gs && thr_x > 0) thr_x else NULL

      contrast_names <- names(res$DEGs)
      n_cont <- length(contrast_names)
      n_gs   <- if (!is.null(genes_arg)) length(genes_arg) else 1L
      n_rows_vol <- n_gs; n_cols_vol <- n_cont
      gsea_vol_h(as.integer(min(10000L, max(550L, n_rows_vol * 550L + 80L))))
      gsea_vol_w(if (n_cols_vol >= 3L) "100%"
                 else paste0(min(100L, n_cols_vol * 50L), "%"))

      tryCatch({
        if (!use_gs) {
          # --- Threshold mode: manual ggplot2 (plotVolcano has issues in this mode) ---
          per_plots <- lapply(contrast_names, function(cn) {
            df <- as.data.frame(res$DEGs[[cn]])
            padj_col <- if ("adj.P.Val" %in% colnames(df)) "adj.P.Val" else
                        if ("padj"      %in% colnames(df)) "padj"      else NULL
            lfc_col  <- if ("logFC"     %in% colnames(df)) "logFC"     else
                        if ("log2FoldChange" %in% colnames(df)) "log2FoldChange" else NULL
            shiny::validate(shiny::need(!is.null(padj_col) && !is.null(lfc_col),
                                        "Cannot find logFC / adj.P.Val columns."))
            df$lfc_   <- df[[lfc_col]]
            df$padj_  <- df[[padj_col]]
            df$logp_  <- -log10(pmax(df$padj_, 1e-300))
            sig <- rep(TRUE, nrow(df))
            if (!is.null(thr_y_raw)) sig <- sig & (df$padj_ <= thr_y_raw)
            if (!is.null(thr_x_arg)) sig <- sig & (abs(df$lfc_) >= thr_x_arg)
            if (is.null(thr_y_raw) && is.null(thr_x_arg)) sig <- rep(FALSE, nrow(df))
            df$.sig <- sig
            gg <- ggplot2::ggplot(df, ggplot2::aes(x = lfc_, y = logp_)) +
              ggplot2::geom_point(data = df[!df$.sig, , drop = FALSE],
                                  color = "#B7B7B7", size = pts, alpha = 0.4) +
              ggplot2::geom_point(data = df[df$.sig,  , drop = FALSE],
                                  color = "#6489B4", size = pts, alpha = 0.7) +
              ggplot2::labs(title = cn, x = lfc_col,
                            y = paste0("-log10(", padj_col, ")")) +
              ggplot2::theme_bw(base_size = lsz)
            if (!is.null(thr_y_raw))
              gg <- gg + ggplot2::geom_hline(yintercept = -log10(thr_y_raw),
                                             linetype = "dashed", color = "#555555")
            if (!is.null(thr_x_arg))
              gg <- gg + ggplot2::geom_vline(xintercept = c(-thr_x_arg, thr_x_arg),
                                             linetype = "dashed", color = "#555555")
            gg
          })
        } else if (gs_only && n_top > 0L && !is.null(genes_arg)) {
          # --- Gene set members mode + gs_only N annotation ---
          # plotVolcano returns a ggarrange so we cannot add layers to it.
          # Build each gene-set x contrast plot manually so ggrepel works.
          per_plots <- lapply(contrast_names, function(cn) {
            df <- as.data.frame(res$DEGs[[cn]])
            lfc_col  <- if ("logFC"     %in% colnames(df)) "logFC"
                        else if ("log2FoldChange" %in% colnames(df)) "log2FoldChange"
                        else "logFC"
            padj_col <- if ("adj.P.Val" %in% colnames(df)) "adj.P.Val"
                        else if ("padj" %in% colnames(df)) "padj"
                        else "adj.P.Val"
            df$gene_ <- rownames(df)
            df$lfc_  <- df[[lfc_col]]
            df$logp_ <- -log10(pmax(df[[padj_col]], 1e-300))

            gs_plots <- lapply(names(genes_arg), function(gsn) {
              gs_genes <- .gs_gene_names(genes_arg[[gsn]])
              in_gs <- df$gene_ %in% gs_genes
              gg <- ggplot2::ggplot(df, ggplot2::aes(x = .data$lfc_, y = .data$logp_)) +
                ggplot2::geom_point(data = df[!in_gs, , drop = FALSE],
                                    color = "#B7B7B7", size = pts, alpha = 0.4) +
                ggplot2::geom_point(data = df[in_gs,  , drop = FALSE],
                                    color = "#05254A", size = pts, alpha = 0.8) +
                ggplot2::labs(title    = gsn,
                              subtitle = cn,
                              x = lfc_col,
                              y = paste0("-log10(", padj_col, ")")) +
                ggplot2::theme_bw(base_size = lsz)
              df_gs <- df[in_gs, , drop = FALSE]
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
              gg
            })
            if (length(gs_plots) == 1L) gs_plots[[1]]
            else ggpubr::ggarrange(plotlist = gs_plots,
                                   nrow = length(gs_plots), ncol = 1L, align = "h")
          })
        } else {
          # --- Gene set members mode: per-contrast plotVolcano calls ---
          .vol_args <- list(
            genes                  = genes_arg,
            N                      = if (n_top > 0L) n_top else NULL,
            x                      = "logFC",
            y                      = "-log10(adj.P.Val)",
            pointSize              = pts,
            color                  = "#6489B4",
            highlightcolor         = "#05254A",
            highlightcolor_upreg   = "#038C65",
            highlightcolor_downreg = "#8C0303",
            nointerestcolor        = "#B7B7B7",
            threshold_y            = NULL,
            threshold_x            = NULL,
            labsize                = lsz,
            widthlabs              = 32L,
            invert                 = FALSE
          )
          per_plots <- lapply(contrast_names, function(cn)
            do.call(plotVolcano, c(list(DEResultsList = res$DEGs[cn]), .vol_args)))
        }

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
        gs_only <- isTRUE(input$vol_topn_gs_only)
        use_gs    <- (hl_mode == "gene_sets")
        genes_arg <- if (use_gs && !is.null(res$gene_sets)) res$gene_sets else NULL
        thr_y_raw <- if (!use_gs && thr_y > 0 && thr_y < 1) thr_y else NULL
        thr_x_arg <- if (!use_gs && thr_x > 0) thr_x else NULL
        contrast_names <- names(res$DEGs)
        if (!use_gs) {
          per_plots <- lapply(contrast_names, function(cn) {
            df <- as.data.frame(res$DEGs[[cn]])
            padj_col <- if ("adj.P.Val" %in% colnames(df)) "adj.P.Val" else "padj"
            lfc_col  <- if ("logFC"     %in% colnames(df)) "logFC"     else "log2FoldChange"
            df$lfc_  <- df[[lfc_col]]; df$padj_ <- df[[padj_col]]
            df$logp_ <- -log10(pmax(df$padj_, 1e-300))
            sig <- if (!is.null(thr_y_raw)) df$padj_ <= thr_y_raw else rep(TRUE, nrow(df))
            if (!is.null(thr_x_arg)) sig <- sig & (abs(df$lfc_) >= thr_x_arg)
            df$.sig <- sig
            gg <- ggplot2::ggplot(df, ggplot2::aes(x = lfc_, y = logp_)) +
              ggplot2::geom_point(data = df[!df$.sig, , drop = FALSE],
                                  color = "#B7B7B7", size = pts, alpha = 0.4) +
              ggplot2::geom_point(data = df[df$.sig, , drop = FALSE],
                                  color = "#6489B4", size = pts, alpha = 0.7) +
              ggplot2::labs(title = cn, x = lfc_col, y = paste0("-log10(", padj_col, ")")) +
              ggplot2::theme_bw(base_size = lsz)
            if (!is.null(thr_y_raw))
              gg <- gg + ggplot2::geom_hline(yintercept = -log10(thr_y_raw),
                                             linetype = "dashed", color = "#555555")
            if (!is.null(thr_x_arg))
              gg <- gg + ggplot2::geom_vline(xintercept = c(-thr_x_arg, thr_x_arg),
                                             linetype = "dashed", color = "#555555")
            gg
          })
        } else if (gs_only && n_top > 0L && !is.null(genes_arg)) {
          # gs_only: build manually per gene set so ggrepel works on ggplots
          per_plots <- lapply(contrast_names, function(cn) {
            df <- as.data.frame(res$DEGs[[cn]])
            lfc_col  <- if ("logFC"     %in% colnames(df)) "logFC"     else "log2FoldChange"
            padj_col <- if ("adj.P.Val" %in% colnames(df)) "adj.P.Val" else "padj"
            df$gene_ <- rownames(df)
            df$lfc_  <- df[[lfc_col]]
            df$logp_ <- -log10(pmax(df[[padj_col]], 1e-300))
            gs_plots <- lapply(names(genes_arg), function(gsn) {
              gs_genes <- .gs_gene_names(genes_arg[[gsn]])
              in_gs <- df$gene_ %in% gs_genes
              gg <- ggplot2::ggplot(df, ggplot2::aes(x = .data$lfc_, y = .data$logp_)) +
                ggplot2::geom_point(data = df[!in_gs, , drop = FALSE],
                                    color = "#B7B7B7", size = pts, alpha = 0.4) +
                ggplot2::geom_point(data = df[in_gs,  , drop = FALSE],
                                    color = "#05254A", size = pts, alpha = 0.8) +
                ggplot2::labs(title = gsn, subtitle = cn, x = lfc_col,
                              y = paste0("-log10(", padj_col, ")")) +
                ggplot2::theme_bw(base_size = lsz)
              df_gs <- df[in_gs, , drop = FALSE]
              if (nrow(df_gs) > 0L) {
                df_gs <- df_gs[order(abs(df_gs$lfc_), decreasing = TRUE), ]
                df_top <- head(df_gs, n_top)
                gg <- gg + ggrepel::geom_text_repel(
                  data = df_top,
                  ggplot2::aes(x = .data$lfc_, y = .data$logp_, label = .data$gene_),
                  size = lsz / 3, max.overlaps = n_top * 2L, force = 10,
                  color = "#05254A", fontface = "bold")
              }
              gg
            })
            if (length(gs_plots) == 1L) gs_plots[[1]]
            else ggpubr::ggarrange(plotlist = gs_plots,
                                   nrow = length(gs_plots), ncol = 1L, align = "h")
          })
        } else {
          .vol_args <- list(
            genes = genes_arg, N = if (n_top > 0L) n_top else NULL,
            x = "logFC", y = "-log10(adj.P.Val)", pointSize = pts,
            color = "#6489B4", highlightcolor = "#05254A",
            highlightcolor_upreg = "#038C65", highlightcolor_downreg = "#8C0303",
            nointerestcolor = "#B7B7B7",
            threshold_y = NULL, threshold_x = NULL,
            labsize = lsz, widthlabs = 32L, invert = FALSE
          )
          per_plots <- lapply(contrast_names, function(cn)
            do.call(plotVolcano, c(list(DEResultsList = res$DEGs[cn]), .vol_args)))
        }
        combined_vol <- if (length(per_plots) == 1L) per_plots[[1]] else
          ggpubr::ggarrange(plotlist = per_plots,
                            ncol = length(per_plots), nrow = 1L, align = "h")
        h_in <- gsea_vol_h() / 150
        ggplot2::ggsave(file, plot = combined_vol,
                        width = 10, height = max(4, h_in), dpi = 150, units = "in")
      }
    )

    # ---- GSEA: Enrichment Plots ---------------------------------------------

    output$gsea_enrich_plot <- shiny::renderPlot({
      res <- gsea_result(); shiny::req(!is.null(res))
      font_s  <- input$gsea_title_size  %||% 10L
      w_title <- input$gsea_width_title %||% 32L
      n_plots <- length(res$gsea_res) * length(res$gene_sets)
      enrich_h(as.integer(min(6000L, max(400L, n_plots * 220L))))
      tryCatch(
        print(plotGSEAenrichment(
          GSEA_results = res$gsea_res, DEGList = res$DEGs,
          gene_sets = res$gene_sets, titlesize = font_s, widthTitle = w_title,
          grid = TRUE,
          nrow = length(res$gene_sets), ncol = length(res$gsea_res))),
        error = function(e)
          shiny::showNotification(paste("Enrichment plot error:", e$message),
                                  type = "warning", duration = 8))
    }, height = function() enrich_h(), res = 150)

    # ---- GSEA: NES Lollipop -------------------------------------------------

    output$nes_lollipop <- shiny::renderPlot({
      res <- gsea_result(); shiny::req(!is.null(res))
      sig_thr  <- input$gsea_sig_threshold %||% 0.05
      font_s   <- input$gsea_title_size    %||% 10L
      w_labels <- input$gsea_width_legend  %||% 32L
      tryCatch(
        print(plotNESlollipop(
          GSEA_results = res$gsea_res, sig_threshold = sig_thr,
          titlesize = font_s, widthlabels = w_labels, grid = TRUE)),
        error = function(e)
          shiny::showNotification(paste("Lollipop error:", e$message),
                                  type = "warning", duration = 8))
    }, height = function() nes_h(), res = 150)

    # ---- GSEA: Combined Volcano (interactive) --------------------------------

    # Build a clean single-layer ggplot for the combined GSEA volcano.
    # Avoids plotCombinedGSEA's two geom_point layers (which ggplotly renders
    # as two overlapping traces of different sizes, making fake outlines).
    .build_combined_gsea_gg <- function(gsea_res, sig_thr, pt_size, w_legend) {
      combined_data <- do.call(rbind, lapply(names(gsea_res), function(cn) {
        df <- gsea_res[[cn]]; df$contrast <- cn; df
      }))
      combined_data$logpadj <- -log10(combined_data$padj)
      combined_data$pathway_w <- vapply(combined_data$pathway,
        function(x) wrap_title(x, w_legend), character(1))
      combined_data$tip <- paste0(
        "Pathway: ",  combined_data$pathway,
        "<br>NES: ",  round(combined_data$NES, 4),
        "<br>padj: ", signif(combined_data$padj, 3),
        "<br>Contrast: ", combined_data$contrast)

      n_paths <- length(unique(combined_data$pathway_w))
      pal <- grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(min(12L, max(3L, n_paths)), "Set3")
      )(n_paths)

      n_cont  <- length(unique(combined_data$contrast))
      shapes  <- 15:(15 - 1L + n_cont)

      ggplot2::ggplot(combined_data,
          ggplot2::aes(x = NES, y = logpadj,
                       colour = factor(pathway_w),
                       shape  = contrast,
                       text   = tip)) +
        ggplot2::geom_point(size = pt_size) +
        ggplot2::geom_hline(yintercept = -log10(sig_thr),
                            linetype = "dashed", color = "black", linewidth = 0.5) +
        ggplot2::geom_vline(xintercept = 0,
                            linetype = "dashed", color = "black", linewidth = 0.5) +
        # name = NULL: plotly would render scale titles as phantom legend items,
        # duplicating the legendgrouptitle entries added by .clean_gsea_volcano_plotly
        ggplot2::scale_colour_manual(values = pal, name = NULL) +
        ggplot2::scale_shape_manual(values = shapes, name = NULL) +
        ggplot2::labs(x = "Normalized Enrichment Score (NES)",
                      y = "-log10(adj. p-value)") +
        ggplot2::theme_bw() +
        ggplot2::facet_grid(. ~ stat_used,
          labeller = ggplot2::labeller(
            stat_used = c("t" = "Enriched/Depleted", "B" = "Altered")),
          scales = "free", switch = "y") +
        ggplot2::theme(strip.background = ggplot2::element_rect(fill = "white"),
                       legend.position  = "right")
    }

    output$gsea_combined_plot <- plotly::renderPlotly({
      res <- gsea_result(); shiny::req(!is.null(res), !is.null(res$gsea_res))
      sig_thr  <- input$gsea_sig_threshold %||% 0.05
      pt_size  <- input$gsea_point_size    %||% 7L
      w_legend <- input$gsea_width_legend  %||% 32L
      tryCatch({
        gg <- .build_combined_gsea_gg(res$gsea_res, sig_thr, pt_size, w_legend)
        known_pathways  <- unique(do.call(c, lapply(res$gsea_res, `[[`, "pathway")))
        known_contrasts <- names(res$gsea_res)
        plotly::ggplotly(gg, tooltip = "text") |>
          .clean_gsea_volcano_plotly(known_pathways, known_contrasts)
      }, error = function(e)
        shiny::showNotification(paste("Combined GSEA error:", e$message),
                                type = "warning", duration = 8))
    })

    # ---- GSEA: Download handlers --------------------------------------------

    output$dl_gsea_enrich <- shiny::downloadHandler(
      filename = function() paste0("gsea_enrichment_", Sys.Date(), ".png"),
      content  = function(file) {
        res <- gsea_result(); shiny::req(!is.null(res))
        font_s  <- input$gsea_title_size  %||% 10L
        w_title <- input$gsea_width_title %||% 32L
        h_in    <- enrich_h() / 96
        ggplot2::ggsave(file,
          plot = plotGSEAenrichment(
            GSEA_results = res$gsea_res, DEGList = res$DEGs,
            gene_sets = res$gene_sets, titlesize = font_s, widthTitle = w_title,
            grid = TRUE,
            nrow = length(res$gene_sets), ncol = length(res$gsea_res)),
          width = 14, height = max(5, h_in), dpi = 150, units = "in")
      }
    )

    output$dl_nes_lollipop <- shiny::downloadHandler(
      filename = function() paste0("nes_lollipop_", Sys.Date(), ".png"),
      content  = function(file) {
        res <- gsea_result(); shiny::req(!is.null(res))
        sig_thr  <- input$gsea_sig_threshold %||% 0.05
        font_s   <- input$gsea_title_size    %||% 10L
        w_labels <- input$gsea_width_legend  %||% 32L
        h_in <- nes_h() / 96
        ggplot2::ggsave(file,
          plot = plotNESlollipop(res$gsea_res, sig_threshold = sig_thr,
                                  titlesize = font_s, widthlabels = w_labels, grid = TRUE),
          width = 12, height = max(4, h_in), dpi = 150, units = "in")
      }
    )

    output$dl_gsea_combined <- shiny::downloadHandler(
      filename = function() paste0("gsea_combined_", Sys.Date(), ".png"),
      content  = function(file) {
        res <- gsea_result(); shiny::req(!is.null(res), !is.null(res$gsea_res))
        sig_thr  <- input$gsea_sig_threshold %||% 0.05
        pt_size  <- input$gsea_point_size    %||% 7L
        w_legend <- input$gsea_width_legend  %||% 32L
        ggplot2::ggsave(file,
          plot = .build_combined_gsea_gg(res$gsea_res, sig_thr, pt_size, w_legend),
          width = 12, height = 6, dpi = 150, units = "in")
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
        shiny::tags$span(class = "badge bg-secondary me-1", paste0("Padj: ", p$padj_m)),
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

    # Helper: rebuild FPR violin/jitter plot for ONE contrast from cached data.
    # data_list: named list of per-signature data frames (from FPR_Simulation$data).
    # Each df has: signature, contrast, method, cohen, padj, type, FPR columns.
    # Rebuild the FPR plot manually from cached data_list (no re-run needed).
    # Mirrors the original FPR_Simulation plot structure exactly:
    #   - One ggplot per gene set, titled with the gene set name
    #   - facet_wrap(. ~ contrast, ncol=1) — contrasts stacked vertically
    #   - x = method (3 violins: ssGSEA / logmedian / ranking)
    #   - Violin = null distribution (simulated); green dot = original gene set value
    #   - Red dashed line = 95th percentile of simulated distribution
    #   - FPR label above each violin
    # Gene set plots are combined in a grid with ggarrange.
    .build_fpr_all_plots <- function(data_list, fontsize = 10L, cohentype = "d",
                                      ncol_g    = NULL,
                                      color_vals = c(Original  = "#68B393",
                                                     Simulated = "#666666")) {
      method_lvls <- c("ssGSEA", "logmedian", "ranking")

      plot_list <- lapply(names(data_list), function(sig) {
        df <- data_list[[sig]]
        if (is.null(df) || nrow(df) == 0L) return(NULL)

        ml <- method_lvls[method_lvls %in% unique(df$method)]
        df$method <- factor(df$method, levels = ml)
        contrasts  <- unique(df$contrast)

        # 95th-percentile dashed lines per (contrast × method)
        q_data <- do.call(rbind, lapply(contrasts, function(ct) {
          do.call(rbind, lapply(seq_along(ml), function(mi) {
            sub <- df[df$contrast == ct & df$method == ml[mi], ]
            if (nrow(sub) == 0L) return(NULL)
            data.frame(contrast = ct,
                       method   = factor(ml[mi], levels = ml),
                       q_high   = stats::quantile(sub$cohen, 0.95, na.rm = TRUE),
                       xmin = mi - 0.35, xmax = mi + 0.35,
                       stringsAsFactors = FALSE)
          }))
        }))

        # Max cohen + FPR label per (contrast × method)
        orig_df <- df[df$type == "Original", ]
        all_max <- stats::aggregate(cohen ~ method + contrast, data = df, FUN = max)
        all_max$FPR <- NA_real_
        for (i in seq_len(nrow(all_max))) {
          idx <- which(orig_df$method   == all_max$method[i] &
                         orig_df$contrast == all_max$contrast[i])
          if (length(idx) > 0L) all_max$FPR[i] <- orig_df$FPR[idx[1L]]
        }
        all_max$label   <- sprintf("FPR=%.2f", all_max$FPR)
        all_max$y_label <- all_max$cohen + 0.3
        all_max$method  <- factor(all_max$method, levels = ml)

        ggplot2::ggplot() +
          ggplot2::geom_jitter(
            data = df[df$type == "Simulated", ],
            ggplot2::aes(y = cohen, x = method, color = type),
            width = 0.3, height = 0, size = 2, alpha = 0.5) +
          ggplot2::geom_violin(
            data = df,
            ggplot2::aes(y = cohen, x = method),
            fill = "#F0F0F0", color = "black", alpha = 0.5) +
          ggplot2::geom_jitter(
            data = df[df$type == "Original", ],
            ggplot2::aes(y = cohen, x = method, color = type),
            width = 0, height = 0, size = 2, alpha = 1) +
          ggplot2::geom_text(
            data = all_max,
            ggplot2::aes(x = method, y = y_label, label = label),
            size = 3, inherit.aes = FALSE) +
          {if (!is.null(q_data) && nrow(q_data) > 0L)
             ggplot2::geom_segment(
               data = q_data,
               ggplot2::aes(x = xmin, xend = xmax, y = q_high, yend = q_high),
               linetype = "dashed", color = "red", inherit.aes = FALSE)
           else NULL} +
          ggplot2::labs(
            title = wrap_title(sig, 32L),
            y     = if (cohentype == "d") "|Cohen's d|" else "|Cohen's f|",
            x     = "Method", color = "") +
          ggplot2::theme_classic(base_size = fontsize) +
          ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, size = fontsize),
            strip.text = ggplot2::element_text(size = max(7L, fontsize - 1L))) +
          ggplot2::facet_wrap(. ~ contrast, scales = "free", ncol = 1L,
                              strip.position = "left") +
          ggplot2::scale_color_manual(values = color_vals)
      })

      plot_list <- Filter(Negate(is.null), plot_list)
      n <- length(plot_list)
      if (n == 0L) return(NULL)
      if (is.null(ncol_g)) ncol_g <- min(4L, ceiling(sqrt(n)))
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
      ncol_g    <- 3L
      # Reactively filter cached FPR data by selected contrasts (no re-run needed).
      # sel_conts is a character vector of contrast strings (e.g. "A-B", "A-C").
      data_list  <- res$fpr_res$data
      sel_conts  <- input$score_selected_groups
      if (!is.null(sel_conts) && length(sel_conts) > 0L) {
        data_list <- lapply(data_list, function(df)
          df[df$contrast %in% sel_conts, , drop = FALSE])
        data_list <- Filter(function(df) nrow(df) > 0L, data_list)
      }
      shiny::req(length(data_list) > 0L)
      tryCatch({
        p <- .build_fpr_all_plots(data_list, fontsize, cohentype, ncol_g)
        if (!is.null(p)) print(p)
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
        ncol_g    <- 3L
        p    <- .build_fpr_all_plots(res$fpr_res$data, fontsize, cohentype, ncol_g)
        h_in <- fpr_h() / 150
        ggplot2::ggsave(file, plot = p,
                        width = 12, height = max(4, h_in),
                        dpi = 150, units = "in")
      }
    )

  })
}
