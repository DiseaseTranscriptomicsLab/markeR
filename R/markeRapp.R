#' @importFrom bslib bs_theme nav_panel page_navbar sidebar layout_sidebar accordion accordion_panel navset_card_tab card card_header nav_spacer nav_item
#' @import shiny
#' @import GEOquery
#' @importFrom DT DTOutput renderDT datatable
#' @importFrom htmltools findDependencies
#' @importFrom readxl read_xls read_xlsx
#' @importFrom R.utils gunzip
#' @importFrom data.table fread
#' @importFrom Biobase exprs pData annotation experimentData
#' @importFrom shinyWidgets pickerInput pickerOptions
#' @importFrom utils download.file read.csv read.delim read.table untar unzip packageVersion
#' @importFrom tools file_ext file_path_sans_ext
#' @importFrom ComplexHeatmap draw

# Session log preamble (R version + package versions the Shiny app runs
# on): this information is identical for every session and known the
# moment the app process starts - it doesn't depend on anything
# per-session. Computing it once here, at package/app-build time, means
# it can be baked directly into the STATIC UI (see the global_log_entries
# placeholder below), so it's genuinely present in the page's HTML from
# the very first byte, with no dependency on Shiny's reactive flush timing
# at all. dataServer() (R/markeR_tab_data.R) calls this same function again
# per session, purely so its "Session started" timestamp reflects that
# session's own start time; the package/version lines themselves are
# always identical.
#
# NOTE: reads DESCRIPTION straight off disk with read.dcf() rather than
# via utils::packageDescription(): the latter's single-field return value
# is a bare character string rather than a named/list-like object (so
# `desc$Imports`-style access silently errors), and its lookup can come up
# empty under devtools::load_all()-style dev workflows where the package
# isn't "installed" in the usual sense. find.package()/system.file() are
# both load_all()-aware, so this path resolves correctly either way.
.markeR_session_preamble <- function() {
  tryCatch({
    .parse_dep_field <- function(field) {
      if (is.null(field) || length(field) == 0L || is.na(field)) return(character(0))
      pkgs <- strsplit(field, ",")[[1]]
      pkgs <- trimws(sub("\\s*\\(.*?\\)", "", pkgs))  # strip version specs like (>= 4.0)
      pkgs[nzchar(pkgs) & pkgs != "R"]
    }
    desc_path <- system.file("DESCRIPTION", package = "markeR")
    if (!nzchar(desc_path)) {
      desc_path <- file.path(find.package("markeR", quiet = TRUE), "DESCRIPTION")
    }
    dcf <- read.dcf(desc_path)
    imports_field <- if ("Imports" %in% colnames(dcf)) dcf[1, "Imports"] else NA_character_
    # These Suggests are genuinely used by the app at runtime (loaded
    # conditionally via requireNamespace()): biomaRt/AnnotationDbi/
    # org.Hs.eg.db/org.Mm.eg.db for the Preprocessing tab's gene-ID
    # conversion, stringdist/tximport for GEO import, BiocManager for
    # both - so a "not installed" against one of them is meaningful, not
    # noise. Everything else in Suggests (devtools, testthat, roxygen2,
    # knitr, rmarkdown, covr, mockery, magick, BiocStyle, sva, ...) is
    # dev/test/doc tooling the app never touches, so those stay excluded.
    app_suggests <- c("BiocManager", "biomaRt", "AnnotationDbi",
                       "org.Hs.eg.db", "org.Mm.eg.db", "stringdist", "tximport")
    pkgs <- sort(unique(c(
      "markeR",
      .parse_dep_field(imports_field),
      app_suggests
    )))
    pkg_lines <- vapply(pkgs, function(p) {
      v <- tryCatch(as.character(utils::packageVersion(p)), error = function(e) "not installed")
      paste0("[Session]   ", p, ": ", v)   # note: ALL lines prefixed with [Session]
    }, character(1L))
    c(
      paste0("[Session] Session started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      paste0("[Session] R version: ", R.version$major, ".", R.version$minor,
             " (", R.version$`svn rev`, ") | Platform: ", R.version$platform),
      "[Session] Package versions:",
      pkg_lines
    )
  }, error = function(e) {
    # Never let this come up completely empty - a visible error line is far
    # easier to diagnose than a session log that silently shows nothing at
    # all, which is what a swallowed error here used to look like.
    c(
      paste0("[Session] Session started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      paste0("[Session] Package info unavailable: ", conditionMessage(e))
    )
  })
}

# Computed once, at package/app-build time (NOT per-session) - shared by
# every session, and used below to bake real session/package info directly
# into the static UI. See .markeR_session_preamble()'s comment above.
.markeR_initial_session_lines <- .markeR_session_preamble()

ui <- bslib::page_navbar(

  # fillable = FALSE at the root: page_navbar() defaults to a flexbox "fill"
  # layout that constrains all tab content to the viewport height, clipping
  # anything taller instead of letting the page scroll. Individual tabs'
  # layout_sidebar(fillable = FALSE) calls only override this within their
  # own card, not the root container establishing that height constraint in
  # the first place - setting it here too removes the constraint at its
  # source for every tab (About included), which is more robust across
  # different rendering contexts (e.g. plain browser vs. an iframe-embedded
  # deployment) than relying on a single nested override.
  fillable = FALSE,

  # ---- Main NavBar & Title ----
  title = shiny::tags$span(
    #shiny::a("markeR "),
    shiny::a(         
      shiny::HTML(
        paste0("<span style='color:#306F1D;'>m</span>",
               "<span style='color:#020201;'>a</span>",
               "<span style='color:#801F4F;'>r</span>",
               "<span style='color:#020201;'>k</span>",
               "<span style='color:#EBB43E;'>e</span>",
               "<span style='color:#1E497B;'>R </span>" 
        )    ),
      style = "font-size: 1.5em;"),
    
    # shiny::a(
    #   paste0("v", utils::packageVersion("markeR")),
    #   href = paste0(
    #     "https://github.com/DiseaseTranscriptomicsLab/markeR/releases/tag/v",
    #     utils::packageVersion("markeR")
    #   ),
    #   target = "_blank",
    #   style = "font-size: 0.85em; color: #949494; text-decoration: none;"
    # )
    shiny::a("Web App 0.99.10",
             style = "font-size: 0.7em; color: #949494; text-decoration: none;")
  ),
  
  theme = bslib::bs_theme(bootswatch = "yeti",
                          primary = "rgb(186, 89, 0)"),
  
  
  # ---- Tabs ----
  
  ##### ABOUT #####
  
  bslib::nav_panel(
    "About",
    # fillable = FALSE: otherwise this tab participates in bslib's flexbox
    # "fill" chain (inherited from page_navbar()), which constrains content
    # height to the viewport and clips overflow instead of letting the page
    # scroll normally - the same class of bug already worked around in the
    # other tabs' layout_sidebar() calls.
    fillable = FALSE,
    shiny::tags$div(
      style = "
    display: flex;
    flex-wrap: wrap;          /* allows drop on very narrow screens */
    align-items: center;      /* vertically center */
    justify-content: center;
    max-width: 1000px;        /* increase to fit larger logo */
    margin: 0 auto;
  ",
      
      # Logo on the left
      shiny::tags$div(
        style = "flex: 0 0 200px; margin-right: 20px;",  # fixed 200px width
        shiny::img(
          src = "logo/logo.png",
          style = "width: 100%; height: auto;"
        )
      ),
      
      # Text on the right
      shiny::tags$div(
        style = "flex: 1 1 600px;",
        #shiny::h3("Welcome to markeR!"),
        shiny::h2(
          shiny::HTML(
            paste0("Welcome to ",
                   "<span style='color:#306F1D;'>m</span>",
                   "<span style='color:#020201;'>a</span>",
                   "<span style='color:#801F4F;'>r</span>",
                   "<span style='color:#020201;'>k</span>",
                   "<span style='color:#EBB43E;'>e</span>",
                   "<span style='color:#1E497B;'>R</span>",
                   "!"
            )
          )
        ),
        
        shiny::p(
          
          shiny::HTML(
            paste0("You may have a collection of genes that putatively mark a phenotype but be unsure how to combine them into ",
                   "<span style='color:#306F1D;'> meaningful metrics for phenotype quantification</span>",
                   ". markeR helps you evaluate gene sets as phenotype markers and explore their information. Below is a quick guide to each
                   section of the app:"
            )
          ),
          
          style = "font-size: 1.2em;"
        )
      )
    ),

    # ── Version + citation banners (side by side) ───────────────────────────
    # The version quoted here (v1.3.1) should be kept in sync with the
    # "built with markeR v1.3.1" text below and with DESCRIPTION's Version
    # field - it names the markeR release this app was developed against,
    # not necessarily whatever markeR version happens to be installed at
    # runtime. The citation banner cites the NAR Genomics and Bioinformatics
    # paper (same one linked as "our paper" below) rather than the
    # Bioconductor package DOI, since that paper is markeR's citable
    # reference.
    shiny::tags$div(
      style = paste(
        "display:flex; flex-wrap:wrap; gap:16px; align-items:stretch;",
        "max-width: 1000px; margin: 0 auto 20px;"
      ),

      # Version banner (narrower)
      shiny::tags$div(
        style = paste(
          "flex: 1 1 240px; max-width: 300px; padding: 14px 18px;",
          "background:#f8f9fa; border:1px solid #dee2e6; border-radius:8px;",
          "font-size: 0.9em; color:#495057; text-align:center;",
          "display:flex; flex-direction:column; align-items:center;",
          "justify-content:center; gap:8px;"
        ),
        shiny::icon("box-open", style = "color:#306F1D; font-size:2em;"),
        shiny::tags$div(
          "This app was built with the ",
          shiny::tags$b("markeR"), " Bioconductor package ",
          shiny::a(shiny::tags$b("v1.3.1"),
            href = "https://github.com/DiseaseTranscriptomicsLab/markeR/releases/tag/v1.3.1",
            target = "_blank", rel = "noopener noreferrer"),
          ". For more on the methodology behind it, see ",
          shiny::a("our paper",
            href = "https://doi.org/10.1093/nargab/lqag057",
            target = "_blank", rel = "noopener noreferrer"),
          "."
        )
      ),

      # Citation banner (wider)
      shiny::tags$div(
        style = paste(
          "flex: 2 1 400px; padding: 14px 18px;",
          "background:#f8f9fa; border:1px solid #dee2e6; border-radius:8px;",
          "font-size: 0.9em; color:#495057; text-align:center;",
          "display:flex; flex-direction:column; align-items:center;",
          "justify-content:center; gap:8px;"
        ),
        shiny::icon("quote-left", style = "color:#306F1D; font-size:2em;"),
        shiny::tags$div(
          shiny::tags$b("To cite this app, please cite:"),
          shiny::tags$br(),
          shiny::tags$span(style = "font-style:italic;",
            "Martins-Silva R, Kaizeler A, Barbosa-Morais NL (2026). Exploring molecular ",
            "signatures of senescence with markeR, an R toolkit for evaluating gene sets as ",
            "phenotypic markers. NAR Genomics and Bioinformatics, 8(2), lqag057. ",
            shiny::a("doi:10.1093/nargab/lqag057",
              href = "https://doi.org/10.1093/nargab/lqag057",
              target = "_blank", rel = "noopener noreferrer"),
            "."
          )
        )
      )
    ),

    shiny::tags$div(
      style = "
      display: flex;
      flex-wrap: wrap;
      align-items: flex-start;
      margin-left: 50px;
    ",
      
      # Left: text
      shiny::tags$div(
        style = "
        flex: 1 1 500px;    /* grow/shrink, minimum width 400px */
        padding-right: 15px;
      ",
        shiny::h5("1. Data", style = "color: #306F1D;"),         # blue
        shiny::p("Load your gene expression datasets and associated metadata,
           as well as the gene sets of interest that are hypothesised to mark
           a given phenotype. This is the first step to get your data into the app.
           Alternatively, import data directly from GEO or explore example data on senescence."),
        
        shiny::h5("2. Preprocessing", style = "color: #EBB43E;"), # orange
        shiny::p("Filter lowly expressed genes, normalise your data, and explore 
           it with Principal Component Analysis (PCA). This step is optional 
           if your input data is already preprocessed."),
        
        shiny::h5("3. Gene Sets", style = "color: #801F4F;"),      # green
        shiny::p("Gene sets can be compared using the Jaccard index and odds 
           ratio to assess similarity with reference collections like 
           MsigDB Hallmarks. Results are displayed as a heatmap, and 
           the expression of individual genes can also be visualised to 
           see how each contributes to phenotype quantification."),
        
        shiny::h5("4. Benchmarking Mode", style = "color: #020201;"), # red
        shiny::p(
          "Test the robustness of (one or more) gene sets in distinguishing phenotypic groups using two modules:
     Enrichment and Score. Enrichment uses the Gene Set Enrichment Analysis (GSEA) method based on differential gene expression,
     while Score offers three methods for sample-wise quantification: ranking, logmedian, and ssGSEA.
     More details are available in ",
          shiny::a("markeR's paper", href = "https://doi.org/10.1093/nargab/lqag057", target = "_blank"),
          " or ",
          shiny::a("vignette", href = "https://diseasetranscriptomicslab.github.io/markeR/", target = "_blank"),
          "."
        ),
        
        shiny::h5("5. Discovery Mode", style = "color: #1E497B;"), # purple
        shiny::p(
          "Explores all metrics (scores, enrichment, etc.)
     for every user-selected variable to find potential associations
     with the phenotype marked by a gene set. It is intended for
     hypothesis generation and assumes the gene set is robust in
     marking the phenotype. This differs from benchmarking mode,
     where multiple gene sets can be tested to identify the
     best-performing one."
        )
      ),
      
      # Right: image
      shiny::tags$div(
        style = "
        flex: 1 1 700px;   /* grow/shrink, minimum width 300px */
        display: flex;
        justify-content: center;
        align-items: flex-start;
      ",
        shiny::img(
          src = "methods/methods.png",
          style = "max-width: 85%; height: auto; margin-top: 10px;"
        )
      )
    ),

    shiny::hr(),
    shiny::uiOutput("data_log_ui_about"),

    shiny::div(
      # Extra bottom margin so this clears the app's fixed, always-on-top
      # log bar at the very bottom of the viewport (#global-log-bar in
      # markeRapp.R's <head> styles) - without it, the fixed bar sits over
      # this text and hides it regardless of how short the About page is.
      style = "text-align:right; color:#8a8f96; font-size:0.75em; margin:10px 6px 60px;",
      "markeR's Bioconductor package was converted into this Shiny web app",
      " with the help of Claude Sonnet 5 (Anthropic)."
    )

  ),

  ##### DATA #####
  bslib::nav_panel(
    "Data",
    dataUI()
  ),

  ##### PREPROCESSING #####
  bslib::nav_panel(
    "Preprocessing",
    preprocessingUI("preprocessing")
  ),

  ##### GENE SETS #####
  bslib::nav_panel(
    "Gene Sets",
    geneSetsUI("gene_sets")
  ),

  ##### BENCHMARKING MODE #####
  bslib::nav_panel(
    "Benchmarking Mode",
    benchmarkingUI("benchmarking")
  ),

  ##### DISCOVERY MODE #####
  bslib::nav_panel(
    "Discovery Mode",
    discoveryUI("discovery")
  ),

  bslib::nav_spacer(),

  ##### TUTORIAL (right-aligned navbar link) #####
  bslib::nav_item(
    shiny::tags$a(
      shiny::icon("circle-info"), " Tutorial",
      href = "tutorial/Tutorial_markeRShinyWebApp_shinyv0_99_7.pdf",
      target = "_blank", rel = "noopener noreferrer",
      class = "nav-link",
      title = "Open the markeR user tutorial (PDF)"
    )
  ),

  ##### REPORT BUG (right-aligned navbar link) #####
  bslib::nav_item(
    shiny::tags$a(
      shiny::icon("bug"), " Report bug",
      href = "https://github.com/DiseaseTranscriptomicsLab/markeR/issues",
      target = "_blank", rel = "noopener noreferrer",
      class = "nav-link"
    )
  ),

  ##### SOURCE CODE (right-aligned navbar link) #####
  bslib::nav_item(
    shiny::tags$a(
      shiny::icon("github"), " Source Code",
      href = "https://github.com/DiseaseTranscriptomicsLab/markeR/tree/ShinyApp",
      target = "_blank", rel = "noopener noreferrer",
      class = "nav-link",
      title = "View the markeR source code on GitHub"
    )
  )



)

# Then wrap whole page with head() + global fixed-bottom log bar
ui <- tagList(
  tags$head(
    tags$style(HTML("
      /* Belt-and-braces on top of fillable = FALSE above: force the page
         itself to always be scrollable and never height-clipped, since
         bslib/htmltools can rebuild this constraint in ways that vary by
         rendering context (plain browser vs. an iframe-embedded deployment
         like ShinyProxy). */
      html, body {
        height: auto !important;
        min-height: 100% !important;
        overflow-y: auto !important;
      }
      /* Reserve enough space at the bottom of every page so its own content
         never ends up underneath the fixed log bar below (#global-log-bar
         is position:fixed, so it always covers this much of the viewport's
         bottom edge regardless of scroll position). */
      body { padding-bottom: 60px; }
      #shiny-notification-panel {
        top: auto !important;
        bottom: 56px !important;
        right: 30px !important;
        left: auto !important;
        width: 320px !important;
      }
      /* withProgress bar color */
      .shiny-progress .progress-bar {
        background-color: rgb(186, 89, 0) !important;
      }
      /* Warning notifications */
      #shiny-notification-panel .shiny-notification-warning {
        border-left: 4px solid #ffc107 !important;
      }
      /* ── Global log bar (subtle dark grey, matching the app's neutral
         palette instead of standing out) ───────────── */
      #global-log-bar {
        position: fixed; bottom: 0; left: 0; right: 0;
        z-index: 1050;
        background: #495057;       /* medium-dark grey background */
        color: #f1f3f5;
        font-family: monospace; font-size: 0.79em;
        box-shadow: 0 -3px 10px rgba(0,0,0,0.2);
      }
      #global-log-header {
        display: flex; align-items: center;
        justify-content: space-between;
        padding: 5px 16px; cursor: pointer;
        border-top: 2px solid #6c757d;
        user-select: none;
      }
      #global-log-header:hover { background: #565e64; }
      #global-log-body {
        max-height: 0px; overflow: hidden;
        transition: max-height 0.22s ease;
        background: #3d4348;
      }
      #global-log-entries {
        padding: 6px 16px 10px 16px;
        max-height: 200px; overflow-y: auto;
        line-height: 1.9; color: #dee2e6;
      }
      /* Most recent entry highlighted, subtly, in off-white */
      #global-log-entries .glog-entry:first-child { color: #f8f9fa; font-weight: 600; }
      /* Download button in log bar */
      #global-log-entries .btn { color: #dee2e6; border-color: #6c757d; }

      /* ── Locked nav tabs ─────────────────────────────── */
      .nav-tab-locked {
        color: #adb5bd !important;
        pointer-events: none !important;
        cursor: not-allowed !important;
        opacity: 0.5 !important;
        position: relative;
      }
      /* Wrapper to re-enable pointer-events so tooltip fires */
      .nav-tab-locked-wrap {
        pointer-events: auto !important;
        cursor: not-allowed !important;
        display: inline-block;
        position: relative;
      }
      /* tooltip rendered by JS into <body> to avoid overflow clipping */
      #lock-tip-floating {
        position: fixed;
        background: rgba(30,30,30,0.93);
        color: #fff;
        font-size: 0.83em;
        line-height: 1.45;
        max-width: 300px;
        padding: 7px 11px;
        border-radius: 5px;
        border: 1px solid rgba(255,255,255,0.14);
        box-shadow: 0 3px 10px rgba(0,0,0,0.4);
        pointer-events: none;
        z-index: 99999;
        transition: opacity 0.15s;
      }
    "))
  ),
  # DT outputs in this app are built dynamically inside shiny::renderUI()
  # (only appearing after the user clicks a "Run" button), and none exist in
  # the app's static UI otherwise. With no DT widget present in the very
  # first page load, the browser can end up never being told to load
  # DataTables' JS at all. This injects DataTables' JS/CSS dependencies into
  # the initial HTML <head> directly, regardless of which tab loads first.
  # NOTE: this used to be a hidden, always-on DT::DTOutput()/renderDT() pair
  # instantiating a real (empty) table off-screen purely to trigger asset
  # loading - that live widget instance was the source of a client-side
  # "Cannot read properties of null (reading 'lazyRender')" error thrown
  # from DT's own htmlwidgets binding code (datatables.js) on every page
  # load, which could disrupt whatever other Shiny output updates happened
  # to be dispatched in the same batch (e.g. the global log bar rendering
  # blank until an unrelated later reactive flush). Declaring the
  # dependencies directly - with no live DT widget ever created - preloads
  # the same JS/CSS without that failure mode.
  tags$head(
    htmltools::findDependencies(DT::datatable(data.frame(x = 1)))
  ),
  ui,
  # Global log bar (fixed to bottom of viewport, visible on all tabs)
  tags$div(
    id = "global-log-bar",
    tags$div(
      id = "global-log-header",
      onclick = "toggleMarkeRLog()",
      # NOTE: deliberately NOT inline = TRUE here (see markeR_tab_data.R,
      # near outputOptions(..., suspendWhenHidden = FALSE), for the full
      # story on why). Separately: uiOutput()'s `...` is passed straight
      # through to its container tag, so anything given here shows up as
      # placeholder content in the raw, pre-Shiny HTML - and gets replaced
      # once the real renderUI() value arrives. Even with
      # suspendWhenHidden = FALSE, this output's very first computed value
      # does not appear to land on the literal first flush (only from the
      # next one onward, e.g. the first time you interact with the app) -
      # this placeholder means the bar never looks empty in the meantime.
      shiny::uiOutput("global_log_header", inline = FALSE, container = shiny::tags$div, fill = FALSE,
        shiny::tags$span(style = "color:#e2e8f0; font-size:0.88em;", "📋 Session log")
      ),
      tags$span(id = "glog-arrow", style = "font-size:0.85em;", "▲")
    ),
    tags$div(
      id = "global-log-body",
      tags$div(
        id = "global-log-entries",
        # Placeholder mirrors real global_log_entries output's markup
        # exactly (see R/markeR_tab_data.R), built from the same session
        # preamble computed once above - so real session/package info is
        # visible from the very first byte of HTML, with no dependency on
        # when (or whether) Shiny's own render happens to catch up. Once
        # the user does anything, the real renderUI() output replaces this
        # with the identical "Session info" block plus the action log.
        shiny::uiOutput("global_log_entries", inline = FALSE, container = shiny::tags$div, fill = FALSE,
          shiny::tags$details(
            style = "margin-bottom:6px;",
            shiny::tags$summary(
              style = "font-size:0.8em; color:#94a3b8; cursor:pointer; user-select:none;",
              "Session info"
            ),
            shiny::tags$pre(
              style = paste0(
                "font-size:0.75em; color:#94a3b8; padding:4px 0; margin:4px 0 0;",
                "background:transparent; border:none; white-space:pre-wrap; word-break:break-all;",
                "max-height:160px; overflow-y:auto;"
              ),
              paste(sub("^\\[Session\\]\\s*", "", .markeR_initial_session_lines), collapse = "\n")
            )
          ),
          shiny::tags$div(style = "color:#64748b; padding:4px 0;", "No actions logged yet.")
        )
      )
    )
  ),
  tags$script(HTML("
    function toggleMarkeRLog() {
      var body  = document.getElementById('global-log-body');
      var arrow = document.getElementById('glog-arrow');
      if (!body.style.maxHeight || body.style.maxHeight === '0px') {
        body.style.maxHeight = '220px';
        arrow.textContent = '▼';
      } else {
        body.style.maxHeight = '0px';
        arrow.textContent = '▲';
      }
    }

    /* ── Nav tab locking ─────────────────────────────── */
    // Find a navbar <a> by its visible text label (trim & case-insensitive)
    function _findNavLink(label) {
      var links = document.querySelectorAll('.navbar-nav .nav-link, .nav-tabs .nav-link');
      for (var i = 0; i < links.length; i++) {
        if (links[i].textContent.trim().toLowerCase() === label.trim().toLowerCase())
          return links[i];
      }
      return null;
    }

    // Floating tooltip helper - appended to <body> so it is never clipped by overflow:hidden
    (function() {
      var _tip = null;
      function _showTip(text, rect) {
        if (!_tip) {
          _tip = document.createElement('div');
          _tip.id = 'lock-tip-floating';
          document.body.appendChild(_tip);
        }
        _tip.textContent = text;
        _tip.style.opacity = '0';
        _tip.style.display = 'block';
        // Position above the hovered element, centred
        var tipW = _tip.offsetWidth || 220;
        var left = rect.left + rect.width / 2 - tipW / 2;
        var top  = rect.top - _tip.offsetHeight - 8;
        if (left < 6) left = 6;
        if (top  < 6) top  = rect.bottom + 6;   // flip below if no room above
        _tip.style.left = left + 'px';
        _tip.style.top  = top  + 'px';
        _tip.style.opacity = '1';
      }
      function _hideTip() {
        if (_tip) { _tip.style.opacity = '0'; _tip.style.display = 'none'; }
      }
      window._lockTipShow = _showTip;
      window._lockTipHide = _hideTip;
    })();

    Shiny.addCustomMessageHandler('lockNavTabs', function(msg) {
      // msg: { tabs: ['Preprocessing','Gene Sets',...], locked: true/false, tip: '...' }
      var tabs = msg.tabs || [];
      var locked = msg.locked !== false;
      var tip = msg.tip || 'Load and validate data first';
      tabs.forEach(function(label) {
        var a = _findNavLink(label);
        if (!a) return;
        if (locked) {
          a.classList.add('nav-tab-locked');
          // Wrap in span if not already, attach mouse listeners
          var parent = a.parentNode;
          if (!parent.classList.contains('nav-tab-locked-wrap')) {
            var wrap = document.createElement('span');
            wrap.className = 'nav-tab-locked-wrap';
            wrap.setAttribute('data-lock-tip', tip);
            parent.insertBefore(wrap, a);
            wrap.appendChild(a);
            wrap.addEventListener('mouseenter', function() {
              window._lockTipShow(this.getAttribute('data-lock-tip'),
                                  this.getBoundingClientRect());
            });
            wrap.addEventListener('mouseleave', window._lockTipHide);
          } else {
            parent.setAttribute('data-lock-tip', tip);
          }
        } else {
          a.classList.remove('nav-tab-locked');
          // Unwrap tooltip span if present
          var wrap = a.closest('.nav-tab-locked-wrap');
          if (wrap) {
            wrap.parentNode.insertBefore(a, wrap);
            wrap.parentNode.removeChild(wrap);
          }
          window._lockTipHide();
        }
      });
    });
  "))
)

# ---- Server ----

server <- function(input, output, session) {

  #bs_themer()

  ###### DATA TAB ######
  # All expression/metadata/gene-set loading, validation, and preview logic
  # lives in R/markeR_tab_data.R (dataUI()/dataServer()) - see that file's
  # header comment for why it isn't a namespaced shiny module like the other
  # tabs. Everything downstream (preprocessing, gene sets, benchmarking,
  # discovery, and the nav-tab locking below) consumes the reactives it
  # returns, exactly as it consumed the equivalent local reactiveVals before.
  data_mod          <- dataServer(input, output, session)
  expr_data         <- data_mod$expr_data
  meta_data         <- data_mod$meta_data
  gene_sets         <- data_mod$gene_sets
  log_step          <- data_mod$log_step
  data_log          <- data_mod$data_log
  expr_quality_warn <- data_mod$expr_quality_warn

  ###### PREPROCESSING MODULE ######

  pp_module <- preprocessingServer(
    id       = "preprocessing",
    get_expr = expr_data,
    get_meta = meta_data,
    log_fn   = log_step,
    get_log  = data_log
  )

  # When the user finalises preprocessing, push the processed data back into
  # the app's main expr_data/meta_data reactives so downstream tabs pick it up.
  # Metadata must be pushed back too (not just expr): if samples were removed
  # during preprocessing, meta_data() otherwise keeps its original, larger row
  # count forever - Data tab and downstream tabs would then see fewer
  # expression columns than metadata rows.
  shiny::observeEvent(pp_module$finalized(), {
    shiny::req(pp_module$finalized() > 0L, pp_module$final_data())
    expr_data(pp_module$final_data())
    if (!is.null(pp_module$final_meta())) meta_data(pp_module$final_meta())
    expr_quality_warn(NULL)   # clear any earlier quality warning
  }, ignoreInit = TRUE)

  ###### ANALYSIS MODULES ######

  geneSetsServer(
    id            = "gene_sets",
    get_expr      = expr_data,
    get_meta      = meta_data,
    get_gene_sets = gene_sets
  )

  benchmarkingServer(
    id            = "benchmarking",
    get_expr      = expr_data,
    get_meta      = meta_data,
    get_gene_sets = gene_sets
  )

  discoveryServer(
    id            = "discovery",
    get_expr      = expr_data,
    get_meta      = meta_data,
    get_gene_sets = gene_sets
  )

  ###### NAV TAB LOCKING ######

  # Tabs that require valid, matched data AND finalised preprocessing.
  .analysis_tabs <- c("Benchmarking Mode", "Discovery Mode")

  # Mirror exactly the sample_mismatch_banner logic.
  # Tabs are locked ONLY when both expr and meta are loaded AND sample names/counts mismatch.
  # If either is missing, tabs are unlocked (user is still configuring - don't block navigation).
  # .tab_lock_state() is defined at top level in R/markeR_tab_data.R (shared
  # with that tab's own "proceed" banner).
  shiny::observe({
    expr      <- expr_data()
    meta      <- meta_data()
    finalized <- pp_module$finalized()   # reactive dependency

    # ── Preprocessing tab: locked only when data is not loaded / mismatched.
    pp_state <- .tab_lock_state(expr, meta)
    session$sendCustomMessage("lockNavTabs", list(
      tabs   = list("Preprocessing"),
      locked = pp_state$locked,
      tip    = if (!is.null(pp_state$tip)) pp_state$tip else ""
    ))

    # ── Gene Sets / Benchmarking / Discovery: all require preprocessing to be finalised.
    analysis_state <- if (pp_state$locked) {
      pp_state   # same "fix your data first" message
    } else if (finalized == 0L) {
      list(locked = TRUE,
           tip = paste0("Complete preprocessing first - click ",
                        "‘Finalise’ or ‘Use Data As-Is’ in the Preprocessing tab."))
    } else {
      list(locked = FALSE, tip = NULL)
    }

    session$sendCustomMessage("lockNavTabs", list(
      tabs   = as.list(c("Gene Sets", .analysis_tabs)),
      locked = analysis_state$locked,
      tip    = if (!is.null(analysis_state$tip)) analysis_state$tip else ""
    ))
  })


}

# ---- Launcher ----

#' Launch the markeR Shiny Application
#'
#' Starts the markeR interactive Shiny app (Benchmarking and Discovery Mode
#' tabs) in the default web browser, or as a standalone server when called
#' with `host`/`port` (e.g. from a Docker container).
#'
#' @param ... Additional arguments passed to \code{\link[shiny]{runApp}}
#'   (e.g. \code{host}, \code{port}, \code{launch.browser}).
#'
#' @return Does not return; runs the Shiny application until interrupted.
#'
#' @examples
#' \dontrun{
#' markeRapp()
#'
#' # As typically run inside a Docker container:
#' markeRapp(host = "0.0.0.0", port = 3838, launch.browser = FALSE)
#' }
#'
#' @export
markeRapp <- function(...){

  # Register a URL path "/logo" that points to your package's inst/figures folder
  shiny::addResourcePath(
    "logo",
    system.file("figures", package = "markeR")
  )
  
  # Workflow schematic
  shiny::addResourcePath(
    "methods",
    system.file("shiny/www", package = "markeR")   # points to inst/shiny/www/
  ) 

  # User tutorial (PDF), served the same way as the workflow schematic above -
  # a separate alias so the URL is self-explanatory (/tutorial/...) even
  # though it currently resolves to the same inst/shiny/www/ folder. Works
  # identically when deployed on a web server, since addResourcePath() (not a
  # local file:// link) is what's linked to in the navbar below.
  shiny::addResourcePath(
    "tutorial",
    system.file("shiny/www", package = "markeR")
  )

  options(shiny.maxRequestSize = 1000 * 1024^2)
  
  app <- shiny::shinyApp(ui, server)
  shiny::runApp(app, ...)
}