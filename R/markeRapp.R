#' @importFrom bslib bs_theme nav_panel page_navbar sidebar layout_sidebar accordion accordion_panel navset_card_tab card card_header nav_spacer nav_item
#' @import shiny
#' @import GEOquery
#' @importFrom DT DTOutput renderDT datatable
#' @importFrom readxl read_xls read_xlsx
#' @importFrom R.utils gunzip
#' @importFrom data.table fread
#' @importFrom Biobase exprs pData annotation experimentData
#' @importFrom shinyWidgets pickerInput pickerOptions
#' @importFrom utils download.file read.csv read.delim read.table untar unzip packageVersion
#' @importFrom tools file_ext file_path_sans_ext
#' @importFrom ComplexHeatmap draw

ui <- bslib::page_navbar(
  
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
    shiny::a("Web App 0.99.1",
             style = "font-size: 0.7em; color: #949494; text-decoration: none;")
  ),
  
  theme = bslib::bs_theme(bootswatch = "yeti",
                          primary = "rgb(186, 89, 0)"),
  
  
  # ---- Tabs ----
  
  ##### ABOUT #####
  
  bslib::nav_panel(
    "About", 
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
            paste0("Have a collection of genes that putatively mark a phenotype, but not sure how to combine them and turn them into  ",
                   "<span style='color:#306F1D;'> meaningful metrics for phenotype quantification</span>", 
                   "? markeR helps you evaluate gene sets as phenotype markers and explore their information. Here’s a quick guide to each 
                   section of the app:"
            )
          ),
          
          style = "font-size: 1.2em;"
        )
      )
    ),

    shiny::tags$div(
      style = paste(
        "max-width: 1000px; margin: 0 auto 20px; padding: 14px 18px;",
        "background:#f8f9fa; border:1px solid #dee2e6; border-radius:8px;",
        "font-size: 0.9em; color:#495057; text-align:center;"
      ),
      shiny::icon("dna", style = "color:#306F1D; margin-right:6px;"),
      "This app was built with the ",
      shiny::tags$b("markeR"), " Bioconductor package ",
      shiny::a(shiny::tags$b("v1.2"),
        href = "https://github.com/DiseaseTranscriptomicsLab/markeR/releases/tag/v1.2",
        target = "_blank", rel = "noopener noreferrer"),
      ". For more on the methodology behind it, see ",
      shiny::a("our paper",
        href = "https://doi.org/10.1093/nargab/lqag057",
        target = "_blank", rel = "noopener noreferrer"),
      "."
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
           Alternatively, you can import data directly from GEO or explore example data on senescence."),
        
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
     Enrichment and Score. Enrichment uses the Gene Set Enrichment Analyses (GSEA) method based on differential gene expression, 
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
     marking the phenotype, in contrast to benchmarking mode, 
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
    
    bslib::nav_panel(
      "Data",
      
      bslib::layout_sidebar(
        
        ###### SIDEBAR ######
        
        sidebar = bslib::sidebar(
          width = 350,
          
          shiny::div(
            style = "padding-bottom:15px;",
            shiny::h4("Data", style = "font-weight:700; color:#306F1D; margin-bottom:4px;"),
            shiny::p(
              "Load your gene expression matrix, sample metadata, and gene sets.",
              " Everything downstream depends on what you configure here.",
              style = "color:#6c757d; font-size:0.85em; margin-bottom:6px;"
            ),
            shiny::hr()
          ),
          
          bslib::accordion(
            
            ####### > EXPRESSION DATA #######
            
            bslib::accordion_panel(
              "Gene Expression Data",
              
              shiny::radioButtons(
                "expr_source",
                "Source:",
                choices = c(
                  "Use Example Data \n(Unprocessed)" = "example_raw",
                  "Use Example Data \n(Processed)"   = "example_proc",
                  "Upload File"                    = "upload",
                  "Retrieve from GEO"              = "geo"
                )
              ),
              
              # Help text for unprocessed example
              shiny::conditionalPanel(
                "input.expr_source == 'example_raw'",
                shiny::helpText(
                  shiny::tags$strong("Example Unprocessed Data:"),
                  "This dataset is a manual compilation of RNA-seq experiments on senescence in human cell lines, treated with different senescence inducers and their respective proliferative and quiescent controls. ",
                  "It has been used in Martins-Silva et al., 2026 (",
                  shiny::tags$a(href="https://doi.org/10.1093/nargab/lqag057", "NAR Genomics and Bioinformatics", target="_blank"),
                  "). The raw read counts for all samples are available for download through ",
                  shiny::tags$a(href="https://zenodo.org/records/18714122", "Zenodo", target="_blank"),
                  ". Genes are in rows, samples in columns."
                )
              ),
              
              # Help text for processed example
              shiny::conditionalPanel(
                "input.expr_source == 'example_proc'",
                shiny::helpText(
                  shiny::tags$strong("Example Processed Data:"),
                  "This dataset contains the same RNA-seq experiments on senescence in human cell lines, as described for the unprocessed data above, but filtered for lowly expressed genes, normalised, and batch-corrected. ",
                  "Processing followed the methods described in Martins-Silva et al., 2026 (",
                  shiny::tags$a(href="https://doi.org/10.1093/nargab/lqag057", "NAR Genomics and Bioinformatics", target="_blank"),
                  "). The processed data is available via ",
                  shiny::tags$a(href="https://zenodo.org/records/18714122", "Zenodo", target="_blank"),
                  ". Genes are in rows, samples in columns."
                )
              ),
              shiny::conditionalPanel(
                "input.expr_source == 'upload'",
                
                shiny::fileInput(
                  "expr_file",
                  "Upload expression matrix (.csv, .txt, .rds)",
                  accept = c(".csv", ".rds", ".txt")
                ),
                
                shiny::helpText(
                  
                  shiny::tags$p(shiny::tags$b("💡 Expected structure:")),
                  shiny::tags$ul(
                    shiny::tags$li("Genes in rows"),
                    shiny::tags$li("Samples in columns")
                  ),
                  
                  shiny::tags$p(shiny::tags$b("Text or CSV files (.csv, .txt)")),
                  shiny::tags$ul(
                    shiny::tags$li("First column: gene names"),
                    shiny::tags$li("First row: sample names"),
                    shiny::tags$li("First cell can be empty or contain 'Gene'"),
                    shiny::tags$li("Separator can be comma, semicolon, or tab (the app will try to detect it automatically)"),
                    shiny::tags$li("Numeric values can use '.' or ',' as decimal separator"),
                    shiny::tags$li("Trailing separators in a row will be ignored"),
                    shiny::tags$li("Do not use quotation marks around values")
                  ),
                  
                  shiny::tags$p(shiny::tags$b(".rds")),
                  shiny::tags$ul(
                    shiny::tags$li("Must contain a matrix or data frame"),
                    shiny::tags$li("Gene names as row names"),
                    shiny::tags$li("Sample names as column names")
                  ),
                ),

                # Sample management controls (shown after upload)
                shiny::uiOutput("upload_expr_manage_ui")
              ),

              shiny::conditionalPanel(
                "input.expr_source == 'geo'",
                # ── New adaptive GEO import module ──────────────────────
                geoImportUI("geo_import")
              )
            ),
            
            
            ####### > METADATA #######
            
            
            bslib::accordion_panel(
              "Metadata",
              
              shiny::radioButtons(
                "meta_source",
                "Source:",
                choices = c(
                  "Use Example Metadata" = "example",
                  "Upload Metadata"      = "upload",
                  "From GEO Selection"   = "geo"
                )
              ),
              
              shiny::conditionalPanel(
                "input.meta_source == 'upload'",
                shiny::fileInput(
                  "meta_file",
                  "Upload metadata (.csv, .txt, .tsv, .xlsx, .rds)",
                  accept = c(".csv", ".txt", ".tsv", ".xls", ".xlsx", ".rds")
                ),
                shiny::helpText(
                  shiny::tags$p(shiny::tags$b("💡 Expected structure:")),
                  shiny::tags$ul(
                    shiny::tags$li("One row per sample"),
                    shiny::tags$li(
                      shiny::tags$b("First column must contain the Sample IDs"),
                      " - these will be matched to expression matrix column names"
                    ),
                    shiny::tags$li("Remaining columns: any sample-level variables (group, condition, batch, …)")
                  ),
                  shiny::tags$p(shiny::tags$b("Text / CSV / TSV files:")),
                  shiny::tags$ul(
                    shiny::tags$li("Separator detected automatically (comma, semicolon, or tab)"),
                    shiny::tags$li("First row must be a header with column names"),
                    shiny::tags$li("Quoted values are handled automatically")
                  ),
                  shiny::tags$p(shiny::tags$b(".rds / .xlsx:")),
                  shiny::tags$ul(
                    shiny::tags$li("Must contain a data frame with samples as rows"),
                    shiny::tags$li("Column names used as variable labels")
                  )
                ),
                # Match info + alignment controls (shown after upload)
                shiny::uiOutput("upload_meta_manage_ui")
              ),
              # Show GEO SampleID column selector when using GEO as metadata source
              shiny::conditionalPanel(
                "input.meta_source == 'geo'",
                shiny::uiOutput("geo_sampleid_selector"),
                shiny::uiOutput("geo_meta_manage_ui")
              )
            ),
            
            
            ####### > GENE SETS #######
            
            
            bslib::accordion_panel(
              "Gene Sets",
              
              shiny::radioButtons(
                "geneset_source",
                "Source:",
                choices = c(
                  "Use Example Gene Sets" = "example",
                  "Upload Gene Sets"      = "upload"
                )
              ),
              
              shiny::conditionalPanel(
                "input.geneset_source == 'upload'",
                shiny::fileInput(
                  "geneset_files",
                  "Upload one or more files (.csv, .rds, .rda)",
                  accept = c(".csv", ".rds", ".rda"),
                  multiple = TRUE
                ),
                shiny::helpText(
                  # string of text explaining that if the input is a csv it should have either one column with gene names, or two columns the first with gene names and the second with expected direction of change (+1/-1). If a .rds, can be either a vector of gene names, a data frame where the first column is gene names and the second the direction of change, or a named list with one entry per gene set (following the rules for vectors and data frames)
                  
                  shiny::helpText(
                    shiny::tags$strong("Accepted input formats for gene sets:"),
                    shiny::tags$ul(
                      shiny::tags$li(
                        shiny::strong(".csv files:"),
                        " must contain either:",
                        shiny::tags$ul(
                          shiny::tags$li("One column with gene names"),
                          shiny::tags$li("Two columns where the first column contains gene names and the second column contains the expected direction of enrichment (+1 or -1)")
                        ),
                        shiny::tags$em("The file name (without extension) will be used as the gene set name.")
                      ),
                      shiny::tags$li(
                        shiny::strong(".rds/.rda files:"),
                        " may contain:",
                        shiny::tags$ul(
                          shiny::tags$li("A character vector of gene names"),
                          shiny::tags$li("A data frame (first column = gene names, second column = direction of change)"),
                          shiny::tags$li("A named list of gene sets, where each element follows one of the formats above")
                        ),
                        shiny::tags$em("If uploading individual files, the file name will be used as the gene set name.")
                      )
                    )
                  )
                  
                )
              )
            )
          )
        ),
        
        
        ###### > MAIN PREVIEW AREA ###### 
        tags$style(
          HTML("
    /* Make row names column wider ONLY for expr_preview */
    #expr_preview table.dataTable th:first-child,
    #expr_preview table.dataTable td:first-child {
      width: 200px !important;   /* desired width */
      max-width: 200px !important;
      white-space: nowrap;
      overflow: hidden;
      text-overflow: ellipsis;
    }
  ")
        ),
        
        bslib::card(
          bslib::card_header("Data Preview"),

          # Proceed / mismatch banner - always visible at the top of the card
          shiny::uiOutput("data_proceed_btn"),

          shiny::div(style = "overflow-y: auto; max-height: 900px;",  # <-- scrollable container
                     shiny::uiOutput("summary_ui"),
                     shiny::hr(),
                     shiny::h5("Expression Preview"),
                     shiny::uiOutput("expr_quality_banner"),
                     DT::DTOutput("expr_preview"),
                     shiny::h5("Metadata Preview"),
                     DT::DTOutput("meta_preview"),
                     shiny::h5("Gene Sets Overview"),
                     shiny::tags$div(
                       style = "margin-bottom: 10px; font-style: italic; color: #555;",
                       "💡 Tip: Click on a gene set row to see its genes and expected direction of regulation change."
                     ),
                     DT::DTOutput("geneset_summary")
          )
        )
        
        
      )
    )
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

  ##### REPORT BUG (right-aligned navbar link) #####
  bslib::nav_item(
    shiny::tags$a(
      shiny::icon("bug"), " Report bug",
      href = "https://github.com/DiseaseTranscriptomicsLab/markeR/issues",
      target = "_blank", rel = "noopener noreferrer",
      class = "nav-link"
    )
  )



)

# Then wrap whole page with head() + global fixed-bottom log bar
ui <- tagList(
  tags$head(
    tags$style(HTML("
      body { padding-bottom: 40px; }
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
      /* ── Global log bar (orange palette) ───────────── */
      #global-log-bar {
        position: fixed; bottom: 0; left: 0; right: 0;
        z-index: 1050;
        background: #7a3800;       /* deep burnt-orange background */
        color: #fde8c8;
        font-family: monospace; font-size: 0.79em;
        box-shadow: 0 -3px 10px rgba(0,0,0,0.28);
      }
      #global-log-header {
        display: flex; align-items: center;
        justify-content: space-between;
        padding: 5px 16px; cursor: pointer;
        border-top: 2px solid #EBB43E;
        user-select: none;
      }
      #global-log-header:hover { background: #8f4200; }
      #global-log-body {
        max-height: 0px; overflow: hidden;
        transition: max-height 0.22s ease;
        background: #5c2a00;
      }
      #global-log-entries {
        padding: 6px 16px 10px 16px;
        max-height: 200px; overflow-y: auto;
        line-height: 1.9; color: #fde8c8;
      }
      /* Most recent entry highlighted in bright orange */
      #global-log-entries .glog-entry:first-child { color: #EBB43E; font-weight: 600; }
      /* Download button in log bar */
      #global-log-entries .btn { color: #fde8c8; border-color: #EBB43E; }

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
  # DataTables' JS at all. This hidden, always-present output forces that
  # dependency into the initial HTML regardless of which tab loads first.
  # NOTE: deliberately NOT `visibility:hidden` / zero-size - some widget
  # libraries treat invisible/zero-area elements as "not really on screen"
  # and skip or defer full initialization, which would silently defeat the
  # whole point of this block. Positioning it off-screen (but still a real,
  # normally-sized, "visible" element as far as layout/rendering is
  # concerned) keeps it out of the way without triggering that shortcut.
  shiny::div(
    style = "position:absolute; top:0; left:-9999px; width:300px; height:200px; overflow:hidden;",
    DT::DTOutput("markeR_dep_dt")
  ),
  ui,
  # Global log bar (fixed to bottom of viewport, visible on all tabs)
  tags$div(
    id = "global-log-bar",
    tags$div(
      id = "global-log-header",
      onclick = "toggleMarkeRLog()",
      shiny::uiOutput("global_log_header", inline = TRUE),
      tags$span(id = "glog-arrow", style = "font-size:0.85em;", "▲")
    ),
    tags$div(
      id = "global-log-body",
      tags$div(
        id = "global-log-entries",
        shiny::uiOutput("global_log_entries")
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

  # Minimal, invisible render matching the hidden DT output declared in the
  # static UI above; see the comment there for why it exists.
  output$markeR_dep_dt <- DT::renderDT({
    DT::datatable(data.frame(x = numeric(0)))
  })

  ###### REACTIVE STORAGE CONTAINERS ######
  expr_data   <- shiny::reactiveVal(NULL)  # Expression matrix/data.frame
  meta_data   <- shiny::reactiveVal(NULL)  # Sample metadata
  gene_sets   <- shiny::reactiveVal(NULL)  # Named list of gene sets
  geo_objects <- shiny::reactiveVal(NULL)  # Raw GEO objects (can be multiple)
  geo_meta_cols <- reactiveVal(NULL) # Will store the original metadata column names (excluding SampleID)
  full_geo_meta      <- reactiveVal(NULL)
  upload_meta_raw    <- reactiveVal(NULL)   # Raw uploaded metadata (for SampleID re-mapping)
  module_meta_active <- reactiveVal(FALSE)  # TRUE while geo_import_module owns meta_data()
  # Session log: initialise with R + all DESCRIPTION dependency versions.
  .session_preamble <- tryCatch({
    # Parse Imports + Suggests + Depends from the installed package DESCRIPTION
    .parse_dep_field <- function(field) {
      if (is.null(field) || is.na(field)) return(character(0))
      pkgs <- strsplit(field, ",")[[1]]
      pkgs <- trimws(sub("\\s*\\(.*?\\)", "", pkgs))  # strip version specs like (>= 4.0)
      pkgs[nzchar(pkgs) & pkgs != "R"]
    }
    desc  <- utils::packageDescription("markeR",
                                        fields = c("Imports", "Suggests", "Depends"))
    pkgs  <- sort(unique(c(
      "markeR",
      .parse_dep_field(desc$Imports),
      .parse_dep_field(desc$Suggests),
      .parse_dep_field(desc$Depends)
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
  }, error = function(e) character(0))
  data_log <- reactiveVal(.session_preamble)  # Session log with timestamped entries
  expr_quality_warn  <- reactiveVal(NULL)           # Warning message for suspicious expression data
  geo_supp_files  <- reactiveVal(NULL)
  expr_candidates_available <- reactiveVal(FALSE) # Flag to indicate multiple candidate files are available after auto detection
  
  # Append a timestamped entry to the processing log.
  # isolate() prevents the calling observer from taking a reactive dependency
  # on data_log - without it, writing data_log() re-invalidates the observer
  # that called log_step(), causing an infinite reactive loop.
  log_step <- function(msg) {
    ts <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
    data_log(c(isolate(data_log()), paste(ts, msg)))
  }

  ###### PREPROCESSING MODULE ######

  pp_module <- preprocessingServer(
    id       = "preprocessing",
    get_expr = expr_data,
    get_meta = meta_data,
    log_fn   = log_step,
    get_log  = data_log
  )

  # When the user finalises preprocessing, push the processed data back into
  # the app's main expr_data reactive so downstream tabs pick it up.
  shiny::observeEvent(pp_module$finalized(), {
    shiny::req(pp_module$finalized() > 0L, pp_module$final_data())
    expr_data(pp_module$final_data())
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
  .tab_lock_state <- function(expr, meta) {
    # No expression data → always lock (includes failed / pending downloads)
    if (is.null(expr))
      return(list(locked = TRUE,
                  tip = "No expression data loaded. Import data in the Data tab first."))
    if (!is.data.frame(expr) || ncol(expr) == 0L)
      return(list(locked = TRUE,
                  tip = "Expression data is empty. Import data in the Data tab first."))
    # No metadata → lock
    if (is.null(meta) || !is.data.frame(meta) || nrow(meta) == 0L)
      return(list(locked = TRUE,
                  tip = "No metadata loaded. Import metadata in the Data tab first."))

    n_expr <- ncol(expr)
    n_meta <- nrow(meta)
    e_nms  <- colnames(expr)

    count_ok <- (n_expr == n_meta)

    # Same multi-candidate check as sample_mismatch_banner:
    # try first column, SampleID column, and rownames as possible ID vectors.
    names_ok <- FALSE
    if (count_ok) {
      candidates <- unique(list(
        as.character(meta[[ colnames(meta)[1] ]]),
        if (!is.null(meta[["SampleID"]])) as.character(meta[["SampleID"]]) else NULL,
        rownames(meta)
      ))
      candidates <- Filter(Negate(is.null), candidates)
      names_ok <- any(vapply(candidates, function(ids) setequal(e_nms, ids), logical(1)))
    }

    if (count_ok && names_ok)
      return(list(locked = FALSE, tip = NULL))

    tip <- if (!count_ok)
      paste0("Sample mismatch: ", n_expr, " expression samples vs ", n_meta,
             " metadata rows. Fix in the Data tab first.")
    else
      paste0("Sample names do not match between expression data and metadata. Fix in the Data tab first.")

    list(locked = TRUE, tip = tip)
  }

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
                        "’Finalize’ or ‘Use Data As-Is’ in the Preprocessing tab."))
    } else {
      list(locked = FALSE, tip = NULL)
    }

    session$sendCustomMessage("lockNavTabs", list(
      tabs   = as.list(c("Gene Sets", .analysis_tabs)),
      locked = analysis_state$locked,
      tip    = if (!is.null(analysis_state$tip)) analysis_state$tip else ""
    ))
  })


  ###### HELPER FUNCTIONS ######

  # Show a centred spinner modal (same style as the GEO import module).
  # Call on.exit(shiny::removeModal()) in the caller so it always closes.
  show_loading_modal <- function(title, subtitle = NULL) {
    shiny::showModal(shiny::modalDialog(
      title     = NULL,
      footer    = NULL,
      easyClose = FALSE,
      shiny::tags$div(
        style = "text-align:center;padding:40px 20px;",
        shiny::tags$div(
          class = "spinner-border text-primary",
          style = "width:3rem;height:3rem;",
          role  = "status",
          shiny::tags$span(class = "visually-hidden", "Loading…")
        ),
        shiny::tags$h5(title, style = "margin-top:18px;color:#333;"),
        if (!is.null(subtitle))
          shiny::tags$p(subtitle,
                        style = "color:#777;font-size:0.9em;margin-top:6px;")
      )
    ))
  }

  # Show a Zenodo download-failure modal with a "Try Again" button.
  # retry_btn_id: inputId of the actionButton placed in the footer.
  show_dl_error_modal <- function(msg, retry_btn_id) {
    shiny::showModal(shiny::modalDialog(
      title = shiny::tags$span(
        shiny::icon("circle-exclamation", style = "color:#dc3545;margin-right:6px;"),
        "Example Data Download Failed"
      ),
      shiny::tags$div(
        shiny::tags$p("Could not download the example data from Zenodo."),
        shiny::tags$p(style = "font-size:0.84em;color:#888;font-family:monospace;word-break:break-all;", msg),
        shiny::tags$p("Check your internet connection and click ", shiny::tags$strong("Try Again"), "."),
        shiny::tags$p(
          style = "background:#f0fdf4;border:1px solid #bbf7d0;border-radius:5px;padding:8px 12px;font-size:0.88em;color:#166534;margin-top:8px;",
          shiny::icon("circle-check", style="margin-right:4px;"),
          shiny::tags$strong("Note:"), " Uploading your own files and importing from GEO ",
          "still work normally - only the example data requires an internet connection. ",
          "Click ", shiny::tags$strong("Close"), " to proceed with another data source."
        )
      ),
      footer = shiny::tagList(
        shiny::actionButton(
          retry_btn_id,
          shiny::tagList(shiny::icon("rotate-right"), " Try Again"),
          class = "btn-primary"
        ),
        shiny::modalButton("Close")
      ),
      easyClose = TRUE
    ))
  }

  clean_name <- function(x) tools::file_path_sans_ext(base::basename(x))
  
  safe_read_table <- function(path) {

    ext <- tolower(tools::file_ext(path))

    # ── Excel ──────────────────────────────────────────────────────────────────
    if (ext %in% c("xls", "xlsx")) {
      df <- if (ext == "xls") readxl::read_xls(path) else readxl::read_xlsx(path)
      df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
      # Auto-detect gene column: if first column is non-numeric, use as rownames
      first_col_num <- suppressWarnings(
        all(!is.na(as.numeric(as.character(df[[1]]))))
      )
      if (!first_col_num) {
        rownames(df) <- as.character(df[[1]])
        df <- df[, -1, drop = FALSE]
      }
      # Coerce remaining columns to numeric
      df[] <- lapply(df, function(x) suppressWarnings(as.numeric(as.character(x))))
      return(df)
    }

    # ── Text / delimited ───────────────────────────────────────────────────────
    if (ext %in% c("txt", "csv", "tsv")) {

      # --- peek at raw bytes to detect separator ---
      first_line <- readLines(path, n = 1, warn = FALSE, encoding = "UTF-8")
      # Strip BOM if present
      first_line <- sub("^\xef\xbb\xbf", "", first_line)

      sep <- if (grepl("\t", first_line)) {
        "\t"
      } else if (grepl(";", first_line)) {
        ";"
      } else if (grepl(",", first_line)) {
        ","
      } else {
        " "
      }

      # --- read everything as character (maximally flexible) ---
      df <- utils::read.table(
        path,
        header           = TRUE,
        sep              = sep,
        check.names      = FALSE,
        stringsAsFactors = FALSE,
        comment.char     = "",
        quote            = "\"’",   # handle both quote styles
        fill             = TRUE,
        colClasses       = "character"
      )

      # --- strip any remaining surrounding quotes from all cells ---
      df[] <- lapply(df, function(x) gsub("^[\"’]|[\"’]$", "", trimws(x)))

      # --- remove fully-empty trailing columns ---
      df <- df[, colSums(!is.na(df) & df != "") > 0, drop = FALSE]

      if (ncol(df) == 0) stop("File appears to be empty after parsing.")

      # --- auto-detect gene name column ---
      # If the first column can’t be coerced to numeric it’s the gene column.
      # If ALL columns are numeric the row numbers serve as gene identifiers
      # (unusual, but we leave rownames as-is).
      first_col_num <- suppressWarnings(
        all(!is.na(as.numeric(gsub(",", ".", df[[1]]))))
      )
      if (!first_col_num) {
        rownames(df) <- df[[1]]
        df <- df[, -1, drop = FALSE]
      }

      if (ncol(df) == 0) stop("No data columns remain after removing the gene name column.")

      # --- normalize decimals (comma → dot) and coerce to numeric ---
      df[] <- lapply(df, function(x) suppressWarnings(as.numeric(gsub(",", ".", x))))

      if (nrow(df) == 0) stop("File has no data rows.")

      return(base::as.data.frame(df))
    }

    # ── RDS ────────────────────────────────────────────────────────────────────
    if (ext == "rds") {
      df <- readRDS(path)
      if (inherits(df, "matrix")) df <- as.data.frame(df, stringsAsFactors = FALSE)
      if (!inherits(df, "data.frame")) stop("Uploaded RDS must contain a data.frame or matrix.")
      if (is.null(rownames(df)) || nrow(df) == 0) stop("Uploaded RDS must have row names (gene names).")
      if (is.null(colnames(df)) || ncol(df) == 0) stop("Uploaded RDS must have column names (sample names).")
      # If any column is character, try to coerce (handles quoted numerics in RDS)
      df[] <- lapply(df, function(x) if (is.character(x)) suppressWarnings(as.numeric(x)) else x)
      return(df)
    }

    stop("Unsupported file format: .", ext,
         ". Supported formats: .txt, .csv, .tsv, .xls, .xlsx, .rds")
  }
  
  # ── Gene-set–specific reader ─────────────────────────────────────────────────
  # Unlike safe_read_table (designed for numeric expression matrices), this reader
  # preserves character data and handles the many formats gene sets come in:
  #   • RDS: character vector, named list, or data frame with gene names in col 1
  #   • RDA: same as RDS but stores an object by variable name (first object used)
  #   • CSV / TXT / TSV: gene names in first column; optional direction in second
  safe_read_geneset <- function(path) {
    ext <- tolower(tools::file_ext(path))

    if (ext == "rds") {
      obj <- readRDS(path)
      # Already the right type - return as-is for parse_geneset_object to handle
      return(obj)
    }

    if (ext == "rda") {
      e <- new.env(parent = emptyenv())
      load(path, envir = e)
      nms <- ls(e)
      if (length(nms) == 0) stop("RDA file is empty.")
      return(get(nms[1], envir = e))
    }

    if (ext %in% c("xls", "xlsx")) {
      df <- if (ext == "xls") readxl::read_xls(path) else readxl::read_xlsx(path)
      return(as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE))
    }

    if (ext %in% c("txt", "csv", "tsv")) {
      first_line <- readLines(path, n = 1, warn = FALSE, encoding = "UTF-8")
      first_line <- sub("^\xef\xbb\xbf", "", first_line)   # strip BOM
      sep <- if (grepl("\t", first_line)) "\t" else
             if (grepl(";",  first_line)) ";" else
             if (grepl(",",  first_line)) "," else " "
      df <- utils::read.table(
        path, header = TRUE, sep = sep, check.names = FALSE,
        stringsAsFactors = FALSE, comment.char = "", fill = TRUE,
        colClasses = "character"
      )
      df[] <- lapply(df, function(x) gsub("^[\"']|[\"']$", "", trimws(x)))
      df   <- df[, colSums(!is.na(df) & df != "") > 0, drop = FALSE]
      if (ncol(df) == 0) stop("File appears empty after parsing.")
      return(df)
    }

    stop("Unsupported gene set format: .", ext,
         ". Supported: .csv, .txt, .tsv, .rds, .rda, .xls, .xlsx")
  }

  # ── Metadata-specific reader: preserves character columns ───────────────────
  safe_read_meta <- function(path) {
    ext <- tolower(tools::file_ext(path))

    if (ext %in% c("xls", "xlsx")) {
      df <- if (ext == "xls") readxl::read_xls(path) else readxl::read_xlsx(path)
      return(as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE))
    }

    if (ext %in% c("txt", "csv", "tsv")) {
      first_line <- readLines(path, n = 1, warn = FALSE)
      first_line <- sub("^\xef\xbb\xbf", "", first_line)  # strip BOM
      sep <- if (grepl("\t", first_line)) "\t"
             else if (grepl(";", first_line)) ";"
             else if (grepl(",", first_line)) ","
             else " "
      df <- utils::read.table(
        path,
        header           = TRUE,
        sep              = sep,
        check.names      = FALSE,
        stringsAsFactors = FALSE,
        comment.char     = "",
        quote            = "\"'",
        fill             = TRUE
      )
      # Strip surrounding quotes
      df[] <- lapply(df, function(x) gsub("^[\"']|[\"']$", "", trimws(x)))
      # Remove fully-empty columns
      df   <- df[, colSums(!is.na(df) & df != "") > 0, drop = FALSE]
      if (nrow(df) == 0 || ncol(df) == 0) stop("Metadata file appears empty.")
      return(df)
    }

    if (ext == "rds") {
      df <- readRDS(path)
      if (inherits(df, "matrix")) df <- as.data.frame(df, stringsAsFactors = FALSE)
      if (!inherits(df, "data.frame")) stop("Uploaded RDS must contain a data.frame.")
      return(df)
    }

    stop("Unsupported metadata format: .", ext,
         ". Supported: .txt, .csv, .tsv, .xls, .xlsx, .rds")
  }

  # ── Gene set format rules ────────────────────────────────────────────────────
  # Accepted input formats (single gene set):
  #   • character vector       → gene names
  #   • data frame / matrix    → col 1 = gene names (character), col 2 = direction (optional)
  # Accepted input formats (multiple gene sets):
  #   • named list of the above
  # After parsing each element is either a character vector or a 2-col data frame
  # with colnames c("gene","direction").
  parse_geneset_object <- function(obj, name_hint) {
    # Matrix → coerce to data frame (keeps colnames/rownames)
    if (is.matrix(obj)) obj <- as.data.frame(obj, stringsAsFactors = FALSE)

    # Character vector → single gene set
    if (is.character(obj)) return(stats::setNames(list(obj), name_hint))

    # Data frame → single gene set
    if (is.data.frame(obj)) {
      if (ncol(obj) == 0) stop("Gene set data frame has no columns.")
      # If first column is numeric (e.g. expression matrix passed by mistake) flag it
      col1 <- obj[[1]]
      if (!is.character(col1) && !is.factor(col1)) {
        stop("The first column of the gene set must contain gene names (character), not numeric values. ",
             "Check your file format: col 1 = gene name, col 2 = direction (+1 / -1, optional).")
      }
      if (ncol(obj) == 1) {
        # Just gene names
        return(stats::setNames(list(as.character(col1)), name_hint))
      }
      # Take first two columns as gene + direction
      df <- obj[, 1:2, drop = FALSE]
      colnames(df) <- c("gene", "direction")
      df$gene      <- as.character(df$gene)
      df$direction <- suppressWarnings(as.numeric(df$direction))
      return(stats::setNames(list(df), name_hint))
    }

    # Named list → each element is its own gene set; recurse
    if (is.list(obj)) {
      nms <- names(obj)
      if (is.null(nms) || any(!nzchar(nms))) {
        # Give unnamed elements generic names
        nms <- ifelse(is.null(nms) | !nzchar(nms),
                      paste0(name_hint, "_", seq_along(obj)), nms)
      }
      parsed <- unlist(mapply(parse_geneset_object, obj, nms,
                               SIMPLIFY = FALSE, USE.NAMES = FALSE), recursive = FALSE)
      return(parsed)
    }

    stop("Invalid gene set format. Expected a character vector, data frame/matrix ",
         "(col 1 = gene names, col 2 = direction), or a named list of those.")
  }
  
  # ── Helper: reorder metadata rows to match expression column order ──────────
  # Called whenever meta_data() is set. If the SampleID column is a perfect
  # set-match to colnames(expr), rows are silently reordered to match.
  align_meta_to_expr <- function(meta, expr) {
    if (is.null(expr) || is.null(meta) || ncol(expr) == 0 || nrow(meta) == 0) return(meta)
    expr_cols <- colnames(expr)
    id_col    <- if ("SampleID" %in% colnames(meta)) "SampleID" else colnames(meta)[1]
    meta_ids  <- as.character(meta[[id_col]])
    if (setequal(expr_cols, meta_ids) && !identical(expr_cols, meta_ids)) {
      idx  <- match(expr_cols, meta_ids)
      meta <- meta[idx, , drop = FALSE]
      rownames(meta) <- expr_cols
      log_step("Metadata rows reordered to match expression column order.")
    }
    meta
  }

  ###### DATA ######

  # ---- Zenodo example download helpers (closures - capture server-scope reactives) ----
  # Each function manages its own modals: shows spinner, wraps download in tryCatch,
  # removes spinner on success, replaces it with an error modal (with retry button) on failure.
  # Called from both the main radio-button observer AND the retry observeEvent so the
  # same logic runs on first attempt and on any number of retries without code duplication.

  # R's utils::download.file() default timeout is only 60 seconds, which is
  # comfortably enough on the app's server (fast, close connection to
  # Zenodo) but can cut off the larger example files (the raw counts file is
  # ~70MB) partway through on a slower local connection - producing a
  # truncated download and a "Timeout was reached" warning rather than a
  # clean failure. Temporarily raising it just for the download call (and
  # restoring it on exit) avoids that without changing timeout behaviour
  # anywhere else in the app.
  .zenodo_download <- function(url, destfile, timeout_s = 300) {
    old_timeout <- getOption("timeout")
    on.exit(options(timeout = old_timeout), add = TRUE)
    options(timeout = max(timeout_s, old_timeout))
    utils::download.file(url = url, destfile = destfile, mode = "wb", quiet = TRUE)
  }

  .dl_expr_raw <- function() {
    show_loading_modal("Loading example data…",
                       "Downloading unprocessed counts from Zenodo.")
    dl_ok <- tryCatch({
      tmp <- tempfile(fileext = ".rds")
      status <- .zenodo_download(
        url      = "https://zenodo.org/records/18714122/files/counts.rds?download=1",
        destfile = tmp
      )
      if (status != 0L) stop("Download failed (status ", status, "). Check your internet connection.")
      e <- readRDS(tmp)
      expr_data(e)
      log_step(paste0("Expression data loaded: example (unprocessed) - ",
                      nrow(e), " genes × ", ncol(e), " samples."))
      TRUE
    }, error = function(err) {
      shiny::removeModal()
      show_dl_error_modal(conditionMessage(err), "retry_expr_raw")
      FALSE
    })
    if (isTRUE(dl_ok)) shiny::removeModal()
  }

  .dl_expr_proc <- function() {
    show_loading_modal("Loading example data…",
                       "Downloading processed counts from Zenodo.")
    dl_ok <- tryCatch({
      tmp <- tempfile(fileext = ".rds")
      status <- .zenodo_download(
        url      = "https://zenodo.org/records/18714122/files/corrcounts.rds?download=1",
        destfile = tmp
      )
      if (status != 0L) stop("Download failed (status ", status, "). Check your internet connection.")
      e <- readRDS(tmp)
      expr_data(e)
      log_step(paste0("Expression data loaded: example (processed) - ",
                      nrow(e), " genes × ", ncol(e), " samples."))
      TRUE
    }, error = function(err) {
      shiny::removeModal()
      show_dl_error_modal(conditionMessage(err), "retry_expr_proc")
      FALSE
    })
    if (isTRUE(dl_ok)) shiny::removeModal()
  }

  .dl_meta_example <- function() {
    show_loading_modal("Loading example metadata…",
                       "Downloading from Zenodo.")
    dl_ok <- tryCatch({
      tmp <- tempfile(fileext = ".rds")
      status <- .zenodo_download(
        url      = "https://zenodo.org/records/18714122/files/metadata.rds?download=1",
        destfile = tmp
      )
      if (status != 0L) stop("Download failed (status ", status, "). Check your internet connection.")
      m <- readRDS(tmp)
      m <- align_meta_to_expr(m, expr_data())
      meta_data(m)
      log_step(paste0("Metadata loaded: example - ",
                      nrow(m), " samples × ", ncol(m), " variables."))
      TRUE
    }, error = function(err) {
      shiny::removeModal()
      show_dl_error_modal(conditionMessage(err), "retry_meta_example")
      FALSE
    })
    if (isTRUE(dl_ok)) shiny::removeModal()
  }

  .dl_gs_example <- function() {
    show_loading_modal("Loading example gene sets…",
                       "Downloading from Zenodo.")
    dl_ok <- tryCatch({
      tmp <- tempfile(fileext = ".rds")
      status <- .zenodo_download(
        url      = "https://zenodo.org/records/18714122/files/SenescenceGeneSets.rds?download=1",
        destfile = tmp
      )
      if (status != 0L) stop("Download failed (status ", status, "). Check your internet connection.")
      gs <- readRDS(tmp)
      gene_sets(gs)
      log_step(paste0("Gene sets loaded: example - ", length(gs), " set(s), ",
                      sum(sapply(gs, function(x) if (is.data.frame(x)) nrow(x) else length(x))),
                      " genes total."))
      TRUE
    }, error = function(err) {
      shiny::removeModal()
      show_dl_error_modal(conditionMessage(err), "retry_gs_example")
      FALSE
    })
    if (isTRUE(dl_ok)) shiny::removeModal()
  }

  ####### > Example Data #######

  observeEvent(input$expr_source, {
    src_label <- switch(input$expr_source,
      example_raw  = "example (unprocessed)",
      example_proc = "example (processed)",
      upload       = "user upload",
      geo          = "GEO import",
      input$expr_source
    )
    log_step(paste0("Expression source changed to: ", src_label, "."))
    if (input$expr_source == "example_raw")  .dl_expr_raw()
    if (input$expr_source == "example_proc") .dl_expr_proc()
    # Auto-switch metadata to the example set when example expression data is selected
    if (input$expr_source %in% c("example_raw", "example_proc")) {
      shiny::updateRadioButtons(session, "meta_source", selected = "example")
    }
  })

  # "Try Again" retry handlers - call the same helper function directly,
  # no need to toggle the radio button (avoids unreliable updateRadioButtons tricks).
  observeEvent(input$retry_expr_raw,  { .dl_expr_raw()  }, ignoreInit = TRUE)
  observeEvent(input$retry_expr_proc, { .dl_expr_proc() }, ignoreInit = TRUE)
  
  ###### > User Input Data ######
  
  observeEvent(input$expr_file, {
    req(input$expr_file)
    show_loading_modal(paste0("Reading ", input$expr_file$name, "…"),
                       "Parsing expression matrix.")
    on.exit(shiny::removeModal(), add = TRUE)
    result <- tryCatch(
      safe_read_table(input$expr_file$datapath),
      error = function(err) {
        shiny::showNotification(
          paste0("Could not parse file: ", conditionMessage(err)),
          type = "error", duration = 10
        )
        NULL
      }
    )
    if (is.null(result)) return()

    # Check data quality: flag if >20 % of values are NA (coercion failures)
    # or if any column is non-numeric (should not happen after safe_read_table,
    # but catches RDS files with character columns).
    na_frac      <- mean(is.na(unlist(result)))
    has_non_num  <- !all(sapply(result, is.numeric))
    warn_msg <- NULL
    if (has_non_num) {
      warn_msg <- "Some columns are non-numeric - this file may not be an expression matrix. Check that genes are rows and samples are columns."
    } else if (na_frac > 0.20) {
      warn_msg <- sprintf(
        "%.0f%% of values are missing after parsing - the file may not be a numeric expression matrix, or the separator / decimal character was misdetected.",
        na_frac * 100
      )
    }
    expr_quality_warn(warn_msg)

    expr_data(result)
    module_meta_active(FALSE)
    log_step(paste0("Expression data loaded: uploaded file '", input$expr_file$name, "' - ",
                    nrow(result), " genes × ", ncol(result), " samples."))
  })
  
  ###### > GEO data - delegated to geoImportModule ######

  # ── Instantiate the module ─────────────────────────────────────────────────
  geo_module <- geoImportServer("geo_import", log_fn = log_step)

  # ── Wire module outputs into the app's reactive storage ───────────────────
  observe({
    req(input$expr_source == "geo")

    mod_expr <- geo_module$expr_data()
    mod_meta <- geo_module$meta_data()

    if (!is.null(mod_expr)) {
      e <- as.data.frame(mod_expr)
      expr_data(e)
      log_step(paste0("Expression data loaded: GEO import - ",
                      nrow(e), " genes × ", ncol(e), " samples."))
    }

    if (!is.null(mod_meta)) {
      # Set flag FIRST so the app-level column-selector observers skip overwriting
      module_meta_active(TRUE)
      full_geo_meta(mod_meta)
      meta_data(mod_meta)
      log_step(paste0("Metadata loaded: GEO - ",
                      nrow(mod_meta), " samples × ", ncol(mod_meta), " variables."))
      # Automatically switch metadata source to geo so column selectors appear
      updateRadioButtons(session, "meta_source", selected = "geo")
    }
  })

  # ── Legacy helpers kept for the non-GEO upload path ─────────────────────
  
  # [read_geo_supp_matrix replaced by geo_import_module.R::parse_expression_file]
  read_geo_supp_matrix_LEGACY <- function(path, original_name = NULL, geo_accession = NULL) {
    fname <- if (!is.null(original_name)) original_name else basename(path)
    ext <- tolower(tools::file_ext(fname))
    
    #showNotification(paste("Processing file:", fname), type = "message", duration = 2)
    
    # --- Handle gz ---
    if (ext == "gz" || grepl("\\.gz$", fname, ignore.case = TRUE)) {
      tmp <- tempfile(fileext = sub("\\.gz$", "", fname))
      R.utils::gunzip(path, destname = tmp, overwrite = TRUE, remove = FALSE)
      path <- tmp
      ext <- tolower(tools::file_ext(tmp))
    }
    
    # --- Handle TAR archives ---
    if (ext == "tar") {
      extract_dir <- tempfile()
      dir.create(extract_dir)
      utils::untar(path, exdir = extract_dir)
      files <- list.files(extract_dir, full.names = TRUE, recursive = TRUE)
      
      # detect gz gene count files
      gz_files <- files[grepl("_geneCOUNT\\.txt\\.gz$", files, ignore.case = TRUE)]
      
      if (length(gz_files) == 0) {
        geo_link <- if (!is.null(geo_accession)) paste0(
          "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", geo_accession
        ) else "#"
        showNotification(
          HTML(paste0(
            "⚠️ TAR archive has no recognized gene count files.<br>",
            "Check GEO entry: <a href='", geo_link, "' target='_blank'>", geo_accession, "</a>"
          )),
          type = "warning", duration = 10
        )
        updateSelectInput(session, "geo_selected_supp", selected = "")
        return(NULL)
      }
      
      # Warn about nonstandard files
      nonstandard_files <- gz_files[!grepl("_S\\d+_geneCOUNT", basename(gz_files))]
      if (length(nonstandard_files) > 0) {
        geo_link <- if (!is.null(geo_accession)) paste0(
          "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", geo_accession
        ) else "#"
        showNotification(
          HTML(paste0(
            "⚠️ Some files in TAR archive do not match expected RNA-seq naming.<br>",
            "Check GEO entry: <a href='", geo_link, "' target='_blank'>", geo_accession, "</a>"
          )),
          type = "warning", duration = NULL
        )
        updateSelectInput(session, "geo_selected_supp", selected = "")
      }
      
      # --- Load all gene count files ---
      count_list <- list()
      for (f in gz_files) {
        tmp_txt <- tempfile(fileext = ".txt")
        R.utils::gunzip(f, destname = tmp_txt, overwrite = TRUE, remove = FALSE)
        df <- read.delim(tmp_txt, header = TRUE, stringsAsFactors = FALSE)
        gene_ids <- df[[1]]
        counts   <- as.numeric(df[[2]])
        sample_id <- sub("_.*", "", basename(f))
        count_list[[sample_id]] <- data.frame(Gene = gene_ids, Count = counts, stringsAsFactors = FALSE)
      }
      
      # Rename counts column to sample names
      count_list <- Map(function(df, sn) { colnames(df)[2] <- sn; df }, count_list, names(count_list))
      merged <- Reduce(function(x, y) merge(x, y, by = "Gene", all = FALSE), count_list)
      rownames(merged) <- merged$Gene
      merged$Gene <- NULL
      expr_matrix <- as.matrix(merged)
      storage.mode(expr_matrix) <- "numeric"
      
      if (!is.numeric(expr_matrix) || nrow(expr_matrix) == 0 || ncol(expr_matrix) == 0) {
        geo_link <- if (!is.null(geo_accession)) paste0(
          "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", geo_accession
        ) else "#"
        showNotification(
          HTML(paste0(
            "⚠️ Matrix is empty or non-numeric.<br>",
            "Check GEO entry: <a href='", geo_link, "' target='_blank'>", geo_accession, "</a>"
          )),
          type = "warning", duration = NULL
        )
        return()
      }
      return(expr_matrix)
    }
    
    # --- Read regular files ---
    df <- switch(ext,
                 "xls"  = readxl::read_xls(path),
                 "xlsx" = readxl::read_xlsx(path),
                 "csv"  = utils::read.csv(path, check.names = FALSE),
                 "tsv"  = utils::read.delim(path, check.names = FALSE),
                 "txt"  = utils::read.delim(path, check.names = FALSE),
                 stop("Unsupported supplementary file format: ", ext)
    )
    
    # First column as rownames if needed
    if (!all(sapply(df, is.numeric))) {
      rownames(df) <- df[[1]]
      df <- df[, -1, drop = FALSE]
    }
    
    mat <- as.matrix(df)
    if (!is.numeric(mat)) {
      geo_link <- if (!is.null(geo_accession)) paste0(
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", geo_accession
      ) else "#"
      showNotification(
        HTML(paste0(
          "⚠️ Extracted file not numeric.<br>",
          "Check GEO entry: <a href='", geo_link, "' target='_blank'>", geo_accession, "</a>"
        )),
        type = "warning", duration = NULL
      )
      return()
    }
    
    return(mat)
  }
  
  # load_supplementary <- function(file_path, geo_accession = NULL) {
  #   tryCatch({
  #     mat <- read_geo_supp_matrix(file_path, basename(file_path), geo_accession)
  #     if (!is.null(mat) && is.numeric(mat) && nrow(mat) > 0 && ncol(mat) > 0) return(mat)
  #     stop("File could not be parsed as numeric matrix")
  #   }, error = function(e) {
  #     showNotification(paste("Supplementary file failed:", e$message),
  #                      type = "warning", duration = 8)
  #     return(NULL)
  #   })
  # }
  
  # --- Extraction helpers for archives (legacy - superseded by geo_import_module) ---
  extract_files <- function(file_path, depth = 0, max_depth = 5) {
    
    if (depth > max_depth) return(character(0))
    if (!file.exists(file_path)) return(character(0))
    
    ext <- tolower(tools::file_ext(file_path))
    
    # Create extraction directory
    tmp_dir <- tempfile(pattern = "geo_extract_")
    dir.create(tmp_dir)
    
    extracted_files <- character(0)
    
    if (ext == "gz") {
      
      # If .tar.gz, untar directly
      if (grepl("\\.tar\\.gz$", file_path, ignore.case = TRUE)) {
        utils::untar(file_path, exdir = tmp_dir)
        extracted_files <- list.files(tmp_dir, full.names = TRUE, recursive = TRUE)
      } else {
        # plain .gz
        out_file <- file.path(tmp_dir, tools::file_path_sans_ext(basename(file_path)))
        R.utils::gunzip(file_path, destname = out_file, overwrite = TRUE, remove = FALSE)
        extracted_files <- out_file
      }
      
    } else if (ext == "tar") {
      
      utils::untar(file_path, exdir = tmp_dir)
      extracted_files <- list.files(tmp_dir, full.names = TRUE, recursive = TRUE)
      
    } else if (ext == "zip") {
      
      utils::unzip(file_path, exdir = tmp_dir)
      extracted_files <- list.files(tmp_dir, full.names = TRUE, recursive = TRUE)
      
    } else {
      
      # Not archive → return directly
      extracted_files <- file_path
    }
    
    # --- Recursive extraction of nested archives ---
    nested_archives <- extracted_files[
      grepl("\\.(gz|tar|zip)$", extracted_files, ignore.case = TRUE)
    ]
    
    if (length(nested_archives) > 0) {
      nested_results <- unlist(
        lapply(nested_archives, extract_files, depth = depth + 1)
      )
      extracted_files <- c(extracted_files, nested_results)
    }
    
    # --- Keep only table-like files ---
    extracted_files[
      grepl("\\.(txt|tsv|csv|xls|xlsx)$",
            extracted_files,
            ignore.case = TRUE)
    ]
  }
  
  detect_expr_candidates <- function(file_paths) {
    all_files <- unlist(lapply(file_paths, extract_files))
    candidates <- data.frame(
      file = all_files,
      prob_expr = 0,
      split_prob = 0,
      stringsAsFactors = FALSE
    )
    
    for (i in seq_along(all_files)) {
      f <- all_files[i]
      df <- tryCatch(read.delim(f, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
      if (is.null(df)) next
      
      df <- df[, colSums(is.na(df)) < nrow(df), drop = FALSE]
      numeric_fraction <- sapply(df, function(col) mean(suppressWarnings(!is.na(as.numeric(col)))))
      
      # Probability of being expression data
      candidates$prob_expr[i] <- min(1, max(numeric_fraction) * ncol(df) / 10)
      
      # Probability that file is single-sample (split)
      candidates$split_prob[i] <- ifelse(length(which(numeric_fraction > 0.8)) == 1, 0.9, 0.1)
    }
    
    # Filter low probability
    candidates <- candidates[candidates$prob_expr > 0.1, ]
    candidates <- candidates[order(-candidates$prob_expr), ]
    candidates
  }
  
  # --- Read single candidate expression file ---
  read_expr_candidate <- function(file, split_prob = 0.1) {
    df <- tryCatch(read.delim(file, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(df)) return(NULL)
    
    # Remove fully NA columns
    df <- df[, colSums(is.na(df)) < nrow(df), drop = FALSE]
    
    # Identify numeric vs gene columns
    numeric_fraction <- sapply(df, function(col) mean(suppressWarnings(!is.na(as.numeric(col)))))
    numeric_cols <- which(numeric_fraction > 0.8)
    gene_col <- which(numeric_fraction < 0.2)[1]
    
    if (is.na(gene_col) || length(numeric_cols) < 1) return(NULL)
    
    # Extract genes
    genes <- as.character(df[[gene_col]])
    
    # Subset numeric expression data
    expr <- as.matrix(df[, numeric_cols, drop = FALSE])
    storage.mode(expr) <- "numeric"
    
    # Set rownames as genes
    rownames(expr) <- genes
    
    # --- Determine sample names ---
    # Try to detect GSM code in file name
    fname <- basename(file)
    gsm_match <- regmatches(fname, regexpr("GSM[0-9]+", fname))
    
    if (length(gsm_match) == 1 && ncol(expr) == 1) {
      sample_name <- gsm_match
    } else if (length(gsm_match) >= 1) {
      sample_name <- gsm_match
      # If more columns than matches, append _1, _2...
      if (length(sample_name) < ncol(expr)) {
        sample_name <- paste0(sample_name[1], "_", seq_len(ncol(expr)))
      }
    } else {
      # fallback: use file name + _1, _2 ...
      sample_name <- if (ncol(expr) == 1) fname else paste0(fname, "_", seq_len(ncol(expr)))
    }
    
    colnames(expr) <- sample_name
    
    expr
  }
  
  # --- Combine multiple split single-sample files ---
  combine_split_files <- function(files) {
    expr_list <- lapply(files, read_expr_candidate)
    expr_list <- expr_list[!sapply(expr_list, is.null)]
    if (length(expr_list) == 0) return(NULL)
    
    # Align genes across files
    common_genes <- Reduce(intersect, lapply(expr_list, rownames))
    expr_list <- lapply(expr_list, function(x) x[common_genes, , drop = FALSE])
    
    # Column-bind
    expr <- do.call(cbind, expr_list)
    
    expr
  }
  
  
  
  # --- Deterministic / LLM Rescue ---
  try_deterministic_rescue <- function(file_paths) {
    
    showNotification("Attempting deterministic rescue...",
                     type = "warning", duration = 5)
    
    # Step 1: Extract all table-like files
    all_files <- unlist(lapply(file_paths, extract_files))
    if (length(all_files) == 0)
      return(NULL)
    
    # Step 2: Detect candidate expression files
    candidates <- detect_expr_candidates(all_files)
    if (is.null(candidates) || nrow(candidates) == 0)
      return(NULL)
    
    split_files <- candidates$file[candidates$split_prob > 0.5]
    full_files  <- candidates$file[candidates$split_prob <= 0.5]
    
    matrices <- list()
    
    # -----------------------------
    # Step 3a: Process full matrices
    # -----------------------------
    for (f in full_files) {
      
      df <- tryCatch(
        read.delim(f, check.names = FALSE, stringsAsFactors = FALSE),
        error = function(e) NULL
      )
      if (is.null(df)) next
      
      df <- df[, colSums(is.na(df)) < nrow(df), drop = FALSE]
      
      numeric_fraction <- sapply(df, function(col)
        mean(suppressWarnings(!is.na(as.numeric(col)))))
      
      numeric_cols <- which(numeric_fraction > 0.8)
      gene_col     <- which(numeric_fraction < 0.2)[1]
      
      if (!is.na(gene_col) && length(numeric_cols) >= 2) {
        
        genes <- trimws(as.character(df[[gene_col]]))
        valid <- !is.na(genes) & genes != ""
        
        expr <- as.matrix(df[valid, numeric_cols, drop = FALSE])
        storage.mode(expr) <- "numeric"
        
        genes <- genes[valid]
        
        # remove duplicates
        dup <- duplicated(genes)
        if (any(dup)) {
          genes <- genes[!dup]
          expr  <- expr[!dup, , drop = FALSE]
        }
        
        rownames(expr) <- genes
        
        if (nrow(expr) > 1 && ncol(expr) > 1)
          matrices[[basename(f)]] <- expr
      }
      
      # Transpose fallback
      if (ncol(df) > nrow(df)) {
        
        df_t <- as.data.frame(t(df))
        
        numeric_fraction <- sapply(df_t, function(col)
          mean(suppressWarnings(!is.na(as.numeric(col)))))
        
        numeric_cols <- which(numeric_fraction > 0.8)
        gene_col     <- which(numeric_fraction < 0.2)[1]
        
        if (!is.na(gene_col) && length(numeric_cols) >= 2) {
          
          genes <- trimws(as.character(df_t[[gene_col]]))
          valid <- !is.na(genes) & genes != ""
          
          expr <- as.matrix(df_t[valid, numeric_cols, drop = FALSE])
          storage.mode(expr) <- "numeric"
          
          genes <- genes[valid]
          
          dup <- duplicated(genes)
          if (any(dup)) {
            genes <- genes[!dup]
            expr  <- expr[!dup, , drop = FALSE]
          }
          
          rownames(expr) <- genes
          
          if (nrow(expr) > 1 && ncol(expr) > 1)
            matrices[[paste0(basename(f), "_transposed")]] <- expr
        }
      }
    }
    
    # -----------------------------
    # Step 3b: Combine split files
    # -----------------------------
    if (length(split_files) > 1) {
      
      combined <- combine_split_files(split_files)
      
      if (!is.null(combined) &&
          nrow(combined) > 1 &&
          ncol(combined) > 1) {
        
        matrices[["Split_combined"]] <- combined
      }
    }
    
    if (length(matrices) == 0)
      return(NULL)
    
    return(matrices)
  }
  
  load_llm_fallback <- function(file_paths) {
    showNotification("Attempting LLM-based parsing...", type = "warning", duration = 5)
    files <- unlist(lapply(file_paths, extract_files))
    candidate_summaries <- list()
    for (f in files) {
      df <- tryCatch(read.delim(f, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
      if (is.null(df) || nrow(df) < 2 || ncol(df) < 2) next
      df <- df[, colSums(!is.na(df)) > 0, drop = FALSE]
      summary_text <- paste0("File: ", basename(f), " | Rows: ", nrow(df), " | Columns: ", ncol(df),
                             " | First 3 rows:\n", paste(capture.output(print(head(df, 3))), collapse = "\n"))
      candidate_summaries[[f]] <- summary_text
    }
    if (length(candidate_summaries) == 0) return(NULL)
    prompt <- paste(
      "Parse gene expression matrices from the following files:\n\n",
      paste(candidate_summaries, collapse = "\n\n"),
      "\nReturn only numeric matrices with genes as rows and samples as columns."
    )
    parsed <- llm_call_to_parse_table(NULL, prompt)
    if (!is.null(parsed) && is.numeric(parsed) && nrow(parsed) > 0 && ncol(parsed) > 0) return(parsed)
    NULL
  }
  
  load_auto_deterministic <- function(files_path) { try_deterministic_rescue(files_path) }
  
 
  # NOTE: GEO fetch/parse/rescue logic has been replaced by geoImportServer("geo_import").
  # The block below is DEAD CODE - kept only for reference. It is never called.
  .LEGACY_observeEvent_geo_fetch <- function() {
  observeEvent(input$geo_fetch, {
    req(input$geo_accession)
    
    withProgress(message = paste("Fetching GEO accession", input$geo_accession), value = 0, {
      
      tmp_dir <- tempdir()
      
      incProgress(0.3, "Downloading GEO objects...")
      objs <- GEOquery::getGEO(input$geo_accession, GSEMatrix = TRUE)
      geo_objects(objs)
      
      # ---- CHECK IF ANY OBJECT HAS EXPRESSION 
      has_expr <- FALSE
      for (obj in objs) {
        expr <- tryCatch(Biobase::exprs(obj), error = function(e) NULL)
        if (!is.null(expr) && nrow(expr) > 0 && ncol(expr) > 0) {
          has_expr <- TRUE
          break
        }
      }
      
      incProgress(0.3, "Fetching supplementary files...")
      supp <- GEOquery::getGEOSuppFiles(
        input$geo_accession,
        makeDirectory = FALSE,
        baseDir = tmp_dir
      )
      geo_supp_files(supp)
      
      incProgress(0.4, "Finalizing...")
      
      if (!has_expr) {
        showNotification(
          "No expression matrix found in GEO object. Please select a supplementary file.",
          type = "warning",
          duration = 12
        )
      } else {
        showNotification("GEO accession loaded.", type = "default")
      }
    })
  })
  
  
  # --- 2. GEO object selection with descriptive labels -----
  output$geo_object_selector <- renderUI({
    req(geo_objects())
    objs <- geo_objects()
    
    obj_labels <- sapply(objs, function(obj) {
      title <- tryCatch(Biobase::experimentData(obj)@title, error = function(e) "")
      platform <- tryCatch(Biobase::annotation(obj), error = function(e) "")
      paste0(title, " [", platform, "]")
    })
    
    selectInput(
      "geo_selected_object",
      "Select GEO object (platform):",
      choices = setNames(seq_along(objs), obj_labels)
    )
  })
  
  # --- 3. Load expression & metadata, handle empty expression matrices -----
  observeEvent(input$geo_selected_object, {
    
    
    # Reset supplementary selection whenever object changes
    updateSelectInput(
      session,
      "geo_selected_supp",
      selected = ""
    )
    
    req(geo_objects())
    idx <- as.numeric(input$geo_selected_object)
    obj <- geo_objects()[[idx]]
    
    shiny::withProgress(
      message = "Loading expression and metadata...",
      value = 0,
      {
        # Metadata
        incProgress(0.3, detail = "Extracting sample metadata...")
        Sys.sleep(0.1)
        meta <- Biobase::pData(obj)
        full_geo_meta(meta)
        meta_data(cbind(SampleID = rownames(meta),
                        meta[, setdiff(colnames(meta), "SampleID"), drop = FALSE]))
        
        # Expression
        incProgress(0.3, detail = "Extracting expression matrix...")
        Sys.sleep(0.1)
        expr <- Biobase::exprs(obj)
        
        if (nrow(expr) == 0) {
          # Expression matrix empty -> show supplementary file dropdown
          # showNotification(
          #   "Expression matrix is empty! Please select a supplementary file.",
          #   type = "warning",
          #   duration = 8
          # )
          expr_data(NULL)
        } else {
          expr_data(expr)
          geo_supp_files(NULL)  # no need for supplementary selection
        }
        
        incProgress(0.4, detail = "Finalizing...")
        Sys.sleep(0.1)
        
        # Automatically switch metadata source to GEO
        updateRadioButtons(
          session,
          inputId = "meta_source",
          selected = "geo"
        )
        
        showNotification(
          "Metadata loaded successfully.",
          type = "default",
          duration = 4
        )
      }
    )
  })
  
  # --- 4. Render supplementary file selector if expression matrix is empty -----
  output$geo_supp_selector <- renderUI({
    supp <- geo_supp_files()
    if (is.null(supp) || nrow(supp) == 0) return(NULL)  # ensure files exist
    
    files <- rownames(supp)
    
    tagList(
      # --- Help text ---

      # --- File selector ---
      selectInput(
        "geo_selected_supp",
        "Select a supplementary file containing expression data:",
        choices = c("- Please choose a file -" = "",
                    setNames(files, basename(files))),
        selected = ""
      ),
      helpText(
        "Sometimes the GEO data doesn't appear in a standard format, and, as such, the usual workflow for reading supplementary files may fail. ",
        "If it fails, please try the automatic detector, where we attempt different approaches (including LLM-based parsing) ",
        "to infer the file organization and check if any files contain gene expression data."
      ),
      
      # --- Automatic detection button ---
      actionButton("try_auto_supp", "Try automatic detection/rescue")
    )
  })
  
  # --- 5. Load selected supplementary file -----
  observeEvent(input$geo_selected_supp, {
    req(input$geo_selected_supp != "")
    req(input$geo_selected_supp)
    
    file_path <- input$geo_selected_supp  # full path
    expr <- NULL
    
    withProgress(message = "Processing supplementary file...", value = 0, {
      
      incProgress(0.2, detail = "Reading file...")
      Sys.sleep(0.3)  # fake delay for user experience
      
      tryCatch({
        
        # --- read gzipped or plain text files ---
        expr <- read_geo_supp_matrix(
          path = file_path,
          original_name = basename(file_path),
          geo_accession = input$geo_accession
        )
        
        incProgress(0.3, detail = "Validating matrix...")
        Sys.sleep(0.3)
        
        # --- check numeric matrix ---
        if (!is.numeric(expr) || nrow(expr) == 0 || ncol(expr) == 0) {
          
          # Reset dropdown selection
          updateSelectInput(
            session,
            "geo_selected_supp",
            selected = ""
          )
          
          stop("Matrix not numeric or empty")
        }
        
        
        incProgress(0.3, detail = "Finalizing matrix...")
        Sys.sleep(0.3)
        
        if (is.null(expr)) return()
        expr_data(expr)
        # Auto-align metadata to the supplementary file's sample columns
        meta <- meta_data()
        if (!is.null(meta) && nrow(meta) > 0L) {
          expr_cols <- colnames(expr)
          id_col    <- if ("SampleID" %in% colnames(meta)) "SampleID" else colnames(meta)[1]
          meta_ids  <- as.character(meta[[id_col]])
          common    <- intersect(expr_cols, meta_ids)
          if (length(common) == length(expr_cols)) {
            # Perfect match (possibly reordered) - just reorder
            meta_data(align_meta_to_expr(meta, expr))
            showNotification("Expression data loaded from supplementary file.", type = "default")
          } else if (length(common) > 0L) {
            # Partial match - restrict to common samples automatically
            idx       <- match(common, meta_ids)
            meta_sub  <- meta[idx, , drop = FALSE]
            # Also restrict expression to matched samples
            expr_data(expr[, common, drop = FALSE])
            meta_data(align_meta_to_expr(meta_sub, expr[, common, drop = FALSE]))
            showNotification(
              shiny::HTML(paste0("Expression data loaded. ",
                "<b>", length(common), "/", ncol(expr), " samples</b> matched to metadata - ",
                "restricted automatically. Check the Data tab if needed.")),
              type = "warning", duration = 10)
          } else {
            showNotification(
              shiny::HTML(paste0("Expression data loaded but <b>no sample names matched</b> metadata. ",
                "Use the sample management controls in the sidebar to align them.")),
              type = "warning", duration = 12)
          }
        } else {
          showNotification("Expression data loaded from supplementary file.", type = "default")
        }
        
      }, error = function(e) {
        expr_data(NULL)
        
        # Reset dropdown selection
        updateSelectInput(
          session,
          "geo_selected_supp",
          selected = ""
        )
        
        showNotification(
          paste("Cannot load selected supplementary file:", e$message),
          type = "error",
          duration = 10
        )
      })
      
      incProgress(0.2, detail = "Done.")
      Sys.sleep(0.2)
    })
  })
    
  
  # --- 6. Automatic deterministic / LLM rescue -----
  
  # --- 0. Initialize reactive values ---
  expr_candidates_available <- reactiveVal(FALSE)
  geo_candidates <- reactiveVal(NULL)
  
  # --- 1. Automatic deterministic / LLM rescue ---
  observeEvent(input$try_auto_supp, {
    
    supp <- geo_supp_files()
    req(supp)
    
    files_path <- rownames(supp)
    
    # Reset reactive flags
    expr_candidates_available(FALSE)
    geo_candidates(NULL)
    
    # ----------------------------------
    # Step 1: Deterministic rescue
    # ----------------------------------
    matrices <- try_deterministic_rescue(files_path)
 
        # ----------------------------------
    # Step 2: LLM fallback (only if needed)
    # ----------------------------------
    if (is.null(matrices)) {
      
      matrices <- load_llm_fallback(files_path)
      
      if (is.null(matrices)) {
        
        showNotification(
          HTML(paste0(
            "<b>⚠️ Could not convert GEO data into matrix.</b><br>",
            "Check your GEO entry: <a href='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",
            input$geo_accession,
            "' target='_blank'>",
            input$geo_accession,
            "</a>.<br>",
            "If the data exists in genes × samples format, format manually and upload."
          )),
          type = "error",
          duration = NULL
        )
        
        return()
      }
    }
    
    # ----------------------------------
    # Step 3: Operate at MATRIX level
    # ----------------------------------
    if (length(matrices) == 1) {
      
      expr_data(matrices[[1]])
      
      showNotification(
        "Expression matrix loaded automatically.",
        type = "default"
      )
      
    } else {
      
      # Multiple valid processed matrices
      geo_candidates(matrices)
      expr_candidates_available(TRUE)
      
      showNotification(
        "Multiple expression matrices detected. Please select one.",
        type = "message"
      )
    }
  })
  
  # --- 2. Candidate dropdown after automatic detection ---
  output$geo_candidate_selector <- renderUI({
    # Only build UI if multiple candidates are available
    if (!expr_candidates_available()) return(NULL)
    
    candidates <- geo_candidates()
    req(candidates)
    
    # Build choices with probability display
    choices <- paste0(
      basename(candidates$file),
      " (prob expr: ", round(candidates$prob_expr, 2),
      ", split: ", round(candidates$split_prob, 2), ")"
    )
    names(choices) <- candidates$file
    
    tagList(
      p("Multiple candidate expression files detected. You may choose one to load manually:"),
      selectInput("geo_selected_candidate",
                  "Candidate expression files:",
                  choices = choices,
                  selected = NULL)
    )
  })
  
  # --- 3. Load selected candidate file ---
  observeEvent(input$geo_selected_candidate, {
    req(input$geo_selected_candidate)
    f <- input$geo_selected_candidate
    expr <- read_expr_candidate(f)
    
    if (!is.null(expr)) {
      expr_data(expr)
      showNotification("Expression matrix loaded from selected candidate file.", type = "default")
    } else {
      showNotification("Selected file does not appear to be a valid expression matrix.", type = "error")
    }
  }) } # end .LEGACY_observeEvent_geo_fetch

  ##### METADATA  ######################################### 
  
  
  observeEvent(input$meta_source, {
    src_label <- switch(input$meta_source,
      example = "example",
      upload  = "user upload",
      geo     = "GEO import",
      input$meta_source
    )
    log_step(paste0("Metadata source changed to: ", src_label, "."))
    if (input$meta_source == "example") .dl_meta_example()
  })

  observeEvent(input$retry_meta_example, { .dl_meta_example() }, ignoreInit = TRUE)
  
  
  
  observeEvent(input$meta_file, {
    req(input$meta_file)
    show_loading_modal(paste0("Reading ", input$meta_file$name, "…"),
                       "Parsing metadata file.")
    on.exit(shiny::removeModal(), add = TRUE)
    raw <- safe_read_meta(input$meta_file$datapath)
    upload_meta_raw(raw)
    module_meta_active(FALSE)
    meta_data(align_meta_to_expr(raw, expr_data()))
    # Log basic info + match summary against expression data
    expr <- expr_data()
    match_info <- if (!is.null(expr)) {
      ids     <- as.character(raw[[colnames(raw)[1]]])
      n_match <- length(intersect(ids, colnames(expr)))
      paste0(", ", n_match, "/", ncol(expr), " samples matched to expression columns")
    } else ""
    log_step(paste0("Metadata loaded: '", input$meta_file$name, "' - ",
                    nrow(raw), " samples × ", ncol(raw), " variables", match_info, "."))
  })
  
  
  
  # --- GEO metadata handling with SampleID mapping and subsetting ---
  observe({
    req(full_geo_meta(), input$meta_source)

    if (input$meta_source == "geo") {
      df <- full_geo_meta()
      
      # Determine default SampleID column.
      # When the geo_import_module provided the data it already built a "SampleID"
      # column - prefer that unconditionally. For the legacy pData path, fall back
      # to geo_accession (the canonical GEO sample identifier).
      default_sampleid <- if (isTRUE(module_meta_active()) && "SampleID" %in% colnames(df)) {
        "SampleID"
      } else if ("geo_accession" %in% colnames(df)) {
        "geo_accession"
      } else if ("SampleID" %in% colnames(df)) {
        "SampleID"
      } else {
        colnames(df)[1]
      }
      
      # Render UI for SampleID column selection
      output$geo_sampleid_selector <- renderUI({
        selectInput(
          inputId = "geo_sampleid_col",
          label   = "Select column containing Sample IDs:",
          choices = colnames(df),
          selected = default_sampleid
        )
      })
      
    }
  })

  # --- Update meta_data() based on SampleID column (legacy non-module GEO path) ---
  observe({
    req(full_geo_meta(), input$geo_sampleid_col, input$meta_source)
    if (isTRUE(module_meta_active())) return(NULL)

    if (input$meta_source == "geo") {
      df <- full_geo_meta()

      sample_ids <- if (input$geo_sampleid_col %in% colnames(df)) {
        df[[input$geo_sampleid_col]]
      } else {
        rownames(df)
      }

      # Keep all columns except the chosen ID col and any duplicate "SampleID"
      subset_cols <- setdiff(colnames(df), c(input$geo_sampleid_col, "SampleID"))
      subset_meta <- df[, subset_cols, drop = FALSE]

      combined <- cbind(SampleID = sample_ids, subset_meta)
      meta_data(align_meta_to_expr(combined, expr_data()))
    }
  })
  
  # ── GEO metadata mismatch management UI ───────────────────────────────────
  output$geo_meta_manage_ui <- shiny::renderUI({
    shiny::req(input$meta_source == "geo")
    meta <- meta_data(); expr <- expr_data()
    if (is.null(meta) || is.null(expr)) return(NULL)
    id_col  <- if ("SampleID" %in% colnames(meta)) "SampleID" else colnames(meta)[1]
    meta_ids <- as.character(meta[[id_col]])
    e_nms    <- colnames(expr)
    common   <- intersect(e_nms, meta_ids)
    if (setequal(e_nms, meta_ids)) return(NULL)  # already matched - nothing to show
    n_common <- length(common)
    n_expr   <- length(e_nms)
    shiny::tagList(
      shiny::hr(style="margin:8px 0;"),
      shiny::div(
        style="background:#fff3cd;border:1px solid #ffc107;border-radius:6px;padding:9px 12px;",
        shiny::tags$strong(style="font-size:0.88em;",
          paste0("⚠️ ", n_common, "/", n_expr, " samples matched")),
        shiny::tags$p(style="font-size:0.83em;margin:4px 0 8px;",
          paste0("Expression columns don't fully match the SampleID column. ",
                 if (n_common > 0)
                   paste0(n_common, " common samples found.")
                 else "No overlap found.")),
        if (n_common > 0 && n_common < n_expr)
          shiny::actionButton("geo_restrict_common", paste0("Use ", n_common, " matched samples"),
            class="btn-sm btn-warning"),
        shiny::tags$p(style="font-size:0.8em;color:#856404;margin:6px 0 0;",
          "Or change the SampleID column above to one that matches expression column names.")
      )
    )
  })

  shiny::observeEvent(input$geo_restrict_common, {
    meta <- meta_data(); expr <- expr_data()
    shiny::req(meta, expr)
    id_col   <- if ("SampleID" %in% colnames(meta)) "SampleID" else colnames(meta)[1]
    meta_ids <- as.character(meta[[id_col]])
    e_nms    <- colnames(expr)
    common   <- intersect(e_nms, meta_ids)
    if (length(common) == 0L) {
      shiny::showNotification("No common samples found.", type="error"); return()
    }
    idx  <- match(common, meta_ids)
    meta_sub <- meta[idx, , drop=FALSE]
    expr_sub <- expr[, common, drop=FALSE]
    meta_data(align_meta_to_expr(meta_sub, expr_sub))
    expr_data(expr_sub)
    log_step(paste0("GEO metadata restricted to ", length(common), " matched samples."))
    shiny::showNotification(
      paste0("Restricted to ", length(common), " matched samples."),
      type="default", duration=5)
  })

  # ── Persistent sample mismatch banner ─────────────────────────────────────
  # ── Upload: expression sample management ───────────────────────────────────
  output$upload_expr_manage_ui <- shiny::renderUI({
    req(input$expr_source == "upload", expr_data())
    nms <- colnames(expr_data())
    if (is.null(nms) || length(nms) == 0L) return(NULL)

    shiny::tagList(
      shiny::hr(style = "margin:8px 0;"),
      shiny::tags$strong(style = "font-size:0.87em;", "Remove samples"),
      shinyWidgets::pickerInput(
        "upload_expr_keep",
        label = NULL,
        choices  = nms,
        selected = nms,
        multiple = TRUE,
        options  = list(
          `actions-box`            = TRUE,
          `live-search`            = TRUE,
          `selected-text-format`   = "count > 3",
          `count-selected-text`    = "{0} of {1} kept"
        )
      ),
      shiny::actionButton(
        "apply_upload_expr_remove", "Apply",
        class = "btn-sm btn-outline-primary",
        style = "margin-top:2px;"
      )
    )
  })

  observeEvent(input$apply_upload_expr_remove, {
    req(expr_data(), input$upload_expr_keep)
    old_n <- ncol(expr_data())
    keep  <- intersect(input$upload_expr_keep, colnames(expr_data()))
    if (length(keep) > 0L) {
      expr_data(expr_data()[, keep, drop = FALSE])
      removed <- old_n - length(keep)
      if (removed > 0L) {
        log_step(paste0("Expression samples removed: ", removed,
                        " removed, ", length(keep), " kept."))
        shiny::showNotification(
          paste0(removed, " expression sample(s) removed - ", length(keep), " kept."),
          type = "warning", duration = 5
        )
      }
    }
  })

  # ── Upload: metadata sample management ─────────────────────────────────────
  output$upload_meta_manage_ui <- shiny::renderUI({
    req(input$meta_source == "upload", upload_meta_raw())
    raw  <- upload_meta_raw()
    cols <- colnames(raw)
    if (is.null(cols) || length(cols) == 0L) return(NULL)

    id_col  <- input$upload_meta_sampleid_col
    if (is.null(id_col) || !id_col %in% cols) id_col <- cols[1]
    row_ids <- as.character(raw[[id_col]])

    # Compute match against expression data (if available)
    expr    <- expr_data()
    matched <- if (!is.null(expr)) intersect(row_ids, colnames(expr)) else character(0)
    n_meta  <- length(row_ids)
    n_expr  <- if (!is.null(expr)) ncol(expr) else 0L
    n_match <- length(matched)
    fully_matched <- !is.null(expr) && n_match == n_expr && n_match == n_meta

    # ── Helper: bordered option card ──────────────────────────────────────────
    opt_card <- function(num, title, subtitle, body, border_color = "#dee2e6") {
      shiny::tags$div(
        style = paste0("border:1px solid ", border_color,
                       ";border-radius:6px;margin-bottom:10px;"),
        shiny::tags$div(
          style = paste0("display:flex;align-items:center;gap:8px;",
                         "padding:8px 10px;background:#f8f9fa;",
                         "border-bottom:1px solid ", border_color, ";"),
          shiny::tags$span(
            style = paste0("display:inline-flex;align-items:center;",
                           "justify-content:center;width:22px;height:22px;",
                           "border-radius:50%;background:#495057;color:#fff;",
                           "font-size:0.75em;font-weight:700;flex-shrink:0;"),
            num
          ),
          shiny::tags$div(
            shiny::tags$div(style = "font-weight:600;font-size:0.9em;line-height:1.2;", title),
            shiny::tags$div(style = "font-size:0.8em;color:#666;margin-top:2px;", subtitle)
          )
        ),
        shiny::tags$div(style = "padding:10px;", body)
      )
    }

    # ── SampleID selector always shown ───────────────────────────────────────
    id_selector <- shiny::tagList(
      shiny::hr(style = "margin:10px 0 6px;"),
      shiny::tags$strong(style = "font-size:0.87em;", "Sample ID column"),
      shiny::selectInput("upload_meta_sampleid_col", label = NULL,
                         choices = cols, selected = id_col)
    )

    # ── Match status banner ───────────────────────────────────────────────────
    match_banner <- if (is.null(expr)) {
      NULL
    } else if (fully_matched) {
      shiny::tags$div(
        style = paste0("background:#d4edda;color:#155724;border-radius:6px;",
                       "padding:8px 12px;margin-bottom:8px;font-size:0.87em;"),
        paste0("✅ All ", n_match, "/", n_expr, " samples matched.")
      )
    } else {
      dropped <- n_expr - n_match
      shiny::tags$div(
        style = paste0("background:#fff3cd;color:#856404;border-radius:6px;",
                       "padding:10px 12px;margin-bottom:8px;"),
        shiny::tags$strong(paste0("⚠️ ", n_match, "/", n_expr, " samples matched")),
        if (dropped > 0)
          shiny::tags$p(style = "margin:4px 0 0;font-size:0.85em;",
                        paste0(dropped, " expression column(s) not found in metadata. ",
                               "Use the options below to resolve this."))
      )
    }

    # ── Alignment options (only shown when there is a mismatch) ──────────────
    alignment_opts <- if (!is.null(expr) && !fully_matched) {
      shiny::tagList(

        opt_card("1", "Remove samples",
          "Deselect samples to remove from the metadata. Apply when done.",
          shiny::tagList(
            shinyWidgets::pickerInput(
              "upload_meta_keep", label = NULL,
              choices = row_ids, selected = row_ids, multiple = TRUE,
              options = list(`actions-box` = TRUE, `live-search` = TRUE,
                             `selected-text-format` = "count > 3",
                             `count-selected-text` = "{0} of {1} kept")
            ),
            shiny::actionButton("apply_upload_meta_manage", "Apply",
                                class = "btn-sm btn-outline-primary",
                                style = "margin-top:4px;")
          ),
          border_color = "#f5c6cb"
        ),

        opt_card("2", "Match by keyword in a metadata column",
          paste0("Use when a metadata column contains text that partially overlaps ",
                 "with expression column names (e.g. 'Library: sampleA' vs 'sampleA')."),
          shiny::tagList(
            shiny::selectInput(
              "upload_meta_match_col", "Select column to search:",
              choices = cols, selected = cols[1]
            ),
            shiny::actionButton("apply_upload_meta_match_col",
                                shiny::tagList(shiny::icon("magnifying-glass"), "Try keyword match"),
                                class = "btn-warning btn-sm")
          ),
          border_color = "#ffc107"
        ),

        opt_card("3", "Rename Sample IDs manually",
          paste0("Paste ", n_meta, " new IDs (one per line) to replace the current Sample IDs."),
          shiny::tagList(
            shiny::textAreaInput(
              "upload_meta_rename_ids", label = NULL,
              value = paste(row_ids, collapse = "\n"),
              rows  = min(6, n_meta), resize = "vertical", width = "100%"
            ),
            shiny::actionButton("apply_upload_meta_rename", "Rename & apply",
                                class = "btn-sm btn-outline-secondary",
                                style = "margin-top:4px;")
          ),
          border_color = "#dee2e6"
        )
      )
    } else {
      # Fully matched - still show simple remove picker
      shiny::tagList(
        shiny::tags$strong(style = "font-size:0.87em;", "Remove samples"),
        shinyWidgets::pickerInput(
          "upload_meta_keep", label = NULL,
          choices = row_ids, selected = row_ids, multiple = TRUE,
          options = list(`actions-box` = TRUE, `live-search` = TRUE,
                         `selected-text-format` = "count > 3",
                         `count-selected-text` = "{0} of {1} kept")
        ),
        shiny::actionButton("apply_upload_meta_manage", "Apply",
                            class = "btn-sm btn-outline-primary",
                            style = "margin-top:4px;")
      )
    }

    shiny::tagList(id_selector, match_banner, alignment_opts)
  })

  observeEvent(input$apply_upload_meta_manage, {
    req(upload_meta_raw(), input$upload_meta_sampleid_col)
    raw    <- upload_meta_raw()
    id_col <- input$upload_meta_sampleid_col
    if (!id_col %in% colnames(raw)) return()

    keep <- if (!is.null(input$upload_meta_keep)) {
      input$upload_meta_keep
    } else {
      as.character(raw[[id_col]])
    }

    # Filter rows, then rebuild with SampleID as first column
    mask     <- as.character(raw[[id_col]]) %in% keep
    filtered <- raw[mask, , drop = FALSE]
    ids      <- as.character(filtered[[id_col]])
    other    <- setdiff(colnames(filtered), c(id_col, "SampleID"))
    result   <- data.frame(
      SampleID = ids,
      filtered[, other, drop = FALSE],
      check.names = FALSE, stringsAsFactors = FALSE
    )
    rownames(result) <- ids
    meta_data(align_meta_to_expr(result, expr_data()))
    removed <- nrow(raw) - nrow(result)
    log_step(paste0("Metadata updated: SampleID column = '", id_col, "'",
                    if (removed > 0L) paste0(", ", removed, " sample(s) removed") else "",
                    " - ", nrow(result), " samples kept."))
    if (removed > 0L)
      shiny::showNotification(
        paste0(removed, " metadata sample(s) removed - ", nrow(result), " kept."),
        type = "warning", duration = 5
      )
  })

  # ── Upload metadata: rename Sample IDs ────────────────────────────────────
  observeEvent(input$apply_upload_meta_rename, {
    req(upload_meta_raw(), input$upload_meta_sampleid_col, input$upload_meta_rename_ids)
    raw    <- upload_meta_raw()
    id_col <- input$upload_meta_sampleid_col
    if (!id_col %in% colnames(raw)) return()

    new_ids <- trimws(strsplit(input$upload_meta_rename_ids, "\n")[[1]])
    new_ids <- new_ids[nchar(new_ids) > 0]
    old_ids <- as.character(raw[[id_col]])

    if (length(new_ids) != length(old_ids)) {
      shiny::showNotification(
        paste0("Number of new IDs (", length(new_ids), ") must match number of metadata rows (",
               length(old_ids), ")."),
        type = "error", duration = 8
      )
      return()
    }

    raw[[id_col]] <- new_ids
    other  <- setdiff(colnames(raw), c(id_col, "SampleID"))
    result <- data.frame(SampleID = new_ids, raw[, other, drop = FALSE],
                         check.names = FALSE, stringsAsFactors = FALSE)
    rownames(result) <- new_ids
    upload_meta_raw(result)  # update raw so future Apply uses new IDs
    meta_data(align_meta_to_expr(result, expr_data()))

    matched <- intersect(new_ids, colnames(expr_data()))
    log_step(paste0("Metadata Sample IDs renamed - ", length(matched), "/",
                    ncol(expr_data()), " now matched to expression columns."))
    shiny::showNotification(
      paste0(length(matched), "/", ncol(expr_data()),
             " samples matched after renaming."),
      type = if (length(matched) == ncol(expr_data())) "default" else "warning",
      duration = 6
    )
  }, ignoreInit = TRUE)

  # ── Upload metadata: match by keyword in a column ─────────────────────────
  observeEvent(input$apply_upload_meta_match_col, {
    req(upload_meta_raw(), input$upload_meta_sampleid_col,
        input$upload_meta_match_col, expr_data())
    raw      <- upload_meta_raw()
    id_col   <- input$upload_meta_sampleid_col
    match_col <- input$upload_meta_match_col
    if (!match_col %in% colnames(raw)) return()

    expr_cols  <- colnames(expr_data())
    meta_vals  <- as.character(raw[[match_col]])

    # For each expression column, find a metadata row whose match_col value
    # contains the expression column name as a substring (or vice versa).
    new_ids <- vapply(meta_vals, function(mv) {
      hit <- expr_cols[sapply(expr_cols, function(ec) {
        grepl(ec, mv, fixed = TRUE) || grepl(mv, ec, fixed = TRUE)
      })]
      if (length(hit) == 1L) hit else NA_character_
    }, character(1))

    n_matched <- sum(!is.na(new_ids))

    if (n_matched == 0L) {
      shiny::showNotification(
        paste0("No matches found between expression columns and '", match_col, "'. ",
               "Try a different column or use the rename option."),
        type = "warning", duration = 8
      )
      return()
    }

    # Replace SampleID with the matched expression column name (NA for unmatched)
    other  <- setdiff(colnames(raw), c(id_col, "SampleID"))
    result <- data.frame(SampleID = new_ids, raw[, other, drop = FALSE],
                         check.names = FALSE, stringsAsFactors = FALSE)
    # Keep only matched rows
    result <- result[!is.na(result$SampleID), , drop = FALSE]
    rownames(result) <- result$SampleID

    upload_meta_raw(result)
    meta_data(align_meta_to_expr(result, expr_data()))

    dropped <- nrow(raw) - nrow(result)
    log_step(paste0("Metadata matched by keyword in column '", match_col, "' - ",
                    n_matched, "/", ncol(expr_data()), " samples matched",
                    if (dropped > 0) paste0(", ", dropped, " unmatched rows dropped") else "", "."))
    shiny::showNotification(
      paste0(n_matched, "/", ncol(expr_data()), " samples matched via '", match_col, "'",
             if (dropped > 0) paste0(" - ", dropped, " unmatched rows removed") else "", "."),
      type = if (n_matched == ncol(expr_data())) "default" else "warning",
      duration = 7
    )
  }, ignoreInit = TRUE)

  # ── Session log UI + download ─────────────────────────────────────────────
  output$data_log_ui <- shiny::renderUI({
    log <- data_log()
    # Separate session preamble ([Session] lines) from action entries
    is_session <- startsWith(log, "[Session]")
    session_lines <- log[is_session]
    action_lines  <- log[!is_session]

    shiny::tagList(
      shiny::tags$div(
        style = "display:flex;align-items:center;justify-content:space-between;margin-bottom:6px;",
        shiny::tags$p(
          style = "font-size:0.75em;font-weight:600;color:#6b7280;text-transform:uppercase;letter-spacing:.05em;margin:0;",
          "Session log"
        ),
        shiny::downloadButton(
          "download_log", "Download log",
          class = "btn-sm btn-outline-secondary",
          style = "font-size:0.78em;padding:2px 10px;"
        )
      ),
      # Session preamble (collapsible, closed by default)
      if (length(session_lines) > 0)
        shiny::tags$details(
          style = "margin-bottom:6px;",
          shiny::tags$summary(
            style = "font-size:0.75em;color:#9ca3af;cursor:pointer;user-select:none;",
            "Session info"
          ),
          shiny::tags$pre(
            style = paste0(
              "background:#f1f5f9;border:1px solid #e2e8f0;border-radius:4px;",
              "padding:6px 10px;margin-top:4px;max-height:180px;overflow-y:auto;",
              "font-family:monospace;font-size:0.74em;line-height:1.6;color:#6b7280;",
              "white-space:pre-wrap;word-break:break-all;"
            ),
            paste(sub("^\\[Session\\]\\s*", "", session_lines), collapse = "\n")
          )
        ),
      # Action log entries (newest first)
      shiny::tags$div(
        style = paste0(
          "background:#f8fafc;border:1px solid #e2e8f0;border-radius:6px;",
          "padding:8px 12px;max-height:200px;overflow-y:auto;",
          "font-family:monospace;font-size:0.78em;line-height:1.7;color:#374151;"
        ),
        if (length(action_lines) == 0)
          shiny::tags$div(style = "color:#9ca3af;", "No actions logged yet.")
        else
          lapply(rev(action_lines), function(entry) {
            shiny::tags$div(entry)
          })
      )
    )
  })

  # ── Proceed / mismatch banner (shown at the TOP of the Data card) ────────────
  # Combines: persistent mismatch warning (same logic as sample_mismatch_banner)
  # with the "proceed to preprocessing" button when everything is OK.
  output$data_proceed_btn <- shiny::renderUI({
    expr <- tryCatch(expr_data(), error = function(x) NULL)
    meta <- tryCatch(meta_data(), error = function(x) NULL)

    # Nothing loaded yet - show nothing
    if (is.null(expr)) return(NULL)

    n_g <- nrow(expr); n_s <- ncol(expr)

    # Reuse the same mismatch check as .tab_lock_state / sample_mismatch_banner
    lock <- .tab_lock_state(expr, meta)

    if (lock$locked) {
      # ── Mismatch state: show the yellow warning banner ──────────────────────
      n_expr <- ncol(expr)
      n_meta <- if (!is.null(meta)) nrow(meta) else 0L
      e_nms  <- colnames(expr)
      msg <- if (!is.null(meta) && n_expr != n_meta) {
        paste0(n_expr, " expression sample", if (n_expr != 1) "s" else "",
               " vs ", n_meta, " metadata row", if (n_meta != 1) "s" else "", ".")
      } else if (!is.null(meta)) {
        m_ids <- as.character(meta[[ colnames(meta)[1] ]])
        n_miss <- length(setdiff(e_nms, m_ids))
        paste0(n_expr, " samples in both, but names differ - ",
               n_miss, " expression column", if (n_miss != 1) "s" else "",
               " not found in metadata.")
      } else {
        "Metadata not loaded."
      }
      shiny::tags$div(
        style = paste0(
          "background:#fff3cd;border:1px solid #ffc107;border-left:4px solid #ffc107;",
          "border-radius:4px;padding:9px 14px;margin:8px 8px 0 8px;"
        ),
        shiny::tags$div(
          style = "display:flex;align-items:center;gap:7px;",
          shiny::tags$span(style = "font-size:1em;", "⚠️"),
          shiny::tags$strong(style = "font-size:0.9em;", "Sample mismatch - "),
          shiny::tags$span(style = "font-size:0.87em;", msg)
        ),
        shiny::tags$div(
          style = "font-size:0.8em;color:#856404;margin-top:3px;",
          "Use the sample management controls in the sidebar to fix this before proceeding."
        )
      )
    } else if (!is.null(expr)) {
      # ── All good: show the green proceed banner ─────────────────────────────
      shiny::div(
        style = paste0(
          "background:#f0fdf4;border:1px solid #bbf7d0;border-left:4px solid #22c55e;",
          "border-radius:4px;padding:9px 14px;margin:8px 8px 0 8px;",
          "display:flex;align-items:center;justify-content:space-between;flex-wrap:wrap;gap:8px;"
        ),
        shiny::tags$div(
          shiny::tags$span(style = "font-weight:600;color:#166534;font-size:0.9em;",
            paste0("✅ Data ready: ",
                   format(n_g, big.mark = ","), " genes × ", n_s, " samples")),
          shiny::tags$div(style = "font-size:0.82em;color:#4b7a5e;margin-top:1px;",
            "Proceed to Preprocessing to filter, normalise, and review your data.")
        ),
        shiny::tags$button(
          onclick = paste0(
            "var lnk=Array.from(document.querySelectorAll('.nav-link'))",
            ".find(function(el){return el.textContent.trim()==='Preprocessing';});",
            "if(lnk&&!lnk.classList.contains('nav-tab-locked'))lnk.click();"
          ),
          class = "btn btn-success btn-sm",
          style = "white-space:nowrap;",
          shiny::icon("arrow-right"), " Proceed to Preprocessing"
        )
      )
    }
  })

  # ── Global log bar (fixed bottom - visible across all tabs) ───────────────
  output$global_log_header <- shiny::renderUI({
    log <- data_log()
    n   <- sum(!startsWith(log, "[Session]"))  # count only action entries
    shiny::tags$span(
      style = "color:#e2e8f0; font-size:0.88em;",
      if (n == 0L) "📋 Session log"
      else paste0("📋 Session log (", n,
                   if (n == 1L) " action" else " actions",
                   " - click to expand)")
    )
  })

  output$global_log_entries <- shiny::renderUI({
    log <- data_log()
    action_lines  <- log[!startsWith(log, "[Session]")]
    session_lines <- log[startsWith(log, "[Session]")]
    shiny::tagList(
      # Session info collapsed at top
      if (length(session_lines) > 0)
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
            paste(sub("^\\[Session\\]\\s*", "", session_lines), collapse = "\n")
          )
        ),
      if (length(action_lines) == 0)
        shiny::tags$div(style = "color:#64748b; padding:4px 0;", "No actions logged yet.")
      else
        lapply(rev(action_lines), function(e)
          shiny::tags$div(class = "glog-entry", e)
        ),
      shiny::downloadButton(
        "download_log", "Download log",
        class = "btn-sm btn-outline-light",
        style = "margin-top:8px; font-size:0.78em; padding:2px 10px;"
      )
    )
  })

  output$download_log <- shiny::downloadHandler(
    filename = function() paste0("markeR_session_log_",
                                 format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"),
    content  = function(file) {
      log <- data_log()
      # Write session info first, then a blank line, then action entries
      session_lines <- log[startsWith(log, "[Session]")]
      action_lines  <- log[!startsWith(log, "[Session]")]
      writeLines(c(session_lines, "", action_lines), file)
    }
  )

  ##### GENE SETS  #########################################

  observeEvent(input$geneset_source, {
    src_label <- switch(input$geneset_source,
      example = "example",
      upload  = "user upload",
      input$geneset_source
    )
    log_step(paste0("Gene set source changed to: ", src_label, "."))
    if (input$geneset_source == "example") .dl_gs_example()
  })

  observeEvent(input$retry_gs_example, { .dl_gs_example() }, ignoreInit = TRUE)



  observeEvent(input$geneset_files, {
    req(input$geneset_files)
    show_loading_modal(
      paste0("Reading ", nrow(input$geneset_files), " gene set file(s)…"),
      "Parsing and validating gene sets."
    )
    on.exit(shiny::removeModal(), add = TRUE)
    compiled <- list()
    for (i in seq_len(nrow(input$geneset_files))) {
      path <- input$geneset_files$datapath[i]
      name <- clean_name(input$geneset_files$name[i])
      obj <- tryCatch(
        safe_read_geneset(path),
        error = function(e) {
          shiny::showNotification(
            paste0("Could not read '", input$geneset_files$name[i], "': ", conditionMessage(e)),
            type = "error", duration = 10)
          NULL
        }
      )
      if (is.null(obj)) next
      parsed <- tryCatch(
        parse_geneset_object(obj, name),
        error = function(e) {
          shiny::showNotification(
            paste0("Could not parse '", input$geneset_files$name[i], "': ", conditionMessage(e)),
            type = "error", duration = 10)
          NULL
        }
      )
      if (!is.null(parsed)) compiled <- c(compiled, parsed)
    }
    gene_sets(compiled)
    log_step(paste0("Gene sets loaded: ",
                    nrow(input$geneset_files), " file(s) uploaded - ",
                    length(compiled), " set(s) parsed."))
  })
  
  observe({
    shiny::req(gene_sets())       # stop if no gene sets
    gs <- gene_sets()
    
    normalized <- lapply(names(gs), function(n) {
      obj <- gs[[n]]
      
      if (is.data.frame(obj)) {
        if (ncol(obj) > 2) {
          shiny::showNotification(
            paste0("Gene set '", n, "' has more than 2 columns. Only the first two will be used."),
            type = "warning"
          )
          obj <- obj[, 1:2]
        }
        if (ncol(obj) == 2) colnames(obj) <- c("Gene", "Direction")
        if (ncol(obj) == 1) obj <- as.vector(obj[[1]])
      }
      
      obj
    })
    
    names(normalized) <- names(gs)
    gene_sets(normalized)   # overwrite reactiveVal
  }) 
  
  
  
  ##### DATA SUMMARY (REACTIVE OUTPUT) ###############################
  output$summary_ui <- renderUI({
    expr <- expr_data()
    meta <- meta_data()
    gs   <- gene_sets()

    if (is.null(expr) && is.null(meta) && is.null(gs)) {
      return(shiny::tags$p(
        style = "color:#9ca3af;font-size:0.88em;padding:6px 0;",
        "No data loaded yet - use the sidebar to load expression data, metadata, and gene sets."
      ))
    }

    # Helper: one stat tile
    stat_tile <- function(icon_name, value, label, bg, fg, border) {
      shiny::tags$div(
        style = paste0(
          "flex:1 1 130px;min-width:110px;",
          "background:", bg, ";border:1px solid ", border, ";",
          "border-radius:8px;padding:10px 14px;display:flex;",
          "flex-direction:column;gap:2px;"
        ),
        shiny::tags$div(
          style = paste0("font-size:1.5em;font-weight:700;color:", fg, ";line-height:1.1;"),
          value
        ),
        shiny::tags$div(
          style = "font-size:0.77em;color:#6b7280;font-weight:500;",
          shiny::icon(icon_name, style = paste0("color:", fg, ";margin-right:4px;")),
          label
        )
      )
    }

    tiles <- list()
    if (!is.null(expr)) {
      tiles <- c(tiles, list(
        stat_tile("dna", format(nrow(expr), big.mark=","), "genes",
                  "#eff6ff", "#1d4ed8", "#bfdbfe"),
        stat_tile("vials", format(ncol(expr), big.mark=","), "expression samples",
                  "#f0fdf4", "#15803d", "#bbf7d0")
      ))
    }
    if (!is.null(meta)) {
      n_meta_vars <- ncol(meta) - ("SampleID" %in% colnames(meta))
      tiles <- c(tiles, list(
        stat_tile("table", format(nrow(meta), big.mark=","),
                  paste0("metadata rows × ", n_meta_vars, " vars"),
                  "#faf5ff", "#7e22ce", "#e9d5ff")
      ))
    }
    if (!is.null(gs)) {
      tiles <- c(tiles, list(
        stat_tile("layer-group", length(gs), "gene set(s)",
                  "#fff7ed", "#c2410c", "#fed7aa")
      ))
    }

    shiny::tags$div(
      style = "margin-bottom:4px;",
      shiny::tags$p(
        style = "font-size:0.75em;font-weight:600;color:#6b7280;text-transform:uppercase;letter-spacing:.05em;margin-bottom:6px;",
        "Loaded data summary"
      ),
      shiny::tags$div(
        style = "display:flex;flex-wrap:wrap;gap:8px;",
        tiles
      )
    )
  })
  
  ##### PREVIEW TABLES ################################################
  
  
  output$expr_quality_banner <- shiny::renderUI({
    msg <- expr_quality_warn()
    if (is.null(msg)) return(NULL)
    shiny::tags$div(
      style = paste0(
        "background:#fff3cd;border:1px solid #ffc107;border-left:4px solid #ffc107;",
        "border-radius:4px;padding:9px 14px;margin:8px 8px 4px 8px;"
      ),
      shiny::tags$div(
        style = "display:flex;align-items:center;gap:7px;",
        shiny::tags$span(style = "font-size:1em;", "⚠️"),
        shiny::tags$strong(style = "font-size:0.9em;", "Suspicious data - "),
        shiny::tags$span(style = "font-size:0.87em;", msg)
      )
    )
  })

  output$expr_preview <- DT::renderDT({
    shiny::req(expr_data())
    
    DT::datatable(
      round(expr_data(),2),  # round for better display
      options = list(
        pageLength = 5,
        lengthMenu = c(5, 10, 20, 50),
        scrollX = TRUE,
        scrollY = "250px",      # <-- added vertical scroll
        paging = TRUE
      ),
      rownames = TRUE
    )
  })
  
  output$meta_preview <- DT::renderDT({
    shiny::req(meta_data())
    
    DT::datatable(
      meta_data(),
      options = list(
        pageLength = 5,
        lengthMenu = c(5, 10, 20, 50),
        scrollX = TRUE,
        scrollY = "250px",      # <-- added vertical scroll
        paging = TRUE
      ),
      rownames = FALSE
    )
  })
  
  output$geneset_summary <- DT::renderDT({
    shiny::req(gene_sets())  # ensures gene sets are ready and normalized
    
    gs <- gene_sets()        # always safe now: cleaned objects
    
    df <- do.call(rbind, lapply(names(gs), function(n) {
      obj <- gs[[n]]
      
      if (is.data.frame(obj) && "Direction" %in% colnames(obj)) {
        data.frame(
          `Gene Set` = n,
          Genes   = nrow(obj),
          `Positive (+1)`     = sum(obj$Direction == 1, na.rm = TRUE),
          `Negative (-1)`     = sum(obj$Direction == -1, na.rm = TRUE),
          check.names = FALSE
        )
      } else {
        data.frame(
          `Gene Set` = n,
          Genes   = length(obj),
          `Positive (+1)`     = NA,
          `Negative (-1)`     = NA,
          check.names = FALSE
        )
      }
    }))
    
    DT::datatable(
      df,
      selection = "single",
      options = list(
        pageLength = 5,
        scrollX = TRUE,
        scrollY = "250px",
        lengthMenu = c(5, 10, 20, 50),
        escape = FALSE  # allows + and parentheses to display properly
        
      ),
      rownames = FALSE
    )
  })
  
  # pop up window with gene set composition
  shiny::observeEvent(input$geneset_summary_rows_selected, {
    shiny::req(input$geneset_summary_rows_selected) #input automatically created by DT
    gs <- gene_sets()
    
    # Get selected gene set name
    selected_row <- input$geneset_summary_rows_selected
    selected_name <- names(gs)[selected_row]
    geneset_obj <- gs[[selected_name]]
    
    # Prepare table to show in modal
    if (is.data.frame(geneset_obj) && "Direction" %in% colnames(geneset_obj)) {
      modal_df <- geneset_obj
    } else {
      modal_df <- data.frame(Gene = geneset_obj)
    }
    
    # Show modal
    shiny::showModal(
      shiny::modalDialog(
        title = paste("Genes in set:", selected_name),
        
        # Wrap DT in a div to control width
        shiny::tags$div(
          style = "max-width: 700px; margin: auto;",  # narrow table, centered in modal
          DT::renderDT({
            DT::datatable(
              modal_df,
              options = list(
                pageLength = 10,
                scrollX = TRUE,
                scrollY = "300px",        #  height of scrollable area
                scrollCollapse = TRUE,    #  collapse if fewer rows
                paging = FALSE,           #  optional, better for modal scroll
                columnDefs = list(
                  list(className = 'dt-center', targets = "_all")
                )
              ),
              rownames = FALSE
            )
          })
        ),
        
        easyClose = TRUE,
        size = "l"  # modal is large
      )
    )
    
    
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
  
  options(shiny.maxRequestSize = 1000 * 1024^2)
  
  app <- shiny::shinyApp(ui, server)
  shiny::runApp(app, ...)
}