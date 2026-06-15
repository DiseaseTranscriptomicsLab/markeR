#' @importFrom bslib bs_theme
#' @import shiny
#' @import GEOquery
#' @import DT
#' @importFrom readxl read_xls read_xlsx
#' @importFrom R.utils gunzip
#' @importFrom data.table fread

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
    shiny::a("dev",
             style = "font-size: 0.7em; color: #949494; text-decoration: none;")
  ),
  
  theme = bslib::bs_theme(bootswatch = "sketchy"),
  
  
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
          shiny::a("markeR's paper", href = "https://www.biorxiv.org/content/10.64898/2025.12.05.692517", target = "_blank"),
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
            shiny::h4("Data Configuration"),
            shiny::p(
              "Load gene expression, metadata and gene sets.",
              style = "color:#6c757d;"
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
                  shiny::tags$a(href="https://www.biorxiv.org/content/10.64898/2025.12.05.692517", "bioRxiv", target="_blank"),
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
                  shiny::tags$a(href="https://www.biorxiv.org/content/10.64898/2025.12.05.692517", "bioRxiv", target="_blank"),
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
                  "Upload (.csv, .rds)",
                  accept = c(".csv", ".rds")
                ),
                # Sample management controls (shown after upload)
                shiny::uiOutput("upload_meta_manage_ui")
              ),
              # Show GEO metadata column selection only if 'geo' is selected
              # Only show if using GEO as metadata source
              shiny::conditionalPanel(
                "input.meta_source == 'geo'",
                
                # SampleID column selection
                shiny::uiOutput("geo_sampleid_selector"),
                
                # Metadata columns selection
                shinyWidgets::pickerInput(
                  inputId = "selected_meta_cols",
                  label   = "Select metadata columns to keep:",
                  choices = NULL,   # updated server-side
                  selected = NULL,
                  multiple = TRUE,
                  options = list(`actions-box` = TRUE)
                )
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

          # Persistent mismatch banner — stays until sample counts/names align
          shiny::uiOutput("sample_mismatch_banner"),

          shiny::div(style = "overflow-y: auto; max-height: 900px;",  # <-- scrollable container
                     shiny::uiOutput("summary_ui"),
                     shiny::hr(),
                     shiny::h5("Expression Preview"),
                     DT::DTOutput("expr_preview"),
                     shiny::h5("Metadata Preview"),
                     DT::DTOutput("meta_preview"),
                     shiny::h5("Gene Sets Overview"),
                     shiny::tags$div(
                       style = "margin-bottom: 10px; font-style: italic; color: #555;",
                       "💡 Tip: Click on a gene set row to see its genes and expected direction of regulation change."
                     ),
                     DT::DTOutput("geneset_summary"),

                     # ── Processing log ────────────────────────────────────
                     shiny::hr(),
                     shiny::uiOutput("data_log_ui")
          )
        )
        
        
      )
    )
  ),
  
  ##### PREPROCESSING #####
  bslib::nav_panel(
    "Preprocessing",
    shiny::h4("Coming soon..."),
    shiny::h1("👁️👄👁")
  ),
  
  ##### GENE SETS #####
  bslib::nav_panel(
    "Gene Sets",
    shiny::h4("Coming soon..."),
    shiny::h1("👁️👄👁")
  ),
  
  ##### BENCHMARKING MODE #####
  bslib::nav_panel(
    "Benchmarking Mode",
    shiny::h4("Coming soon..."),
    shiny::h1("👁️👄👁")
  ),
  
  ##### DISCOVERY MODE #####
  bslib::nav_panel(
    "Discovery Mode",
    shiny::h4("Coming soon..."),
    shiny::h1("👁️👄👁")
  ) 
  
  
  
)

# Then wrap whole page with head()
ui <- tagList(
  tags$head(
    tags$style(HTML("
      #shiny-notification-panel {
        top: auto !important;
        bottom: 50px !important;
        right: 30px !important;
        left: auto !important;
      }
    "))
  ),
  ui
)

# ---- Server ----

server <- function(input, output, session) {
  
  ###### REACTIVE STORAGE CONTAINERS ######
  expr_data   <- shiny::reactiveVal(NULL)  # Expression matrix/data.frame
  meta_data   <- shiny::reactiveVal(NULL)  # Sample metadata
  gene_sets   <- shiny::reactiveVal(NULL)  # Named list of gene sets
  geo_objects <- shiny::reactiveVal(NULL)  # Raw GEO objects (can be multiple)
  geo_meta_cols <- reactiveVal(NULL) # Will store the original metadata column names (excluding SampleID)
  full_geo_meta      <- reactiveVal(NULL)
  upload_meta_raw    <- reactiveVal(NULL)   # Raw uploaded metadata (for SampleID re-mapping)
  module_meta_active <- reactiveVal(FALSE)  # TRUE while geo_import_module owns meta_data()
  data_log           <- reactiveVal(character(0))  # Timestamped processing log
  geo_supp_files  <- reactiveVal(NULL)
  expr_candidates_available <- reactiveVal(FALSE) # Flag to indicate multiple candidate files are available after auto detection
  
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

  # Append a timestamped entry to the processing log.
  # isolate() prevents the calling observer from taking a reactive dependency
  # on data_log — without it, writing data_log() re-invalidates the observer
  # that called log_step(), causing an infinite reactive loop.
  log_step <- function(msg) {
    ts <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
    data_log(c(isolate(data_log()), paste(ts, msg)))
  }

  # Show a Zenodo download-failure modal with a "Try Again" button.
  # retry_btn_id: inputId of the actionButton placed in the footer.
  show_dl_error_modal <- function(msg, retry_btn_id) {
    shiny::showModal(shiny::modalDialog(
      title = shiny::tags$span(
        shiny::icon("circle-exclamation", style = "color:#dc3545;margin-right:6px;"),
        "Download Failed"
      ),
      shiny::tags$div(
        shiny::tags$p("Could not download the file from Zenodo."),
        shiny::tags$p(style = "font-size:0.84em;color:#888;font-family:monospace;word-break:break-all;", msg),
        shiny::tags$p("Check your internet connection and click “Try Again”, or try a different data source.")
      ),
      footer = shiny::tagList(
        shiny::actionButton(
          retry_btn_id,
          shiny::tagList(shiny::icon("rotate-right"), " Try Again"),
          class = "btn-primary"
        ),
        shiny::modalButton("Cancel")
      ),
      easyClose = TRUE
    ))
  }

  clean_name <- function(x) tools::file_path_sans_ext(base::basename(x))
  
  safe_read_table <- function(path) {
    
    ext <- tools::file_ext(path)
    
    if (ext %in% c("txt", "csv")) {
      
      # --- guess separator from first line ---
      first_line <- readLines(path, n = 1)
      sep <- if (grepl("\t", first_line)) {
        "\t"
      } else if (grepl(";", first_line)) {
        ";"
      } else if (grepl(",", first_line)) {
        ","
      } else {
        ""  # fallback: any whitespace
      }
      
      # --- read everything as character ---
      df <- utils::read.table(
        path,
        header = TRUE,
        row.names = 1,      # first column = gene names
        check.names = FALSE,
        stringsAsFactors = FALSE,
        sep = sep,
        comment.char = "",
        quote = "\"",
        fill = TRUE,
        colClasses = "character"
      )
      
      # --- remove empty trailing columns (from trailing separators) ---
      df <- df[, colSums(!is.na(df) & df != "") > 0, drop = FALSE]
      
      # --- normalize decimals and coerce to numeric ---
      df[] <- lapply(df, function(x) {
        x <- gsub(",", ".", x)
        as.numeric(x)
      })
      
      # --- validation ---
      if (!all(sapply(df, is.numeric))) {
        stop(
          "Uploaded file could not be parsed as numeric expression matrix. ",
          "Check the first row, separator, and decimal symbols."
        )
      }
      
      if (is.null(rownames(df))) {
        stop("First column must contain gene names.")
      }
      
      return(base::as.data.frame(df))
      
    } else if (ext == "rds") {
      
      df <- readRDS(path)
      
      # if it’s a matrix, coerce to data.frame
      if (inherits(df, "matrix")) {
        df <- as.data.frame(df, stringsAsFactors = FALSE)
      }
      
      # now check type: must be data.frame
      if (!inherits(df, "data.frame")) {
        stop("Uploaded RDS must contain a data.frame or matrix.")
      }
      
      # validate row and column names
      if (is.null(rownames(df)) || nrow(df) == 0) {
        stop("Uploaded RDS must have row names (gene names).")
      }
      
      if (is.null(colnames(df)) || ncol(df) == 0) {
        stop("Uploaded RDS must have column names (sample names).")
      }
      
      return(df)
      
    } else {
      
      stop("Unsupported file format.")
      
    }
  }
  
  parse_geneset_object <- function(obj, name_hint) {
    if (is.character(obj)) return(stats::setNames(list(obj), name_hint))
    if (is.data.frame(obj)) {
      if (ncol(obj) == 1) return(stats::setNames(list(obj[[1]]), name_hint))
      if (ncol(obj) >= 2) {
        df <- obj[, 1:2]
        colnames(df) <- c("gene", "direction")
        return(stats::setNames(list(df), name_hint))
      }
    }
    if (is.list(obj)) return(obj)
    stop("Invalid gene set format.")
  }
  
  ###### DATA ######

  # ---- Zenodo example download helpers (closures — capture server-scope reactives) ----
  # Each function manages its own modals: shows spinner, wraps download in tryCatch,
  # removes spinner on success, replaces it with an error modal (with retry button) on failure.
  # Called from both the main radio-button observer AND the retry observeEvent so the
  # same logic runs on first attempt and on any number of retries without code duplication.

  .dl_expr_raw <- function() {
    show_loading_modal("Loading example data…",
                       "Downloading unprocessed counts from Zenodo.")
    dl_ok <- tryCatch({
      shiny::withProgress(message = "Downloading unprocessed example data...", value = 0, {
        incProgress(0.1, detail = "Preparing download...")
        Sys.sleep(0.2)
        tmp <- tempfile(fileext = ".rds")
        incProgress(0.4, detail = "Downloading from Zenodo...")
        utils::download.file(
          url      = "https://zenodo.org/records/18714122/files/counts.rds?download=1",
          destfile = tmp,
          mode     = "wb"
        )
        incProgress(0.3, detail = "Loading expression matrix...")
        e <- readRDS(tmp)
        expr_data(e)
        log_step(paste0("Expression data loaded: example (unprocessed) — ",
                        nrow(e), " genes × ", ncol(e), " samples."))
        incProgress(0.2, detail = "Finalizing...")
      })
      TRUE
    }, error = function(err) {
      shiny::removeModal()
      show_dl_error_modal(conditionMessage(err), "retry_expr_raw")
      FALSE
    })
    if (isTRUE(dl_ok)) {
      shiny::removeModal()
      shiny::showNotification("Unprocessed expression data loaded successfully.",
                              type = "default", duration = 3)
    }
  }

  .dl_expr_proc <- function() {
    show_loading_modal("Loading example data…",
                       "Downloading processed counts from Zenodo.")
    dl_ok <- tryCatch({
      shiny::withProgress(message = "Downloading processed example data...", value = 0, {
        incProgress(0.1, detail = "Preparing download...")
        Sys.sleep(0.2)
        tmp <- tempfile(fileext = ".rds")
        incProgress(0.4, detail = "Downloading from Zenodo...")
        utils::download.file(
          url      = "https://zenodo.org/records/18714122/files/corrcounts.rds?download=1",
          destfile = tmp,
          mode     = "wb"
        )
        incProgress(0.3, detail = "Loading processed matrix...")
        e <- readRDS(tmp)
        expr_data(e)
        log_step(paste0("Expression data loaded: example (processed) — ",
                        nrow(e), " genes × ", ncol(e), " samples."))
        incProgress(0.2, detail = "Finalizing...")
      })
      TRUE
    }, error = function(err) {
      shiny::removeModal()
      show_dl_error_modal(conditionMessage(err), "retry_expr_proc")
      FALSE
    })
    if (isTRUE(dl_ok)) {
      shiny::removeModal()
      shiny::showNotification("Processed expression data loaded successfully.",
                              type = "default", duration = 3)
    }
  }

  .dl_meta_example <- function() {
    show_loading_modal("Loading example metadata…",
                       "Downloading from Zenodo.")
    dl_ok <- tryCatch({
      shiny::withProgress(message = "Downloading example metadata...", value = 0, {
        incProgress(0.2, detail = "Preparing download...")
        Sys.sleep(0.2)
        tmp <- tempfile(fileext = ".rds")
        incProgress(0.5, detail = "Downloading from Zenodo...")
        utils::download.file(
          url      = "https://zenodo.org/records/18714122/files/metadata.rds?download=1",
          destfile = tmp,
          mode     = "wb"
        )
        incProgress(0.2, detail = "Loading metadata...")
        m <- readRDS(tmp)
        meta_data(m)
        log_step(paste0("Metadata loaded: example — ",
                        nrow(m), " samples × ", ncol(m), " variables."))
        incProgress(0.1, detail = "Finalizing...")
      })
      TRUE
    }, error = function(err) {
      shiny::removeModal()
      show_dl_error_modal(conditionMessage(err), "retry_meta_example")
      FALSE
    })
    if (isTRUE(dl_ok)) {
      shiny::removeModal()
      shiny::showNotification("Metadata loaded successfully.", type = "default", duration = 3)
    }
  }

  .dl_gs_example <- function() {
    show_loading_modal("Loading example gene sets…",
                       "Downloading from Zenodo.")
    dl_ok <- tryCatch({
      shiny::withProgress(message = "Downloading example gene sets...", value = 0, {
        incProgress(0.2, detail = "Preparing download...")
        Sys.sleep(0.2)
        tmp <- tempfile(fileext = ".rds")
        incProgress(0.5, detail = "Downloading from Zenodo...")
        utils::download.file(
          url      = "https://zenodo.org/records/18714122/files/SenescenceGeneSets.rds?download=1",
          destfile = tmp,
          mode     = "wb"
        )
        incProgress(0.2, detail = "Loading gene sets...")
        gs <- readRDS(tmp)
        gene_sets(gs)
        log_step(paste0("Gene sets loaded: example — ", length(gs), " set(s), ",
                        sum(sapply(gs, function(x) if (is.data.frame(x)) nrow(x) else length(x))),
                        " genes total."))
        incProgress(0.1, detail = "Finalizing...")
      })
      TRUE
    }, error = function(err) {
      shiny::removeModal()
      show_dl_error_modal(conditionMessage(err), "retry_gs_example")
      FALSE
    })
    if (isTRUE(dl_ok)) {
      shiny::removeModal()
      shiny::showNotification("Gene sets loaded successfully.", type = "default", duration = 3)
    }
  }

  ####### > Example Data #######

  observeEvent(input$expr_source, {
    updateRadioButtons(session, "meta_source", selected = "example")
    if (input$expr_source == "example_raw")  .dl_expr_raw()
    if (input$expr_source == "example_proc") .dl_expr_proc()
  })

  # "Try Again" retry handlers — call the same helper function directly,
  # no need to toggle the radio button (avoids unreliable updateRadioButtons tricks).
  observeEvent(input$retry_expr_raw,  { shiny::removeModal(); .dl_expr_raw()  }, ignoreInit = TRUE)
  observeEvent(input$retry_expr_proc, { shiny::removeModal(); .dl_expr_proc() }, ignoreInit = TRUE)
  
  ###### > User Input Data ######
  
  observeEvent(input$expr_file, {
    req(input$expr_file)
    show_loading_modal(paste0("Reading ", input$expr_file$name, "…"),
                       "Parsing expression matrix.")
    on.exit(shiny::removeModal(), add = TRUE)
    e <- safe_read_table(input$expr_file$datapath)
    expr_data(e)
    module_meta_active(FALSE)
    log_step(paste0("Expression data loaded: uploaded file '", input$expr_file$name, "' — ",
                    nrow(e), " genes × ", ncol(e), " samples."))
  })
  
  ###### > GEO data — delegated to geoImportModule ######

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
      log_step(paste0("Expression data loaded: GEO import — ",
                      nrow(e), " genes × ", ncol(e), " samples."))
    }

    if (!is.null(mod_meta)) {
      # Set flag FIRST so the app-level column-selector observers skip overwriting
      module_meta_active(TRUE)
      full_geo_meta(mod_meta)
      meta_data(mod_meta)
      log_step(paste0("Metadata loaded: GEO — ",
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
  
  # --- Extraction helpers for archives (legacy — superseded by geo_import_module) ---
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
  # The block below is DEAD CODE — kept only for reference. It is never called.
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
        choices = c("— Please choose a file —" = "",
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
        showNotification("Expression data loaded from supplementary file.", type = "default")
        
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
    if (input$meta_source == "example") .dl_meta_example()
  })

  observeEvent(input$retry_meta_example, { shiny::removeModal(); .dl_meta_example() }, ignoreInit = TRUE)
  
  
  
  observeEvent(input$meta_file, {
    req(input$meta_file)
    show_loading_modal(paste0("Reading ", input$meta_file$name, "…"),
                       "Parsing metadata file.")
    on.exit(shiny::removeModal(), add = TRUE)
    raw <- safe_read_table(input$meta_file$datapath)
    upload_meta_raw(raw)   # keep unmodified copy so SampleID re-mapping is non-destructive
    module_meta_active(FALSE)
    meta_data(raw)
    log_step(paste0("Metadata loaded: uploaded file '", input$meta_file$name, "' — ",
                    nrow(raw), " samples × ", ncol(raw), " variables."))
  })
  
  
  
  # --- GEO metadata handling with SampleID mapping and subsetting ---
  observe({
    req(full_geo_meta(), input$meta_source)

    if (input$meta_source == "geo") {
      df <- full_geo_meta()
      
      # Determine default SampleID column.
      # When the geo_import_module provided the data it already built a "SampleID"
      # column — prefer that unconditionally. For the legacy pData path, fall back
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
      
      # Update pickerInput for metadata column selection (exclude SampleID)
      cols <- setdiff(colnames(df), default_sampleid)
      shinyWidgets::updatePickerInput(
        session,
        "selected_meta_cols",
        choices = cols,
        selected = cols
      )
    }
  })
  
  # --- Update meta_data() based on SampleID and selected columns ---
  observe({
    req(full_geo_meta(), input$geo_sampleid_col, input$selected_meta_cols, input$meta_source)
    # Skip when geo_import_module owns meta_data() — module already built the final metadata
    if (isTRUE(module_meta_active())) return(NULL)

    if (input$meta_source == "geo") {
      df <- full_geo_meta()
      
      # Map SampleID column
      sample_ids <- if (input$geo_sampleid_col %in% colnames(df)) {
        df[[input$geo_sampleid_col]]
      } else {
        rownames(df)
      }
      
      # Subset remaining metadata columns, ensuring SampleID is first.
      # Exclude both the chosen id_col AND any pre-existing "SampleID" column so
      # we never end up with two columns of the same name (can happen when the
      # geo_import_module already added a SampleID column to full_geo_meta).
      subset_cols <- setdiff(input$selected_meta_cols,
                             c(input$geo_sampleid_col, "SampleID"))
      subset_meta <- df[, intersect(subset_cols, colnames(df)), drop = FALSE]

      meta_data(cbind(SampleID = sample_ids, subset_meta))
    }
  })
  
  # ── Persistent sample mismatch banner ─────────────────────────────────────
  output$sample_mismatch_banner <- shiny::renderUI({
    expr <- expr_data()
    meta <- meta_data()
    if (is.null(expr) || is.null(meta)) return(NULL)

    n_expr <- ncol(expr)
    n_meta <- nrow(meta)
    e_nms  <- colnames(expr)

    count_ok <- (n_expr == n_meta)

    # Primary check: first column of metadata vs expression column names.
    # Also try rownames and any "SampleID" column as fallbacks so we never
    # produce a false positive when data is internally consistent.
    names_ok <- FALSE
    if (count_ok) {
      candidates <- unique(list(
        as.character(meta[[ colnames(meta)[1] ]]),          # first column  ← primary
        if (!is.null(meta[["SampleID"]])) as.character(meta[["SampleID"]]) else NULL,
        rownames(meta)
      ))
      candidates <- Filter(Negate(is.null), candidates)
      names_ok <- any(vapply(candidates, function(ids) setequal(e_nms, ids), logical(1)))
    }

    if (count_ok && names_ok) return(NULL)

    # Build the best-guess ID vector for the message (prefer first column)
    m_ids <- as.character(meta[[ colnames(meta)[1] ]])

    msg <- if (!count_ok) {
      paste0(n_expr, " expression sample", if (n_expr != 1) "s" else "",
             " vs ", n_meta, " metadata row", if (n_meta != 1) "s" else "", ".")
    } else {
      n_missing <- length(setdiff(e_nms, m_ids))
      paste0(n_expr, " samples in both, but names differ — ",
             n_missing, " expression column", if (n_missing != 1) "s" else "",
             " not found in metadata.")
    }

    shiny::tags$div(
      style = paste0(
        "background:#fff3cd;border:1px solid #ffc107;border-left:4px solid #ffc107;",
        "border-radius:4px;padding:9px 14px;margin:8px 8px 0 8px;"
      ),
      shiny::tags$div(
        style = "display:flex;align-items:center;gap:7px;",
        shiny::tags$span(style = "font-size:1em;", "⚠️"),
        shiny::tags$strong(style = "font-size:0.9em;", "Sample mismatch — "),
        shiny::tags$span(style = "font-size:0.87em;", msg)
      ),
      shiny::tags$div(
        style = "font-size:0.8em;color:#856404;margin-top:3px;",
        "Use the sample management controls in the sidebar to fix this before running analysis."
      )
    )
  })

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
      if (removed > 0L)
        log_step(paste0("Expression samples removed: ", removed,
                        " removed, ", length(keep), " kept."))
    }
  })

  # ── Upload: metadata sample management ─────────────────────────────────────
  output$upload_meta_manage_ui <- shiny::renderUI({
    req(input$meta_source == "upload", upload_meta_raw())
    raw  <- upload_meta_raw()
    cols <- colnames(raw)
    if (is.null(cols) || length(cols) == 0L) return(NULL)

    # Re-render reactively when user changes the SampleID column
    id_col  <- input$upload_meta_sampleid_col
    if (is.null(id_col) || !id_col %in% cols) id_col <- cols[1]
    row_ids <- as.character(raw[[id_col]])

    shiny::tagList(
      shiny::hr(style = "margin:8px 0;"),
      shiny::tags$strong(style = "font-size:0.87em;", "Sample ID column"),
      shiny::selectInput(
        "upload_meta_sampleid_col",
        label    = NULL,
        choices  = cols,
        selected = id_col
      ),
      shiny::tags$strong(style = "font-size:0.87em;", "Remove samples"),
      shinyWidgets::pickerInput(
        "upload_meta_keep",
        label    = NULL,
        choices  = row_ids,
        selected = row_ids,
        multiple = TRUE,
        options  = list(
          `actions-box`            = TRUE,
          `live-search`            = TRUE,
          `selected-text-format`   = "count > 3",
          `count-selected-text`    = "{0} of {1} kept"
        )
      ),
      shiny::actionButton(
        "apply_upload_meta_manage", "Apply",
        class = "btn-sm btn-outline-primary",
        style = "margin-top:2px;"
      )
    )
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
    meta_data(result)
    removed <- nrow(raw) - nrow(result)
    log_step(paste0("Metadata updated: SampleID column = '", id_col, "'",
                    if (removed > 0L) paste0(", ", removed, " sample(s) removed") else "",
                    " — ", nrow(result), " samples kept."))
  })

  # ── Processing log UI + download ───────────────────────────────────────────
  output$data_log_ui <- shiny::renderUI({
    log <- data_log()
    if (length(log) == 0L) return(NULL)

    shiny::tagList(
      shiny::tags$div(
        style = "display:flex;align-items:center;justify-content:space-between;margin-bottom:6px;",
        shiny::tags$p(
          style = "font-size:0.75em;font-weight:600;color:#6b7280;text-transform:uppercase;letter-spacing:.05em;margin:0;",
          "Processing log"
        ),
        shiny::downloadButton(
          "download_log", "Download log",
          class = "btn-sm btn-outline-secondary",
          style = "font-size:0.78em;padding:2px 10px;"
        )
      ),
      shiny::tags$div(
        style = paste0(
          "background:#f8fafc;border:1px solid #e2e8f0;border-radius:6px;",
          "padding:8px 12px;max-height:200px;overflow-y:auto;",
          "font-family:monospace;font-size:0.78em;line-height:1.7;color:#374151;"
        ),
        lapply(rev(log), function(entry) {  # newest first
          shiny::tags$div(entry)
        })
      )
    )
  })

  output$download_log <- shiny::downloadHandler(
    filename = function() paste0("markeR_processing_log_",
                                 format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"),
    content  = function(file) {
      writeLines(data_log(), file)
    }
  )

  ##### GENE SETS  #########################################

  observeEvent(input$geneset_source, {
    if (input$geneset_source == "example") .dl_gs_example()
  })

  observeEvent(input$retry_gs_example, { shiny::removeModal(); .dl_gs_example() }, ignoreInit = TRUE)



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
      obj <- safe_read_table(path)
      parsed <- parse_geneset_object(obj, name)
      compiled <- c(compiled, parsed)
    }
    gene_sets(compiled)
    log_step(paste0("Gene sets loaded: ",
                    nrow(input$geneset_files), " file(s) uploaded — ",
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
        "No data loaded yet — use the sidebar to load expression data, metadata, and gene sets."
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