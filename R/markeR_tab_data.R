# =============================================================================
# Data Tab
# =============================================================================
# Handles loading and validating everything the rest of the app depends on:
# the expression matrix, sample metadata, and gene sets (example data, file
# upload, or GEO import). Organised the same way as the other tabs (own file,
# a UI builder and a server function), but NOT wrapped in shiny::moduleServer /
# NS(): the existing UI relies on plain (non-namespaced) conditionalPanel JS
# conditions (e.g. "input.expr_source == 'upload'") and on a session log bar
# that lives outside this tab in the app's global UI, both of which would
# need reworking to use module namespacing safely. dataServer() is called
# once from the main server with the app's own input/output/session, exactly
# as the code behaved before this file existed - only the location changed.
#
# dataServer() returns the reactives the rest of the app depends on:
#   $expr_data, $meta_data, $gene_sets  - the loaded data (reactiveVal; call
#     with an argument to set, e.g. expr_data(new_df), same as before)
#   $log_step, $data_log               - session log (passed into
#     preprocessingServer() as log_fn / get_log)
#   $expr_quality_warn                 - reactiveVal flagging suspicious
#     expression data (set here, read by this tab's own quality banner)
# =============================================================================

if (!exists("%||%", mode = "function"))
  `%||%` <- function(a, b) if (!is.null(a)) a else b

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

# ---- UI --------------------------------------------------------------------

dataUI <- function() {
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
                  shiny::HTML(paste0(
                    "This dataset is a manual compilation of RNA-seq experiments on senescence in human cell lines, treated with different senescence inducers and their respective proliferative and quiescent controls. ",
                    "It has been used in Martins-Silva et al., 2026 (",
                    "<a href=\"https://doi.org/10.1093/nargab/lqag057\" target=\"_blank\">NAR Genomics and Bioinformatics</a>",
                    "). Loads instantly from the data installed with the app. The raw read counts for all samples are also permanently archived on ",
                    "<a href=\"https://zenodo.org/records/18714122\" target=\"_blank\">Zenodo</a>",
                    " for independent download or citation. Genes are in rows, samples in columns."
                  ))
                )
              ),
              
              # Help text for processed example
              shiny::conditionalPanel(
                "input.expr_source == 'example_proc'",
                shiny::helpText(
                  shiny::tags$strong("Example Processed Data:"),
                  shiny::HTML(paste0(
                    "This dataset contains the same RNA-seq experiments on senescence in human cell lines, as described for the unprocessed data above, but filtered for lowly expressed genes, normalised, and batch-corrected. ",
                    "Processing followed the methods described in Martins-Silva et al., 2026 (",
                    "<a href=\"https://doi.org/10.1093/nargab/lqag057\" target=\"_blank\">NAR Genomics and Bioinformatics</a>",
                    "). Loads instantly from the data installed with the app. The processed data is also permanently archived on ",
                    "<a href=\"https://zenodo.org/records/18714122\" target=\"_blank\">Zenodo</a>",
                    " for independent download or citation. Genes are in rows, samples in columns."
                  ))
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
                  "Upload Gene Sets"      = "upload",
                  "Paste Gene Set(s)"     = "paste"
                )
              ),
              
              shiny::conditionalPanel(
                "input.geneset_source == 'upload'",
                shiny::fileInput(
                  "geneset_files",
                  "Upload one or more files (.csv, .txt, .tsv, .rds, .rda)",
                  accept = c(".csv", ".txt", ".tsv", ".rds", ".rda"),
                  multiple = TRUE
                ),
                shiny::helpText(
                  shiny::tags$strong("Accepted input formats for gene sets:"),
                  shiny::tags$ul(
                    shiny::tags$li(
                      shiny::strong(".csv / .txt / .tsv files:"),
                      " must contain either:",
                      shiny::tags$ul(
                        shiny::tags$li("One column with gene names (undirected gene set)"),
                        shiny::tags$li("Two columns where the first column contains gene names and the second column contains the expected direction of enrichment (+1 or -1)")
                      ),
                      shiny::tags$em(
                        "Do not include a header row. Every line is read as data, starting from ",
                        "the first line (gene names only, or gene plus direction). ",
                        "The separator (comma, semicolon, tab, or plain whitespace) is detected automatically."
                      ),
                      shiny::tags$br(),
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
                  ),
                  shiny::tags$p(
                    style = "margin-top:6px;",
                    shiny::tags$strong("Multiple files: "),
                    "you can select more than one file at once. Each file becomes its own gene set, and all are collected into the same list."
                  )
                )
              ),

              shiny::conditionalPanel(
                "input.geneset_source == 'paste'",
                shiny::helpText(
                  shiny::tags$strong("Paste gene sets directly - no file needed:"),
                  shiny::tags$ul(
                    shiny::tags$li("Give the set a name, then paste one gene per line into the box."),
                    shiny::tags$li(
                      "For a ", shiny::strong("directional"), " gene set, add the expected ",
                      "direction after a semicolon on ", shiny::strong("every"), " line, e.g."
                    ),
                    shiny::tags$li(style = "list-style:none;margin-left:-4px;",
                      shiny::tags$code("TP53;up"), shiny::tags$br(), shiny::tags$code("MKI67;down"),
                      " (", shiny::tags$code("+1"), "/", shiny::tags$code("-1"), " also work)."
                    ),
                    shiny::tags$li(
                      "For an ", shiny::strong("undirected"), " gene set, list genes with no ",
                      "semicolon at all, e.g. ", shiny::tags$code("ACTB")
                    ),
                    shiny::tags$li(shiny::strong("A single gene set cannot mix the two:"),
                      " either every line has a direction, or none of them do."),
                    shiny::tags$li("Add more boxes below for additional gene sets - each is collected into the same list.")
                  )
                ),
                shiny::uiOutput("paste_gs_boxes"),
                shiny::div(style = "display:flex;gap:8px;margin:8px 0;",
                  shiny::actionButton("add_paste_gs_box", shiny::tagList(shiny::icon("plus"), " Add another gene set"),
                    class = "btn-sm btn-outline-secondary"),
                  shiny::actionButton("apply_paste_gs", shiny::tagList(shiny::icon("check"), " Use these gene sets"),
                    class = "btn-sm btn-primary")
                ),
                shiny::uiOutput("paste_gs_status")
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
}

# ---- Server -----------------------------------------------------------------

dataServer <- function(input, output, session) {

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

  ###### HELPER FUNCTIONS ######

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
        shiny::tags$p("Could not find the example data installed with the app, and the fallback download from Zenodo also failed."),
        shiny::tags$p(style = "font-size:0.84em;color:#888;font-family:monospace;word-break:break-all;", msg),
        shiny::tags$p("Check your internet connection and click ", shiny::tags$strong("Try Again"), "."),
        shiny::tags$p(
          style = "background:#f0fdf4;border:1px solid #bbf7d0;border-radius:5px;padding:8px 12px;font-size:0.88em;color:#166534;margin-top:8px;",
          shiny::icon("circle-check", style="margin-right:4px;"),
          shiny::tags$strong("Note:"), " Uploading your own files and importing from GEO ",
          "still work normally. The example data is installed together with the app and does not normally require an internet connection; ",
          "this message appears because that copy could not be found. ",
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
      return(.read_geneset_delimited(path))
    }

    stop("Unsupported gene set format: .", ext,
         ". Supported: .csv, .txt, .tsv, .rds, .rda, .xls, .xlsx")
  }

  # ── Delimited gene-set reader (.csv / .txt / .tsv) ──────────────────────────
  # These files must NOT have a header row: every line is data, starting from
  # the very first line. (Previously header=TRUE was always assumed, which
  # silently dropped the first gene whenever there was no header.)
  #   • Separator auto-detected: tab, semicolon, comma, else any whitespace
  #     (so plain space/whitespace-separated .txt files work too).
  #   • One column  → gene names only (undirected gene set).
  #   • Two columns → gene name + direction (+1 / -1). Extra columns ignored.
  .read_geneset_delimited <- function(path) {
    first_line <- readLines(path, n = 1, warn = FALSE, encoding = "UTF-8")
    first_line <- sub("^\xef\xbb\xbf", "", first_line)  # strip BOM

    sep <- if (grepl("\t", first_line)) "\t" else
           if (grepl(";",  first_line)) ";"  else
           if (grepl(",",  first_line)) ","  else ""  # "" -> any whitespace

    raw <- utils::read.table(
      path, header = FALSE, sep = sep, check.names = FALSE,
      stringsAsFactors = FALSE, comment.char = "", fill = TRUE,
      quote = "\"'", colClasses = "character", strip.white = TRUE
    )

    if (nrow(raw) == 0) stop("File appears to be empty after parsing.")

    # Strip stray quotes, drop fully-empty columns, keep at most 2 columns.
    raw[] <- lapply(raw, function(x) gsub("^[\"']|[\"']$", "", trimws(x)))
    raw   <- raw[, colSums(!is.na(raw) & raw != "") > 0, drop = FALSE]
    if (ncol(raw) == 0) stop("File appears to be empty after parsing.")
    raw <- raw[, seq_len(min(2, ncol(raw))), drop = FALSE]

    genes <- as.character(raw[[1]])
    if (ncol(raw) == 1) return(genes)

    data.frame(
      gene      = genes,
      direction = suppressWarnings(as.numeric(raw[[2]])),
      stringsAsFactors = FALSE
    )
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

  # Example data ships bundled with the package under inst/appdata (installed
  # alongside the app, so system.file() finds it in any deployment - notably
  # the Docker image). When present, it's read directly - no network request,
  # no modal delay. If a given file is missing (e.g. a dev checkout that
  # hasn't pulled it, or a build that intentionally strips inst/appdata),
  # we transparently fall back to the original Zenodo download so the app
  # still works, just slower on first load.
  .local_appdata_path <- function(filename) {
    p <- system.file("appdata", filename, package = "markeR")
    if (nzchar(p) && file.exists(p)) p else NULL
  }

  .dl_expr_raw <- function() {
    local_path <- .local_appdata_path("counts.rds")
    if (!is.null(local_path)) {
      e <- readRDS(local_path)
      expr_data(e)
      log_step(paste0("Expression data loaded: example (unprocessed) - ",
                      nrow(e), " genes × ", ncol(e), " samples."))
      return(invisible(TRUE))
    }
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
    local_path <- .local_appdata_path("corrcounts.rds")
    if (!is.null(local_path)) {
      e <- readRDS(local_path)
      expr_data(e)
      log_step(paste0("Expression data loaded: example (processed) - ",
                      nrow(e), " genes × ", ncol(e), " samples."))
      return(invisible(TRUE))
    }
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
    local_path <- .local_appdata_path("metadata.rds")
    if (!is.null(local_path)) {
      m <- readRDS(local_path)
      m <- align_meta_to_expr(m, expr_data())
      meta_data(m)
      log_step(paste0("Metadata loaded: example - ",
                      nrow(m), " samples × ", ncol(m), " variables."))
      return(invisible(TRUE))
    }
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
    local_path <- .local_appdata_path("SenescenceGeneSets.rds")
    if (!is.null(local_path)) {
      gs <- readRDS(local_path)
      gene_sets(gs)
      log_step(paste0("Gene sets loaded: example - ", length(gs), " set(s), ",
                      sum(sapply(gs, function(x) if (is.data.frame(x)) nrow(x) else length(x))),
                      " genes total."))
      return(invisible(TRUE))
    }
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
      warn_msg <- "Some columns are non-numeric. This file may not be an expression matrix. Check that genes are rows and samples are columns."
    } else if (na_frac > 0.20) {
      warn_msg <- sprintf(
        "%.0f%% of values are missing after parsing. The file may not be a numeric expression matrix, or the separator or decimal character may have been misdetected.",
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
          paste0("Expression columns do not fully match the SampleID column. ",
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
        paste0(n_expr, " samples in both, but names differ. ",
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

  ##### GENE SETS - paste (text box) input #####

  # One gene per line; "gene;direction" for directional sets (direction is
  # optional per line and accepts up/down/+1/-1/pos/neg). Reuses
  # parse_geneset_object() so pasted sets go through the exact same
  # validation as uploaded files.
  .parse_pasted_geneset_lines <- function(txt) {
    if (is.null(txt) || !nzchar(trimws(txt))) return(NULL)
    lines <- trimws(strsplit(txt, "\n")[[1]])
    lines <- lines[nzchar(lines)]
    if (length(lines) == 0L) return(NULL)

    .dir_to_num <- function(x) {
      xl <- tolower(trimws(x))
      if (!nzchar(xl)) return(NA_character_)
      if (xl %in% c("up", "+", "+1", "1", "pos", "positive")) return("1")
      if (xl %in% c("down", "-", "-1", "neg", "negative")) return("-1")
      x  # leave as-is (e.g. already numeric-ish) for as.numeric() downstream
    }

    parts   <- strsplit(lines, ";", fixed = TRUE)
    genes   <- vapply(parts, function(p) trimws(p[1]), character(1))
    has_dir <- vapply(parts, length, integer(1)) >= 2L
    keep    <- nzchar(genes)
    if (!any(keep)) return(NULL)
    genes   <- genes[keep]; has_dir <- has_dir[keep]; parts <- parts[keep]

    if (!any(has_dir)) return(genes)  # plain undirected gene set

    # A gene set is either fully directional or fully undirected - never a
    # mix, since downstream directional methods need a direction for every
    # gene. Flag it rather than silently treating the direction-less lines
    # as NA.
    if (!all(has_dir)) {
      missing_genes <- genes[!has_dir]
      stop("has both directional and non-directional lines. Add a direction to ",
           if (length(missing_genes) <= 5L) paste(missing_genes, collapse = ", ")
           else paste0(paste(missing_genes[1:5], collapse = ", "), ", … (",
                       length(missing_genes), " genes missing a direction)"),
           ", or remove directions from all lines to make it undirected.")
    }

    dirs <- vapply(parts, function(p) .dir_to_num(p[2]), character(1))
    data.frame(gene = genes, direction = dirs, stringsAsFactors = FALSE)
  }

  paste_gs_n <- shiny::reactiveVal(1L)
  observeEvent(input$add_paste_gs_box, { paste_gs_n(paste_gs_n() + 1L) })

  output$paste_gs_boxes <- shiny::renderUI({
    n <- paste_gs_n()
    shiny::tagList(lapply(seq_len(n), function(i) {
      shiny::div(style = "border:1px solid #e5e7eb;border-radius:6px;padding:8px 10px;margin-bottom:8px;",
        shiny::textInput(paste0("paste_gs_name_", i), "Gene set name:",
          value = shiny::isolate(input[[paste0("paste_gs_name_", i)]]) %||% "",
          placeholder = "e.g. MyUpregulatedSet", width = "100%"),
        shiny::textAreaInput(paste0("paste_gs_genes_", i), "Genes (one per line):",
          value = shiny::isolate(input[[paste0("paste_gs_genes_", i)]]) %||% "",
          rows = 6, placeholder = "Directional:\nTP53;up\nMKI67;down\n\nor undirected:\nACTB\nGAPDH",
          width = "100%")
      )
    }))
  })

  observeEvent(input$apply_paste_gs, {
    n <- paste_gs_n()
    compiled <- list()
    problems <- character(0)
    for (i in seq_len(n)) {
      nm  <- trimws(input[[paste0("paste_gs_name_", i)]] %||% "")
      txt <- input[[paste0("paste_gs_genes_", i)]] %||% ""
      if (!nzchar(nm) && !nzchar(trimws(txt))) next  # skip fully empty boxes
      if (!nzchar(nm)) { problems <- c(problems, paste0("Gene set #", i, " needs a name.")); next }
      had_parse_error <- FALSE
      obj <- tryCatch(
        .parse_pasted_geneset_lines(txt),
        error = function(e) {
          had_parse_error <<- TRUE
          problems <<- c(problems, paste0("'", nm, "' ", conditionMessage(e)))
          NULL
        }
      )
      if (is.null(obj)) {
        if (!had_parse_error) problems <- c(problems, paste0("'", nm, "' has no genes."))
        next
      }
      parsed <- tryCatch(
        parse_geneset_object(obj, clean_name(nm)),
        error = function(e) {
          problems <<- c(problems, paste0("'", nm, "': ", conditionMessage(e)))
          NULL
        }
      )
      if (!is.null(parsed)) compiled <- c(compiled, parsed)
    }
    if (length(problems) > 0L)
      shiny::showNotification(paste(problems, collapse = "\n"), type = "error", duration = 10)
    if (length(compiled) == 0L) return()
    gene_sets(compiled)
    log_step(paste0("Gene sets loaded: ", length(compiled), " pasted set(s) parsed."))
    output$paste_gs_status <- shiny::renderUI({
      shiny::tags$small(style = "color:#28a745;display:block;margin-top:4px;",
        shiny::icon("circle-check"),
        sprintf(" %d gene set(s) loaded: %s", length(compiled), paste(names(compiled), collapse = ", ")))
    })
  })



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
        "No data loaded yet. Use the sidebar to load expression data, metadata, and gene sets."
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
  

  ###### PUBLIC INTERFACE ######
  list(
    expr_data         = expr_data,
    meta_data         = meta_data,
    gene_sets         = gene_sets,
    log_step          = log_step,
    data_log          = data_log,
    expr_quality_warn = expr_quality_warn
  )
}
