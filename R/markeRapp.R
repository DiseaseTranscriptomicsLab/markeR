#' @importFrom bslib bs_theme
#' @import shiny
#' @import GEOquery
#' @import DT

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
                  "Upload (.csv, .rds, .rda)",
                  accept=c(".csv", ".rds", ".rda")
                ),
                shiny::helpText(
                  "Genes in rows, samples in columns."
                )
              ),
              
              shiny::conditionalPanel( # id data comes from geo, user needs to have a section to then choose which file is the expression data and (optionally) which one is the metadata
                "input.expr_source == 'geo'",
                shiny::textInput(
                  "geo_accession",
                  "Enter GEO Accession (e.g., GSE130727)", #GSE130727
                  placeholder = "GSE130727"
                ),
                shiny::actionButton(
                  "geo_fetch",
                  "Download GEO"
                ),
                shiny::hr(),
                shiny::uiOutput("geo_object_selector")
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
                     DT::DTOutput("geneset_summary")
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
  
  ##### REACTIVE STORAGE CONTAINERS ##################################
  expr_data   <- shiny::reactiveVal(NULL)  # Expression matrix/data.frame
  meta_data   <- shiny::reactiveVal(NULL)  # Sample metadata
  gene_sets   <- shiny::reactiveVal(NULL)  # Named list of gene sets
  geo_objects <- shiny::reactiveVal(NULL)  # Raw GEO objects (can be multiple)
  
  ##### HELPER FUNCTIONS (NOT REACTIVE) ##############################
  clean_name <- function(x) tools::file_path_sans_ext(base::basename(x))
  
  safe_read_table <- function(path) {
    ext <- tolower(tools::file_ext(path))
    if (ext == "csv") return(utils::read.csv(path, row.names = 1, check.names = FALSE))
    if (ext == "rds") return(base::readRDS(path))
    if (ext == "rda") {
      e <- new.env()
      base::load(path, envir = e)
      return(e[[base::ls(e)[1]]])
    }
    stop("Unsupported file format.")
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
  
  ##### EXAMPLE DATA LOADING ######################################### 
  
  # Solution with "fake" progress bar for better UI experience 
  
  observeEvent(input$expr_source, {
    
    # ---- UNPROCESSED DATA ----
    if (input$expr_source == "example_raw") {
      
      shiny::withProgress(
        message = "Downloading unprocessed example data...",
        value = 0,
        {
          
          incProgress(0.1, detail = "Preparing download...")
          Sys.sleep(0.2)
          
          tmp <- tempfile(fileext = ".rds")
          
          incProgress(0.4, detail = "Downloading from Zenodo...")
          utils::download.file(
            url  = "https://zenodo.org/records/18714122/files/counts.rds?download=1",
            destfile = tmp,
            mode = "wb"
          )
          
          incProgress(0.3, detail = "Loading expression matrix...")
          expr_data(readRDS(tmp))
          
          incProgress(0.2, detail = "Finalizing...")
        }
      )
      
      showNotification(
        "Unprocessed expression data loaded successfully.",
        type = "default",
        duration = 3
      )
    }
    
    
    # ---- PROCESSED DATA ----
    if (input$expr_source == "example_proc") {
      
      shiny::withProgress(
        message = "Downloading processed example data...",
        value = 0,
        {
          
          incProgress(0.1, detail = "Preparing download...")
          Sys.sleep(0.2)
          
          tmp <- tempfile(fileext = ".rds")
          
          incProgress(0.4, detail = "Downloading from Zenodo...")
          utils::download.file(
            url  = "https://zenodo.org/records/18714122/files/corrcounts.rds?download=1",
            destfile = tmp,
            mode = "wb"
          )
          
          incProgress(0.3, detail = "Loading processed matrix...")
          expr_data(readRDS(tmp))
          
          incProgress(0.2, detail = "Finalizing...")
        }
      )
      
      showNotification(
        "Processed expression data loaded successfully.",
        type = "default",
        duration = 3
      )
    }
  })
   
  observeEvent(input$meta_source, {
    
    if (input$meta_source == "example") {
      
      shiny::withProgress(
        message = "Downloading example metadata...",
        value = 0,
        {
          
          incProgress(0.2, detail = "Preparing download...")
          Sys.sleep(0.2)
          
          tmp <- tempfile(fileext = ".rds")
          
          incProgress(0.5, detail = "Downloading from Zenodo...")
          utils::download.file(
            url  = "https://zenodo.org/records/18714122/files/metadata.rds?download=1",
            destfile = tmp,
            mode = "wb"
          )
          
          incProgress(0.2, detail = "Loading metadata...")
          meta_data(readRDS(tmp))
          
          incProgress(0.1, detail = "Finalizing...")
        }
      )
      
      showNotification(
        "Metadata loaded successfully.",
        type = "default",
        duration = 3
      )
    }
  })
   
  observeEvent(input$geneset_source, {
    
    if (input$geneset_source == "example") {
      
      shiny::withProgress(
        message = "Downloading example gene sets...",
        value = 0,
        {
          
          incProgress(0.2, detail = "Preparing download...")
          Sys.sleep(0.2)
          
          tmp <- tempfile(fileext = ".rds")
          
          incProgress(0.5, detail = "Downloading from Zenodo...")
          utils::download.file(
            url  = "https://zenodo.org/records/18714122/files/SenescenceGeneSets.rds?download=1",
            destfile = tmp,
            mode = "wb"
          )
          
          incProgress(0.2, detail = "Loading gene sets...")
          gene_sets(readRDS(tmp))
          
          incProgress(0.1, detail = "Finalizing...")
        }
      )
      
      showNotification(
        "Gene sets loaded successfully.",
        type = "default",
        duration = 3
      )
    }
  })
  # observeEvent(input$expr_source, {
  #   if (input$expr_source == "example_raw") {
  #     shiny::withProgress(message = "Downloading unprocessed example data...", value = 0.7, {
  #       tmp <- tempfile(fileext = ".rds")
  #       utils::download.file(
  #         url  = "https://zenodo.org/records/18714122/files/counts.rds?download=1",
  #         destfile = tmp, mode = "wb"
  #       )
  #       expr_data(readRDS(tmp))
  #     })
  #   }
  #   
  #   if (input$expr_source == "example_proc") {
  #     shiny::withProgress(message = "Downloading processed example data...", value = 0.7, {
  #       tmp <- tempfile(fileext = ".rds")
  #       utils::download.file(
  #         url  = "https://zenodo.org/records/18714122/files/corrcounts.rds?download=1",
  #         destfile = tmp, mode = "wb"
  #       )
  #       expr_data(readRDS(tmp))
  #     })
  #   }
  # })
  # observeEvent(input$meta_source, {
  #   if (input$meta_source == "example") {
  #     shiny::withProgress(message = "Downloading example metadata...", value = 0.7, {
  #       tmp <- tempfile(fileext = ".rds")
  #       utils::download.file(
  #         url = "https://zenodo.org/records/18714122/files/metadata.rds?download=1",
  #         destfile = tmp, mode = "wb"
  #       )
  #       meta_data(readRDS(tmp))
  #     })
  #   }
  # })
  # observeEvent(input$geneset_source, {
  #   if (input$geneset_source == "example") {
  #     shiny::withProgress(message = "Downloading example gene sets...", value = 0.7, {
  #       tmp <- tempfile(fileext = ".rds")
  #       utils::download.file(
  #         url = "https://zenodo.org/records/18714122/files/SenescenceGeneSets.rds?download=1",
  #         destfile = tmp, mode = "wb"
  #       )
  #       gene_sets(readRDS(tmp))
  #     })
  #   }
  # })
  
  ##### FILE UPLOAD HANDLERS #########################################
  
  observeEvent(input$expr_file, {
    req(input$expr_file)
    expr_data(safe_read_table(input$expr_file$datapath))
  })
  
  observeEvent(input$meta_file, {
    req(input$meta_file)
    meta_data(safe_read_table(input$meta_file$datapath))
  })
  
  observeEvent(input$geneset_files, {
    req(input$geneset_files)
    compiled <- list()
    for (i in seq_len(nrow(input$geneset_files))) {
      path <- input$geneset_files$datapath[i]
      name <- clean_name(input$geneset_files$name[i])
      obj <- safe_read_table(path)
      parsed <- parse_geneset_object(obj, name)
      compiled <- c(compiled, parsed)
    }
    gene_sets(compiled)
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
  
  ##### GEO DOWNLOAD AND SELECTION ###################################
  observeEvent(input$geo_fetch, {
    req(input$geo_accession)
    tryCatch({
      g <- GEOquery::getGEO(input$geo_accession, GSEMatrix = TRUE)
      geo_objects(g)
    }, error = function(e) showNotification("Failed to download GEO dataset.", type = "error"))
  })
  
  output$geo_object_selector <- renderUI({
    req(geo_objects())
    selectInput("geo_selected_object", "Select GEO object:", choices = seq_along(geo_objects()))
  })
  
  observeEvent(input$geo_selected_object, {
    req(geo_objects())
    idx <- as.numeric(input$geo_selected_object)
    obj <- geo_objects()[[idx]]
    expr_data(Biobase::exprs(obj))
    meta_data(Biobase::pData(obj))
  })
  
  ##### DATA SUMMARY (REACTIVE OUTPUT) ###############################
  output$summary_ui <- renderUI({
    expr <- expr_data()
    meta <- meta_data()
    gs   <- gene_sets()
    tagList(
      if (!is.null(expr)) p(strong("Expression data: "), nrow(expr), " genes × ", ncol(expr), " samples"),
      if (!is.null(meta)) p(strong("Metadata: "), nrow(meta), " samples × ", ncol(meta), " variables"),
      if (!is.null(gs))   p(strong("Gene sets loaded: "), length(gs))
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

  
 shiny::observeEvent(input$geneset_summary_rows_selected, {
   shiny::req(input$geneset_summary_rows_selected)
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
               columnDefs = list(list(className = 'dt-center', targets = "_all"))  # center all columns
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