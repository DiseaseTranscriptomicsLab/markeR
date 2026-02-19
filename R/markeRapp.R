#' @importFrom bslib bs_theme
#' @import shiny
#' @import GEOquery

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
                  "Use Example Data" = "example",
                  "Upload File"      = "upload",
                  "Retrieve from GEO"= "geo"
                )
              ),
              
              shiny::conditionalPanel(
                "input.expr_source == 'upload'",
                shiny::fileInput(
                  "expr_file",
                  "Upload (.csv, .rds, .rda)"
                ),
                shiny::helpText(
                  "CSV: genes in rows, samples in columns."
                )
              ),
              
              shiny::conditionalPanel(
                "input.expr_source == 'geo'",
                shiny::textInput(
                  "geo_accession",
                  "Enter GEO Accession (e.g., GSE12345)"
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
                  "Upload (.csv, .rds)"
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
                  "Upload one or more files",
                  multiple = TRUE
                ),
                shiny::helpText(
                  "Vector OR data.frame: col1=gene, col2=direction (+1/-1)"
                )
              )
            )
          )
        ),
        
 
        ###### MAIN PREVIEW AREA ###### 
 
        
        bslib::card(
          bslib::card_header("Data Preview"),
          
          shiny::uiOutput("summary_ui"),
          
          shiny::hr(),
          
          shiny::h5("Expression Preview"),
          shiny::tableOutput("expr_preview"),
          
          shiny::h5("Metadata Preview"),
          shiny::tableOutput("meta_preview"),
          
          shiny::h5("Gene Sets Overview"),
          shiny::uiOutput("geneset_summary")
        )
      )
    )
  ),
  
  ##### PREPROCESSING #####
  bslib::nav_panel(
    "Preprocessing",
    shiny::h4("About markeR"),
    shiny::p("This app provides a graphical interface to evaluate gene sets as phenotype markers in R.")
  ),
  
  ##### GENE SETS #####
  bslib::nav_panel(
    "Gene Sets",
    shiny::h4("About markeR"),
    shiny::p("This app provides a graphical interface to evaluate gene sets as phenotype markers in R.")
  ),
  
  ##### BENCHMARKING MODE #####
  bslib::nav_panel(
    "Benchmarking Mode",
    shiny::h4("About markeR"),
    shiny::p("This app provides a graphical interface to evaluate gene sets as phenotype markers in R.")
  ),
  
  ##### DISCOVERY MODE #####
  bslib::nav_panel(
    "Discovery Mode",
    shiny::h4("About markeR"),
    shiny::p("This app provides a graphical interface to evaluate gene sets as phenotype markers in R.")
  ) 
  
  
  
)

# ---- Server ----
server <- function(input, output, session) {
  
  ####################################################
  # EXAMPLE DATA
  ####################################################
  
  example_expr <- matrix(
    rnorm(1000),
    nrow = 100,
    dimnames = list(
      paste0("Gene", 1:100),
      paste0("Sample", 1:10)
    )
  )
  
  example_meta <- data.frame(
    Sample = paste0("Sample", 1:10),
    Group  = rep(c("Control", "Senescent"), each = 5)
  )
  
  example_genesets <- list(
    Senescence_Up = data.frame(
      gene = paste0("Gene", 1:10),
      direction = 1
    ),
    Senescence_Down = data.frame(
      gene = paste0("Gene", 20:30),
      direction = -1
    )
  )
  
  ####################################################
  # GEO STORAGE (SIMULATED)
  ####################################################
  geo_objects <- shiny::reactiveVal(NULL)
  
  shiny::observeEvent(input$geo_fetch, {
    req(input$geo_accession)
    
    shiny::withProgress(message = "Downloading GEO data...", {
      gse_raw <- GEOquery::getGEO(input$geo_accession, GSEMatrix = TRUE, getGPL = FALSE)
      
      # Wrap single ExpressionSet into a list if needed
      gse_list <- if (inherits(gse_raw, "ExpressionSet")) list(gse_raw) else gse_raw
      
      geo_objects(gse_list)
    })
  })
  
  
  ####################################################
  # EXPRESSION REACTIVE
  ####################################################
  
  expr_data <- shiny::reactive({
    
    if (input$expr_source == "example")
      return(example_expr)
    
    if (input$expr_source == "upload" &&
        !is.null(input$expr_file)) {
      
      ext <- tools::file_ext(input$expr_file$name)
      
      if (ext == "csv")
        return(as.matrix(utils::read.csv(
          input$expr_file$datapath,
          row.names = 1
        )))
      
      if (ext == "rds")
        return(readRDS(input$expr_file$datapath))
    }
    
    if (input$expr_source == "geo") {
      req(geo_objects())
      return(geo_objects()[[input$geo_expr_choice]])
    }
    
    NULL
  })
  
  ####################################################
  # METADATA REACTIVE
  ####################################################
  
  meta_data <- shiny::reactive({
    
    if (input$meta_source == "example")
      return(example_meta)
    
    if (input$meta_source == "upload" &&
        !is.null(input$meta_file)) {
      
      ext <- tools::file_ext(input$meta_file$name)
      
      if (ext == "csv")
        return(utils::read.csv(input$meta_file$datapath))
      
      if (ext == "rds")
        return(readRDS(input$meta_file$datapath))
    }
    
    if (input$meta_source == "geo") {
      req(geo_objects())
      return(geo_objects()[[input$geo_meta_choice]])
    }
    
    NULL
  })
  
  ####################################################
  # GENE SETS REACTIVE
  ####################################################
  
  gene_sets <- shiny::reactive({
    
    if (input$geneset_source == "example")
      return(example_genesets)
    
    if (input$geneset_source == "upload" &&
        !is.null(input$geneset_files)) {
      
      files <- input$geneset_files
      out <- list()
      
      for (i in seq_len(nrow(files))) {
        
        ext <- tools::file_ext(files$name[i])
        name <- tools::file_path_sans_ext(files$name[i])
        
        if (ext == "csv")
          df <- utils::read.csv(files$datapath[i])
        else if (ext == "rds")
          df <- readRDS(files$datapath[i])
        else next
        
        out[[name]] <- df
      }
      
      return(out)
    }
    
    NULL
  })
  
  ####################################################
  # SUMMARY
  ####################################################
  
  output$summary_ui <- shiny::renderUI({
    
    shiny::tagList(
      
      if (!is.null(expr_data()))
        shiny::p(
          shiny::strong("Expression: "),
          nrow(expr_data()), " genes | ",
          ncol(expr_data()), " samples"
        ),
      
      if (!is.null(meta_data()))
        shiny::p(
          shiny::strong("Metadata: "),
          nrow(meta_data()), " samples | ",
          ncol(meta_data()), " variables"
        ),
      
      if (!is.null(gene_sets()))
        shiny::p(
          shiny::strong("Gene Sets: "),
          length(gene_sets()), " sets"
        )
    )
  })
  
  ####################################################
  # PREVIEWS
  ####################################################
  
  output$expr_preview <- shiny::renderTable({
    req(expr_data())
    head(expr_data()[, 1:min(5, ncol(expr_data()))])
  }, rownames = TRUE)
  
  output$meta_preview <- shiny::renderTable({
    req(meta_data())
    head(meta_data())
  })
  
  output$geneset_summary <- shiny::renderUI({
    
    gs <- gene_sets()
    if (is.null(gs)) return(NULL)
    
    shiny::tagList(
      lapply(names(gs), function(name) {
        
        df <- gs[[name]]
        total <- nrow(df)
        pos <- sum(df[,2] == 1, na.rm = TRUE)
        neg <- sum(df[,2] == -1, na.rm = TRUE)
        
        shiny::div(
          style = "margin-bottom:15px;",
          shiny::h6(name),
          shiny::p("Genes: ", total),
          shiny::p("+1: ", pos, " | -1: ", neg),
          shiny::tableOutput(NULL)
        )
      })
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
  
  
  app <- shiny::shinyApp(ui, server)
  shiny::runApp(app, ...)
}