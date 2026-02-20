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
                  "Upload (.csv, .rds, .rda)",
                  accept=c(".csv", ".rds", ".rda")
                ),
                shiny::helpText(
                  "Genes in rows, samples in columns."
                )
              ),
              
              shiny::conditionalPanel(
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
  
  ##### DATA IMPORT #####
  
  ######   EXAMPLE DATA  ######  
  
  # import example data here
  
  ######  GEO ###### 
  
  # import geo data here
  
  ######   INPUT DATA  ######  
  
  # input data
  
  
  ###### SUMMARY  ###### 
  
  # Number of Samples
  # Number of genes
  # Number of gene sets
  
  ####### PREVIEWS  ###### 
  
  
  
  
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