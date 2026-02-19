#' @importFrom bslib bs_theme
#' @import shiny

ui <- bslib::page_navbar(
  
  # ---- Title with inline version ----
  title = shiny::tags$span(
    shiny::a("markeR "),
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
             style = "font-size: 0.6em; color: #949494; text-decoration: none;")
  ),
  
  theme = bslib::bs_theme(bootswatch = "sketchy"),
  
  # ---- Tabs ----
  
  bslib::nav_panel(
    "About",
    
    # shiny::tags$div(
    #   style = "text-align: center; max-width: 800px; margin: 0 auto;",
    #   shiny::h3("Welcome to markeR!"),
    #   shiny::p(
    #     "Have a collection of genes that putatively mark a phenotype, but not sure how to combine them and turn them into meaningful metrics for phenotype quantification? 
    #  markeR helps you evaluate gene sets as phenotype markers and explore their information.
    #  Here’s a quick guide to each section of the app:",
    #     style = "font-size: 1.2em;"
    #   )
    # ),
    
    
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
        style = "flex: 0 0 200px; margin-right: 15px;",  # fixed 200px width
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
          #"Have a collection of genes that putatively mark a phenotype, but not sure how to combine them and turn them into meaningful metrics for phenotype quantification? 
       #markeR helps you evaluate gene sets as phenotype markers and explore their information.
       #Here’s a quick guide to each section of the app:",
          
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
        flex: 1 1 400px;    /* grow/shrink, minimum width 400px */
        padding-right: 20px;
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
        flex: 1 1 800px;   /* grow/shrink, minimum width 300px */
        display: flex;
        justify-content: center;
        align-items: flex-start;
      ",
        shiny::img(
          src = "methods/methods.png",
          style = "max-width: 100%; height: auto; margin-top: 10px;"
        )
      )
    )
    
  ),
  
  
  bslib::nav_panel(
    "Data",
    shiny::h4("About markeR"),
    shiny::p("This app provides a graphical interface to evaluate gene sets as phenotype markers in R.")
  ),
  
  bslib::nav_panel(
    "Preprocessing Data",
    shiny::h4("About markeR"),
    shiny::p("This app provides a graphical interface to evaluate gene sets as phenotype markers in R.")
  ),
  
  bslib::nav_panel(
    "Gene Sets",
    shiny::h4("About markeR"),
    shiny::p("This app provides a graphical interface to evaluate gene sets as phenotype markers in R.")
  ),
  
  bslib::nav_panel(
    "Benchmarking Mode",
    shiny::h4("About markeR"),
    shiny::p("This app provides a graphical interface to evaluate gene sets as phenotype markers in R.")
  ),
  
  
  bslib::nav_panel(
    "Discovery Mode",
    shiny::h4("About markeR"),
    shiny::p("This app provides a graphical interface to evaluate gene sets as phenotype markers in R.")
  ) 
  
  
  
)

# ---- Server ----
server <- function(input, output) {
  
  output$distPlot <- shiny::renderPlot({
    x    <- faithful[, 2]
    bins <- seq(min(x), max(x), length.out = input$bins + 1)
    
    hist(x, breaks = bins, col = 'darkgray', border = 'white',
         xlab = 'Waiting time to next eruption (in mins)',
         main = 'Histogram of waiting times')
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