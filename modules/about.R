# ENABL Biomarker/Target Analysis Tool - About Module
# Provides documentation, help, and application information

# =============================================================================
# UI Function
# =============================================================================

aboutUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    div(
      class = "about-module",
      style = "padding: 20px;",
      
      # Header
      div(
        class = "module-header text-center",
        style = "background: linear-gradient(135deg, #2196F3 0%, #1976D2 100%); color: white; padding: 30px; border-radius: 10px; margin-bottom: 30px;",
        h1("ENABL Biomarker / Target Analysis Tool (EB-TAT)", style = "margin: 0 0 10px 0; font-weight: bold;"),
        h4("Version 2.0 - Modularized Architecture", style = "margin: 0; opacity: 0.9;"),
        p("A comprehensive Shiny application for transcriptomics data analysis", style = "margin: 10px 0 0 0; font-size: 16px;")
      ),
      
      # Main content
      fluidRow(
        # Quick info column
        column(
          4,
          # Application info card
          div(
            class = "card",
            style = "margin-bottom: 20px;",
            div(
              class = "card-header bg-primary text-white",
              h4("Application Information", style = "margin: 0;")
            ),
            div(
              class = "card-body",
              tags$table(
                class = "table table-borderless",
                tags$tr(
                  tags$td(strong("Version:")),
                  tags$td("2.0.0")
                ),
                tags$tr(
                  tags$td(strong("Release Date:")),
                  tags$td("2024")
                ),
                tags$tr(
                  tags$td(strong("Created By:")),
                  tags$td("Vijay Singh (GIS)")
                ),
                tags$tr(
                  tags$td(strong("First Version:")),
                  tags$td("2021")
                ),
                tags$tr(
                  tags$td(strong("Maintainer:")),
                  tags$td("ENABL Bioinformatics Team")
                )
              )
            )
          ),
          
          # Quick links card
          div(
            class = "card",
            style = "margin-bottom: 20px;",
            div(
              class = "card-header bg-info text-white",
              h4("Quick Links", style = "margin: 0;")
            ),
            div(
              class = "card-body",
              tags$ul(
                class = "list-unstyled",
                tags$li(style = "margin-bottom: 10px;", 
                       tags$a(href = "#data-input", "Data Input Guide", 
                              class = "btn btn-outline-primary btn-sm btn-block")),
                tags$li(style = "margin-bottom: 10px;", 
                       tags$a(href = "#volcano-plot", "Volcano Plot Tutorial", 
                              class = "btn btn-outline-primary btn-sm btn-block")),
                tags$li(style = "margin-bottom: 10px;", 
                       tags$a(href = "#enrichment", "Enrichment Analysis Help", 
                              class = "btn btn-outline-primary btn-sm btn-block")),
                tags$li(style = "margin-bottom: 10px;", 
                       actionButton(ns("show_shortcuts"), "Keyboard Shortcuts", 
                                   class = "btn btn-outline-primary btn-sm btn-block"))
              )
            )
          ),
          
          # System info card
          div(
            class = "card",
            div(
              class = "card-header bg-success text-white",
              h4("System Information", style = "margin: 0;")
            ),
            div(
              class = "card-body",
              verbatimTextOutput(ns("system_info")),
              actionButton(ns("refresh_system_info"), "Refresh", 
                          class = "btn btn-outline-success btn-sm")
            )
          )
        ),
        
        # Main content column
        column(
          8,
          # Features card
          div(
            class = "card",
            style = "margin-bottom: 20px;",
            div(
              class = "card-header bg-warning text-dark",
              h4("Key Features", style = "margin: 0;")
            ),
            div(
              class = "card-body",
              fluidRow(
                column(
                  6,
                  tags$ul(
                    class = "feature-list",
                    tags$li(tags$i(class = "fas fa-chart-line text-primary"), 
                           "Interactive volcano plots"),
                    tags$li(tags$i(class = "fas fa-database text-primary"), 
                           "Multiple data input methods"),
                    tags$li(tags$i(class = "fas fa-project-diagram text-primary"), 
                           "Gene set enrichment analysis"),
                    tags$li(tags$i(class = "fas fa-braille text-primary"), 
                           "Expression visualization")
                  )
                ),
                column(
                  6,
                  tags$ul(
                    class = "feature-list",
                    tags$li(tags$i(class = "fas fa-network-wired text-primary"), 
                           "Protein-protein interaction networks"),
                    tags$li(tags$i(class = "fas fa-filter text-primary"), 
                           "Advanced filtering options"),
                    tags$li(tags$i(class = "fas fa-download text-primary"), 
                           "Export functionality"),
                    tags$li(tags$i(class = "fas fa-cogs text-primary"), 
                           "Modular architecture")
                  )
                )
              )
            )
          ),
          
          # Data sources card
          div(
            class = "card",
            style = "margin-bottom: 20px;",
            div(
              class = "card-header bg-danger text-white",
              h4("Data Sources & Integration", style = "margin: 0;")
            ),
            div(
              class = "card-body",
              tags$ul(
                class = "list-group",
                tags$li(
                  class = "list-group-item d-flex justify-content-between align-items-center",
                  "ENABL Transcriptomics Data",
                  tags$span(class = "badge bg-primary rounded-pill", "Primary")
                ),
                tags$li(
                  class = "list-group-item d-flex justify-content-between align-items-center",
                  "Human Protein Atlas",
                  tags$span(class = "badge bg-success rounded-pill", "Integrated")
                ),
                tags$li(
                  class = "list-group-item d-flex justify-content-between align-items-center",
                  "STRING Database",
                  tags$span(class = "badge bg-success rounded-pill", "Integrated")
                ),
                tags$li(
                  class = "list-group-item d-flex justify-content-between align-items-center",
                  "Gene Ontology (GO) Databases",
                  tags$span(class = "badge bg-success rounded-pill", "Integrated")
                ),
                tags$li(
                  class = "list-group-item d-flex justify-content-between align-items-center",
                  "KEGG Pathway Databases",
                  tags$span(class = "badge bg-success rounded-pill", "Integrated")
                ),
                tags$li(
                  class = "list-group-item d-flex justify-content-between align-items-center",
                  "Reactome Pathways",
                  tags$span(class = "badge bg-success rounded-pill", "Integrated")
                )
              )
            )
          ),
          
          # Getting started card
          div(
            class = "card",
            style = "margin-bottom: 20px;",
            div(
              class = "card-header bg-secondary text-white",
              h4("Getting Started", style = "margin: 0;")
            ),
            div(
              class = "card-body",
              tags$ol(
                style = "padding-left: 20px;",
                tags$li(strong("Data Input:"), "Load your differential expression data using the Data tab"),
                tags$li(strong("Visualization:"), "Explore patterns using Volcano Plot and Expression Plot tabs"),
                tags$li(strong("Analysis:"), "Use Analytics and Enrichment tabs for deeper insights"),
                tags$li(strong("Export:"), "Download plots and results for publications")
              ),
              hr(),
              p("For detailed documentation, please refer to the module descriptions and tooltips throughout the application.")
            )
          ),
          
          # Technical information card
          div(
            class = "card",
            div(
              class = "card-header bg-dark text-white",
              h4("Technical Information", style = "margin: 0;")
            ),
            div(
              class = "card-body",
              p("EB-TAT is built using R and Shiny with a modular architecture for maintainability and extensibility. The application utilizes several specialized R packages for bioinformatics analysis and visualization."),
              tags$h5("Core Technologies:"),
              tags$ul(
                tags$li(strong("Frontend:"), "Shiny, HTML5, CSS3, JavaScript"),
                tags$li(strong("Backend:"), "R with specialized bioinformatics packages"),
                tags$li(strong("Visualization:"), "ggplot2, pheatmap, enrichplot"),
                tags$li(strong("Analysis:"), "DESeq2, clusterProfiler, DOSE")
              ),
              tags$h5("Architecture:"),
              p("The application follows a modular design pattern with separate UI and server components for each functional area, enabling easy maintenance and future enhancements.")
            )
          )
        )
      ),
      
      # Footer
      div(
        class = "text-center",
        style = "margin-top: 30px; padding: 20px; background-color: #f8f9fa; border-radius: 10px;",
        p(
          "For questions, support, or feature requests, please contact:",
          tags$br(),
          tags$strong("Vijay Singh"), 
          tags$br(),
          tags$a(href = "mailto:vijay.s.gautam@gmail.com", "vijay.s.gautam@gmail.com"),
          style = "margin: 0;"
        ),
        p(
          "ENABL Consortium © 2024",
          style = "margin: 10px 0 0 0; font-size: 12px; color: #666;"
        )
      )
    )
  )
}

# =============================================================================
# Server Function
# =============================================================================

aboutServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Reactive values for this module
    module_rv <- reactiveValues(
      system_info = NULL,
      last_update = Sys.time()
    )
    
    # Get system information
    get_system_info <- function() {
      info <- list()
      
      # R version and platform
      info$r_version <- R.version.string
      info$platform <- R.version$platform
      
      # Shiny version
      info$shiny_version <- as.character(packageVersion("shiny"))
      
      # Key package versions
      key_packages <- c("ggplot2", "dplyr", "DT", "clusterProfiler", "DESeq2")
      package_versions <- sapply(key_packages, function(pkg) {
        if (requireNamespace(pkg, quietly = TRUE)) {
          as.character(packageVersion(pkg))
        } else {
          "Not available"
        }
      })
      info$package_versions <- package_versions
      
      # System locale
      info$locale <- Sys.getlocale()
      
      # Working directory
      info$working_dir <- getwd()
      
      # Memory usage
      info$memory <- paste0(round(pryr::mem_used() / 1024^2, 1), " MB")
      
      # Session info
      info$session_time <- Sys.time()
      
      return(info)
    }
    
    # Format system information for display
    format_system_info <- function(info) {
      if (is.null(info)) return("Loading system information...")
      
      text <- c()
      text <- c(text, "R VERSION:")
      text <- c(text, paste(" ", info$r_version))
      text <- c(text, paste(" ", info$platform))
      text <- c(text, "")
      
      text <- c(text, "KEY PACKAGES:")
      for (pkg in names(info$package_versions)) {
        text <- c(text, paste(" ", pkg, ":", info$package_versions[pkg]))
      }
      text <- c(text, "")
      
      text <- c(text, "SYSTEM INFO:")
      text <- c(text, paste(" Shiny:", info$shiny_version))
      text <- c(text, paste(" Locale:", info$locale))
      text <- c(text, paste(" Memory:", info$memory))
      text <- c(text, paste(" Working dir:", info$working_dir))
      text <- c(text, paste(" Last update:", format(info$session_time)))
      
      return(paste(text, collapse = "\n"))
    }
    
    # Update system information
    update_system_info <- function() {
      module_rv$system_info <- get_system_info()
      module_rv$last_update <- Sys.time()
    }
    
    # Initial system info load
    update_system_info()
    
    # Event observers
    observeEvent(input$refresh_system_info, {
      showNotification("Refreshing system information...", type = "message")
      update_system_info()
      showNotification("System information updated", type = "message")
    })
    
    observeEvent(input$show_shortcuts, {
      showModal(modalDialog(
        title = "Keyboard Shortcuts & Tips",
        size = "l",
        easyClose = TRUE,
        footer = modalButton("Close"),
        tags$div(
          class = "container-fluid",
          tags$h4("Navigation Shortcuts"),
          tags$ul(
            tags$li(strong("Ctrl + 1:"), "Data Input tab"),
            tags$li(strong("Ctrl + 2:"), "Volcano Plot tab"),
            tags$li(strong("Ctrl + 3:"), "Analytics tab"),
            tags$li(strong("Ctrl + 4:"), "Enrichment tab"),
            tags$li(strong("Ctrl + 5:"), "Expression Plot tab")
          ),
          tags$hr(),
          tags$h4("General Tips"),
          tags$ul(
            tags$li("Hover over controls to see tooltips with additional information"),
            tags$li("Click on points in plots to select genes"),
            tags$li("Use the filter controls to focus on significant results"),
            tags$li("Download buttons are available in each module for results export"),
            tags$li("The application automatically saves your settings during the session")
          ),
          tags$hr(),
          tags$h4("Troubleshooting"),
          tags$ul(
            tags$li("If plots don't render, check that data is loaded in the Data Input tab"),
            tags$li("For enrichment analysis, ensure you have at least 5 significant genes"),
            tags$li("Large datasets may take longer to process - watch for loading indicators"),
            tags$li("Check the System Information panel for package compatibility issues")
          )
        )
      ))
    })
    
    # Outputs
    output$system_info <- renderText({
      format_system_info(module_rv$system_info)
    })
    
    # Return reactive values (though this module doesn't need to return much)
    return(reactive({
      list(
        system_info = module_rv$system_info,
        last_update = module_rv$last_update
      )
    }))
  })
}

# =============================================================================
# CSS Styles (included in the app.R header)
# =============================================================================

# Note: These styles should be added to the app.R header
about_module_styles <- HTML("
<style>
.about-module .feature-list {
  list-style: none;
  padding-left: 0;
}

.about-module .feature-list li {
  margin-bottom: 12px;
  padding-left: 30px;
  position: relative;
}

.about-module .feature-list li i {
  position: absolute;
  left: 0;
  top: 2px;
  font-size: 16px;
}

.about-module .card-header {
  border-radius: 5px 5px 0 0 !important;
}

.about-module .module-header {
  box-shadow: 0 4px 6px rgba(0,0,0,0.1);
}

.about-module .btn-block {
  text-align: left;
  padding: 8px 12px;
}

.about-module .list-group-item {
  border: 1px solid rgba(0,0,0,.125);
  margin-bottom: 5px;
  border-radius: 5px;
}

.about-module .badge {
  font-size: 10px;
}
</style>
")