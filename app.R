# ENABL Biomarker/Target Analysis Tool (EB-TAT)
# Main application file - Updated with Literature Trends module

# Load required packages and utilities
if (!require("shiny")) install.packages("shiny")
if (!require("shinythemes")) install.packages("shinythemes")

library(shiny)
library(shinythemes)

# Check if global.R exists and source it
if (file.exists("global.R")) {
  source("global.R", local = TRUE)
} else {
  stop("global.R file not found. Please run install.sh to set up the application.")
}

# Improved source_module function with better error handling and flexible naming
source_module <- function(file) {
  if (!file.exists(file)) {
    warning("Module file not found: ", file)
    # Return a placeholder function
    return(function(id) {
      ns <- NS(id)
      div(
        class = "module-error",
        h3("Module Not Available"),
        p("The ", basename(file), " module could not be loaded."),
        p("Please check that all required files are present.")
      )
    })
  }

  tryCatch({
    source(file, local = TRUE)
    
    # Try multiple possible UI function names
    base_name <- tools::file_path_sans_ext(basename(file))
    
    # Possible UI function name patterns
    possible_ui_names <- c(
      # CamelCase version (data_input -> dataInputUI)
      paste0(gsub("_(.)", "\\U\\1", base_name, perl = TRUE), "UI"),
      # Original snake_case with UI suffix
      paste0(base_name, "UI"),
      # Direct mapping for consistency
      "dataInputUI", "volcanoPlotUI", "analyticsUI", "enrichmentUI", 
      "expressionPlotUI", "literatureTrendsUI", "aboutUI"
    )
    
    # Find which UI function exists
    ui_function <- NULL
    for (ui_name in possible_ui_names) {
      if (exists(ui_name)) {
        ui_function <- get(ui_name)
        message("Found UI function: ", ui_name)
        break
      }
    }
    
    if (!is.null(ui_function)) {
      return(ui_function)
    } else {
      warning("No UI function found for module: ", base_name)
      warning("Tried: ", paste(possible_ui_names, collapse = ", "))
      return(function(id) {
        ns <- NS(id)
        div(
          class = "module-error",
          h3("Module UI Function Missing"),
          p("The UI function for ", basename(file), " was not found."),
          p("Tried: ", paste(possible_ui_names, collapse = ", "))
        )
      })
    }
    
  }, error = function(e) {
    warning("Error loading module ", file, ": ", e$message)
    # Return a placeholder function
    return(function(id) {
      ns <- NS(id)
      div(
        class = "module-error",
        h3("Module Load Error"),
        p("Error loading ", basename(file), " module:"),
        p(style = "font-family: monospace; font-size: 12px;", e$message),
        p("Please check the file for syntax errors.")
      )
    })
  })
}

# Source utility files with error handling
source_utility <- function(file) {
  if (file.exists(file)) {
    tryCatch({
      source(file, local = TRUE)
      message("Successfully sourced: ", file)
    }, error = function(e) {
      warning("Error sourcing utility file ", file, ": ", e$message)
    })
  } else {
    warning("Utility file not found: ", file)
  }
}

# Source utility files
utility_files <- c(
  "utils/helpers.R",
  "utils/data_loading.R",
  "utils/plotting.R"
)

for (file in utility_files) {
  source_utility(file)
}

# Load modules
module_files <- c(
  "modules/data_input.R",
  "modules/volcano_plot.R",
  "modules/analytics.R",
  "modules/enrichment.R",
  "modules/expression_plot.R",
  "modules/literature_trends.R",  # New literature trends module
  "modules/about.R"
)

# Load module UI functions
module_ui_functions <- list()
for (i in seq_along(module_files)) {
  module_name <- tools::file_path_sans_ext(basename(module_files[i]))
  module_ui_functions[[module_name]] <- source_module(module_files[i])
}

# Define UI
ui <- navbarPage(
  title = div(
    style = "display: flex; align-items: center;",
    img(src = "data:image/svg+xml;base64,PHN2ZyB3aWR0aD0iMzAiIGhlaWdodD0iMzAiIHZpZXdCb3g9IjAgMCAzMCAzMCIgZmlsbD0ibm9uZSIgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnIj4KPHJlY3Qgd2lkdGg9IjMwIiBoZWlnaHQ9IjMwIiBmaWxsPSIjMjE5NkYzIi8+Cjx0ZXh0IHg9IjgiIHk9IjIwIiBmaWxsPSJ3aGl0ZSIgZm9udC1mYW1pbHk9IkFyaWFsIiBmb250LXNpemU9IjE0Ij5FQlRBVDwvdGV4dD4KPC9zdmc+",
        height = "30px", style = "margin-right:10px;"),
    span("ENABL Biomarker / Target Analysis Tool (EB-TAT)", style = "font-weight: bold;")
  ),
  id = "main_nav",
  theme = shinythemes::shinytheme("flatly"),
  collapsible = TRUE,
  header = tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0/css/all.min.css"),
    tags$style(HTML("
      /* Global background and text colors */
      body { 
        background-color: white !important; 
        color: #333333;
      }
      
      .navbar { 
        font-weight: bold; 
        background-color: #2196F3 !important;
      }
      
      .shiny-notification { 
        position: fixed; 
        top: 20px; 
        right: 20px; 
        width: 350px; 
      }
      
      .loading-spinner { 
        position: absolute; 
        top: 50%; 
        left: 50%; 
        transform: translate(-50%, -50%); 
      }
      
      .module-error { 
        background-color: #ffebee; 
        padding: 20px; 
        border-radius: 5px; 
        margin: 20px; 
      }
      
      .success-message { 
        color: green; 
        font-weight: bold; 
      }
      
      .error-message { 
        color: red; 
        font-weight: bold; 
      }
      
      .warning-message { 
        color: orange; 
        font-weight: bold; 
      }
      
      .feature-list { 
        list-style: none; 
        padding-left: 0; 
      }
      
      .feature-list li { 
        margin-bottom: 12px; 
        padding-left: 30px; 
        position: relative; 
      }
      
      .feature-list li i { 
        position: absolute; 
        left: 0; 
        top: 2px; 
        font-size: 16px; 
      }
      
      /* Card styling */
      .card {
        background-color: white;
        border: 1px solid #dee2e6;
        border-radius: 0.375rem;
        box-shadow: 0 0.125rem 0.25rem rgba(0, 0, 0, 0.075);
      }
      
      .card-header {
        background-color: #f8f9fa;
        border-bottom: 1px solid #dee2e6;
        padding: 0.75rem 1.25rem;
      }
      
      .card-body {
        padding: 1.25rem;
        background-color: white;
      }
      
      /* Ensure all content areas have white background */
      .main-content, 
      .tab-content, 
      .shiny-tab-content,
      .dataTables_wrapper {
        background-color: white !important;
      }
      
      /* Data input module specific styling */
      .data-input-module {
        background-color: white;
        min-height: 100vh;
      }
      
      /* Button styling */
      .btn-block {
        width: 100%;
      }
      
      /* Tab styling */
      .nav-tabs .nav-link.active {
        background-color: white;
        border-color: #dee2e6 #dee2e6 white;
      }
    "))
  ),
  footer = div(
    style = "text-align: center; padding: 15px; font-size: 12px; color: #666; background-color: #f8f9fa; border-top: 1px solid #dee2e6;",
    " 2024 ENABL Consortium | ",
    a("Contact Support", href = "mailto:support@enabl.org", style = "color: #007bff; text-decoration: none;"),
    " | ",
    a("Documentation", href = "#", onclick = "window.open('README.md', '_blank'); return false;", style = "color: #007bff; text-decoration: none;")
  ),

  # Data input tab
  tabPanel(
    title = span(icon("database"), "Data"),
    value = "data_tab",
    module_ui_functions$data_input("data_input")
  ),

  # Volcano plot tab
  tabPanel(
    title = span(icon("chart-line"), "Volcano Plot"),
    value = "volcano_tab",
    module_ui_functions$volcano_plot("volcano_plot")
  ),

  # Analytics tab
  tabPanel(
    title = span(icon("chart-bar"), "Analytics"),
    value = "analytics_tab",
    module_ui_functions$analytics("analytics")
  ),

  # Enrichment tab
  tabPanel(
    title = span(icon("project-diagram"), "Enrichment"),
    value = "enrichment_tab",
    module_ui_functions$enrichment("enrichment")
  ),

  # Expression plot tab
  tabPanel(
    title = span(icon("braille"), "Expression Plot"),
    value = "expression_tab",
    module_ui_functions$expression_plot("expression_plot")
  ),

  # Literature Trends tab - NEW
  tabPanel(
    title = span(icon("book"), "Literature Trends"),
    value = "literature_trends_tab",
    module_ui_functions$literature_trends("literature_trends")
  ),

  # About tab
  tabPanel(
    title = span(icon("info-circle"), "About"),
    value = "about_tab",
    module_ui_functions$about("about")
  )
)

# Define server logic
server <- function(input, output, session) {
  # Initialize reactive values with defaults
  rv <- reactiveValues(
    data = NULL,
    filtered_data = NULL,
    user_selection = NULL,
    data_loaded = FALSE,
    last_error = NULL
    # Removed module_status from here to fix the reactive context error
  )
  
  # Initialize module-specific reactive values
  rv_analytics <- reactiveValues(
    fc_cut = c(-0.5, 0.5),
    p_cut = 1.301,
    top_hit = 20,
    user_gene_list1 = NULL,
    plot_type = "reg",
    venntype = "combined",
    boxp_ind = "D",
    indication = "D",
    indication_p = "F",
    pub_tr = "F",
    ord_of = "F",
    ppiscore = "900",
    hse = FALSE,
    rec_lig = FALSE,
    sec = FALSE,
    sec2 = FALSE,
    nessen = FALSE,
    bld = FALSE,
    tf = FALSE,
    immun = FALSE,
    show_selected_table = FALSE
  )

  rv_enrichment <- reactiveValues(
    fc_cut_e = c(-1.0, 1.0),
    p_cut_e = 1.31,
    e_plot_type = "kpe",
    olap = FALSE
  )

  rv_expression <- reactiveValues(
    user_gene_list2 = NULL,
    ptype = "ap"
  )

  # Show startup message
  showNotification(
    HTML("<div style='font-weight: bold;'>EB-TAT is loading. Please wait...</div>"),
    type = "message",
    duration = 5
  )

  # Safe module calling function with flexible naming - FIXED VERSION
  safe_call_module <- function(module_name, display_name, id, ...) {
    tryCatch({
      # Try multiple possible server function names
      possible_server_names <- c(
        # CamelCase version
        paste0(gsub("_(.)", "\\U\\1", module_name, perl = TRUE), "Server"),
        # Original snake_case with Server suffix
        paste0(module_name, "Server"),
        # Direct mapping for consistency
        "dataInputServer", "volcanoPlotServer", "analyticsServer", "enrichmentServer", 
        "expressionPlotServer", "literatureTrendsServer", "aboutServer"
      )
      
      server_function <- NULL
      for (server_name in possible_server_names) {
        if (exists(server_name)) {
          server_function <- get(server_name)
          message("Found server function: ", server_name)
          break
        }
      }
      
      if (!is.null(server_function)) {
        result <- do.call(server_function, list(id, ...))
        message("Successfully loaded module: ", display_name)
        return(result)
      } else {
        error_msg <- paste("Module server function not found. Tried:", paste(possible_server_names, collapse = ", "))
        warning(error_msg)
        return(NULL)
      }
    }, error = function(e) {
      error_msg <- paste("Error in", display_name, "module:", e$message)
      rv$last_error <- error_msg
      warning(error_msg)
      return(NULL)
    })
  }

  # Call modules with error handling - FIXED: Removed module_status tracking
  data_input <- safe_call_module("data_input", "Data Input", "data_input", rv)
  volcano_plot <- safe_call_module("volcano_plot", "Volcano Plot", "volcano_plot", rv)
  analytics <- safe_call_module("analytics", "Analytics", "analytics", rv, rv_analytics)
  enrichment <- safe_call_module("enrichment", "Enrichment", "enrichment", rv, rv_enrichment)
  expression_plot <- safe_call_module("expression_plot", "Expression Plot", "expression_plot", rv, rv_expression)
  literature_trends <- safe_call_module("literature_trends", "Literature Trends", "literature_trends", rv)
  about <- safe_call_module("about", "About", "about")

  # Observe changes from data input module
  observe({
    if (!is.null(data_input)) {
      data_input_result <- data_input()
      if (!is.null(data_input_result)) {
        rv$data <- data_input_result$data
        rv$filtered_data <- data_input_result$filtered_data
        rv$data_loaded <- data_input_result$data_loaded

        # Show success message
        if (!is.null(rv$data) && nrow(rv$data) > 0 && rv$data_loaded) {
          showNotification(
            HTML(paste0("<div class='success-message'>Data loaded successfully: ",
                       nrow(rv$data), " rows, ", ncol(rv$data), " columns</div>")),
            type = "message",
            duration = 3
          )
        }
      }
    }
  })

  # Add debug logging - FIXED VERSION (no reactive value access issues)
  observe({
    # Create a safe version that won't cause reactive errors
    tryCatch({
      # Use reactiveValuesToList to safely access all values
      current_state <- reactiveValuesToList(rv)
      
      cat("DEBUG: Current reactive values:\n")
      cat("  data_loaded:", current_state$data_loaded, "\n")
      cat("  data rows:", ifelse(!is.null(current_state$data), nrow(current_state$data), 0), "\n")
      cat("  last_error:", ifelse(!is.null(current_state$last_error), current_state$last_error, "None"), "\n")
      
    }, error = function(e) {
      cat("DEBUG: Error accessing reactive values:", e$message, "\n")
    })
    
    invalidateLater(10000) # Update every 10 seconds
  })

  # Build session info function - FIXED VERSION (no module_status)
  build_session_info <- function() {
    tryCatch({
      # Safely access reactive values
      current_state <- reactiveValuesToList(rv)
      
      txt <- c(
        paste("Current Time:", format(Sys.time())),
        paste("Session Start:", ifelse(!is.null(session$startTime), format(session$startTime), "Not available")),
        paste("Working Directory:", getwd()),
        paste("R Version:", R.version.string),
        paste("Data loaded:", current_state$data_loaded),
        paste("Data rows:", ifelse(!is.null(current_state$data), nrow(current_state$data), 0)),
        paste("Data columns:", ifelse(!is.null(current_state$data), ncol(current_state$data), 0)),
        paste("Last error:", ifelse(!is.null(current_state$last_error), current_state$last_error, "None"))
      )

      paste(txt, collapse = "\n")
    }, error = function(e) {
      return(paste("Error generating session info:", e$message))
    })
  }

  # Observe URL parameters for pre-loading settings
  observe({
    query <- parseQueryString(session$clientData$url_search)

    # Load data from URL if specified
    if (!is.null(query[['url']])) {
      tryCatch({
        # This would need to be implemented in the data input module
        showNotification("URL parameter detected. Please use the Data tab to load from URL.", type = "default")
      }, error = function(e) {
        message("Error handling URL parameter: ", e$message)
      })
    }
  })

  # Global error handling
  observe({
    if (!is.null(rv$last_error)) {
      showNotification(
        HTML(paste0("<div class='error-message'>Error: ", rv$last_error, "</div>")),
        type = "error",
        duration = 10
      )
    }

    if (!is.null(rv$data) && nrow(rv$data) == 0) {
      showNotification(
        HTML("<div class='warning-message'>Warning: No data available or data format is incorrect.</div>"),
        type = "warning",
        duration = 5
      )
    }
  })

  # Session info for debugging
  output$session_info <- renderPrint({
    cat(build_session_info())
  })

  # Help system - add tooltips
  observe({
    tryCatch({
      # Add tooltips as needed
    }, error = function(e) {
      # Silently fail if tooltip can't be added
    })
  })

  # Handle browser refresh/reload
  observeEvent(input$reload, {
    session$reload()
  })

  # Download handlers
  output$download_session <- downloadHandler(
    filename = function() {
      paste0("session_info_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
    },
    content = function(file) {
      txt <- build_session_info()
      writeLines(txt, file)
    }
  )

  output$download_session_info <- downloadHandler(
    filename = function() {
      paste("ebtat_session_info_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt", sep = "")
    },
    content = function(file) {
      txt <- build_session_info()
      writeLines(txt, file)
    }
  )

  # Handle session end
  session$onSessionEnded(function() {
    message("EB-TAT session ended: ", Sys.time())
    if (!is.null(session$startTime)) {
      message("Session duration: ", round(difftime(Sys.time(), session$startTime, units = "mins"), 1), " minutes")
    }
  })

  # Debug panel (hidden by default)
  output$debug_panel <- renderUI({
    if (isTRUE(getOption("shiny.debug", FALSE))) {
      tagList(
        h4("Debug Information"),
        verbatimTextOutput("session_info"),
        downloadButton("download_session_info", "Download Session Info"),
        actionButton("reload", "Reload App", icon = icon("refresh"))
      )
    }
  })

  # Status monitoring
  observe({
    invalidateLater(30000) # Update every 30 seconds
    message("EB-TAT Status Check - Time: ", Sys.time())
    message("Data loaded: ", rv$data_loaded)
    if (rv$data_loaded) {
      message("Data dimensions: ", nrow(rv$data), " x ", ncol(rv$data))
    }
  })
}

# Enable bookmarking
enableBookmarking(store = "url")

# Safe application runner with comprehensive error handling
safe_runApp <- function() {
  tryCatch({
    # Check if required directories exist
    required_dirs <- c("modules", "utils", "www")
    for (dir in required_dirs) {
      if (!dir.exists(dir)) {
        dir.create(dir, recursive = TRUE)
        warning("Created missing directory: ", dir)
      }
    }

    # Check if required files exist
    required_files <- c("global.R")
    missing_files <- c()

    for (file in required_files) {
      if (!file.exists(file)) {
        missing_files <- c(missing_files, file)
      }
    }

    if (length(missing_files) > 0) {
      warning("Missing required files: ", paste(missing_files, collapse = ", "))
    }

    # Check utility files
    utility_files <- c("utils/helpers.R", "utils/data_loading.R", "utils/plotting.R")
    for (file in utility_files) {
      if (!file.exists(file)) {
        warning("Utility file not found: ", file)
      }
    }

    # Check module files
    module_files <- c(
      "modules/data_input.R",
      "modules/volcano_plot.R",
      "modules/analytics.R",
      "modules/enrichment.R",
      "modules/expression_plot.R",
      "modules/literature_trends.R",
      "modules/about.R"
    )

    for (file in module_files) {
      if (!file.exists(file)) {
        warning("Module file not found: ", file)
      }
    }

    # Run the app
    message("Starting EB-TAT application...")
    message("Server started at: ", Sys.time())
    message("Working directory: ", getwd())
    message("Available files: ", paste(list.files(), collapse = ", "))

    shinyApp(ui = ui, server = server, options = list(
      host = "0.0.0.0",
      port = 3838,
      launch.browser = TRUE,
      display.mode = "normal"
    ))

  }, error = function(e) {
    message("Error starting EB-TAT application: ", e$message)
    message("Working directory: ", getwd())
    message("Files in directory: ", paste(list.files(), collapse = ", "))

    # Create a minimal error UI
    error_ui <- fluidPage(
      theme = shinythemes::shinytheme("flatly"),
      titlePanel(
        div(
          style = "display: flex; align-items: center;",
          img(src = "data:image/svg+xml;base64,PHN2ZyB3aWR0aD0iMzAiIGhlaWdodD0iMzAiIHZpZXdCb3g9IjAgMCAzMCAzMCIgZmlsbD0ibm9uZSIgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnIj4KPHJlY3Qgd2lkdGg9IjMwIiBoZWlnaHQ9IjMwIiBmaWxsPSIjZGMzNTQ1Ii8+Cjx0ZXh0IHg9IjgiIHk9IjIwIiBmaWxsPSJ3aGl0ZSIgZm9udC1mYW1pbHk9IkFyaWFsIiBmb250LXNpemU9IjE0Ij5FQlRBVDwvdGV4dD4KPC9zdmc+",
              height = "30px", style = "margin-right:10px;"),
          span("EB-TAT - Application Error", style = "font-weight: bold; color: #dc3545;")
        )
      ),
      div(
        style = "text-align: center; padding: 50px;",
        h3("Application Failed to Start", style = "color: #dc3545;"),
        div(
          style = "background-color: #f8d7da; color: #721c24; padding: 20px; border-radius: 5px; margin: 20px;",
          h4("Error Details:"),
          p(style = "font-family: monospace;", e$message),
          p("Working directory: ", getwd())
        ),
        p("Please check that all required files are present and run install.sh to set up the application."),
        actionButton("reload", "Reload Application", icon = icon("refresh"), class = "btn-primary")
      )
    )

    error_server <- function(input, output, session) {
      observeEvent(input$reload, {
        session$reload()
      })
    }

    shinyApp(ui = error_ui, server = error_server)
  })
}

# Run the application
if (interactive()) {
  safe_runApp()
} else {
  # For deployment
  shinyApp(ui = ui, server = server)
}