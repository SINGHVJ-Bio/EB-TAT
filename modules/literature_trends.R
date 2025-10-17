# ENABL Biomarker/Target Analysis Tool - Literature Trends Module
# Dynamic PubMed literature trend analysis

# =============================================================================
# UI Function
# =============================================================================

literatureTrendsUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    div(
      class = "literature-trends-module",
      style = "padding: 20px;",
      
      # Header
      div(
        class = "module-header",
        h3("Literature Trends", style = "color: #2196F3;"),
        p("Dynamic PubMed literature trend analysis for genes and diseases")
      ),
      
      fluidRow(
        # Controls column
        column(
          4,
          # Gene selection
          div(
            class = "card",
            style = "margin-bottom: 20px;",
            div(
              class = "card-header",
              h4("Gene Selection", style = "margin: 0;")
            ),
            div(
              class = "card-body",
              textAreaInput(
                ns("genes_input"),
                "Enter gene symbols (one per line):",
                rows = 6,
                placeholder = "TP53\nEGFR\nBRCA1\n..."
              ),
              actionButton(
                ns("load_top_genes"),
                "Load Top Significant Genes",
                class = "btn-outline-primary"
              ),
              helpText("Select up to 10 genes for analysis")
            )
          ),
          
          # Disease terms configuration
          div(
            class = "card",
            style = "margin-bottom: 20px;",
            div(
              class = "card-header",
              h4("Disease Terms", style = "margin: 0;")
            ),
            div(
              class = "card-body",
              textAreaInput(
                ns("disease_terms"),
                "Disease search terms (one per line, combined with AND):",
                rows = 4,
                placeholder = "fibrosis\nliver"
              ),
              helpText("Multiple terms will be combined with AND operator"),
              numericInput(
                ns("years_back"),
                "Years to analyze:",
                value = 10,
                min = 5,
                max = 20,
                step = 1
              )
            )
          ),
          
          # Analysis controls
          div(
            class = "card",
            style = "margin-bottom: 20px;",
            div(
              class = "card-header",
              h4("Analysis Controls", style = "margin: 0;")
            ),
            div(
              class = "card-body",
              actionButton(
                ns("run_analysis"),
                "Run Literature Analysis",
                class = "btn-success",
                icon = icon("search")
              ),
              actionButton(
                ns("clear_results"),
                "Clear Results",
                class = "btn-warning",
                icon = icon("trash")
              ),
              hr(),
              p(strong("Note:"), "PubMed API has rate limits. Please be patient."),
              verbatimTextOutput(ns("api_status"))
            )
          ),
          
          # Export options
          div(
            class = "card",
            div(
              class = "card-header",
              h4("Export Results", style = "margin: 0;")
            ),
            div(
              class = "card-body",
              downloadButton(
                ns("download_trends"),
                "Download Trends Data",
                class = "btn-primary"
              ),
              downloadButton(
                ns("download_plot"),
                "Download Plot",
                class = "btn-info"
              )
            )
          )
        ),
        
        # Results column
        column(
          8,
          # Status and summary
          div(
            class = "card",
            style = "margin-bottom: 20px;",
            div(
              class = "card-header",
              h4("Analysis Status", style = "margin: 0;")
            ),
            div(
              class = "card-body",
              conditionalPanel(
                condition = paste0("output['", ns("analysis_complete"), "']"),
                div(
                  class = "alert alert-success",
                  strong("Analysis complete!"),
                  textOutput(ns("analysis_summary"))
                )
              ),
              verbatimTextOutput(ns("analysis_info"))
            )
          ),
          
          # Main trends plot
          div(
            class = "card",
            style = "margin-bottom: 20px;",
            div(
              class = "card-header",
              h4("Literature Trends", style = "margin: 0;")
            ),
            div(
              class = "card-body",
              withSpinner(
                plotOutput(ns("trends_plot"), height = "500px"),
                type = 4,
                color = "#2196F3"
              )
            )
          ),
          
          # Results tables
          div(
            class = "card",
            div(
              class = "card-header",
              h4("Analysis Results", style = "margin: 0;")
            ),
            div(
              class = "card-body",
              tabsetPanel(
                type = "tabs",
                tabPanel(
                  "Trends Data",
                  DT::dataTableOutput(ns("trends_table"))
                ),
                tabPanel(
                  "Gene Summary",
                  DT::dataTableOutput(ns("summary_table"))
                )
              )
            )
          )
        )
      )
    )
  )
}

# =============================================================================
# Server Function
# =============================================================================

literatureTrendsServer <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Reactive values for this module
    module_rv <- reactiveValues(
      trends_data = NULL,
      analysis_complete = FALSE,
      current_plot = NULL,
      last_analysis = NULL
    )
    
    # Load top significant genes
    observeEvent(input$load_top_genes, {
      if (!rv$data_loaded || is.null(rv$filtered_data)) {
        showNotification("No data available. Please load data first.", type = "warning")
        return()
      }
      
      # Get top 10 significant genes
      top_genes <- get_top_genes(rv$filtered_data, n = 10, by = "padj")
      gene_symbols <- top_genes$symbol
      
      # Update gene input
      updateTextAreaInput(
        session, 
        "genes_input", 
        value = paste(gene_symbols, collapse = "\n")
      )
      
      showNotification(paste("Loaded", length(gene_symbols), "top significant genes"), 
                      type = "message")
    })
    
    # Set default disease terms from config
    observe({
      if (exists("literature_config") && !is.null(literature_config$disease_terms)) {
        default_terms <- paste(literature_config$disease_terms, collapse = "\n")
        updateTextAreaInput(session, "disease_terms", value = default_terms)
      }
      
      if (exists("literature_config") && !is.null(literature_config$years_back)) {
        updateNumericInput(session, "years_back", value = literature_config$years_back)
      }
    })
    
    # Run literature analysis
    run_literature_analysis <- function() {
      # Get genes from input
      gene_text <- input$genes_input
      if (is.null(gene_text) || gene_text == "") {
        showNotification("Please enter gene symbols", type = "warning")
        return(NULL)
      }
      
      genes <- parse_gene_input(gene_text)
      if (length(genes) == 0) {
        showNotification("No valid gene symbols found", type = "error")
        return(NULL)
      }
      
      # Limit to 10 genes for performance
      if (length(genes) > 10) {
        genes <- genes[1:10]
        showNotification("Limited to first 10 genes", type = "warning")
      }
      
      # Get disease terms
      disease_text <- input$disease_terms
      if (is.null(disease_text) || disease_text == "") {
        showNotification("Please enter disease terms", type = "warning")
        return(NULL)
      }
      
      disease_terms <- unlist(strsplit(disease_text, "[\r\n]+"))
      disease_terms <- trimws(disease_terms)
      disease_terms <- disease_terms[disease_terms != ""]
      
      if (length(disease_terms) == 0) {
        showNotification("No valid disease terms found", type = "error")
        return(NULL)
      }
      
      # Get email from config
      email <- if (exists("literature_config")) literature_config$email else NULL
      
      if (is.null(email)) {
        showNotification("Email not configured for PubMed API. Please check config.", 
                        type = "error")
        return(NULL)
      }
      
      # Show progress
      showNotification("Searching PubMed... This may take a few minutes.", 
                      type = "message", duration = 10)
      
      # Run analysis
      trends_data <- get_literature_trends(
        gene_symbols = genes,
        disease_terms = disease_terms,
        years_back = input$years_back,
        email = email
      )
      
      if (is.null(trends_data)) {
        showNotification("No literature trends data found or error occurred", 
                        type = "warning")
        return(NULL)
      }
      
      module_rv$trends_data <- trends_data
      module_rv$analysis_complete <- TRUE
      module_rv$last_analysis <- Sys.time()
      
      showNotification("Literature analysis completed successfully!", 
                      type = "message")
      
      return(trends_data)
    }
    
    # Create trends plot
    create_trends_plot <- function() {
      trends_data <- module_rv$trends_data
      if (is.null(trends_data)) {
        return(create_empty_plot("Run analysis to see literature trends"))
      }
      
      disease_terms <- unlist(strsplit(input$disease_terms, "[\r\n]+"))
      disease_terms <- trimws(disease_terms)
      disease_terms <- disease_terms[disease_terms != ""]
      
      title <- paste("Literature Trends:", paste(disease_terms, collapse = " + "))
      
      p <- create_literature_trends_plot(trends_data, title)
      module_rv$current_plot <- p
      
      return(p)
    }
    
    # Create summary table
    create_summary_table <- function() {
      trends_data <- module_rv$trends_data
      if (is.null(trends_data)) return(NULL)
      
      summary_data <- trends_data %>%
        group_by(gene) %>%
        summarise(
          total_publications = sum(count, na.rm = TRUE),
          max_year = year[which.max(count)],
          max_count = max(count, na.rm = TRUE),
          .groups = 'drop'
        ) %>%
        arrange(desc(total_publications))
      
      return(summary_data)
    }
    
    # Event observers
    observeEvent(input$run_analysis, {
      run_literature_analysis()
    })
    
    observeEvent(input$clear_results, {
      module_rv$trends_data <- NULL
      module_rv$analysis_complete <- FALSE
      module_rv$current_plot <- NULL
      updateTextAreaInput(session, "genes_input", value = "")
      showNotification("Results cleared", type = "warning")
    })
    
    # Outputs
    output$analysis_complete <- reactive({
      module_rv$analysis_complete
    })
    outputOptions(output, "analysis_complete", suspendWhenHidden = FALSE)
    
    output$api_status <- renderText({
      if (exists("literature_config") && !is.null(literature_config$email)) {
        paste("PubMed API: Configured\nEmail:", literature_config$email)
      } else {
        "PubMed API: Not configured\nPlease set email in config/data_paths.yml"
      }
    })
    
    output$analysis_summary <- renderText({
      if (!module_rv$analysis_complete) return("")
      
      trends_data <- module_rv$trends_data
      genes <- unique(trends_data$gene)
      years <- unique(trends_data$year)
      
      paste(
        "Analyzed", length(genes), "genes over", length(years), "years",
        "| Total publications:", sum(trends_data$count, na.rm = TRUE)
      )
    })
    
    output$analysis_info <- renderPrint({
      if (!module_rv$analysis_complete) {
        cat("Literature analysis not run yet.\n")
        cat("Enter genes and disease terms, then click 'Run Literature Analysis'.\n")
        return()
      }
      
      trends_data <- module_rv$trends_data
      genes <- unique(trends_data$gene)
      disease_terms <- unlist(strsplit(input$disease_terms, "[\r\n]+"))
      disease_terms <- trimws(disease_terms)
      
      cat("LITERATURE ANALYSIS SUMMARY\n")
      cat("===========================\n")
      cat("Analysis time:", format(module_rv$last_analysis), "\n")
      cat("Genes analyzed:", length(genes), "\n")
      cat("Disease terms:", paste(disease_terms, collapse = ", "), "\n")
      cat("Years analyzed:", input$years_back, "\n\n")
      
      # Top genes by publication count
      summary_data <- create_summary_table()
      if (!is.null(summary_data)) {
        cat("TOP GENES BY PUBLICATION COUNT:\n")
        for (i in 1:min(5, nrow(summary_data))) {
          cat(i, ". ", summary_data$gene[i], 
              " (", summary_data$total_publications[i], " publications)\n", sep = "")
        }
      }
    })
    
    output$trends_plot <- renderPlot({
      create_trends_plot()
    })
    
    output$trends_table <- DT::renderDataTable({
      trends_data <- module_rv$trends_data
      if (is.null(trends_data)) return(NULL)
      
      DT::datatable(
        trends_data,
        options = list(
          scrollX = TRUE,
          pageLength = 10,
          dom = 'tip'
        ),
        rownames = FALSE,
        caption = "Literature Trends Data"
      )
    })
    
    output$summary_table <- DT::renderDataTable({
      summary_data <- create_summary_table()
      if (is.null(summary_data)) return(NULL)
      
      DT::datatable(
        summary_data,
        options = list(
          scrollX = TRUE,
          pageLength = 10,
          dom = 'tip'
        ),
        rownames = FALSE,
        caption = "Gene Publication Summary"
      )
    })
    
    # Download handlers
    output$download_trends <- downloadHandler(
      filename = function() {
        paste0("literature_trends_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
      },
      content = function(file) {
        trends_data <- module_rv$trends_data
        if (is.null(trends_data)) {
          showNotification("No trends data to download", type = "error")
          return()
        }
        
        write.csv(trends_data, file, row.names = FALSE)
        showNotification("Trends data downloaded successfully!", type = "message")
      }
    )
    
    output$download_plot <- downloadHandler(
      filename = function() {
        paste0("literature_trends_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
      },
      content = function(file) {
        plot_obj <- module_rv$current_plot
        if (is.null(plot_obj)) {
          showNotification("No plot to download", type = "error")
          return()
        }
        
        success <- save_plot(plot_obj, file, width = 12, height = 8, dpi = 300)
        
        if (success) {
          showNotification("Plot downloaded successfully!", type = "message")
        } else {
          showNotification("Error downloading plot", type = "error")
        }
      }
    )
    
    # Return reactive values for other modules
    return(reactive({
      list(
        trends_data = module_rv$trends_data,
        analysis_complete = module_rv$analysis_complete
      )
    }))
  })
}