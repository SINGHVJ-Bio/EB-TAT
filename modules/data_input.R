# ENABL Biomarker/Target Analysis Tool - Data Input Module
# Handles data loading from various sources

# =============================================================================
# UI Function
# =============================================================================

dataInputUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    div(
      class = "data-input-module",
      style = "padding: 20px; background-color: white;",
      
      # Header
      div(
        class = "module-header",
        h3("Data Input", style = "color: #2196F3;"),
        p("Load differential expression data from various sources")
      ),
      
      fluidRow(
        # Left Column - Controls and Filters (30%)
        column(
          4,
          # Data source selection
          div(
            class = "card",
            style = "margin-bottom: 20px;",
            div(
              class = "card-header bg-primary text-white",
              h4("Data Source", style = "margin: 0;")
            ),
            div(
              class = "card-body",
              radioButtons(
                ns("data_source"),
                label = "Choose data source:",
                choices = list(
                  "Pre-loaded dataset" = "preloaded",
                  "Upload CSV/TSV file" = "upload",
                  "Enter gene list" = "manual",
                  "From URL" = "url"
                ),
                selected = "preloaded",
                inline = FALSE
              )
            )
          ),
          
          # Pre-loaded data options
          conditionalPanel(
            condition = paste0("input['", ns("data_source"), "'] == 'preloaded'"),
            div(
              class = "card",
              style = "margin-bottom: 20px;",
              div(
                class = "card-header bg-info text-white",
                h4("Pre-loaded Data", style = "margin: 0;")
              ),
              div(
                class = "card-body",
                p("Using the pre-loaded differential expression dataset."),
                verbatimTextOutput(ns("preloaded_info"))
              )
            )
          ),
          
          # File upload options
          conditionalPanel(
            condition = paste0("input['", ns("data_source"), "'] == 'upload'"),
            div(
              class = "card",
              style = "margin-bottom: 20px;",
              div(
                class = "card-header bg-info text-white",
                h4("Upload Data File", style = "margin: 0;")
              ),
              div(
                class = "card-body",
                fileInput(
                  ns("file_upload"),
                  "Choose CSV/TSV file:",
                  accept = c(".csv", ".tsv", ".txt"),
                  placeholder = "No file selected"
                ),
                helpText("File should contain columns: gene symbols (or ENSEMBL IDs), fold changes, and p-values"),
                helpText("ENSEMBL IDs will be automatically converted to gene symbols"),
                numericInput(
                  ns("skip_rows"),
                  "Skip rows:",
                  value = 0,
                  min = 0,
                  max = 100
                ),
                checkboxInput(
                  ns("header"),
                  "File has header row",
                  value = TRUE
                )
              )
            )
          ),
          
          # Manual input options
          conditionalPanel(
            condition = paste0("input['", ns("data_source"), "'] == 'manual'"),
            div(
              class = "card",
              style = "margin-bottom: 20px;",
              div(
                class = "card-header bg-info text-white",
                h4("Manual Gene Input", style = "margin: 0;")
              ),
              div(
                class = "card-body",
                textAreaInput(
                  ns("manual_genes"),
                  "Enter gene symbols (one per line or comma-separated):",
                  rows = 6,
                  placeholder = "TP53\nEGFR\nBRCA1\n..."
                ),
                helpText("Enter gene symbols. Fold changes and p-values will be simulated for demonstration.")
              )
            )
          ),
          
          # URL input options
          conditionalPanel(
            condition = paste0("input['", ns("data_source"), "'] == 'url'"),
            div(
              class = "card",
              style = "margin-bottom: 20px;",
              div(
                class = "card-header bg-info text-white",
                h4("Load from URL", style = "margin: 0;")
              ),
              div(
                class = "card-body",
                textInput(
                  ns("data_url"),
                  "Data URL:",
                  placeholder = "https://example.com/data.csv"
                ),
                actionButton(
                  ns("load_url"),
                  "Load from URL",
                  class = "btn-primary"
                )
              )
            )
          ),
          
          # Data filtering options - SEPARATE ROWS
          div(
            class = "card",
            style = "margin-bottom: 20px;",
            div(
              class = "card-header bg-warning text-dark",
              h4("Data Filtering", style = "margin: 0;")
            ),
            div(
              class = "card-body",
              # Fold Change in separate row
              div(
                style = "margin-bottom: 20px;",
                h5("Fold Change Range:"),
                sliderInput(
                  ns("fc_filter"),
                  label = NULL,
                  min = -10,
                  max = 10,
                  value = c(-2, 2),
                  step = 0.1
                )
              ),
              
              # P-value in separate row
              div(
                style = "margin-bottom: 20px;",
                h5("P-value Threshold:"),
                numericInput(
                  ns("pval_filter"),
                  label = NULL,
                  value = 0.05,
                  min = 0,
                  max = 1,
                  step = 0.01,
                  width = "100%"
                )
              ),
              
              # Apply filters button
              actionButton(
                ns("apply_filters"),
                "Apply Filters",
                class = "btn-info btn-block"
              )
            )
          ),
          
          # Data controls
          div(
            class = "card",
            div(
              class = "card-header bg-success text-white",
              h4("Data Controls", style = "margin: 0;")
            ),
            div(
              class = "card-body",
              fluidRow(
                column(
                  6,
                  actionButton(
                    ns("load_data"),
                    "Load Data",
                    class = "btn-success btn-block",
                    icon = icon("database")
                  )
                ),
                column(
                  6,
                  actionButton(
                    ns("reset_data"),
                    "Reset",
                    class = "btn-warning btn-block",
                    icon = icon("undo")
                  )
                )
              ),
              br(),
              conditionalPanel(
                condition = paste0("output['", ns("data_loaded"), "']"),
                div(
                  style = "background-color: #d4edda; padding: 10px; border-radius: 5px; text-align: center;",
                  icon("check-circle", style = "color: #28a745;"),
                  strong(" Data loaded successfully!"),
                  textOutput(ns("data_summary"))
                )
              )
            )
          )
        ),
        
        # Right Column - Preview and Information (70%)
        column(
          8,
          # Data preview and information
          div(
            class = "card",
            style = "margin-bottom: 20px; min-height: 600px;",
            div(
              class = "card-header bg-secondary text-white",
              h4("Data Preview & Information", style = "margin: 0;")
            ),
            div(
              class = "card-body",
              style = "min-height: 550px;",
              tabsetPanel(
                type = "tabs",
                tabPanel(
                  "Preview",
                  br(),
                  div(
                    style = "max-height: 500px; overflow-y: auto;",
                    DT::dataTableOutput(ns("data_preview"))
                  )
                ),
                tabPanel(
                  "Summary",
                  br(),
                  div(
                    style = "max-height: 500px; overflow-y: auto;",
                    verbatimTextOutput(ns("data_stats"))
                  )
                ),
                tabPanel(
                  "Column Info",
                  br(),
                  div(
                    style = "max-height: 500px; overflow-y: auto;",
                    verbatimTextOutput(ns("column_info"))
                  )
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

dataInputServer <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Reactive values for this module
    module_rv <- reactiveValues(
      raw_data = NULL,
      filtered_data = NULL,
      data_loaded = FALSE,
      last_error = NULL
    )
    
    # Load pre-loaded data - MODIFIED FOR GRACEFUL HANDLING
    load_preloaded_data <- function() {
      tryCatch({
        if (exists("diffgenes") && !is.null(diffgenes)) {
          message("Loading pre-loaded differential expression data")
          return(diffgenes)
        } else {
          # Create example data if pre-loaded doesn't exist
          message("Creating example data for demonstration")
          set.seed(123)
          example_data <- data.frame(
            symbol = paste0("GENE", 1:1000),
            log2FoldChange = rnorm(1000, 0, 2),
            padj = runif(1000, 0, 1),
            pvalue = runif(1000, 0, 1),
            baseMean = runif(1000, 10, 1000),
            stringsAsFactors = FALSE
          )
          example_data$minus_log10_padj <- -log10(example_data$padj)
          return(example_data)
        }
      }, error = function(e) {
        warning("Error loading pre-loaded data: ", e$message)
        # Create minimal data as fallback
        set.seed(123)
        example_data <- data.frame(
          symbol = paste0("GENE", 1:500),
          log2FoldChange = rnorm(500, 0, 1.5),
          padj = runif(500, 0, 0.5),
          pvalue = runif(500, 0, 0.3),
          baseMean = runif(500, 50, 500),
          stringsAsFactors = FALSE
        )
        example_data$minus_log10_padj <- -log10(example_data$padj)
        return(example_data)
      })
    }
    
    # Load uploaded file - UPDATED FOR ENSEMBL ID CONVERSION
    load_uploaded_data <- function(file_info) {
      if (is.null(file_info)) return(NULL)
      
      tryCatch({
        ext <- tools::file_ext(file_info$name)
        
        if (ext %in% c("csv", "txt")) {
          data <- read.csv(
            file_info$datapath,
            header = input$header,
            skip = input$skip_rows,
            stringsAsFactors = FALSE
          )
        } else if (ext == "tsv") {
          data <- read.delim(
            file_info$datapath,
            header = input$header,
            skip = input$skip_rows,
            stringsAsFactors = FALSE
          )
        } else {
          stop("Unsupported file format: ", ext)
        }
        
        # Add symbol column if missing using ENSEMBL ID conversion
        if (!"symbol" %in% colnames(data)) {
          if ("gene" %in% colnames(data) && all(grepl("^ENSG", data$gene))) {
            message("Converting ENSEMBL IDs in uploaded file to gene symbols")
            data <- add_symbol_column(data, "gene")
          } else if (!is.null(rownames(data)) && all(grepl("^ENSG", rownames(data)))) {
            message("Converting row names from ENSEMBL IDs to gene symbols")
            symbol_names <- convert_ensembl_to_symbol(rownames(data))
            data$symbol <- symbol_names
          } else {
            # Use first column as symbol
            data$symbol <- data[[1]]
            message("Using first column as gene symbols")
          }
        }
        
        # Standardize column names
        data <- standardize_de_columns(data)
        
        # Calculate -log10 p-value if needed
        if (!"minus_log10_padj" %in% colnames(data) && "padj" %in% colnames(data)) {
          data$minus_log10_padj <- -log10(data$padj)
        }
        
        return(data)
      }, error = function(e) {
        warning("Error loading uploaded file: ", e$message)
        return(NULL)
      })
    }
    
    # Create data from manual gene input
    create_manual_data <- function(gene_text) {
      if (is.null(gene_text) || gene_text == "") return(NULL)
      
      tryCatch({
        genes <- parse_gene_input(gene_text)
        
        if (length(genes) == 0) return(NULL)
        
        # Create simulated data for manual input
        set.seed(456)
        manual_data <- data.frame(
          symbol = genes,
          log2FoldChange = rnorm(length(genes), 0, 1.5),
          padj = runif(length(genes), 0, 0.1),
          pvalue = runif(length(genes), 0, 0.05),
          baseMean = runif(length(genes), 100, 1000),
          stringsAsFactors = FALSE
        )
        manual_data$minus_log10_padj <- -log10(manual_data$padj)
        
        return(manual_data)
      }, error = function(e) {
        warning("Error creating manual data: ", e$message)
        return(NULL)
      })
    }
    
    # Load data from URL
    load_url_data <- function(url) {
      if (is.null(url) || url == "") return(NULL)
      
      tryCatch({
        # Download the file
        temp_file <- tempfile()
        download.file(url, temp_file, quiet = TRUE)
        
        # Determine file type and read
        if (grepl("\\.csv$", url, ignore.case = TRUE)) {
          data <- read.csv(temp_file, stringsAsFactors = FALSE)
        } else if (grepl("\\.tsv$", url, ignore.case = TRUE)) {
          data <- read.delim(temp_file, stringsAsFactors = FALSE)
        } else {
          # Try CSV as default
          data <- read.csv(temp_file, stringsAsFactors = FALSE)
        }
        
        # Clean up
        unlink(temp_file)
        
        # Add symbol column if missing using ENSEMBL ID conversion
        if (!"symbol" %in% colnames(data)) {
          if ("gene" %in% colnames(data) && all(grepl("^ENSG", data$gene))) {
            message("Converting ENSEMBL IDs in URL file to gene symbols")
            data <- add_symbol_column(data, "gene")
          } else {
            # Use first column as symbol
            data$symbol <- data[[1]]
          }
        }
        
        # Standardize columns
        data <- standardize_de_columns(data)
        
        return(data)
      }, error = function(e) {
        warning("Error loading data from URL: ", e$message)
        return(NULL)
      })
    }
    
    # Main data loading function
    load_data <- function() {
      source_type <- input$data_source
      new_data <- NULL
      
      tryCatch({
        switch(source_type,
          "preloaded" = {
            new_data <- load_preloaded_data()
          },
          "upload" = {
            new_data <- load_uploaded_data(input$file_upload)
          },
          "manual" = {
            new_data <- create_manual_data(input$manual_genes)
          },
          "url" = {
            new_data <- load_url_data(input$data_url)
          }
        )
        
        if (is.null(new_data)) {
          module_rv$last_error <- "Failed to load data"
          return(FALSE)
        }
        
        # Validate the data
        if (!validate_de_data(new_data)) {
          module_rv$last_error <- "Data validation failed"
          return(FALSE)
        }
        
        # Store the data
        module_rv$raw_data <- new_data
        module_rv$filtered_data <- new_data
        module_rv$data_loaded <- TRUE
        module_rv$last_error <- NULL
        
        # Update global reactive values
        rv$data <- new_data
        rv$filtered_data <- new_data
        rv$data_loaded <- TRUE
        
        message("Data loaded successfully: ", nrow(new_data), " rows")
        return(TRUE)
        
      }, error = function(e) {
        module_rv$last_error <- paste("Data loading error:", e$message)
        warning(module_rv$last_error)
        return(FALSE)
      })
    }
    
    # Apply filters to data
    apply_data_filters <- function() {
      if (is.null(module_rv$raw_data)) return(NULL)
      
      tryCatch({
        filtered <- module_rv$raw_data
        
        # Apply fold change filter if column exists
        if ("log2FoldChange" %in% colnames(filtered)) {
          filtered <- filtered[
            filtered$log2FoldChange >= input$fc_filter[1] & 
            filtered$log2FoldChange <= input$fc_filter[2], 
          ]
        }
        
        # Apply p-value filter if column exists
        if ("padj" %in% colnames(filtered)) {
          filtered <- filtered[filtered$padj <= input$pval_filter, ]
        } else if ("pvalue" %in% colnames(filtered)) {
          filtered <- filtered[filtered$pvalue <= input$pval_filter, ]
        }
        
        module_rv$filtered_data <- filtered
        rv$filtered_data <- filtered
        
        message("Filters applied: ", nrow(filtered), " rows remaining")
        return(TRUE)
        
      }, error = function(e) {
        warning("Error applying filters: ", e$message)
        return(FALSE)
      })
    }
    
    # Event observers
    observeEvent(input$load_data, {
      showNotification("Loading data...", type = "message")
      success <- load_data()
      
      if (success) {
        showNotification("Data loaded successfully!", type = "message", duration = 3)
      } else {
        showNotification(
          paste("Error loading data:", module_rv$last_error), 
          type = "error", 
          duration = 5
        )
      }
    })
    
    observeEvent(input$reset_data, {
      module_rv$raw_data <- NULL
      module_rv$filtered_data <- NULL
      module_rv$data_loaded = FALSE
      module_rv$last_error <- NULL
      
      rv$data <- NULL
      rv$filtered_data <- NULL
      rv$data_loaded <- FALSE
      
      showNotification("Data reset to default", type = "warning")
    })
    
    observeEvent(input$apply_filters, {
      success <- apply_data_filters()
      
      if (success) {
        showNotification("Filters applied successfully!", type = "message", duration = 3)
      } else {
        showNotification("Error applying filters", type = "error", duration = 3)
      }
    })
    
    observeEvent(input$load_url, {
      if (!is.null(input$data_url) && input$data_url != "") {
        showNotification("Loading data from URL...", type = "message")
        success <- load_data()
        
        if (success) {
          showNotification("Data loaded from URL successfully!", type = "message", duration = 3)
        } else {
          showNotification(
            paste("Error loading from URL:", module_rv$last_error), 
            type = "error", 
            duration = 5
          )
        }
      }
    })
    
    # Outputs
    output$data_loaded <- reactive({
      module_rv$data_loaded
    })
    outputOptions(output, "data_loaded", suspendWhenHidden = FALSE)
    
    output$preloaded_info <- renderText({
      if (exists("diffgenes") && !is.null(diffgenes)) {
        paste("Pre-loaded data:", nrow(diffgenes), "genes")
      } else {
        "Using example data for demonstration"
      }
    })
    
    output$data_summary <- renderText({
      if (module_rv$data_loaded) {
        data <- module_rv$filtered_data
        paste(
          "Rows:", nrow(data), 
          "| Columns:", ncol(data),
          "| Genes:", length(unique(data$symbol))
        )
      } else {
        "No data loaded"
      }
    })
    
    output$data_preview <- DT::renderDataTable({
      if (module_rv$data_loaded) {
        preview_data <- create_summary_df(module_rv$filtered_data, max_rows = 10)
        DT::datatable(
          preview_data,
          options = list(
            scrollX = TRUE,
            pageLength = 10,
            dom = 'tip'
          ),
          rownames = FALSE
        )
      }
    })
    
    output$data_stats <- renderPrint({
      if (module_rv$data_loaded) {
        data <- module_rv$filtered_data
        cat("DATA SUMMARY\n")
        cat("============\n")
        cat("Total genes:", nrow(data), "\n")
        cat("Total columns:", ncol(data), "\n\n")
        
        if ("log2FoldChange" %in% colnames(data)) {
          fc_stats <- calculate_fc_stats(data$log2FoldChange)
          cat("FOLD CHANGE STATISTICS\n")
          cat("Mean:", round(fc_stats$mean, 3), "\n")
          cat("Median:", round(fc_stats$median, 3), "\n")
          cat("SD:", round(fc_stats$sd, 3), "\n")
          cat("Up-regulated:", fc_stats$up_regulated, "\n")
          cat("Down-regulated:", fc_stats$down_regulated, "\n\n")
        }
        
        if ("padj" %in% colnames(data)) {
          sig_genes <- sum(data$padj < 0.05, na.rm = TRUE)
          cat("Significant genes (padj < 0.05):", sig_genes, "\n")
        }
      } else {
        cat("No data loaded\n")
      }
    })
    
    output$column_info <- renderPrint({
      if (module_rv$data_loaded) {
        data <- module_rv$filtered_data
        cat("COLUMN INFORMATION\n")
        cat("==================\n")
        for (col in colnames(data)) {
          cat(col, ":", class(data[[col]]), "\n")
          if (is.numeric(data[[col]])) {
            cat("  Range: [", round(min(data[[col]], na.rm = TRUE), 3), 
                ", ", round(max(data[[col]], na.rm = TRUE), 3), "]\n", sep = "")
            cat("  NA values:", sum(is.na(data[[col]])), "\n")
          }
          cat("\n")
        }
      } else {
        cat("No data loaded\n")
      }
    })
    
    # Return reactive values for other modules to use
    return(reactive({
      list(
        data = module_rv$raw_data,
        filtered_data = module_rv$filtered_data,
        data_loaded = module_rv$data_loaded,
        last_error = module_rv$last_error
      )
    }))
  })
}