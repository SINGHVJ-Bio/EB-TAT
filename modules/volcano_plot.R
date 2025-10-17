# ENABL Biomarker/Target Analysis Tool - Volcano Plot Module
# Creates interactive volcano plots for differential expression data

# =============================================================================
# UI Function
# =============================================================================

volcanoPlotUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    div(
      class = "volcano-plot-module",
      style = "padding: 20px;",
      
      # Header
      div(
        class = "module-header",
        h3("Volcano Plot", style = "color: #2196F3;"),
        p("Visualize differential expression results with interactive volcano plots")
      ),
      
      fluidRow(
        # Controls column
        column(
          4,
          # Plot configuration
          div(
            class = "card",
            style = "margin-bottom: 20px;",
            div(
              class = "card-header",
              h4("Plot Configuration", style = "margin: 0;")
            ),
            div(
              class = "card-body",
              # Threshold controls
              sliderInput(
                ns("fc_threshold"),
                "Fold Change Threshold:",
                min = 0,
                max = 5,
                value = 1,
                step = 0.1
              ),
              sliderInput(
                ns("p_threshold"),
                "P-value Threshold (-log10):",
                min = 0,
                max = 10,
                value = 1.3,
                step = 0.1
              ),
              numericInput(
                ns("point_size"),
                "Point Size:",
                value = 2,
                min = 1,
                max = 10,
                step = 0.5
              ),
              numericInput(
                ns("point_alpha"),
                "Point Transparency:",
                value = 0.7,
                min = 0.1,
                max = 1,
                step = 0.1
              ),
              selectInput(
                ns("color_scheme"),
                "Color Scheme:",
                choices = c(
                  "Default" = "default",
                  "Colorblind Friendly" = "colorblind",
                  "High Contrast" = "high_contrast"
                ),
                selected = "default"
              )
            )
          ),
          
          # Gene labeling
          div(
            class = "card",
            style = "margin-bottom: 20px;",
            div(
              class = "card-header",
              h4("Gene Labeling", style = "margin: 0;")
            ),
            div(
              class = "card-body",
              radioButtons(
                ns("label_type"),
                "Label genes by:",
                choices = list(
                  "None" = "none",
                  "Top N significant" = "top_n",
                  "Custom list" = "custom",
                  "All significant" = "all_sig"
                ),
                selected = "top_n"
              ),
              conditionalPanel(
                condition = paste0("input['", ns("label_type"), "'] == 'top_n'"),
                numericInput(
                  ns("top_n_genes"),
                  "Number of top genes to label:",
                  value = 10,
                  min = 1,
                  max = 50,
                  step = 1
                )
              ),
              conditionalPanel(
                condition = paste0("input['", ns("label_type"), "'] == 'custom'"),
                textAreaInput(
                  ns("custom_genes"),
                  "Enter gene symbols to label:",
                  rows = 5,
                  placeholder = "TP53\nEGFR\nBRCA1\n..."
                )
              ),
              checkboxInput(
                ns("avoid_overlap"),
                "Avoid label overlap",
                value = TRUE
              )
            )
          ),
          
          # Download options
          div(
            class = "card",
            div(
              class = "card-header",
              h4("Export", style = "margin: 0;")
            ),
            div(
              class = "card-body",
              numericInput(
                ns("plot_width"),
                "Plot Width (inches):",
                value = 10,
                min = 4,
                max = 20,
                step = 1
              ),
              numericInput(
                ns("plot_height"),
                "Plot Height (inches):",
                value = 8,
                min = 4,
                max = 20,
                step = 1
              ),
              numericInput(
                ns("plot_dpi"),
                "Resolution (DPI):",
                value = 300,
                min = 72,
                max = 600,
                step = 50
              ),
              selectInput(
                ns("plot_format"),
                "File Format:",
                choices = c("PNG" = "png", "PDF" = "pdf", "SVG" = "svg"),
                selected = "png"
              ),
              downloadButton(
                ns("download_plot"),
                "Download Plot",
                class = "btn-primary"
              )
            )
          )
        ),
        
        # Plot column
        column(
          8,
          # Plot output
          div(
            class = "card",
            style = "margin-bottom: 20px;",
            div(
              class = "card-header",
              h4("Volcano Plot", style = "margin: 0;")
            ),
            div(
              class = "card-body",
              # Status message
              conditionalPanel(
                condition = paste0("!output['", ns("data_available"), "']"),
                div(
                  class = "alert alert-warning",
                  strong("No data available."),
                  "Please load data in the Data Input tab."
                )
              ),
              # Plot with loading spinner
              withSpinner(
                plotOutput(
                  ns("volcano_plot"),
                  height = "600px",
                  click = ns("plot_click")
                ),
                type = 4,
                color = "#2196F3"
              )
            )
          ),
          
          # Plot information and interactions
          div(
            class = "card",
            div(
              class = "card-header",
              h4("Plot Information", style = "margin: 0;")
            ),
            div(
              class = "card-body",
              fluidRow(
                column(
                  6,
                  h5("Statistics:"),
                  verbatimTextOutput(ns("plot_stats"))
                ),
                column(
                  6,
                  h5("Clicked Point:"),
                  verbatimTextOutput(ns("click_info"))
                )
              ),
              # Selected genes table
              conditionalPanel(
                condition = paste0("output['", ns("show_selected_table"), "']"),
                hr(),
                h5("Selected Genes:"),
                DT::dataTableOutput(ns("selected_genes_table"))
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

volcanoPlotServer <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Reactive values for this module
    module_rv <- reactiveValues(
      plot_data = NULL,
      current_plot = NULL,
      selected_genes = NULL,
      clicked_point = NULL
    )
    
    # Prepare data for plotting
    prepare_plot_data <- reactive({
      req(rv$data_loaded)
      req(rv$filtered_data)
      
      data <- rv$filtered_data
      
      # Ensure we have required columns
      if (!"log2FoldChange" %in% colnames(data) || 
          (!"padj" %in% colnames(data) && !"pvalue" %in% colnames(data))) {
        return(NULL)
      }
      
      # Calculate -log10 p-value if not present
      if (!"minus_log10_padj" %in% colnames(data)) {
        if ("padj" %in% colnames(data)) {
          data$minus_log10_padj <- safe_minus_log10(data$padj)
        } else if ("pvalue" %in% colnames(data)) {
          data$minus_log10_padj <- safe_minus_log10(data$pvalue)
        }
      }
      
      # Add significance classification
      data <- prepare_volcano_data(
        data, 
        input$fc_threshold, 
        10^(-input$p_threshold)  # Convert back from -log10
      )
      
      return(data)
    })
    
    # Get genes to label
    get_genes_to_label <- reactive({
      plot_data <- prepare_plot_data()
      if (is.null(plot_data)) return(character(0))
      
      label_type <- input$label_type
      
      if (label_type == "none") {
        return(character(0))
      } else if (label_type == "top_n") {
        # Get top N significant genes by p-value
        sig_data <- plot_data[plot_data$significance != "Not Significant", ]
        if (nrow(sig_data) == 0) return(character(0))
        
        # Order by significance (p-value)
        sig_data <- sig_data[order(sig_data$padj %||% sig_data$pvalue), ]
        top_genes <- head(sig_data$symbol, input$top_n_genes)
        return(top_genes)
      } else if (label_type == "custom") {
        # Parse custom gene list
        genes <- parse_gene_input(input$custom_genes)
        # Only return genes that are in the data
        return(intersect(genes, plot_data$symbol))
      } else if (label_type == "all_sig") {
        # All significant genes
        sig_data <- plot_data[plot_data$significance != "Not Significant", ]
        return(sig_data$symbol)
      }
    })
    
    # Create volcano plot
    create_volcano_plot_reactive <- reactive({
      plot_data <- prepare_plot_data()
      if (is.null(plot_data)) return(NULL)
      
      # Get genes to label
      genes_to_label <- get_genes_to_label()
      
      # Create the plot
      p <- create_volcano_plot(
        de_data = plot_data,
        fc_threshold = input$fc_threshold,
        p_threshold = 10^(-input$p_threshold),  # Convert back from -log10
        highlight_genes = genes_to_label,
        title = "Volcano Plot",
        point_size = input$point_size,
        alpha = input$point_alpha
      )
      
      # Apply color scheme
      p <- apply_color_scheme(p, input$color_scheme)
      
      # Store the current plot
      module_rv$current_plot <- p
      
      return(p)
    })
    
    # Apply color scheme to plot
    apply_color_scheme <- function(plot_obj, scheme) {
      if (scheme == "colorblind") {
        plot_obj <- plot_obj + 
          scale_color_manual(
            values = c(
              "Not Significant" = "gray70",
              "Up-regulated" = "#E69F00",  # Orange
              "Down-regulated" = "#56B4E9"  # Blue
            )
          )
      } else if (scheme == "high_contrast") {
        plot_obj <- plot_obj + 
          scale_color_manual(
            values = c(
              "Not Significant" = "gray50",
              "Up-regulated" = "#FF0000",  # Red
              "Down-regulated" = "#0000FF"  # Blue
            )
          )
      }
      # Default scheme is already applied in create_volcano_plot
      return(plot_obj)
    }
    
    # Calculate plot statistics
    calculate_plot_statistics <- reactive({
      plot_data <- prepare_plot_data()
      if (is.null(plot_data)) return(NULL)
      
      total_genes <- nrow(plot_data)
      up_regulated <- sum(plot_data$significance == "Up-regulated", na.rm = TRUE)
      down_regulated <- sum(plot_data$significance == "Down-regulated", na.rm = TRUE)
      not_sig <- sum(plot_data$significance == "Not Significant", na.rm = TRUE)
      
      stats <- list(
        total_genes = total_genes,
        up_regulated = up_regulated,
        down_regulated = down_regulated,
        not_significant = not_sig,
        percent_sig = round((up_regulated + down_regulated) / total_genes * 100, 1)
      )
      
      return(stats)
    })
    
    # Handle plot clicks
    observeEvent(input$plot_click, {
      plot_data <- prepare_plot_data()
      if (is.null(plot_data)) return()
      
      click <- input$plot_click
      
      # Find the closest point
      distances <- sqrt(
        (plot_data$log2FoldChange - click$x)^2 + 
        ((plot_data$minus_log10_padj %||% 0) - click$y)^2
      )
      
      closest_idx <- which.min(distances)
      closest_gene <- plot_data[closest_idx, ]
      
      module_rv$clicked_point <- closest_gene
      
      # Add to selected genes
      if (!closest_gene$symbol %in% module_rv$selected_genes) {
        module_rv$selected_genes <- c(module_rv$selected_genes, closest_gene$symbol)
      }
    })
    
    # Outputs
    output$data_available <- reactive({
      rv$data_loaded && !is.null(rv$filtered_data)
    })
    outputOptions(output, "data_available", suspendWhenHidden = FALSE)
    
    output$show_selected_table <- reactive({
      !is.null(module_rv$selected_genes) && length(module_rv$selected_genes) > 0
    })
    outputOptions(output, "show_selected_table", suspendWhenHidden = FALSE)
    
    output$volcano_plot <- renderPlot({
      if (!rv$data_loaded) {
        return(create_empty_plot("Please load data in the Data Input tab"))
      }
      
      plot_obj <- create_volcano_plot_reactive()
      if (is.null(plot_obj)) {
        return(create_empty_plot("Error creating volcano plot"))
      }
      
      return(plot_obj)
    })
    
    output$plot_stats <- renderPrint({
      stats <- calculate_plot_statistics()
      if (is.null(stats)) {
        cat("No data available\n")
        return()
      }
      
      cat("PLOT STATISTICS\n")
      cat("===============\n")
      cat("Total genes:", stats$total_genes, "\n")
      cat("Up-regulated:", stats$up_regulated, "\n")
      cat("Down-regulated:", stats$down_regulated, "\n")
      cat("Not significant:", stats$not_significant, "\n")
      cat("Significant genes:", stats$up_regulated + stats$down_regulated, 
          paste0("(", stats$percent_sig, "%)"), "\n")
      cat("\nTHRESHOLDS\n")
      cat("Fold change: ±", input$fc_threshold, "\n")
      cat("P-value: ", format(10^(-input$p_threshold), scientific = TRUE), "\n")
    })
    
    output$click_info <- renderPrint({
      gene <- module_rv$clicked_point
      if (is.null(gene)) {
        cat("Click on a point to see gene information\n")
        return()
      }
      
      cat("GENE INFORMATION\n")
      cat("================\n")
      cat("Symbol:", gene$symbol, "\n")
      if (!is.null(gene$log2FoldChange)) {
        cat("Log2 FC:", round(gene$log2FoldChange, 3), "\n")
      }
      if (!is.null(gene$padj)) {
        cat("Adj. p-value:", format(gene$padj, scientific = TRUE), "\n")
      } else if (!is.null(gene$pvalue)) {
        cat("P-value:", format(gene$pvalue, scientific = TRUE), "\n")
      }
      if (!is.null(gene$significance)) {
        cat("Significance:", gene$significance, "\n")
      }
    })
    
    output$selected_genes_table <- DT::renderDataTable({
      if (is.null(module_rv$selected_genes) || length(module_rv$selected_genes) == 0) {
        return(NULL)
      }
      
      plot_data <- prepare_plot_data()
      if (is.null(plot_data)) return(NULL)
      
      # Get data for selected genes
      selected_data <- plot_data[plot_data$symbol %in% module_rv$selected_genes, ]
      
      # Create display table
      display_cols <- c("symbol", "log2FoldChange", "padj", "pvalue", "significance")
      display_cols <- display_cols[display_cols %in% colnames(selected_data)]
      
      display_data <- selected_data[, display_cols, drop = FALSE]
      
      # Format numeric columns
      if ("log2FoldChange" %in% colnames(display_data)) {
        display_data$log2FoldChange <- round(display_data$log2FoldChange, 3)
      }
      if ("padj" %in% colnames(display_data)) {
        display_data$padj <- format_pvalue(display_data$padj)
      }
      if ("pvalue" %in% colnames(display_data)) {
        display_data$pvalue <- format_pvalue(display_data$pvalue)
      }
      
      DT::datatable(
        display_data,
        options = list(
          scrollX = TRUE,
          pageLength = 10,
          dom = 'tip'
        ),
        rownames = FALSE,
        selection = 'none'
      )
    })
    
    # Download handler
    output$download_plot <- downloadHandler(
      filename = function() {
        paste0(
          "volcano_plot_",
          format(Sys.time(), "%Y%m%d_%H%M%S"),
          ".",
          input$plot_format
        )
      },
      content = function(file) {
        plot_obj <- module_rv$current_plot
        if (is.null(plot_obj)) {
          showNotification("No plot to download", type = "error")
          return()
        }
        
        success <- save_plot(
          plot_obj,
          file,
          width = input$plot_width,
          height = input$plot_height,
          dpi = input$plot_dpi
        )
        
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
        selected_genes = module_rv$selected_genes,
        plot_data = prepare_plot_data(),
        plot_statistics = calculate_plot_statistics()
      )
    }))
  })
}