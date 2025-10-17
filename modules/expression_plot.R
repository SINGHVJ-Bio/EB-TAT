# ENABL Biomarker/Target Analysis Tool - Expression Plot Module
# Focuses on gene expression visualization and detailed gene-level analysis

# =============================================================================
# UI Function
# =============================================================================

expressionPlotUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    div(
      class = "expression-plot-module",
      style = "padding: 20px;",
      
      # Header
      div(
        class = "module-header",
        h3("Expression Analysis", style = "color: #2196F3;"),
        p("Detailed gene expression visualization and analysis")
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
                ns("user_gene_list2"),
                "Enter gene symbols (one per line or comma-separated):",
                rows = 8,
                placeholder = "TP53\nEGFR\nBRCA1\n..."
              ),
              actionButton(
                ns("load_genes"),
                "Load Genes",
                class = "btn-primary"
              ),
              actionButton(
                ns("clear_genes"),
                "Clear",
                class = "btn-warning"
              ),
              hr(),
              selectInput(
                ns("gene_source"),
                "Or select from:",
                choices = c(
                  "Top significant genes" = "top_sig",
                  "Up-regulated genes" = "top_up",
                  "Down-regulated genes" = "top_down",
                  "Custom list" = "custom"
                ),
                selected = "top_sig"
              ),
              conditionalPanel(
                condition = paste0("input['", ns("gene_source"), "'] != 'custom'"),
                numericInput(
                  ns("gene_count"),
                  "Number of genes:",
                  value = 10,
                  min = 1,
                  max = 50,
                  step = 1
                )
              )
            )
          ),
          
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
              selectInput(
                ns("plot_type"),
                "Plot Type:",
                choices = c(
                  "Bar Plot" = "bar",
                  "Dot Plot" = "dot",
                  "Heatmap" = "heatmap",
                  "Violin Plot" = "violin",
                  "Box Plot" = "box"
                ),
                selected = "bar"
              ),
              selectInput(
                ns("value_column"),
                "Value to Plot:",
                choices = c(
                  "Log2 Fold Change" = "log2FoldChange",
                  "P-value" = "pvalue",
                  "Adjusted P-value" = "padj",
                  "Base Mean" = "baseMean"
                ),
                selected = "log2FoldChange"
              ),
              conditionalPanel(
                condition = paste0("input['", ns("plot_type"), "'] == 'bar' || input['", ns("plot_type"), "'] == 'dot'"),
                checkboxInput(
                  ns("show_pvalues"),
                  "Show p-values",
                  value = TRUE
                )
              ),
              conditionalPanel(
                condition = paste0("input['", ns("plot_type"), "'] == 'heatmap'"),
                checkboxInput(
                  ns("cluster_genes"),
                  "Cluster genes",
                  value = TRUE
                ),
                checkboxInput(
                  ns("cluster_samples"),
                  "Cluster samples",
                  value = FALSE
                )
              ),
              numericInput(
                ns("plot_width_exp"),
                "Plot Width:",
                value = 12,
                min = 6,
                max = 20,
                step = 1
              ),
              numericInput(
                ns("plot_height_exp"),
                "Plot Height:",
                value = 8,
                min = 6,
                max = 20,
                step = 1
              )
            )
          ),
          
          # Display options
          div(
            class = "card",
            style = "margin-bottom: 20px;",
            div(
              class = "card-header",
              h4("Display Options", style = "margin: 0;")
            ),
            div(
              class = "card-body",
              selectInput(
                ns("color_palette"),
                "Color Palette:",
                choices = c(
                  "Default" = "default",
                  "Colorblind Friendly" = "colorblind",
                  "Red-Blue" = "red_blue",
                  "Viridis" = "viridis"
                ),
                selected = "default"
              ),
              checkboxInput(
                ns("show_labels"),
                "Show gene labels",
                value = TRUE
              ),
              checkboxInput(
                ns("sort_genes"),
                "Sort by expression",
                value = TRUE
              ),
              numericInput(
                ns("font_size"),
                "Font Size:",
                value = 12,
                min = 8,
                max = 20,
                step = 1
              )
            )
          ),
          
          # Export options
          div(
            class = "card",
            div(
              class = "card-header",
              h4("Export", style = "margin: 0;")
            ),
            div(
              class = "card-body",
              downloadButton(
                ns("download_plot_exp"),
                "Download Plot",
                class = "btn-success"
              ),
              downloadButton(
                ns("download_data_exp"),
                "Download Data",
                class = "btn-info"
              )
            )
          )
        ),
        
        # Results column
        column(
          8,
          # Status and gene info
          div(
            class = "card",
            style = "margin-bottom: 20px;",
            div(
              class = "card-header",
              h4("Gene Information", style = "margin: 0;")
            ),
            div(
              class = "card-body",
              conditionalPanel(
                condition = paste0("!output['", ns("data_available_exp"), "']"),
                div(
                  class = "alert alert-warning",
                  strong("No data available."),
                  "Please load data in the Data Input tab."
                )
              ),
              conditionalPanel(
                condition = paste0("output['", ns("genes_loaded"), "']"),
                div(
                  class = "alert alert-success",
                  strong("Genes loaded successfully!"),
                  textOutput(ns("gene_summary"))
                )
              ),
              verbatimTextOutput(ns("gene_info"))
            )
          ),
          
          # Main plot
          div(
            class = "card",
            style = "margin-bottom: 20px;",
            div(
              class = "card-header",
              h4("Expression Plot", style = "margin: 0;")
            ),
            div(
              class = "card-body",
              withSpinner(
                plotOutput(ns("expression_plot"), height = "500px"),
                type = 4,
                color = "#2196F3"
              )
            )
          ),
          
          # Gene details table
          div(
            class = "card",
            style = "margin-bottom: 20px;",
            div(
              class = "card-header",
              h4("Gene Details", style = "margin: 0;")
            ),
            div(
              class = "card-body",
              DT::dataTableOutput(ns("gene_details_table"))
            )
          ),
          
          # Additional visualizations
          div(
            class = "card",
            div(
              class = "card-header",
              h4("Additional Views", style = "margin: 0;")
            ),
            div(
              class = "card-body",
              tabsetPanel(
                type = "tabs",
                tabPanel(
                  "Expression Distribution",
                  withSpinner(
                    plotOutput(ns("distribution_plot_exp"), height = "400px"),
                    type = 4,
                    color = "#2196F3"
                  )
                ),
                tabPanel(
                  "Gene Comparison",
                  withSpinner(
                    plotOutput(ns("comparison_plot"), height = "400px"),
                    type = 4,
                    color = "#2196F3"
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

expressionPlotServer <- function(id, rv, rv_expression) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Reactive values for this module
    module_rv <- reactiveValues(
      selected_genes = NULL,
      expression_data = NULL,
      plot_object = NULL,
      genes_loaded = FALSE
    )
    
    # Get genes based on selection method
    get_genes_from_source <- function() {
      req(rv$data_loaded)
      req(rv$filtered_data)
      
      data <- rv$filtered_data
      source_type <- input$gene_source
      gene_count <- input$gene_count
      
      if (source_type == "custom") {
        # Use genes from text input (handled separately)
        return(NULL)
      }
      
      tryCatch({
        if (source_type == "top_sig") {
          # Top significant genes by p-value
          if ("padj" %in% colnames(data)) {
            sig_data <- data[!is.na(data$padj) & data$padj < 0.05, ]
            if (nrow(sig_data) > 0) {
              sig_data <- sig_data[order(sig_data$padj), ]
              return(head(sig_data$symbol, gene_count))
            }
          } else if ("pvalue" %in% colnames(data)) {
            sig_data <- data[!is.na(data$pvalue) & data$pvalue < 0.05, ]
            if (nrow(sig_data) > 0) {
              sig_data <- sig_data[order(sig_data$pvalue), ]
              return(head(sig_data$symbol, gene_count))
            }
          }
        } else if (source_type == "top_up") {
          # Top up-regulated genes
          if ("log2FoldChange" %in% colnames(data)) {
            up_data <- data[data$log2FoldChange > 0, ]
            if (nrow(up_data) > 0) {
              up_data <- up_data[order(up_data$log2FoldChange, decreasing = TRUE), ]
              return(head(up_data$symbol, gene_count))
            }
          }
        } else if (source_type == "top_down") {
          # Top down-regulated genes
          if ("log2FoldChange" %in% colnames(data)) {
            down_data <- data[data$log2FoldChange < 0, ]
            if (nrow(down_data) > 0) {
              down_data <- down_data[order(down_data$log2FoldChange), ]
              return(head(down_data$symbol, gene_count))
            }
          }
        }
        
        # Fallback: return first N genes
        return(head(data$symbol, gene_count))
        
      }, error = function(e) {
        warning("Error getting genes from source: ", e$message)
        return(NULL)
      })
    }
    
    # Load genes from user input
    load_user_genes <- function() {
      gene_text <- input$user_gene_list2
      if (is.null(gene_text) || gene_text == "") return(NULL)
      
      genes <- parse_gene_input(gene_text)
      return(genes)
    }
    
    # Load genes based on current selection
    load_genes <- function() {
      genes <- NULL
      
      if (input$gene_source == "custom") {
        genes <- load_user_genes()
      } else {
        genes <- get_genes_from_source()
      }
      
      if (is.null(genes) || length(genes) == 0) {
        showNotification("No genes to load", type = "warning")
        return(FALSE)
      }
      
      # Filter to genes that exist in the data
      data <- rv$filtered_data
      existing_genes <- genes[genes %in% data$symbol]
      
      if (length(existing_genes) == 0) {
        showNotification("None of the specified genes found in data", type = "error")
        return(FALSE)
      }
      
      module_rv$selected_genes <- existing_genes
      module_rv$genes_loaded <- TRUE
      
      # Get expression data for selected genes
      module_rv$expression_data <- data[data$symbol %in% existing_genes, ]
      
      msg <- paste("Loaded", length(existing_genes), "genes")
      if (length(existing_genes) < length(genes)) {
        msg <- paste(msg, "(", length(genes) - length(existing_genes), "not found)")
      }
      
      showNotification(msg, type = "message", duration = 3)
      return(TRUE)
    }
    
    # Create main expression plot
    create_expression_plot <- function() {
      if (!module_rv$genes_loaded || is.null(module_rv$expression_data)) {
        return(create_empty_plot("Please load genes first"))
      }
      
      data <- module_rv$expression_data
      plot_type <- input$plot_type
      value_col <- input$value_column
      
      if (!value_col %in% colnames(data)) {
        return(create_empty_plot(paste("Column not found:", value_col)))
      }
      
      tryCatch({
        if (plot_type == "bar") {
          p <- create_expression_barplot(data, value_col)
        } else if (plot_type == "dot") {
          p <- create_expression_dotplot(data, input$show_pvalues)
        } else if (plot_type == "heatmap") {
          p <- create_expression_heatmap_simple(data, value_col)
        } else if (plot_type == "violin") {
          p <- create_violin_plot(data, value_col)
        } else if (plot_type == "box") {
          p <- create_box_plot(data, value_col)
        } else {
          p <- create_expression_barplot(data, value_col)
        }
        
        # Apply color palette
        p <- apply_color_palette(p, input$color_palette)
        
        # Apply font size
        p <- p + theme(text = element_text(size = input$font_size))
        
        # Store plot object for downloading
        module_rv$plot_object <- p
        
        return(p)
        
      }, error = function(e) {
        warning("Error creating expression plot: ", e$message)
        return(create_empty_plot("Error creating expression plot"))
      })
    }
    
    # Create expression bar plot
    create_expression_barplot <- function(data, value_col) {
      # Sort data if requested
      if (input$sort_genes) {
        data <- data[order(data[[value_col]]), ]
      }
      data$symbol <- factor(data$symbol, levels = data$symbol)
      
      p <- ggplot(data, aes_string(x = "symbol", y = value_col, fill = value_col)) +
        geom_col() +
        coord_flip() +
        theme_ebtat() +
        labs(
          title = paste("Gene Expression -", value_col),
          x = "Gene Symbol",
          y = value_col
        ) +
        theme(
          axis.text.y = element_text(size = input$font_size - 2)
        )
      
      # Add value labels
      if (input$show_labels) {
        p <- p + geom_text(
          aes(label = round(!!sym(value_col), 3)),
          hjust = -0.1,
          size = input$font_size / 4
        )
      }
      
      return(p)
    }
    
    # Create expression dot plot
    create_expression_dotplot <- function(data, show_pvalues = TRUE) {
      # Sort data
      if (input$sort_genes) {
        data <- data[order(data$log2FoldChange), ]
      }
      data$symbol <- factor(data$symbol, levels = data$symbol)
      
      # Add significance if requested
      if (show_pvalues && "padj" %in% colnames(data)) {
        data$significance <- cut(data$padj,
                                breaks = c(0, 0.001, 0.01, 0.05, 1),
                                labels = c("***", "**", "*", ""),
                                include.lowest = TRUE)
      } else {
        data$significance <- ""
      }
      
      p <- ggplot(data, aes(x = symbol, y = log2FoldChange, color = log2FoldChange)) +
        geom_point(size = 4) +
        coord_flip() +
        theme_ebtat() +
        labs(
          title = "Gene Expression Dot Plot",
          x = "Gene Symbol",
          y = "Log2 Fold Change"
        ) +
        theme(
          axis.text.y = element_text(size = input$font_size - 2)
        )
      
      # Add significance stars
      if (show_pvalues) {
        p <- p + geom_text(
          aes(label = significance),
          vjust = -1.5,
          color = "black",
          size = input$font_size / 2
        )
      }
      
      # Add value labels
      if (input$show_labels) {
        p <- p + geom_text(
          aes(label = round(log2FoldChange, 3)),
          hjust = -0.5,
          size = input$font_size / 4
        )
      }
      
      return(p)
    }
    
    # Create simple heatmap for expression data
    create_expression_heatmap_simple <- function(data, value_col) {
      # Create a matrix-like structure for heatmap
      heatmap_data <- data.frame(
        Gene = data$symbol,
        Value = data[[value_col]]
      )
      
      p <- ggplot(heatmap_data, aes(x = 1, y = Gene, fill = Value)) +
        geom_tile() +
        scale_fill_gradient2(
          low = "blue",
          mid = "white",
          high = "red",
          midpoint = 0,
          name = value_col
        ) +
        theme_ebtat() +
        labs(
          title = paste("Expression Heatmap -", value_col),
          x = "",
          y = "Gene Symbol"
        ) +
        theme(
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = input$font_size - 2)
        )
      
      # Add value labels
      if (input$show_labels) {
        p <- p + geom_text(
          aes(label = round(Value, 2)),
          color = "black",
          size = input$font_size / 4
        )
      }
      
      return(p)
    }
    
    # Create violin plot
    create_violin_plot <- function(data, value_col) {
      # For violin plot, we need multiple values per gene (simulated here)
      # In a real scenario, this would use raw expression data
      plot_data <- data.frame(
        Gene = rep(data$symbol, each = 3),
        Value = rep(data[[value_col]], each = 3) + rnorm(nrow(data) * 3, 0, 0.1)
      )
      
      p <- ggplot(plot_data, aes(x = Gene, y = Value, fill = Gene)) +
        geom_violin(alpha = 0.7) +
        theme_ebtat() +
        labs(
          title = paste("Expression Distribution -", value_col),
          x = "Gene Symbol",
          y = value_col
        ) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = input$font_size - 2),
          legend.position = "none"
        )
      
      return(p)
    }
    
    # Create box plot
    create_box_plot <- function(data, value_col) {
      # Similar to violin plot but with boxplots
      plot_data <- data.frame(
        Gene = rep(data$symbol, each = 3),
        Value = rep(data[[value_col]], each = 3) + rnorm(nrow(data) * 3, 0, 0.1)
      )
      
      p <- ggplot(plot_data, aes(x = Gene, y = Value, fill = Gene)) +
        geom_boxplot(alpha = 0.7) +
        theme_ebtat() +
        labs(
          title = paste("Expression Box Plot -", value_col),
          x = "Gene Symbol",
          y = value_col
        ) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = input$font_size - 2),
          legend.position = "none"
        )
      
      return(p)
    }
    
    # Create distribution plot
    create_distribution_plot_exp <- function() {
      if (!module_rv$genes_loaded || is.null(module_rv$expression_data)) {
        return(create_empty_plot("Please load genes first"))
      }
      
      data <- module_rv$expression_data
      value_col <- input$value_column
      
      if (!value_col %in% colnames(data)) {
        return(create_empty_plot(paste("Column not found:", value_col)))
      }
      
      p <- ggplot(data, aes_string(x = value_col)) +
        geom_histogram(aes(y = ..density..), bins = 20, fill = "lightblue", alpha = 0.7) +
        geom_density(color = "darkblue", size = 1) +
        theme_ebtat() +
        labs(
          title = paste("Distribution of", value_col, "for Selected Genes"),
          x = value_col,
          y = "Density"
        ) +
        theme(text = element_text(size = input$font_size))
      
      return(p)
    }
    
    # Create comparison plot
    create_comparison_plot <- function() {
      if (!module_rv$genes_loaded || is.null(module_rv$expression_data)) {
        return(create_empty_plot("Please load genes first"))
      }
      
      data <- module_rv$expression_data
      
      if (!"log2FoldChange" %in% colnames(data) || !"padj" %in% colnames(data)) {
        return(create_empty_plot("Required columns not found"))
      }
      
      p <- ggplot(data, aes(x = log2FoldChange, y = -log10(padj), label = symbol)) +
        geom_point(aes(color = log2FoldChange, size = -log10(padj)), alpha = 0.7) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "blue") +
        theme_ebtat() +
        labs(
          title = "Selected Genes: Fold Change vs Significance",
          x = "Log2 Fold Change",
          y = "-Log10 Adjusted P-value"
        ) +
        theme(text = element_text(size = input$font_size))
      
      # Add labels if requested
      if (input$show_labels) {
        p <- p + ggrepel::geom_text_repel(
          size = input$font_size / 4,
          box.padding = 0.5,
          max.overlaps = 20
        )
      }
      
      return(p)
    }
    
    # Apply color palette to plot
    apply_color_palette <- function(plot_obj, palette) {
      if (palette == "colorblind") {
        plot_obj <- plot_obj + scale_fill_viridis_c() + scale_color_viridis_c()
      } else if (palette == "red_blue") {
        plot_obj <- plot_obj + 
          scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
          scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
      } else if (palette == "viridis") {
        plot_obj <- plot_obj + scale_fill_viridis_c() + scale_color_viridis_c()
      }
      # Default palette is applied in individual plot functions
      return(plot_obj)
    }
    
    # Update reactive values when inputs change
    observe({
      rv_expression$user_gene_list2 <- input$user_gene_list2
      rv_expression$plot_type <- input$plot_type
    })
    
    # Event observers
    observeEvent(input$load_genes, {
      load_genes()
    })
    
    observeEvent(input$clear_genes, {
      module_rv$selected_genes <- NULL
      module_rv$expression_data <- NULL
      module_rv$genes_loaded <- FALSE
      updateTextAreaInput(session, "user_gene_list2", value = "")
      showNotification("Genes cleared", type = "warning")
    })
    
    observeEvent(input$gene_source, {
      if (input$gene_source != "custom") {
        # Auto-load genes when source changes
        load_genes()
      }
    })
    
    # Outputs
    output$data_available_exp <- reactive({
      rv$data_loaded && !is.null(rv$filtered_data)
    })
    outputOptions(output, "data_available_exp", suspendWhenHidden = FALSE)
    
    output$genes_loaded <- reactive({
      module_rv$genes_loaded
    })
    outputOptions(output, "genes_loaded", suspendWhenHidden = FALSE)
    
    output$gene_summary <- renderText({
      if (!module_rv$genes_loaded) return("")
      
      data <- module_rv$expression_data
      paste(
        "Loaded", length(module_rv$selected_genes), "genes |",
        "Mean fold change:", round(mean(data$log2FoldChange, na.rm = TRUE), 3)
      )
    })
    
    output$gene_info <- renderPrint({
      if (!module_rv$genes_loaded) {
        cat("No genes loaded.\n")
        cat("Enter gene symbols or select from predefined lists.\n")
        return()
      }
      
      data <- module_rv$expression_data
      
      cat("LOADED GENE INFORMATION\n")
      cat("=======================\n")
      cat("Total genes:", nrow(data), "\n")
      
      if ("log2FoldChange" %in% colnames(data)) {
        fc_stats <- calculate_fc_stats(data$log2FoldChange)
        cat("\nFOLD CHANGE STATISTICS:\n")
        cat("Mean:", round(fc_stats$mean, 3), "\n")
        cat("Median:", round(fc_stats$median, 3), "\n")
        cat("Range: [", round(min(data$log2FoldChange), 3), ", ", 
            round(max(data$log2FoldChange), 3), "]\n", sep = "")
        cat("Up-regulated:", fc_stats$up_regulated, "\n")
        cat("Down-regulated:", fc_stats$down_regulated, "\n")
      }
      
      if ("padj" %in% colnames(data)) {
        sig_genes <- sum(data$padj < 0.05, na.rm = TRUE)
        cat("Significant (padj < 0.05):", sig_genes, "\n")
      }
    })
    
    output$expression_plot <- renderPlot({
      create_expression_plot()
    })
    
    output$distribution_plot_exp <- renderPlot({
      create_distribution_plot_exp()
    })
    
    output$comparison_plot <- renderPlot({
      create_comparison_plot()
    })
    
    output$gene_details_table <- DT::renderDataTable({
      if (!module_rv$genes_loaded || is.null(module_rv$expression_data)) {
        return(NULL)
      }
      
      data <- module_rv$expression_data
      display_data <- create_summary_df(data)
      
      DT::datatable(
        display_data,
        options = list(
          scrollX = TRUE,
          pageLength = 10,
          dom = 'tip'
        ),
        rownames = FALSE,
        caption = "Detailed Gene Information"
      )
    })
    
    # Download handlers
    output$download_plot_exp <- downloadHandler(
      filename = function() {
        paste0(
          "expression_plot_",
          input$plot_type, "_",
          format(Sys.time(), "%Y%m%d_%H%M%S"),
          ".png"
        )
      },
      content = function(file) {
        plot_obj <- module_rv$plot_object
        if (is.null(plot_obj)) {
          showNotification("No plot to download", type = "error")
          return()
        }
        
        success <- save_plot(
          plot_obj,
          file,
          width = input$plot_width_exp,
          height = input$plot_height_exp,
          dpi = 300
        )
        
        if (success) {
          showNotification("Plot downloaded successfully!", type = "message")
        } else {
          showNotification("Error downloading plot", type = "error")
        }
      }
    )
    
    output$download_data_exp <- downloadHandler(
      filename = function() {
        paste0(
          "expression_data_",
          format(Sys.time(), "%Y%m%d_%H%M%S"),
          ".csv"
        )
      },
      content = function(file) {
        data <- module_rv$expression_data
        if (is.null(data)) {
          showNotification("No data to download", type = "error")
          return()
        }
        
        write.csv(data, file, row.names = FALSE)
        showNotification("Expression data downloaded successfully!", type = "message")
      }
    )
    
    # Return reactive values for other modules
    return(reactive({
      list(
        selected_genes = module_rv$selected_genes,
        expression_data = module_rv$expression_data,
        genes_loaded = module_rv$genes_loaded
      )
    }))
  })
}