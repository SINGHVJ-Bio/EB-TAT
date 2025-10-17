# ENABL Biomarker/Target Analysis Tool - Analytics Module
# Provides analytical tools and statistics for differential expression data

# =============================================================================
# UI Function
# =============================================================================

analyticsUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    div(
      class = "analytics-module",
      style = "padding: 20px;",
      
      # Header
      div(
        class = "module-header",
        h3("Data Analytics", style = "color: #2196F3;"),
        p("Advanced analytical tools and statistics for differential expression data")
      ),
      
      fluidRow(
        # Controls column
        column(
          4,
          # Analysis configuration
          div(
            class = "card",
            style = "margin-bottom: 20px;",
            div(
              class = "card-header",
              h4("Analysis Settings", style = "margin: 0;")
            ),
            div(
              class = "card-body",
              # Threshold controls
              sliderInput(
                ns("fc_cut"),
                "Fold Change Threshold:",
                min = -5,  # Changed from 0 to -5
                max = 5,
                value = c(-0.5, 0.5),
                step = 0.1
              ),
              numericInput(
                ns("p_cut"),
                "P-value Threshold (-log10):",
                value = 1.301,
                min = 0,
                max = 10,
                step = 0.1
              ),
              numericInput(
                ns("top_hit"),
                "Top N Genes:",
                value = 20,
                min = 5,
                max = 100,
                step = 5
              ),
              selectInput(
                ns("plot_type"),
                "Plot Type:",
                choices = c(
                  "Regression" = "reg",
                  "Distribution" = "dist",
                  "Correlation" = "corr"
                ),
                selected = "reg"
              )
            )
          ),
          
          # Gene list input
          div(
            class = "card",
            style = "margin-bottom: 20px;",
            div(
              class = "card-header",
              h4("Gene List Analysis", style = "margin: 0;")
            ),
            div(
              class = "card-body",
              textAreaInput(
                ns("user_gene_list1"),
                "Enter gene symbols for analysis:",
                rows = 6,
                placeholder = "TP53\nEGFR\nBRCA1\n..."
              ),
              actionButton(
                ns("analyze_genes"),
                "Analyze Gene List",
                class = "btn-primary"
              ),
              hr(),
              selectInput(
                ns("venn_type"),
                "Venn Diagram Type:",
                choices = c(
                  "Combined" = "combined",
                  "Up/Down Regulated" = "up_down",
                  "Custom Lists" = "custom"
                ),
                selected = "combined"
              )
            )
          ),
          
          # Annotation filters
          div(
            class = "card",
            style = "margin-bottom: 20px;",
            div(
              class = "card-header",
              h4("Annotation Filters", style = "margin: 0;")
            ),
            div(
              class = "card-body",
              checkboxInput(
                ns("hse"),
                "Hepatocyte-Specific Expression",
                value = FALSE
              ),
              checkboxInput(
                ns("rec_lig"),
                "Receptor-Ligand Pairs",
                value = FALSE
              ),
              checkboxInput(
                ns("sec"),
                "Secreted Proteins",
                value = FALSE
              ),
              checkboxInput(
                ns("tf"),
                "Transcription Factors",
                value = FALSE
              ),
              checkboxInput(
                ns("immun"),
                "Immune-Related Genes",
                value = FALSE
              ),
              checkboxInput(
                ns("bld"),
                "Blood Plasma Proteins",
                value = FALSE
              )
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
                ns("download_analytics"),
                "Download Analysis",
                class = "btn-success"
              ),
              downloadButton(
                ns("download_plots"),
                "Download Plots",
                class = "btn-info"
              )
            )
          )
        ),
        
        # Results column
        column(
          8,
          # Summary statistics
          div(
            class = "card",
            style = "margin-bottom: 20px;",
            div(
              class = "card-header",
              h4("Summary Statistics", style = "margin: 0;")
            ),
            div(
              class = "card-body",
              fluidRow(
                column(
                  6,
                  wellPanel(
                    h5("Differential Expression"),
                    verbatimTextOutput(ns("de_summary"))
                  )
                ),
                column(
                  6,
                  wellPanel(
                    h5("Significance Overview"),
                    plotOutput(ns("sig_plot"), height = "200px")
                  )
                )
              )
            )
          ),
          
          # Main analysis plots
          div(
            class = "card",
            style = "margin-bottom: 20px;",
            div(
              class = "card-header",
              h4("Analysis Plots", style = "margin: 0;")
            ),
            div(
              class = "card-body",
              tabsetPanel(
                id = ns("plot_tabs"),
                type = "tabs",
                tabPanel(
                  "MA Plot",
                  withSpinner(
                    plotOutput(ns("ma_plot"), height = "400px"),
                    type = 4,
                    color = "#2196F3"
                  )
                ),
                tabPanel(
                  "Distribution",
                  withSpinner(
                    plotOutput(ns("distribution_plot"), height = "400px"),
                    type = 4,
                    color = "#2196F3"
                  )
                ),
                tabPanel(
                  "Venn Diagram",
                  withSpinner(
                    plotOutput(ns("venn_plot"), height = "400px"),
                    type = 4,
                    color = "#2196F3"
                  )
                )
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
                  "Top Genes",
                  DT::dataTableOutput(ns("top_genes_table"))
                ),
                tabPanel(
                  "Filtered Results",
                  DT::dataTableOutput(ns("filtered_table"))
                ),
                tabPanel(
                  "Gene List Results",
                  DT::dataTableOutput(ns("gene_list_table"))
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

analyticsServer <- function(id, rv, rv_analytics) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Reactive values for this module
    module_rv <- reactiveValues(
      filtered_data = NULL,
      analysis_results = NULL,
      gene_list_results = NULL,
      annotation_data = NULL
    )
    
    # Load annotation data
    load_annotation_data <- reactive({
      if (exists("annotation_data") && !is.null(annotation_data)) {
        return(annotation_data)
      }
      return(list())
    })
    
    # Filter data based on thresholds
    filter_data <- reactive({
      req(rv$data_loaded)
      req(rv$filtered_data)
      
      data <- rv$filtered_data
      
      # Apply fold change filter
      if ("log2FoldChange" %in% colnames(data)) {
        data <- data[
          data$log2FoldChange >= input$fc_cut[1] & 
          data$log2FoldChange <= input$fc_cut[2], 
        ]
      }
      
      # Apply p-value filter
      p_threshold <- 10^(-input$p_cut)
      if ("padj" %in% colnames(data)) {
        data <- data[data$padj <= p_threshold, ]
      } else if ("pvalue" %in% colnames(data)) {
        data <- data[data$pvalue <= p_threshold, ]
      }
      
      # Apply annotation filters if any are active
      annotation_filters <- c(
        input$hse, input$rec_lig, input$sec, input$tf, input$immun, input$bld
      )
      
      if (any(annotation_filters)) {
        data <- apply_annotation_filters(data, input)
      }
      
      module_rv$filtered_data <- data
      return(data)
    })
    
    # Apply annotation-based filters
    apply_annotation_filters <- function(data, input) {
      annotation_data <- load_annotation_data()
      if (length(annotation_data) == 0) return(data)
      
      filtered_genes <- character(0)
      
      # Hepatocyte-specific expression
      if (input$hse && !is.null(annotation_data$hepatocytes)) {
        hep_genes <- annotation_data$hepatocytes$symbol
        filtered_genes <- union(filtered_genes, hep_genes)
      }
      
      # Receptor-ligand pairs
      if (input$rec_lig && !is.null(annotation_data$receptor_ligand)) {
        rl_genes <- annotation_data$receptor_ligand$symbol
        filtered_genes <- union(filtered_genes, rl_genes)
      }
      
      # Secreted proteins
      if (input$sec && !is.null(annotation_data$secreted_proteins)) {
        sec_genes <- annotation_data$secreted_proteins$symbol
        filtered_genes <- union(filtered_genes, sec_genes)
      }
      
      # Transcription factors
      if (input$tf && !is.null(annotation_data$human_tf)) {
        tf_genes <- annotation_data$human_tf$symbol
        filtered_genes <- union(filtered_genes, tf_genes)
      }
      
      # Blood plasma proteins
      if (input$bld && exists("plasma_proteins") && !is.null(plasma_proteins)) {
        bld_genes <- plasma_proteins$symbol
        filtered_genes <- union(filtered_genes, bld_genes)
      }
      
      # Apply the gene filter
      if (length(filtered_genes) > 0) {
        data <- data[data$symbol %in% filtered_genes, ]
      }
      
      return(data)
    }
    
    # Analyze user-provided gene list
    analyze_gene_list <- function() {
      req(input$user_gene_list1)
      
      gene_list <- parse_gene_input(input$user_gene_list1)
      if (length(gene_list) == 0) return(NULL)
      
      data <- rv$filtered_data
      if (is.null(data)) return(NULL)
      
      # Find genes in the data
      found_genes <- data[data$symbol %in% gene_list, ]
      not_found <- setdiff(gene_list, data$symbol)
      
      results <- list(
        found_genes = found_genes,
        not_found = not_found,
        total_requested = length(gene_list),
        total_found = nrow(found_genes)
      )
      
      module_rv$gene_list_results <- results
      return(results)
    }
    
    # Create MA plot (Mean vs Fold Change)
    create_ma_plot <- function() {
      data <- filter_data()
      if (is.null(data) || nrow(data) == 0) {
        return(create_empty_plot("No data available for MA plot"))
      }
      
      if (!"baseMean" %in% colnames(data) || !"log2FoldChange" %in% colnames(data)) {
        return(create_empty_plot("Required columns (baseMean, log2FoldChange) not found"))
      }
      
      p <- ggplot(data, aes(x = log10(baseMean + 1), y = log2FoldChange)) +
        geom_point(alpha = 0.6, color = "#2196F3", size = 1.5) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        geom_smooth(method = "lm", color = "darkred", se = FALSE) +
        theme_ebtat() +
        labs(
          title = "MA Plot (Mean Expression vs Fold Change)",
          x = "Log10(Mean Expression + 1)",
          y = "Log2 Fold Change"
        )
      
      return(p)
    }
    
    # Create distribution plot
    create_distribution_plot <- function() {
      data <- filter_data()
      if (is.null(data) || nrow(data) == 0) {
        return(create_empty_plot("No data available for distribution plot"))
      }
      
      if (!"log2FoldChange" %in% colnames(data)) {
        return(create_empty_plot("log2FoldChange column not found"))
      }
      
      p <- ggplot(data, aes(x = log2FoldChange)) +
        geom_histogram(aes(y = ..density..), bins = 50, fill = "#4CAF50", alpha = 0.7) +
        geom_density(color = "darkgreen", size = 1) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
        theme_ebtat() +
        labs(
          title = "Distribution of Fold Changes",
          x = "Log2 Fold Change",
          y = "Density"
        )
      
      return(p)
    }
    
    # Create Venn diagram
    create_venn_diagram_reactive <- function() {
      data <- rv$filtered_data
      if (is.null(data) || nrow(data) == 0) {
        return(create_empty_plot("No data available for Venn diagram"))
      }
      
      venn_type <- input$venn_type
      
      if (venn_type == "combined") {
        # Up vs Down regulated
        if (!"log2FoldChange" %in% colnames(data)) {
          return(create_empty_plot("log2FoldChange column not found"))
        }
        
        up_genes <- data$symbol[data$log2FoldChange > 0]
        down_genes <- data$symbol[data$log2FoldChange < 0]
        
        gene_lists <- list(Up = up_genes, Down = down_genes)
        
      } else if (venn_type == "up_down") {
        # Significant up vs down
        if (!"log2FoldChange" %in% colnames(data) || 
            (!"padj" %in% colnames(data) && !"pvalue" %in% colnames(data))) {
          return(create_empty_plot("Required columns not found"))
        }
        
        p_threshold <- 10^(-input$p_cut)
        fc_threshold <- input$fc_cut[2]
        
        # Up-regulated significant
        up_sig <- data$symbol[
          data$log2FoldChange >= fc_threshold & 
          (data$padj <= p_threshold | data$pvalue <= p_threshold)
        ]
        
        # Down-regulated significant
        down_sig <- data$symbol[
          data$log2FoldChange <= -fc_threshold & 
          (data$padj <= p_threshold | data$pvalue <= p_threshold)
        ]
        
        gene_lists <- list(`Up Sig` = up_sig, `Down Sig` = down_sig)
        
      } else if (venn_type == "custom") {
        # Custom gene lists from user input
        if (is.null(input$user_gene_list1) || input$user_gene_list1 == "") {
          return(create_empty_plot("Please enter gene lists for custom Venn diagram"))
        }
        
        gene_lists <- strsplit(input$user_gene_list1, "\n")[[1]]
        # Simple implementation - could be enhanced for multiple custom lists
        list1_genes <- parse_gene_input(gene_lists[1])
        gene_lists <- list(List1 = list1_genes)
        
        if (length(gene_lists) > 1) {
          list2_genes <- parse_gene_input(gene_lists[2])
          gene_lists$List2 <- list2_genes
        }
      }
      
      # Remove empty lists
      gene_lists <- gene_lists[sapply(gene_lists, length) > 0]
      
      if (length(gene_lists) < 2) {
        return(create_empty_plot("Need at least 2 non-empty gene lists for Venn diagram"))
      }
      
      return(create_venn_diagram(gene_lists))
    }
    
    # Create significance overview plot
    create_sig_overview_plot <- function() {
      data <- rv$filtered_data
      if (is.null(data) || nrow(data) == 0) {
        return(create_empty_plot("No data available"))
      }
      
      # Calculate significance categories
      if ("padj" %in% colnames(data)) {
        sig_cutoff <- 0.05
        sig_data <- data[!is.na(data$padj), ]
        significant <- sum(sig_data$padj < sig_cutoff)
        not_sig <- nrow(sig_data) - significant
      } else if ("pvalue" %in% colnames(data)) {
        sig_cutoff <- 0.05
        sig_data <- data[!is.na(data$pvalue), ]
        significant <- sum(sig_data$pvalue < sig_cutoff)
        not_sig <- nrow(sig_data) - significant
      } else {
        return(create_empty_plot("No p-value data available"))
      }
      
      plot_data <- data.frame(
        Category = c("Significant", "Not Significant"),
        Count = c(significant, not_sig)
      )
      
      p <- ggplot(plot_data, aes(x = Category, y = Count, fill = Category)) +
        geom_col(alpha = 0.8) +
        scale_fill_manual(values = c("Significant" = "#4CAF50", "Not Significant" = "#F44336")) +
        theme_ebtat() +
        labs(title = "Significance Overview") +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none"
        )
      
      return(p)
    }
    
    # Update reactive values when inputs change
    observe({
      rv_analytics$fc_cut <- input$fc_cut
      rv_analytics$p_cut <- input$p_cut
      rv_analytics$top_hit <- input$top_hit
      rv_analytics$plot_type <- input$plot_type
      rv_analytics$venn_type <- input$venn_type
    })
    
    # Event observers
    observeEvent(input$analyze_genes, {
      showNotification("Analyzing gene list...", type = "message")
      results <- analyze_gene_list()
      
      if (!is.null(results)) {
        msg <- paste("Found", results$total_found, "of", results$total_requested, "genes")
        showNotification(msg, type = "message", duration = 3)
      } else {
        showNotification("Error analyzing gene list", type = "error")
      }
    })
    
    # Outputs
    output$de_summary <- renderPrint({
      data <- filter_data()
      if (is.null(data) || nrow(data) == 0) {
        cat("No data available\n")
        return()
      }
      
      cat("FILTERED DATA SUMMARY\n")
      cat("=====================\n")
      cat("Total genes:", nrow(data), "\n")
      
      if ("log2FoldChange" %in% colnames(data)) {
        fc_stats <- calculate_fc_stats(data$log2FoldChange)
        cat("\nFOLD CHANGE:\n")
        cat("Mean:", round(fc_stats$mean, 3), "\n")
        cat("Median:", round(fc_stats$median, 3), "\n")
        cat("Up-regulated:", fc_stats$up_regulated, "\n")
        cat("Down-regulated:", fc_stats$down_regulated, "\n")
      }
      
      if ("padj" %in% colnames(data)) {
        sig_genes <- sum(data$padj < 0.05, na.rm = TRUE)
        cat("Significant (padj < 0.05):", sig_genes, "\n")
      }
    })
    
    output$sig_plot <- renderPlot({
      create_sig_overview_plot()
    })
    
    output$ma_plot <- renderPlot({
      create_ma_plot()
    })
    
    output$distribution_plot <- renderPlot({
      create_distribution_plot()
    })
    
    output$venn_plot <- renderPlot({
      create_venn_diagram_reactive()
    })
    
    output$top_genes_table <- DT::renderDataTable({
      data <- filter_data()
      if (is.null(data) || nrow(data) == 0) return(NULL)
      
      # Get top genes by significance
      top_genes <- get_top_genes(data, n = input$top_hit, by = "padj")
      display_data <- create_summary_df(top_genes)
      
      DT::datatable(
        display_data,
        options = list(
          scrollX = TRUE,
          pageLength = 10,
          dom = 'tip'
        ),
        rownames = FALSE
      )
    })
    
    output$filtered_table <- DT::renderDataTable({
      data <- filter_data()
      if (is.null(data) || nrow(data) == 0) return(NULL)
      
      display_data <- create_summary_df(data, max_rows = 50)
      
      DT::datatable(
        display_data,
        options = list(
          scrollX = TRUE,
          pageLength = 10,
          dom = 'tip'
        ),
        rownames = FALSE
      )
    })
    
    output$gene_list_table <- DT::renderDataTable({
      results <- module_rv$gene_list_results
      if (is.null(results)) return(NULL)
      
      display_data <- create_summary_df(results$found_genes)
      
      DT::datatable(
        display_data,
        options = list(
          scrollX = TRUE,
          pageLength = 10,
          dom = 'tip'
        ),
        rownames = FALSE,
        caption = paste("Found", results$total_found, "of", results$total_requested, "genes")
      )
    })
    
    # Download handlers
    output$download_analytics <- downloadHandler(
      filename = function() {
        paste0("analytics_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
      },
      content = function(file) {
        data <- filter_data()
        if (is.null(data)) {
          showNotification("No data to download", type = "error")
          return()
        }
        
        write.csv(data, file, row.names = FALSE)
        showNotification("Analysis results downloaded successfully!", type = "message")
      }
    )
    
    output$download_plots <- downloadHandler(
      filename = function() {
        paste0("analytics_plots_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip")
      },
      content = function(file) {
        showNotification("Preparing plots for download...", type = "message")
        
        # Create temporary directory for plots
        temp_dir <- tempdir()
        plot_files <- c()
        
        # Generate and save plots
        tryCatch({
          # MA Plot
          ma_plot <- create_ma_plot()
          if (!is.null(ma_plot)) {
            ma_file <- file.path(temp_dir, "ma_plot.png")
            ggsave(ma_file, ma_plot, width = 10, height = 8, dpi = 300)
            plot_files <- c(plot_files, ma_file)
          }
          
          # Distribution Plot
          dist_plot <- create_distribution_plot()
          if (!is.null(dist_plot)) {
            dist_file <- file.path(temp_dir, "distribution_plot.png")
            ggsave(dist_file, dist_plot, width = 10, height = 8, dpi = 300)
            plot_files <- c(plot_files, dist_file)
          }
          
          # Venn Diagram
          venn_plot <- create_venn_diagram_reactive()
          if (!is.null(venn_plot)) {
            venn_file <- file.path(temp_dir, "venn_diagram.png")
            ggsave(venn_file, venn_plot, width = 10, height = 8, dpi = 300)
            plot_files <- c(plot_files, venn_file)
          }
          
          # Create zip file
          if (length(plot_files) > 0) {
            zip(file, plot_files, flags = "-j")
            showNotification("Plots downloaded successfully!", type = "message")
          } else {
            showNotification("No plots available to download", type = "warning")
          }
          
        }, error = function(e) {
          showNotification("Error downloading plots", type = "error")
        })
      }
    )
    
    # Return reactive values for other modules
    return(reactive({
      list(
        filtered_data = module_rv$filtered_data,
        analysis_results = module_rv$analysis_results,
        gene_list_results = module_rv$gene_list_results
      )
    }))
  })
}