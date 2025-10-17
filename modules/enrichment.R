# ENABL Biomarker/Target Analysis Tool - Enrichment Analysis Module
# Handles gene set enrichment analysis and pathway analysis

# =============================================================================
# UI Function
# =============================================================================

enrichmentUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    div(
      class = "enrichment-module",
      style = "padding: 20px;",
      
      # Header
      div(
        class = "module-header",
        h3("Enrichment Analysis", style = "color: #2196F3;"),
        p("Gene set enrichment analysis and pathway analysis using clusterProfiler")
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
              h4("Enrichment Settings", style = "margin: 0;")
            ),
            div(
              class = "card-body",
              # Gene selection
              selectInput(
                ns("gene_selection"),
                "Genes for enrichment:",
                choices = c(
                  "All filtered genes" = "all",
                  "Significant genes only" = "sig",
                  "Up-regulated genes" = "up",
                  "Down-regulated genes" = "down",
                  "Custom gene list" = "custom"
                ),
                selected = "sig"
              ),
              conditionalPanel(
                condition = paste0("input['", ns("gene_selection"), "'] == 'custom'"),
                textAreaInput(
                  ns("custom_genes_enrich"),
                  "Enter gene symbols:",
                  rows = 6,
                  placeholder = "TP53\nEGFR\nBRCA1\n..."
                )
              ),
              
              # Threshold controls
              sliderInput(
                ns("fc_cut_e"),
                "Fold Change Threshold:",
                min = -5,  # Changed from 0 to -5
                max = 5,
                value = c(-1.0, 1.0),
                step = 0.1
              ),
              numericInput(
                ns("p_cut_e"),
                "P-value Threshold (-log10):",
                value = 1.31,
                min = 0,
                max = 10,
                step = 0.1
              ),
              
              # Enrichment parameters
              selectInput(
                ns("enrichment_type"),
                "Enrichment Type:",
                choices = c(
                  "GO Biological Process" = "GOBP",
                  "GO Molecular Function" = "GOMF",
                  "GO Cellular Component" = "GOCC",
                  "KEGG Pathways" = "KEGG",
                  "Reactome Pathways" = "Reactome"
                ),
                selected = "GOBP"
              ),
              numericInput(
                ns("min_geneset"),
                "Min. GeneSet Size:",
                value = 10,
                min = 5,
                max = 100,
                step = 5
              ),
              numericInput(
                ns("max_geneset"),
                "Max. GeneSet Size:",
                value = 500,
                min = 100,
                max = 1000,
                step = 50
              ),
              numericInput(
                ns("top_pathways"),
                "Top Pathways to Show:",
                value = 20,
                min = 5,
                max = 50,
                step = 5
              ),
              checkboxInput(
                ns("olap"),
                "Show overlapping genes",
                value = FALSE
              )
            )
          ),
          
          # Visualization options
          div(
            class = "card",
            style = "margin-bottom: 20px;",
            div(
              class = "card-header",
              h4("Visualization", style = "margin: 0;")
            ),
            div(
              class = "card-body",
              selectInput(
                ns("e_plot_type"),
                "Plot Type:",
                choices = c(
                  "Dot Plot" = "dot",
                  "Bar Plot" = "bar",
                  "Enrichment Map" = "emap",
                  "Category Net" = "cnet",
                  "Gene Concept Net" = "gcn"
                ),
                selected = "dot"
              ),
              numericInput(
                ns("plot_width_e"),
                "Plot Width:",
                value = 12,
                min = 6,
                max = 20,
                step = 1
              ),
              numericInput(
                ns("plot_height_e"),
                "Plot Height:",
                value = 8,
                min = 6,
                max = 20,
                step = 1
              )
            )
          ),
          
          # Action buttons
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
                ns("run_enrichment"),
                "Run Enrichment Analysis",
                class = "btn-success",
                icon = icon("play")
              ),
              actionButton(
                ns("reset_enrichment"),
                "Reset",
                class = "btn-warning",
                icon = icon("refresh")
              ),
              hr(),
              downloadButton(
                ns("download_enrichment"),
                "Download Results",
                class = "btn-primary"
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
                condition = paste0("!output['", ns("data_available_e"), "']"),
                div(
                  class = "alert alert-warning",
                  strong("No data available."),
                  "Please load data in the Data Input tab."
                )
              ),
              conditionalPanel(
                condition = paste0("output['", ns("enrichment_run"), "']"),
                div(
                  class = "alert alert-success",
                  strong("Enrichment analysis completed!"),
                  textOutput(ns("enrichment_summary"))
                )
              ),
              verbatimTextOutput(ns("enrichment_info"))
            )
          ),
          
          # Main enrichment plots
          div(
            class = "card",
            style = "margin-bottom: 20px;",
            div(
              class = "card-header",
              h4("Enrichment Results", style = "margin: 0;")
            ),
            div(
              class = "card-body",
              withSpinner(
                plotOutput(ns("enrichment_plot"), height = "500px"),
                type = 4,
                color = "#2196F3"
              )
            )
          ),
          
          # Results tables
          div(
            class = "card",
            style = "margin-bottom: 20px;",
            div(
              class = "card-header",
              h4("Pathway Details", style = "margin: 0;")
            ),
            div(
              class = "card-body",
              tabsetPanel(
                type = "tabs",
                tabPanel(
                  "Enriched Pathways",
                  DT::dataTableOutput(ns("pathway_table"))
                ),
                tabPanel(
                  "Gene-Pathway Mapping",
                  DT::dataTableOutput(ns("gene_pathway_table"))
                )
              )
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
                  "Volcano with Pathways",
                  withSpinner(
                    plotOutput(ns("volcano_pathways"), height = "400px"),
                    type = 4,
                    color = "#2196F3"
                  )
                ),
                tabPanel(
                  "Pathway Network",
                  withSpinner(
                    plotOutput(ns("pathway_network"), height = "400px"),
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

enrichmentServer <- function(id, rv, rv_enrichment) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Reactive values for this module
    module_rv <- reactiveValues(
      enrichment_results = NULL,
      gene_universe = NULL,
      selected_genes = NULL,
      last_analysis = NULL,
      enrichment_run = FALSE
    )
    
    # Get genes for enrichment analysis
    get_genes_for_enrichment <- reactive({
      req(rv$data_loaded)
      req(rv$filtered_data)
      
      data <- rv$filtered_data
      selection_type <- input$gene_selection
      
      if (selection_type == "custom") {
        # Parse custom gene list
        genes <- parse_gene_input(input$custom_genes_enrich)
        return(genes)
      }
      
      # Apply thresholds
      p_threshold <- 10^(-input$p_cut_e)
      fc_low <- input$fc_cut_e[1]
      fc_high <- input$fc_cut_e[2]
      
      filtered_data <- data
      
      # Apply p-value filter
      if ("padj" %in% colnames(data)) {
        filtered_data <- filtered_data[filtered_data$padj <= p_threshold, ]
      } else if ("pvalue" %in% colnames(data)) {
        filtered_data <- filtered_data[filtered_data$pvalue <= p_threshold, ]
      }
      
      # Apply fold change filter
      if ("log2FoldChange" %in% colnames(filtered_data)) {
        if (selection_type == "up") {
          filtered_data <- filtered_data[filtered_data$log2FoldChange > 0, ]
        } else if (selection_type == "down") {
          filtered_data <- filtered_data[filtered_data$log2FoldChange < 0, ]
        } else if (selection_type == "sig") {
          # Already filtered by p-value, now filter by fold change magnitude
          filtered_data <- filtered_data[
            abs(filtered_data$log2FoldChange) >= max(abs(fc_low), abs(fc_high)), 
          ]
        }
        # For "all", we keep all filtered genes
      }
      
      return(unique(filtered_data$symbol))
    })
    
    # Get gene universe (background genes)
    get_gene_universe <- reactive({
      req(rv$data_loaded)
      req(rv$filtered_data)
      
      # Use all genes in the filtered data as background
      return(unique(rv$filtered_data$symbol))
    })
    
    # Convert gene symbols to ENTREZ IDs
    convert_to_entrez <- function(gene_symbols) {
      if (length(gene_symbols) == 0) return(character(0))
      
      tryCatch({
        # Map gene symbols to ENTREZ IDs
        map <- clusterProfiler::bitr(
          gene_symbols,
          fromType = "SYMBOL",
          toType = "ENTREZID",
          OrgDb = org.Hs.eg.db
        )
        return(map$ENTREZID)
      }, error = function(e) {
        warning("Error converting gene symbols to ENTREZ IDs: ", e$message)
        return(character(0))
      })
    }
    
    # Run enrichment analysis
    run_enrichment_analysis <- function() {
      showNotification("Running enrichment analysis...", type = "message")
      
      genes <- get_genes_for_enrichment()
      universe <- get_gene_universe()
      
      if (length(genes) < 5) {
        showNotification("Too few genes for enrichment analysis (need at least 5)", 
                        type = "error")
        return(NULL)
      }
      
      # Convert to ENTREZ IDs
      gene_entrez <- convert_to_entrez(genes)
      universe_entrez <- convert_to_entrez(universe)
      
      if (length(gene_entrez) == 0) {
        showNotification("Error converting genes to ENTREZ IDs", type = "error")
        return(NULL)
      }
      
      enrichment_type <- input$enrichment_type
      
      tryCatch({
        if (enrichment_type %in% c("GOBP", "GOMF", "GOCC")) {
          # GO enrichment
          result <- clusterProfiler::enrichGO(
            gene = gene_entrez,
            universe = universe_entrez,
            OrgDb = org.Hs.eg.db,
            ont = switch(enrichment_type,
                        "GOBP" = "BP",
                        "GOMF" = "MF", 
                        "GOCC" = "CC"),
            pAdjustMethod = "BH",
            pvalueCutoff = 0.05,
            qvalueCutoff = 0.2,
            minGSSize = input$min_geneset,
            maxGSSize = input$max_geneset,
            readable = TRUE
          )
        } else if (enrichment_type == "KEGG") {
          # KEGG enrichment
          result <- clusterProfiler::enrichKEGG(
            gene = gene_entrez,
            universe = universe_entrez,
            organism = 'hsa',
            pAdjustMethod = "BH",
            pvalueCutoff = 0.05,
            qvalueCutoff = 0.2,
            minGSSize = input$min_geneset,
            maxGSSize = input$max_geneset
          )
        } else if (enrichment_type == "Reactome") {
          # Reactome enrichment
          result <- DOSE::enrichPathway(
            gene = gene_entrez,
            universe = universe_entrez,
            organism = "human",
            pAdjustMethod = "BH",
            pvalueCutoff = 0.05,
            qvalueCutoff = 0.2,
            minGSSize = input$min_geneset,
            maxGSSize = input$max_geneset,
            readable = TRUE
          )
        }
        
        if (is.null(result) || nrow(result) == 0) {
          showNotification("No significant enrichment found", type = "warning")
          return(NULL)
        }
        
        module_rv$enrichment_results <- result
        module_rv$gene_universe <- universe
        module_rv$selected_genes <- genes
        module_rv$last_analysis <- Sys.time()
        module_rv$enrichment_run <- TRUE
        
        showNotification(
          paste("Enrichment analysis completed:", nrow(result), "pathways found"),
          type = "message",
          duration = 5
        )
        
        return(result)
        
      }, error = function(e) {
        showNotification(paste("Enrichment analysis failed:", e$message), 
                        type = "error")
        return(NULL)
      })
    }
    
    # Create enrichment plot
    create_enrichment_plot <- function() {
      results <- module_rv$enrichment_results
      if (is.null(results) || nrow(results) == 0) {
        return(create_empty_plot("No enrichment results available"))
      }
      
      plot_type <- input$e_plot_type
      top_n <- input$top_pathways
      
      tryCatch({
        # Get top pathways
        if (nrow(results) > top_n) {
          results <- results[1:top_n, ]
        }
        
        if (plot_type == "dot") {
          p <- enrichplot::dotplot(results, showCategory = top_n) +
            theme_ebtat() +
            labs(title = "Enrichment Dot Plot")
        } else if (plot_type == "bar") {
          p <- create_enrichment_barplot(results, top_n)
        } else if (plot_type == "emap") {
          # Enrichment map
          if (nrow(results) >= 3) {
            p <- enrichplot::emapplot(results, showCategory = top_n)
          } else {
            p <- create_empty_plot("Need at least 3 pathways for enrichment map")
          }
        } else if (plot_type == "cnet") {
          # Category net plot
          if (nrow(results) >= 2) {
            p <- enrichplot::cnetplot(results, showCategory = min(5, top_n))
          } else {
            p <- create_empty_plot("Need at least 2 pathways for category net")
          }
        } else if (plot_type == "gcn") {
          # Gene concept net
          if (nrow(results) >= 2) {
            p <- enrichplot::cnetplot(results, showCategory = min(5, top_n), 
                                     circular = TRUE, colorEdge = TRUE)
          } else {
            p <- create_empty_plot("Need at least 2 pathways for gene concept net")
          }
        }
        
        return(p)
        
      }, error = function(e) {
        warning("Error creating enrichment plot: ", e$message)
        return(create_empty_plot("Error creating enrichment plot"))
      })
    }
    
    # Create volcano plot with pathway highlights
    create_volcano_pathways_plot <- function() {
      data <- rv$filtered_data
      results <- module_rv$enrichment_results
      
      if (is.null(data) || is.null(results) || nrow(results) == 0) {
        return(create_empty_plot("No data or enrichment results available"))
      }
      
      # Get genes from top pathways
      top_pathways <- results[1:min(5, nrow(results)), ]
      pathway_genes <- unique(unlist(strsplit(top_pathways$geneID, "/")))
      
      # Create volcano plot with pathway genes highlighted
      p <- create_volcano_plot(
        de_data = data,
        fc_threshold = input$fc_cut_e[2],
        p_threshold = 10^(-input$p_cut_e),
        highlight_genes = pathway_genes,
        title = "Volcano Plot - Pathway Genes Highlighted",
        point_size = 2,
        alpha = 0.6
      )
      
      return(p)
    }
    
    # Create pathway network plot
    create_pathway_network_plot <- function() {
      results <- module_rv$enrichment_results
      if (is.null(results) || nrow(results) < 2) {
        return(create_empty_plot("Need at least 2 pathways for network"))
      }
      
      tryCatch({
        # Create a simple network visualization
        pathway_data <- as.data.frame(results)
        if (nrow(pathway_data) > 10) {
          pathway_data <- pathway_data[1:10, ]
        }
        
        # Create a bar plot showing pathway relationships
        p <- ggplot(pathway_data, aes(x = reorder(Description, -pvalue), y = -log10(pvalue))) +
          geom_col(aes(fill = -log10(pvalue)), width = 0.7) +
          scale_fill_gradient(low = "lightblue", high = "darkblue", name = "-Log10 P-value") +
          coord_flip() +
          theme_ebtat() +
          labs(
            title = "Top Enriched Pathways",
            x = "Pathway",
            y = "-Log10 P-value"
          ) +
          theme(
            axis.text.y = element_text(size = 9)
          )
        
        return(p)
        
      }, error = function(e) {
        warning("Error creating pathway network: ", e$message)
        return(create_empty_plot("Error creating pathway network"))
      })
    }
    
    # Update reactive values when inputs change
    observe({
      rv_enrichment$fc_cut_e <- input$fc_cut_e
      rv_enrichment$p_cut_e <- input$p_cut_e
      rv_enrichment$e_plot_type <- input$e_plot_type
      rv_enrichment$olap <- input$olap
    })
    
    # Event observers
    observeEvent(input$run_enrichment, {
      run_enrichment_analysis()
    })
    
    observeEvent(input$reset_enrichment, {
      module_rv$enrichment_results <- NULL
      module_rv$gene_universe <- NULL
      module_rv$selected_genes <- NULL
      module_rv$last_analysis <- NULL
      module_rv$enrichment_run <- FALSE
      
      showNotification("Enrichment analysis reset", type = "warning")
    })
    
    # Outputs
    output$data_available_e <- reactive({
      rv$data_loaded && !is.null(rv$filtered_data)
    })
    outputOptions(output, "data_available_e", suspendWhenHidden = FALSE)
    
    output$enrichment_run <- reactive({
      module_rv$enrichment_run
    })
    outputOptions(output, "enrichment_run", suspendWhenHidden = FALSE)
    
    output$enrichment_summary <- renderText({
      if (!module_rv$enrichment_run) return("")
      
      results <- module_rv$enrichment_results
      genes <- module_rv$selected_genes
      
      paste(
        "Input genes:", length(genes),
        "| Enriched pathways:", ifelse(!is.null(results), nrow(results), 0)
      )
    })
    
    output$enrichment_info <- renderPrint({
      if (!module_rv$enrichment_run) {
        cat("Enrichment analysis not run yet.\n")
        cat("Click 'Run Enrichment Analysis' to start.\n")
        return()
      }
      
      results <- module_rv$enrichment_results
      genes <- module_rv$selected_genes
      
      cat("ENRICHMENT ANALYSIS SUMMARY\n")
      cat("===========================\n")
      cat("Analysis time:", format(module_rv$last_analysis), "\n")
      cat("Input genes:", length(genes), "\n")
      cat("Gene universe:", length(module_rv$gene_universe), "\n")
      cat("Enrichment type:", input$enrichment_type, "\n")
      cat("Significant pathways:", ifelse(!is.null(results), nrow(results), 0), "\n\n")
      
      if (!is.null(results) && nrow(results) > 0) {
        cat("TOP PATHWAYS:\n")
        top_pathways <- results[1:min(5, nrow(results)), ]
        for (i in 1:nrow(top_pathways)) {
          cat(i, ". ", top_pathways$Description[i], 
              " (p=", format(top_pathways$pvalue[i], scientific = TRUE, digits = 3), ")\n", sep = "")
        }
      }
    })
    
    output$enrichment_plot <- renderPlot({
      create_enrichment_plot()
    })
    
    output$volcano_pathways <- renderPlot({
      create_volcano_pathways_plot()
    })
    
    output$pathway_network <- renderPlot({
      create_pathway_network_plot()
    })
    
    output$pathway_table <- DT::renderDataTable({
      results <- module_rv$enrichment_results
      if (is.null(results) || nrow(results) == 0) return(NULL)
      
      # Convert to data frame and format
      pathway_data <- as.data.frame(results)
      
      # Select and format columns
      display_cols <- c("Description", "pvalue", "p.adjust", "qvalue", "Count", "geneID")
      display_cols <- display_cols[display_cols %in% colnames(pathway_data)]
      
      display_data <- pathway_data[, display_cols, drop = FALSE]
      
      # Format numeric columns
      if ("pvalue" %in% colnames(display_data)) {
        display_data$pvalue <- format_pvalue(display_data$pvalue)
      }
      if ("p.adjust" %in% colnames(display_data)) {
        display_data$p.adjust <- format_pvalue(display_data$p.adjust)
      }
      if ("qvalue" %in% colnames(display_data)) {
        display_data$qvalue <- format_pvalue(display_data$qvalue)
      }
      
      DT::datatable(
        display_data,
        options = list(
          scrollX = TRUE,
          pageLength = 10,
          dom = 'tip',
          order = list(2, 'asc')  # Sort by p-value
        ),
        rownames = FALSE,
        caption = "Enriched Pathways (sorted by p-value)"
      )
    })
    
    output$gene_pathway_table <- DT::renderDataTable({
      results <- module_rv$enrichment_results
      if (is.null(results) || nrow(results) == 0) return(NULL)
      
      # Create gene-pathway mapping
      pathway_data <- as.data.frame(results)
      gene_pathway_list <- list()
      
      for (i in 1:nrow(pathway_data)) {
        pathway <- pathway_data$Description[i]
        genes <- unlist(strsplit(pathway_data$geneID[i], "/"))
        
        for (gene in genes) {
          gene_pathway_list[[length(gene_pathway_list) + 1]] <- data.frame(
            Gene = gene,
            Pathway = pathway,
            PValue = pathway_data$pvalue[i],
            stringsAsFactors = FALSE
          )
        }
      }
      
      if (length(gene_pathway_list) == 0) return(NULL)
      
      gene_pathway_df <- do.call(rbind, gene_pathway_list)
      
      # Format p-values
      gene_pathway_df$PValue <- format_pvalue(gene_pathway_df$PValue)
      
      DT::datatable(
        gene_pathway_df,
        options = list(
          scrollX = TRUE,
          pageLength = 10,
          dom = 'tip'
        ),
        rownames = FALSE,
        caption = "Gene-Pathway Associations"
      )
    })
    
    # Download handler
    output$download_enrichment <- downloadHandler(
      filename = function() {
        paste0(
          "enrichment_results_",
          input$enrichment_type, "_",
          format(Sys.time(), "%Y%m%d_%H%M%S"),
          ".csv"
        )
      },
      content = function(file) {
        results <- module_rv$enrichment_results
        if (is.null(results)) {
          showNotification("No enrichment results to download", type = "error")
          return()
        }
        
        # Convert to data frame
        results_df <- as.data.frame(results)
        write.csv(results_df, file, row.names = FALSE)
        showNotification("Enrichment results downloaded successfully!", type = "message")
      }
    )
    
    # Return reactive values for other modules
    return(reactive({
      list(
        enrichment_results = module_rv$enrichment_results,
        selected_genes = module_rv$selected_genes,
        enrichment_run = module_rv$enrichment_run
      )
    }))
  })
}