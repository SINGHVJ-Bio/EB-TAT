# ENABL Biomarker/Target Analysis Tool - Plotting Utilities
# Optimized for current data structure without gene_feature.xlsx

# =============================================================================
# Volcano Plot Functions
# =============================================================================

#' Create a volcano plot from differential expression data
create_volcano_plot <- function(de_data, 
                               fc_threshold = 1, 
                               p_threshold = 0.05,
                               highlight_genes = NULL,
                               title = "Volcano Plot",
                               point_size = 2,
                               alpha = 0.6) {
  
  if (is.null(de_data) || nrow(de_data) == 0) {
    return(create_empty_plot("No data available for volcano plot"))
  }
  
  # Check required columns
  if (!"log2FoldChange" %in% colnames(de_data)) {
    return(create_empty_plot("Missing log2FoldChange column"))
  }
  
  if (!"padj" %in% colnames(de_data) && !"pvalue" %in% colnames(de_data)) {
    return(create_empty_plot("Missing p-value columns"))
  }
  
  # Prepare data for plotting
  plot_data <- prepare_volcano_data(de_data, fc_threshold, p_threshold)
  
  # Create base plot
  p <- ggplot(plot_data, aes(x = log2FoldChange, y = minus_log10_padj)) +
    geom_point(aes(color = significance, alpha = alpha), size = point_size) +
    scale_color_manual(
      values = c(
        "Not Significant" = "gray60",
        "Up-regulated" = "#E31A1C",
        "Down-regulated" = "#1F78B4"
      )
    ) +
    scale_alpha_identity() +
    theme_ebtat() +
    labs(
      title = title,
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-value",
      color = "Significance"
    ) +
    guides(alpha = "none")
  
  # Add significance thresholds
  p <- p + 
    geom_vline(xintercept = c(-fc_threshold, fc_threshold), 
               linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = -log10(p_threshold), 
               linetype = "dashed", alpha = 0.5)
  
  # Highlight specific genes if provided
  if (!is.null(highlight_genes) && length(highlight_genes) > 0) {
    highlight_data <- plot_data[plot_data$symbol %in% highlight_genes, ]
    if (nrow(highlight_data) > 0) {
      p <- p +
        geom_point(data = highlight_data, 
                  color = "black", size = point_size * 1.5, shape = 1) +
        ggrepel::geom_text_repel(
          data = highlight_data,
          aes(label = symbol),
          size = 3,
          box.padding = 0.5,
          max.overlaps = 20
        )
    }
  }
  
  return(p)
}

#' Prepare data for volcano plot
prepare_volcano_data <- function(de_data, fc_threshold, p_threshold) {
  plot_data <- de_data
  
  # Ensure we have minus_log10_padj column
  if (!"minus_log10_padj" %in% colnames(plot_data)) {
    if ("padj" %in% colnames(plot_data)) {
      plot_data$minus_log10_padj <- safe_minus_log10(plot_data$padj)
    } else if ("pvalue" %in% colnames(plot_data)) {
      plot_data$minus_log10_padj <- safe_minus_log10(plot_data$pvalue)
    }
  }
  
  # Add significance classification
  plot_data$significance <- "Not Significant"
  
  if ("log2FoldChange" %in% colnames(plot_data) && "minus_log10_padj" %in% colnames(plot_data)) {
    # Up-regulated
    up_mask <- plot_data$log2FoldChange >= fc_threshold & 
      plot_data$minus_log10_padj >= -log10(p_threshold)
    plot_data$significance[up_mask] <- "Up-regulated"
    
    # Down-regulated
    down_mask <- plot_data$log2FoldChange <= -fc_threshold & 
      plot_data$minus_log10_padj >= -log10(p_threshold)
    plot_data$significance[down_mask] <- "Down-regulated"
  }
  
  # Ensure symbol column exists for labeling
  if (!"symbol" %in% colnames(plot_data)) {
    plot_data$symbol <- paste0("Gene_", 1:nrow(plot_data))
  }
  
  return(plot_data)
}

# =============================================================================
# Expression Plot Functions
# =============================================================================

#' Create expression bar plot for selected genes
create_expression_barplot <- function(de_data, 
                                     selected_genes,
                                     value_column = "log2FoldChange",
                                     title = "Gene Expression") {
  
  if (is.null(de_data) || nrow(de_data) == 0) {
    return(create_empty_plot("No data available for expression plot"))
  }
  
  if (is.null(selected_genes) || length(selected_genes) == 0) {
    return(create_empty_plot("No genes selected for expression plot"))
  }
  
  # Filter data for selected genes
  plot_data <- de_data[de_data$symbol %in% selected_genes, ]
  
  if (nrow(plot_data) == 0) {
    return(create_empty_plot("Selected genes not found in data"))
  }
  
  if (!value_column %in% colnames(plot_data)) {
    return(create_empty_plot(paste("Column not found:", value_column)))
  }
  
  # Order by expression value
  plot_data <- plot_data[order(plot_data[[value_column]]), ]
  plot_data$symbol <- factor(plot_data$symbol, levels = plot_data$symbol)
  
  # Create bar plot
  p <- ggplot(plot_data, aes_string(x = "symbol", y = value_column, fill = value_column)) +
    geom_col() +
    scale_fill_gradient2(
      low = "#1F78B4",
      mid = "white",
      high = "#E31A1C",
      midpoint = 0,
      name = "Log2 FC"
    ) +
    coord_flip() +
    theme_ebtat() +
    labs(
      title = title,
      x = "Gene Symbol",
      y = "Log2 Fold Change"
    ) +
    theme(
      axis.text.y = element_text(size = 10),
      legend.position = "right"
    )
  
  return(p)
}

#' Create expression dot plot with p-values
create_expression_dotplot <- function(de_data, 
                                     selected_genes,
                                     title = "Gene Expression with Significance") {
  
  if (is.null(de_data) || nrow(de_data) == 0) {
    return(create_empty_plot("No data available"))
  }
  
  if (is.null(selected_genes) || length(selected_genes) == 0) {
    return(create_empty_plot("No genes selected"))
  }
  
  # Filter and prepare data
  plot_data <- de_data[de_data$symbol %in% selected_genes, ]
  
  if (nrow(plot_data) == 0) {
    return(create_empty_plot("Selected genes not found"))
  }
  
  # Ensure we have required columns
  required_cols <- c("symbol", "log2FoldChange")
  if (!all(required_cols %in% colnames(plot_data))) {
    return(create_empty_plot("Missing required columns for dot plot"))
  }
  
  # Add significance stars if p-value available
  if ("padj" %in% colnames(plot_data)) {
    plot_data$significance <- cut(plot_data$padj,
                                 breaks = c(0, 0.001, 0.01, 0.05, 1),
                                 labels = c("***", "**", "*", ""),
                                 include.lowest = TRUE)
  } else if ("pvalue" %in% colnames(plot_data)) {
    plot_data$significance <- cut(plot_data$pvalue,
                                 breaks = c(0, 0.001, 0.01, 0.05, 1),
                                 labels = c("***", "**", "*", ""),
                                 include.lowest = TRUE)
  } else {
    plot_data$significance <- ""
  }
  
  # Order by fold change
  plot_data <- plot_data[order(plot_data$log2FoldChange), ]
  plot_data$symbol <- factor(plot_data$symbol, levels = plot_data$symbol)
  
  # Create dot plot
  p <- ggplot(plot_data, aes(x = symbol, y = log2FoldChange, color = log2FoldChange)) +
    geom_point(size = 4) +
    geom_text(aes(label = significance), 
              vjust = -1.5, 
              color = "black", 
              size = 5) +
    scale_color_gradient2(
      low = "#1F78B4",
      mid = "gray80",
      high = "#E31A1C",
      midpoint = 0,
      name = "Log2 FC"
    ) +
    coord_flip() +
    theme_ebtat() +
    labs(
      title = title,
      x = "Gene Symbol",
      y = "Log2 Fold Change"
    ) +
    theme(
      axis.text.y = element_text(size = 10),
      legend.position = "right"
    )
  
  return(p)
}

# =============================================================================
# Enrichment Plot Functions
# =============================================================================

#' Create enrichment bar plot
create_enrichment_barplot <- function(enrichment_results, 
                                     top_n = 20,
                                     title = "Top Enriched Pathways") {
  
  if (is.null(enrichment_results) || nrow(enrichment_results) == 0) {
    return(create_empty_plot("No enrichment results available"))
  }
  
  # Prepare data
  plot_data <- enrichment_results[1:min(top_n, nrow(enrichment_results)), ]
  plot_data <- plot_data[order(plot_data$pvalue), ]
  
  # Create bar plot
  p <- ggplot(plot_data, aes(x = reorder(Description, -pvalue), y = -log10(pvalue))) +
    geom_col(aes(fill = -log10(pvalue)), width = 0.7) +
    scale_fill_gradient(low = "lightblue", high = "darkblue", name = "-Log10 P-value") +
    coord_flip() +
    theme_ebtat() +
    labs(
      title = title,
      x = "Pathway",
      y = "-Log10 P-value"
    ) +
    theme(
      axis.text.y = element_text(size = 9),
      legend.position = "right"
    )
  
  return(p)
}

#' Create enrichment dot plot
create_enrichment_dotplot <- function(enrichment_results,
                                     top_n = 20,
                                     title = "Pathway Enrichment") {
  
  if (is.null(enrichment_results) || nrow(enrichment_results) == 0) {
    return(create_empty_plot("No enrichment results available"))
  }
  
  # Prepare data
  plot_data <- enrichment_results[1:min(top_n, nrow(enrichment_results)), ]
  
  # Create dot plot
  p <- ggplot(plot_data, aes(x = -log10(pvalue), y = reorder(Description, pvalue))) +
    geom_point(aes(size = Count, color = -log10(pvalue))) +
    scale_color_gradient(low = "lightblue", high = "darkblue", name = "-Log10 P-value") +
    scale_size_continuous(name = "Gene Count", range = c(2, 8)) +
    theme_ebtat() +
    labs(
      title = title,
      x = "-Log10 P-value",
      y = "Pathway"
    ) +
    theme(
      axis.text.y = element_text(size = 9)
    )
  
  return(p)
}

# =============================================================================
# Analytical Plot Functions
# =============================================================================

#' Create Venn diagram for gene lists
create_venn_diagram <- function(gene_lists, 
                               list_names = NULL,
                               title = "Gene List Overlap") {
  
  if (is.null(gene_lists) || length(gene_lists) == 0) {
    return(create_empty_plot("No gene lists provided for Venn diagram"))
  }
  
  # Set default names if not provided
  if (is.null(list_names)) {
    list_names <- paste("List", 1:length(gene_lists))
  }
  
  # Filter out empty lists
  non_empty_lists <- sapply(gene_lists, function(x) length(x) > 0)
  gene_lists <- gene_lists[non_empty_lists]
  list_names <- list_names[non_empty_lists]
  
  if (length(gene_lists) == 0) {
    return(create_empty_plot("All gene lists are empty"))
  }
  
  # Create Venn diagram
  tryCatch({
    if (length(gene_lists) <= 4) {
      # Use ggVennDiagram for up to 4 sets
      p <- ggVennDiagram::ggVennDiagram(gene_lists, label = "count", category.names = list_names) +
        scale_fill_gradient(low = "white", high = "red") +
        labs(title = title) +
        theme(plot.title = element_text(hjust = 0.5))
    } else {
      # For more than 4 sets, create an upset plot
      p <- create_upset_plot(gene_lists, list_names, title)
    }
    return(p)
  }, error = function(e) {
    warning("Error creating Venn diagram: ", e$message)
    return(create_empty_plot("Error creating Venn diagram"))
  })
}

#' Create UpSet plot for multiple gene lists
create_upset_plot <- function(gene_lists, list_names, title = "Gene List Intersections") {
  tryCatch({
    # Convert to list format for UpSetR
    all_genes <- unique(unlist(gene_lists))
    binary_matrix <- sapply(gene_lists, function(genes) {
      as.integer(all_genes %in% genes)
    })
    rownames(binary_matrix) <- all_genes
    colnames(binary_matrix) <- list_names
    
    # Create UpSet plot
    p <- UpSetR::upset(as.data.frame(binary_matrix), 
                      nsets = length(gene_lists),
                      main.bar.color = "steelblue",
                      sets.bar.color = "darkred",
                      mb.ratio = c(0.7, 0.3),
                      order.by = "freq")
    return(p)
  }, error = function(e) {
    warning("Error creating UpSet plot: ", e$message)
    return(create_empty_plot("Error creating UpSet plot"))
  })
}

# =============================================================================
# Heatmap Functions
# =============================================================================

#' Create expression heatmap
create_expression_heatmap <- function(expression_matrix,
                                     sample_info = NULL,
                                     gene_info = NULL,
                                     title = "Expression Heatmap",
                                     scale = "row",
                                     show_rownames = TRUE,
                                     show_colnames = TRUE) {
  
  if (is.null(expression_matrix) || nrow(expression_matrix) == 0) {
    return(create_empty_plot("No expression data available for heatmap"))
  }
  
  # Limit number of genes for performance
  if (nrow(expression_matrix) > 100) {
    expression_matrix <- expression_matrix[1:100, ]
    message("Limited heatmap to top 100 genes for performance")
  }
  
  # Create annotation data if provided
  annotation_col <- NULL
  annotation_row <- NULL
  
  if (!is.null(sample_info)) {
    annotation_col <- sample_info
  }
  
  if (!is.null(gene_info)) {
    annotation_row <- gene_info
  }
  
  # Create heatmap
  tryCatch({
    pheatmap::pheatmap(
      expression_matrix,
      scale = scale,
      main = title,
      show_rownames = show_rownames,
      show_colnames = show_colnames,
      annotation_col = annotation_col,
      annotation_row = annotation_row,
      color = colorRampPalette(c("blue", "white", "red"))(100),
      fontsize = 8,
      border_color = NA
    )
  }, error = function(e) {
    warning("Error creating heatmap: ", e$message)
    return(create_empty_plot("Error creating heatmap"))
  })
}

# =============================================================================
# Utility Plot Functions
# =============================================================================

#' Create an empty plot with message
create_empty_plot <- function(message = "No data available") {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = message, size = 6) +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

#' Create a distribution plot for numeric data
create_distribution_plot <- function(data, column, title = "Distribution") {
  if (is.null(data) || !column %in% colnames(data)) {
    return(create_empty_plot(paste("Column not found:", column)))
  }
  
  p <- ggplot(data, aes_string(x = column)) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", alpha = 0.7) +
    geom_density(color = "darkblue", size = 1) +
    theme_ebtat() +
    labs(
      title = title,
      x = column,
      y = "Density"
    )
  
  return(p)
}

#' Create a PCA plot
create_pca_plot <- function(expression_matrix, groups = NULL, title = "PCA Plot") {
  if (is.null(expression_matrix) || nrow(expression_matrix) < 2) {
    return(create_empty_plot("Insufficient data for PCA"))
  }
  
  tryCatch({
    # Perform PCA
    pca_result <- prcomp(t(expression_matrix), scale. = TRUE)
    pca_data <- as.data.frame(pca_result$x[, 1:2])
    
    # Add groups if provided
    if (!is.null(groups)) {
      pca_data$Group <- groups
    }
    
    # Create plot
    p <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
      geom_point(aes(color = if (!is.null(groups)) Group else NULL), size = 3) +
      theme_ebtat() +
      labs(
        title = title,
        x = paste0("PC1 (", round(summary(pca_result)$importance[2,1] * 100, 1), "%)"),
        y = paste0("PC2 (", round(summary(pca_result)$importance[2,2] * 100, 1), "%)")
      )
    
    if (!is.null(groups)) {
      p <- p + scale_color_discrete(name = "Group")
    }
    
    return(p)
    
  }, error = function(e) {
    warning("Error creating PCA plot: ", e$message)
    return(create_empty_plot("Error creating PCA plot"))
  })
}

# =============================================================================
# Plot Export Functions
# =============================================================================

#' Save plot to file with automatic format detection
save_plot <- function(plot_obj, filename, width = 10, height = 8, dpi = 300) {
  if (is.null(plot_obj)) {
    warning("No plot to save")
    return(FALSE)
  }
  
  file_ext <- tolower(tools::file_ext(filename))
  
  tryCatch({
    switch(file_ext,
           png = ggsave(filename, plot_obj, width = width, height = height, dpi = dpi),
           jpg =,
           jpeg = ggsave(filename, plot_obj, width = width, height = height, dpi = dpi),
           pdf = ggsave(filename, plot_obj, width = width, height = height),
           tiff = ggsave(filename, plot_obj, width = width, height = height, dpi = dpi),
           {
             # Default to PNG
             ggsave(paste0(tools::file_path_sans_ext(filename), ".png"), 
                    plot_obj, width = width, height = height, dpi = dpi)
           }
    )
    message("Plot saved as: ", filename)
    return(TRUE)
  }, error = function(e) {
    warning("Error saving plot: ", e$message)
    return(FALSE)
  })
}

#' Create downloadable plot handler
create_download_handler <- function(plot_func, filename, ...) {
  downloadHandler(
    filename = filename,
    content = function(file) {
      plot_obj <- plot_func(...)
      save_plot(plot_obj, file)
    }
  )
}

# =============================================================================
# Theme and Aesthetic Functions
# =============================================================================

#' Custom theme for EB-TAT plots
theme_ebtat_custom <- function(base_size = 14, base_family = "") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = base_size * 1.2),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(size = base_size * 0.8),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = base_size * 0.8),
      panel.grid.major = element_line(color = "grey90", size = 0.2),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey70", fill = NA, size = 0.5),
      strip.background = element_rect(fill = "grey90", color = "grey70"),
      strip.text = element_text(face = "bold")
    )
}

#' Color palette for conditions
get_condition_palette <- function() {
  c(
    "Control" = "#1F78B4",
    "Treatment" = "#E31A1C", 
    "Case" = "#33A02C",
    "Normal" = "#FF7F00",
    "Disease" = "#6A3D9A"
  )
}

#' Color palette for significance
get_significance_palette <- function() {
  c(
    "Not Significant" = "gray70",
    "Up-regulated" = "#E31A1C",
    "Down-regulated" = "#1F78B4",
    "Significant" = "#33A02C"
  )
}