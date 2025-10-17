# ENABL Biomarker/Target Analysis Tool - Helper Functions
# General utility functions for data manipulation and analysis

# =============================================================================
# String Manipulation Functions - UPDATED FOR ENSEMBL ID SUPPORT
# =============================================================================

#' Clean and standardize gene symbols - UPDATED FOR ENSEMBL ID SUPPORT
clean_gene_symbols <- function(symbols) {
  if (is.null(symbols)) return(character(0))
  
  cleaned <- toupper(trimws(as.character(symbols)))
  # Allow ENSEMBL IDs (ENSG...) and standard gene symbols
  cleaned <- gsub("[^A-Z0-9.-]", "", cleaned)  # Remove special characters but keep dots and hyphens
  cleaned <- cleaned[cleaned != "" & !is.na(cleaned)]
  
  return(unique(cleaned))
}

#' Extract genes from text input - UPDATED FOR ENSEMBL ID SUPPORT
parse_gene_input <- function(input_text, sep = "[,;\\s\\n]+") {
  if (is.null(input_text) || input_text == "") {
    return(character(0))
  }
  
  genes <- unlist(strsplit(input_text, sep))
  genes <- clean_gene_symbols(genes)
  
  return(genes)
}

#' Convert various gene identifiers to standard symbols
convert_gene_identifiers <- function(identifiers, from_type = "auto") {
  if (is.null(identifiers) || length(identifiers) == 0) {
    return(character(0))
  }
  
  # Simple cleaning for now - could be extended with biomaRt or similar
  converted <- clean_gene_symbols(identifiers)
  return(converted)
}

# =============================================================================
# Statistical Helper Functions
# =============================================================================

#' Calculate -log10 of p-values with safety checks
safe_minus_log10 <- function(p_values, na_value = 0) {
  if (is.null(p_values)) return(numeric(0))
  
  # Handle zeros and NAs
  safe_p <- pmax(p_values, .Machine$double.eps)  # Avoid log(0)
  safe_p[is.na(safe_p)] <- 1  # Set NA to non-significant
  
  result <- -log10(safe_p)
  result[is.infinite(result)] <- na_value  # Handle infinite values
  
  return(result)
}

#' Apply multiple testing correction
apply_multiple_testing_correction <- function(p_values, method = "BH") {
  if (is.null(p_values)) return(numeric(0))
  
  tryCatch({
    adjusted_p <- p.adjust(p_values, method = method)
    return(adjusted_p)
  }, error = function(e) {
    warning("Error in multiple testing correction: ", e$message)
    return(p_values)
  })
}

#' Calculate fold change statistics
calculate_fc_stats <- function(fold_changes) {
  if (is.null(fold_changes) || length(fold_changes) == 0) {
    return(list(
      mean = NA,
      median = NA,
      sd = NA,
      up_regulated = 0,
      down_regulated = 0
    ))
  }
  
  stats <- list(
    mean = mean(fold_changes, na.rm = TRUE),
    median = median(fold_changes, na.rm = TRUE),
    sd = sd(fold_changes, na.rm = TRUE),
    up_regulated = sum(fold_changes > 0, na.rm = TRUE),
    down_regulated = sum(fold_changes < 0, na.rm = TRUE)
  )
  
  return(stats)
}

# =============================================================================
# Data Validation Functions - UPDATED FOR ENSEMBL ID SUPPORT
# =============================================================================

#' Validate numeric input with ranges
validate_numeric_input <- function(value, min_val = -Inf, max_val = Inf, default = NULL) {
  if (is.null(value) || is.na(value)) {
    return(default)
  }
  
  num_value <- as.numeric(value)
  if (is.na(num_value)) {
    return(default)
  }
  
  # Apply bounds
  bounded_value <- max(min(num_value, max_val), min_val)
  return(bounded_value)
}

#' Validate gene list input - UPDATED FOR ENSEMBL ID SUPPORT
validate_gene_list <- function(gene_list, max_genes = 1000) {
  if (is.null(gene_list) || length(gene_list) == 0) {
    return(list(valid = FALSE, message = "No genes provided"))
  }
  
  # Check length
  if (length(gene_list) > max_genes) {
    return(list(
      valid = FALSE,
      message = paste("Too many genes (", length(gene_list), "). Maximum allowed:", max_genes)
    ))
  }
  
  # Check for valid gene symbols (expanded to include ENSEMBL IDs)
  valid_genes <- grepl("^[A-Z0-9]+[A-Z0-9.-]*$", gene_list) | grepl("^ENSG", gene_list)
  if (!all(valid_genes)) {
    invalid_count <- sum(!valid_genes)
    return(list(
      valid = FALSE,
      message = paste("Contains", invalid_count, "invalid gene symbols")
    ))
  }
  
  return(list(valid = TRUE, message = "Valid gene list"))
}

# =============================================================================
# File Handling Functions
# =============================================================================

#' Safe file reading with format detection
safe_read_file <- function(file_path, ...) {
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  
  file_ext <- tolower(tools::file_ext(file_path))
  
  tryCatch({
    switch(file_ext,
           csv = read.csv(file_path, ...),
           tsv = read.delim(file_path, ...),
           txt = read.delim(file_path, ...),
           xls = readxl::read_excel(file_path, ...),
           xlsx = readxl::read_excel(file_path, ...),
           rds = readRDS(file_path),
           {
             # Default to read.csv with check.names = FALSE
             read.csv(file_path, check.names = FALSE, ...)
           }
    )
  }, error = function(e) {
    stop("Error reading file ", file_path, ": ", e$message)
  })
}

#' Generate safe file names
safe_filename <- function(prefix, extension, timestamp = TRUE) {
  if (timestamp) {
    time_str <- format(Sys.time(), "%Y%m%d_%H%M%S")
    filename <- paste0(prefix, "_", time_str, ".", extension)
  } else {
    filename <- paste0(prefix, ".", extension)
  }
  
  # Remove problematic characters
  filename <- gsub("[^a-zA-Z0-9._-]", "_", filename)
  return(filename)
}

# =============================================================================
# UI Helper Functions
# =============================================================================

#' Create a loading spinner
create_loading_spinner <- function(id, text = "Loading...") {
  div(
    id = id,
    class = "loading-spinner",
    style = "text-align: center; padding: 20px;",
    tags$div(
      class = "spinner-border text-primary",
      role = "status",
      style = "width: 3rem; height: 3rem;"
    ),
    tags$br(),
    tags$span(style = "color: #6c757d;", text)
  )
}

#' Create a notification wrapper
create_notification <- function(message, type = "info", duration = 5) {
  list(
    message = message,
    type = type,
    duration = duration * 1000  # Convert to milliseconds
  )
}

#' Create a download handler factory
create_download_handler <- function(filename, content_func, content_type = "text/csv") {
  downloadHandler(
    filename = filename,
    content = content_func,
    contentType = content_type
  )
}

# =============================================================================
# Data Formatting Functions
# =============================================================================

#' Format p-values for display
format_pvalue <- function(p_values, digits = 4) {
  if (is.null(p_values)) return(character(0))
  
  sapply(p_values, function(p) {
    if (is.na(p)) {
      return("NA")
    } else if (p < 10^(-digits)) {
      return(paste0("<", format(10^(-digits), scientific = FALSE)))
    } else {
      return(format(p, digits = digits, scientific = FALSE))
    }
  })
}

#' Format fold changes for display
format_fold_change <- function(fold_changes, digits = 3) {
  if (is.null(fold_changes)) return(character(0))
  
  sapply(fold_changes, function(fc) {
    if (is.na(fc)) {
      return("NA")
    } else {
      return(format(round(fc, digits), nsmall = digits))
    }
  })
}

#' Create a summary data frame for display
create_summary_df <- function(data, max_rows = 100) {
  if (is.null(data) || nrow(data) == 0) {
    return(data.frame(Message = "No data available"))
  }
  
  # Limit rows for display
  display_data <- data[1:min(max_rows, nrow(data)), ]
  
  # Format numeric columns
  for (col in colnames(display_data)) {
    if (is.numeric(display_data[[col]])) {
      if (grepl("pval|padj", col, ignore.case = TRUE)) {
        display_data[[col]] <- format_pvalue(display_data[[col]])
      } else if (grepl("foldchange|log2fc", col, ignore.case = TRUE)) {
        display_data[[col]] <- format_fold_change(display_data[[col]])
      } else {
        display_data[[col]] <- round(display_data[[col]], 4)
      }
    }
  }
  
  return(display_data)
}

# =============================================================================
# Error Handling Functions
# =============================================================================

#' Safe function execution with error capture
safe_execute <- function(expr, error_value = NULL, warning_msg = NULL) {
  tryCatch({
    eval(expr)
  }, error = function(e) {
    if (!is.null(warning_msg)) {
      warning(warning_msg, ": ", e$message)
    }
    return(error_value)
  })
}

#' Validate function parameters
validate_parameters <- function(..., .checks) {
  errors <- character(0)
  
  for (check in .checks) {
    param_name <- check$param
    value <- get(param_name)
    
    if (check$type == "numeric") {
      if (!is.numeric(value) || is.na(value)) {
        errors <- c(errors, paste(param_name, "must be numeric"))
      }
      if (!is.null(check$min) && value < check$min) {
        errors <- c(errors, paste(param_name, "must be >=", check$min))
      }
      if (!is.null(check$max) && value > check$max) {
        errors <- c(errors, paste(param_name, "must be <=", check$max))
      }
    } else if (check$type == "character") {
      if (!is.character(value) || nchar(value) == 0) {
        errors <- c(errors, paste(param_name, "must be non-empty string"))
      }
    } else if (check$type == "logical") {
      if (!is.logical(value)) {
        errors <- c(errors, paste(param_name, "must be logical"))
      }
    }
  }
  
  if (length(errors) > 0) {
    stop(paste(errors, collapse = "; "))
  }
  
  return(TRUE)
}

# =============================================================================
# Color and Theme Helpers
# =============================================================================

#' Get color palette for categories
get_category_palette <- function(n_categories) {
  if (n_categories <= 8) {
    return(c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000"))
  } else {
    return(rainbow(n_categories))
  }
}

#' Create gradient color scale
create_gradient_colors <- function(values, low_color = "blue", high_color = "red") {
  if (length(values) == 0) return(character(0))
  
  # Normalize values to 0-1 range
  normalized <- (values - min(values, na.rm = TRUE)) / 
    (max(values, na.rm = TRUE) - min(values, na.rm = TRUE))
  normalized[is.na(normalized)] <- 0.5
  
  # Create color gradient
  colorRampPalette(c(low_color, high_color))(100)[round(normalized * 99) + 1]
}

# =============================================================================
# Debugging and Logging Functions
# =============================================================================

#' Log function execution
log_execution <- function(function_name, parameters = NULL, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message <- paste(timestamp, level, function_name)
  
  if (!is.null(parameters)) {
    message <- paste0(message, " - Parameters: ", paste(names(parameters), parameters, sep = "=", collapse = ", "))
  }
  
  cat(message, "\n")
}

#' Create a debug information panel
create_debug_panel <- function(session_data, include_session_info = TRUE) {
  debug_info <- list(
    timestamp = Sys.time(),
    data_loaded = !is.null(session_data$data),
    data_rows = ifelse(!is.null(session_data$data), nrow(session_data$data), 0),
    memory_usage = paste0(round(pryr::mem_used() / 1024^2, 1), " MB")
  )
  
  if (include_session_info) {
    debug_info$session <- list(
      r_version = R.version.string,
      platform = R.version$platform,
      locale = Sys.getlocale()
    )
  }
  
  return(debug_info)
}