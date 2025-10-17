# ENABL Biomarker/Target Analysis Tool - Encoding Check Utility
# This script checks file encodings and helps diagnose encoding issues

# Function to check file encoding
check_file_encoding <- function(file_path) {
  if (!file.exists(file_path)) {
    message("File not found: ", file_path)
    return(NULL)
  }
  
  tryCatch({
    # Try to read with different encodings
    encodings <- c("UTF-8", "Latin-1", "ASCII", "Windows-1252")
    
    results <- list()
    for (enc in encodings) {
      success <- FALSE
      tryCatch({
        test_data <- readLines(file_path, encoding = enc, n = 10)
        success <- TRUE
        results[[enc]] <- list(
          success = TRUE,
          preview = head(test_data, 5),
          lines_read = length(test_data)
        )
      }, error = function(e) {
        results[[enc]] <- list(success = FALSE, error = e$message)
      })
    }
    
    # Check which encodings worked
    working_encodings <- names(results)[sapply(results, function(x) x$success)]
    
    message("File: ", basename(file_path))
    message("Working encodings: ", paste(working_encodings, collapse = ", "))
    
    if (length(working_encodings) > 0) {
      message("Preview with ", working_encodings[1], ":")
      print(results[[working_encodings[1]]]$preview)
    }
    
    return(results)
    
  }, error = function(e) {
    message("Error checking encoding for ", basename(file_path), ": ", e$message)
    return(NULL)
  })
}

# Function to detect column separators
detect_separator <- function(file_path, n_lines = 10) {
  if (!file.exists(file_path)) {
    message("File not found: ", file_path)
    return(NULL)
  }
  
  tryCatch({
    lines <- readLines(file_path, n = n_lines)
    
    # Count occurrences of common separators
    separators <- list(
      comma = ",",
      tab = "\t",
      semicolon = ";",
      pipe = "\\|"
    )
    
    counts <- sapply(separators, function(sep) {
      sum(sapply(lines, function(line) {
        length(strsplit(line, sep)[[1]])
      }))
    })
    
    # Find the most likely separator
    likely_sep <- names(which.max(counts))
    max_count <- max(counts)
    
    message("File: ", basename(file_path))
    message("Likely separator: ", likely_sep, " (count: ", max_count, ")")
    message("All separator counts:")
    print(counts)
    
    return(list(
      separator = separators[[likely_sep]],
      counts = counts,
      lines_checked = length(lines)
    ))
    
  }, error = function(e) {
    message("Error detecting separator for ", basename(file_path), ": ", e$message)
    return(NULL)
  })
}

# Function to check common data files
check_data_files_encoding <- function(data_dir = "data") {
  if (!dir.exists(data_dir)) {
    message("Data directory not found: ", data_dir)
    return(NULL)
  }
  
  message("Checking encoding for all data files in: ", data_dir)
  
  # Find all data files
  data_files <- list.files(data_dir, 
                          pattern = "\\.(csv|tsv|txt|xls|xlsx|rds|RDS)$", 
                          recursive = TRUE, 
                          full.names = TRUE)
  
  results <- list()
  
  for (file in data_files) {
    message("\n", rep("-", 50))
    message("Checking: ", file)
    
    # Check encoding
    encoding_result <- check_file_encoding(file)
    
    # Check separator (for text files)
    if (grepl("\\.(csv|tsv|txt)$", file)) {
      separator_result <- detect_separator(file)
    } else {
      separator_result <- "Binary file (RDS/Excel)"
    }
    
    results[[file]] <- list(
      encoding = encoding_result,
      separator = separator_result,
      file_size = file.info(file)$size,
      file_mtime = file.info(file)$mtime
    )
  }
  
  return(results)
}

# Function to fix encoding issues
fix_file_encoding <- function(file_path, output_path = NULL, from_encoding = "Latin-1", to_encoding = "UTF-8") {
  if (is.null(output_path)) {
    output_path <- file.path(dirname(file_path), 
                            paste0("fixed_", basename(file_path)))
  }
  
  tryCatch({
    # Read with source encoding
    content <- readLines(file_path, encoding = from_encoding)
    
    # Write with target encoding
    writeLines(content, output_path, useBytes = FALSE)
    
    message("File converted successfully:")
    message("  Input: ", file_path, " (", from_encoding, ")")
    message("  Output: ", output_path, " (", to_encoding, ")")
    
    return(output_path)
    
  }, error = function(e) {
    message("Error converting file: ", e$message)
    return(NULL)
  })
}

# Main function to run comprehensive check
run_comprehensive_encoding_check <- function() {
  message("=== EB-TAT Comprehensive Encoding Check ===")
  message("Time: ", Sys.time())
  message("Working directory: ", getwd())
  
  # Check R version and locale
  message("\n=== System Information ===")
  message("R version: ", R.version.string)
  message("System: ", R.version$system)
  message("Locale: ", Sys.getlocale())
  
  # Check package availability
  message("\n=== Required Packages ===")
  required_pkgs <- c("readr", "readxl", "tools")
  for (pkg in required_pkgs) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      message(pkg, ": Available")
    } else {
      message(pkg, ": NOT AVAILABLE")
    }
  }
  
  # Check data files
  message("\n=== Data Files Encoding Check ===")
  data_check <- check_data_files_encoding()
  
  # Summary
  message("\n=== Summary ===")
  if (!is.null(data_check)) {
    message("Total files checked: ", length(data_check))
    
    # Count files by type
    file_types <- table(sapply(names(data_check), function(x) {
      tools::file_ext(x)
    }))
    message("File types:")
    print(file_types)
  }
  
  return(data_check)
}

# Run if called directly
if (sys.nframe() == 0) {
  results <- run_comprehensive_encoding_check()
  
  # Save results to file
  output_file <- paste0("encoding_check_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
  sink(output_file)
  print(results)
  sink()
  
  message("\nResults saved to: ", output_file)
}