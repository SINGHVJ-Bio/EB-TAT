# ENABL Biomarker/Target Analysis Tool - Data Loading Utilities
# Updated for specific file formats: metadata with "Library" and ordinal regression with "Estimate", "P_value", "FDR"

# =============================================================================
# Configuration Loading
# =============================================================================

load_data_config <- function() {
  config_file <- "config/data_paths.yml"
  if (file.exists(config_file)) {
    tryCatch({
      config <- yaml::read_yaml(config_file)
      message("Data configuration loaded successfully")
      return(config)
    }, error = function(e) {
      warning("Error loading data configuration: ", e$message)
      return(NULL)
    })
  } else {
    warning("Data configuration file not found: ", config_file)
    return(NULL)
  }
}

# =============================================================================
# Enhanced ENSEMBL to Gene Symbol Conversion Functions
# =============================================================================

#' Load local biomart and HGNC mapping files
load_local_mappings <- function(biomart_path = NULL, hgnc_path = NULL) {
  mappings <- list(biomart = NULL, hgnc = NULL)
  
  # Load biomart data
  if (!is.null(biomart_path) && file.exists(biomart_path)) {
    tryCatch({
      if (grepl("\\.tsv$", biomart_path)) {
        mappings$biomart <- read.delim(biomart_path, stringsAsFactors = FALSE, check.names = FALSE)
      } else if (grepl("\\.csv$", biomart_path)) {
        mappings$biomart <- read.csv(biomart_path, stringsAsFactors = FALSE, check.names = FALSE)
      }
      message("Loaded local biomart mapping: ", nrow(mappings$biomart), " entries")
    }, error = function(e) {
      warning("Error loading biomart file: ", e$message)
    })
  } else {
    message("Biomart file not found or path not specified: ", biomart_path)
  }
  
  # Load HGNC data - FIXED for duplicate row names
  if (!is.null(hgnc_path) && file.exists(hgnc_path)) {
    tryCatch({
      if (grepl("\\.csv$", hgnc_path)) {
        mappings$hgnc <- read.csv(hgnc_path, stringsAsFactors = FALSE, check.names = FALSE, row.names = NULL)
      } else if (grepl("\\.tsv$", hgnc_path)) {
        mappings$hgnc <- read.delim(hgnc_path, stringsAsFactors = FALSE, check.names = FALSE, row.names = NULL)
      } else if (grepl("\\.txt$", hgnc_path)) {
        mappings$hgnc <- read.delim(hgnc_path, stringsAsFactors = FALSE, check.names = FALSE, row.names = NULL)
      }
      message("Loaded local HGNC mapping: ", nrow(mappings$hgnc), " entries")
    }, error = function(e) {
      warning("Error loading HGNC file: ", e$message)
    })
  } else {
    message("HGNC file not found or path not specified: ", hgnc_path)
  }
  
  return(mappings)
}

#' Enhanced ENSEMBL to symbol conversion with HGNC file
convert_ensembl_to_symbol_enhanced <- function(ensembl_ids, local_mappings = NULL) {
  if (is.null(ensembl_ids) || length(ensembl_ids) == 0) {
    return(character(0))
  }
  
  # Filter for valid ENSEMBL IDs
  valid_ensembl <- ensembl_ids[grepl("^ENSG", ensembl_ids)]
  if (length(valid_ensembl) == 0) {
    message("No valid ENSEMBL IDs found for conversion")
    return(ensembl_ids)  # Return original IDs
  }
  
  message("Converting ", length(valid_ensembl), " ENSEMBL IDs to gene symbols...")
  
  # Try HGNC mapping first (this is what we want to use)
  if (!is.null(local_mappings$hgnc)) {
    hgnc_data <- local_mappings$hgnc
    
    # Check for ENSEMBL ID column in HGNC file
    ensembl_col <- NULL
    symbol_col <- NULL
    
    # Look for ENSEMBL ID column
    possible_ensembl_cols <- c("ensembl_gene_id", "Ensembl Gene ID", "ensembl", "ensembl_id")
    for (col in possible_ensembl_cols) {
      if (col %in% colnames(hgnc_data)) {
        ensembl_col <- col
        break
      }
    }
    
    # Look for symbol column  
    possible_symbol_cols <- c("symbol", "Symbol", "hgnc_symbol", "HGNC symbol")
    for (col in possible_symbol_cols) {
      if (col %in% colnames(hgnc_data)) {
        symbol_col <- col
        break
      }
    }
    
    if (!is.null(ensembl_col) && !is.null(symbol_col)) {
      message("Using HGNC file for ENSEMBL ID conversion")
      message("ENSEMBL column: ", ensembl_col)
      message("Symbol column: ", symbol_col)
      
      # Create mapping
      hgnc_mapping <- hgnc_data[, c(ensembl_col, symbol_col)]
      colnames(hgnc_mapping) <- c("ensembl_id", "symbol")
      
      # Remove rows with missing values
      hgnc_mapping <- hgnc_mapping[!is.na(hgnc_mapping$ensembl_id) & 
                                   !is.na(hgnc_mapping$symbol) & 
                                   hgnc_mapping$symbol != "", ]
      
      # Remove duplicates (keep first occurrence)
      hgnc_mapping <- hgnc_mapping[!duplicated(hgnc_mapping$ensembl_id), ]
      
      message("HGNC mapping entries: ", nrow(hgnc_mapping))
      
      # Perform conversion
      matched <- match(valid_ensembl, hgnc_mapping$ensembl_id)
      converted <- hgnc_mapping$symbol[matched]
      
      successful <- sum(!is.na(converted) & converted != "")
      message("HGNC mapping conversion: ", successful, "/", length(valid_ensembl), " successful conversions")
      
      if (successful > 0) {
        final_converted <- sapply(ensembl_ids, function(id) {
          if (id %in% valid_ensembl) {
            idx <- match(id, hgnc_mapping$ensembl_id)
            if (!is.na(idx) && hgnc_mapping$symbol[idx] != "") {
              return(hgnc_mapping$symbol[idx])
            }
          }
          return(id)  # Return original if no conversion
        })
        return(unname(final_converted))
      }
    } else {
      message("Could not find appropriate columns in HGNC file")
      message("Available columns: ", paste(colnames(hgnc_data), collapse = ", "))
    }
  }
  
  # Fallback to biomaRt if HGNC mapping fails
  message("Falling back to biomaRt for ENSEMBL ID conversion...")
  tryCatch({
    if (!requireNamespace("biomaRt", quietly = TRUE)) {
      stop("biomaRt package not available for ENSEMBL ID conversion")
    }
    
    mart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    gene_symbols <- biomaRt::getBM(
      attributes = c("ensembl_gene_id", "hgnc_symbol"),
      filters = "ensembl_gene_id", 
      values = valid_ensembl,
      mart = mart
    )
    
    mapping <- setNames(gene_symbols$hgnc_symbol, gene_symbols$ensembl_gene_id)
    
    final_converted <- sapply(ensembl_ids, function(id) {
      if (id %in% names(mapping) && mapping[id] != "") {
        return(mapping[id])
      } else {
        return(id)
      }
    })
    
    successful <- sum(final_converted != ensembl_ids & final_converted != "")
    message("BiomaRt conversion: ", successful, "/", length(valid_ensembl), " successful conversions")
    
    return(unname(final_converted))
    
  }, error = function(e) {
    warning("Error converting ENSEMBL IDs to gene symbols using biomaRt: ", e$message)
    message("Using original ENSEMBL IDs as symbols")
    return(ensembl_ids)
  })
}

#' Convert ENSEMBL gene IDs to gene symbols (wrapper function)
convert_ensembl_to_symbol <- function(ensembl_ids) {
  # This will be overridden in global.R with the enhanced version
  # This is a fallback implementation
  if (exists("convert_ensembl_to_symbol_enhanced", mode = "function")) {
    return(convert_ensembl_to_symbol_enhanced(ensembl_ids))
  } else {
    # Basic fallback - return original IDs
    if (is.null(ensembl_ids) || length(ensembl_ids) == 0) {
      return(character(0))
    }
    return(ensembl_ids)
  }
}

#' Add symbol column to data using ENSEMBL ID conversion
add_symbol_column <- function(data, id_column = "gene") {
  if (is.null(data) || !id_column %in% colnames(data)) {
    return(data)
  }
  
  # Check if we already have a symbol column
  if ("symbol" %in% colnames(data)) {
    message("Symbol column already exists, skipping conversion")
    return(data)
  }
  
  # Check if we have ENSEMBL IDs
  if (all(grepl("^ENSG", data[[id_column]]))) {
    message("Converting ENSEMBL IDs in '", id_column, "' column to gene symbols")
    data$symbol <- convert_ensembl_to_symbol(data[[id_column]])
  } else {
    message("No ENSEMBL IDs found in '", id_column, "' column, using as symbols")
    data$symbol <- data[[id_column]]
  }
  
  return(data)
}

# =============================================================================
# Core Data Loading Functions
# =============================================================================

#' Safely load a data file with error handling
safe_load_file <- function(file_path, loader_func = read.csv, ...) {
  if (is.null(file_path) || !file.exists(file_path)) {
    warning("File not found: ", file_path)
    return(NULL)
  }
  
  tryCatch({
    data <- loader_func(file_path, ...)
    message("Loaded: ", basename(file_path), " (", nrow(data), " rows, ", ncol(data), " columns)")
    return(data)
  }, error = function(e) {
    warning("Error loading ", basename(file_path), ": ", e$message)
    return(NULL)
  })
}

#' Load differential expression data (supports CSV, TSV, RDS)
load_differential_expression <- function(config) {
  if (is.null(config$diffgenes_path)) {
    warning("Differential expression path not configured")
    return(NULL)
  }
  
  file_path <- config$diffgenes_path
  if (!file.exists(file_path)) {
    warning("Differential expression file not found: ", file_path)
    return(NULL)
  }
  
  tryCatch({
    ext <- tolower(tools::file_ext(file_path))
    
    if (ext == "rds") {
      de_data <- readRDS(file_path)
    } else if (ext %in% c("csv", "txt")) {
      de_data <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
    } else if (ext == "tsv") {
      de_data <- read.delim(file_path, stringsAsFactors = FALSE, check.names = FALSE)
    } else {
      warning("Unsupported differential expression format: ", ext)
      return(NULL)
    }
    
    # Add symbol column if missing
    de_data <- add_symbol_column(de_data, "gene")
    
    # Standardize column names
    de_data <- standardize_de_columns(de_data)
    
    # Calculate -log10 p-value if not present
    if (!"minus_log10_padj" %in% colnames(de_data) && "padj" %in% colnames(de_data)) {
      de_data$minus_log10_padj <- -log10(de_data$padj)
    }
    
    message("Loaded differential expression data: ", basename(file_path), 
            " (", nrow(de_data), " rows, ", ncol(de_data), " columns)")
    return(de_data)
    
  }, error = function(e) {
    warning("Error loading differential expression data: ", e$message)
    return(NULL)
  })
}

#' Load count data (supports both CSV and RDS)
load_count_data <- function(config) {
  if (is.null(config$count_data_path)) {
    message("Count data path not configured - app will use minimal example data")
    return(NULL)
  }
  
  file_path <- config$count_data_path
  if (!file.exists(file_path)) {
    message("Count data file not found: ", file_path, " - app will use minimal example data")
    return(NULL)
  }
  
  tryCatch({
    ext <- tolower(tools::file_ext(file_path))
    
    if (ext == "rds") {
      count_data <- readRDS(file_path)
    } else if (ext %in% c("csv", "txt")) {
      # Try to read with row names first, then without
      count_data <- tryCatch({
        read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)
      }, error = function(e) {
        temp_data <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
        # Check if first column looks like gene names
        if (all(grepl("^ENSG|^[A-Z0-9]+", temp_data[[1]]))) {
          rownames(temp_data) <- temp_data[[1]]
          temp_data[[1]] <- NULL
        }
        return(temp_data)
      })
    } else {
      warning("Unsupported count data format: ", ext)
      return(NULL)
    }
    
    # Convert row names to symbols if they are ENSEMBL IDs
    if (!is.null(rownames(count_data)) && all(grepl("^ENSG", rownames(count_data)))) {
      message("Converting count matrix row names from ENSEMBL IDs to gene symbols")
      symbol_names <- convert_ensembl_to_symbol(rownames(count_data))
      rownames(count_data) <- symbol_names
    } else if (!"symbol" %in% colnames(count_data)) {
      # If no row names and no symbol column, create one
      if (!is.null(rownames(count_data))) {
        count_data$symbol <- rownames(count_data)
      } else {
        count_data$symbol <- paste0("GENE_", 1:nrow(count_data))
      }
    }
    
    message("Loaded count data: ", basename(file_path), 
            " (", ifelse(is.matrix(count_data), 
                         paste(dim(count_data), collapse = " x "),
                         paste(nrow(count_data), "rows,", ncol(count_data), "columns")), ")")
    return(count_data)
    
  }, error = function(e) {
    warning("Error loading count data: ", e$message, " - app will use minimal example data")
    return(NULL)
  })
}

#' Load metadata with "Library" as sample identifier
load_metadata <- function(config) {
  if (is.null(config$meta_data_path)) {
    message("Metadata path not configured - app will use minimal example data")
    return(NULL)
  }
  
  file_path <- config$meta_data_path
  if (!file.exists(file_path)) {
    message("Metadata file not found: ", file_path, " - app will use minimal example data")
    return(NULL)
  }
  
  tryCatch({
    ext <- tolower(tools::file_ext(file_path))
    
    if (ext == "rds") {
      metadata <- readRDS(file_path)
    } else if (ext %in% c("csv", "txt")) {
      metadata <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
    } else {
      warning("Unsupported metadata format: ", ext)
      return(NULL)
    }
    
    # Rename "Library" to "sample_id" for consistency
    if ("Library" %in% colnames(metadata)) {
      metadata$sample_id <- metadata$Library
      message("Renamed 'Library' column to 'sample_id'")
    } else if (!"sample_id" %in% colnames(metadata)) {
      # Try to find alternative column names
      sample_cols <- grep("sample|id|library", colnames(metadata), ignore.case = TRUE, value = TRUE)
      if (length(sample_cols) > 0) {
        metadata$sample_id <- metadata[[sample_cols[1]]]
        message("Using '", sample_cols[1], "' as sample_id")
      } else {
        warning("No sample identifier column found in metadata")
        return(NULL)
      }
    }
    
    # Check for mandatory columns and create if missing
    mandatory_cols <- c("Gender", "Age_as_at_enrolment", "Race", "Body_weight", "Body_height", 
                       "BMI", "Histology_Steatosis_grade", "Histology_Lobular_inflammation", 
                       "Histology_Ballooning", "Histology_Fibrosis", "RINe", "DISEASE_STAGES", 
                       "Cohort", "Batch", "Sex")
    
    missing_mandatory <- setdiff(mandatory_cols, colnames(metadata))
    if (length(missing_mandatory) > 0) {
      message("Missing mandatory columns in metadata: ", paste(missing_mandatory, collapse = ", "))
      # Create placeholder columns for missing mandatory ones
      for (col in missing_mandatory) {
        metadata[[col]] <- NA
        message("Created placeholder column: ", col)
      }
    }
    
    # Create condition column if not present (using Cohort or DISEASE_STAGES)
    if (!"condition" %in% colnames(metadata)) {
      if ("Cohort" %in% colnames(metadata)) {
        metadata$condition <- metadata$Cohort
        message("Using 'Cohort' as condition")
      } else if ("DISEASE_STAGES" %in% colnames(metadata)) {
        metadata$condition <- metadata$DISEASE_STAGES
        message("Using 'DISEASE_STAGES' as condition")
      } else {
        metadata$condition <- "Unknown"
        message("Created default 'condition' column")
      }
    }
    
    message("Loaded metadata: ", basename(file_path), " (", nrow(metadata), " rows, ", ncol(metadata), " columns)")
    return(metadata)
    
  }, error = function(e) {
    warning("Error loading metadata: ", e$message, " - app will use minimal example data")
    return(NULL)
  })
}

#' Load plasma proteins data
load_plasma_proteins <- function(config) {
  if (is.null(config$bld_plasma_path)) {
    message("Plasma proteins path not configured")
    return(NULL)
  }
  
  file_path <- config$bld_plasma_path
  if (!file.exists(file_path)) {
    message("Plasma proteins file not found: ", file_path)
    return(NULL)
  }
  
  tryCatch({
    if (grepl("\\.xlsx?$", file_path, ignore.case = TRUE)) {
      plasma_data <- readxl::read_excel(file_path)
    } else if (grepl("\\.csv$", file_path, ignore.case = TRUE)) {
      plasma_data <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
    } else if (grepl("\\.tsv$", file_path, ignore.case = TRUE)) {
      plasma_data <- read.delim(file_path, stringsAsFactors = FALSE, check.names = FALSE)
    } else {
      plasma_data <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
    }
    
    # Add symbol column if missing
    plasma_data <- add_symbol_column(plasma_data, "gene")
    
    message("Loaded plasma proteins data: ", nrow(plasma_data), " rows")
    return(plasma_data)
  }, error = function(e) {
    warning("Error loading plasma proteins: ", e$message)
    return(NULL)
  })
}

#' Load ordinal regression data with "Estimate", "P_value", "FDR" columns
load_ordinal_regression <- function(config) {
  if (is.null(config$ord_reg_fib_path)) {
    message("Ordinal regression path not configured - app will use minimal example data")
    return(NULL)
  }
  
  file_path <- config$ord_reg_fib_path
  if (!file.exists(file_path)) {
    message("Ordinal regression file not found: ", file_path, " - app will use minimal example data")
    return(NULL)
  }
  
  tryCatch({
    ext <- tolower(tools::file_ext(file_path))
    
    if (ext == "rds") {
      ord_data <- readRDS(file_path)
    } else if (ext %in% c("csv", "txt")) {
      ord_data <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
    } else {
      warning("Unsupported ordinal regression format: ", ext)
      return(NULL)
    }
    
    # Handle the specific ordinal regression format with unnamed first column
    if (ncol(ord_data) >= 4 && colnames(ord_data)[1] == "" || colnames(ord_data)[1] == "X") {
      # First column is gene names, rename it
      colnames(ord_data)[1] <- "symbol"
      message("Detected ordinal regression format: first column as gene symbols")
    }
    
    # Standardize column names for ordinal regression
    ord_data <- standardize_ordinal_columns(ord_data)
    
    # Ensure symbol column exists
    if (!"symbol" %in% colnames(ord_data)) {
      if ("gene" %in% colnames(ord_data)) {
        ord_data$symbol <- ord_data$gene
      } else if (ncol(ord_data) > 0 && all(grepl("^ENSG|^[A-Z0-9]+", ord_data[[1]]))) {
        ord_data$symbol <- ord_data[[1]]
        # Remove the original first column if it's not one of our standard columns
        if (!colnames(ord_data)[1] %in% c("coefficient", "p_value", "fdr", "estimate")) {
          ord_data[[1]] <- NULL
        }
      } else {
        ord_data$symbol <- paste0("GENE_", 1:nrow(ord_data))
        message("Created placeholder gene symbols for ordinal regression")
      }
    }
    
    # Convert ENSEMBL IDs if present
    if (all(grepl("^ENSG", ord_data$symbol))) {
      message("Converting ENSEMBL IDs in ordinal regression to gene symbols")
      ord_data$symbol <- convert_ensembl_to_symbol(ord_data$symbol)
    }
    
    message("Loaded ordinal regression data: ", basename(file_path), " (", nrow(ord_data), " rows)")
    return(ord_data)
    
  }, error = function(e) {
    warning("Error loading ordinal regression data: ", e$message, " - app will use minimal example data")
    return(NULL)
  })
}

#' Load STRING database
load_string_db <- function(config) {
  if (is.null(config$string_db_path)) {
    warning("STRING database path not configured")
    return(NULL)
  }
  
  file_path <- config$string_db_path
  if (!file.exists(file_path)) {
    warning("STRING database file not found: ", file_path)
    return(NULL)
  }
  
  string_data <- safe_load_file(file_path, readRDS)
  return(string_data)
}

#' Load annotation data
load_annotation_data <- function(config) {
  annotation_list <- list()
  
  # Define annotation files to load
  annotation_files <- list(
    human_tf = list(path = config$human_tf_path, loader = read.csv),
    receptor_ligand = list(path = config$rec_lig_path, loader = read.delim),
    secreted_proteins = list(path = config$secreted_path, loader = read.delim),
    sp_tm_data = list(path = config$sp_tm_path, loader = read.csv),
    hgnc_mapping = list(path = config$hgnc_data_path, loader = read.delim),
    hepatocytes = list(path = config$hcc_path, loader = read.delim),
    kupffer_cells = list(path = config$kuff_path, loader = read.delim),
    hsc_cells = list(path = config$hsc_path, loader = read.csv),
    liver_expression = list(path = config$liver_path, loader = read.delim),
    innate_db = list(path = config$innatdb_path, loader = readxl::read_excel)
  )
  
  for (name in names(annotation_files)) {
    file_info <- annotation_files[[name]]
    if (!is.null(file_info$path) && file.exists(file_info$path)) {
      annotation_list[[name]] <- safe_load_file(file_info$path, file_info$loader)
      
      # Add symbol column to annotation data if it contains gene information
      if (!is.null(annotation_list[[name]]) && 
          ("gene" %in% colnames(annotation_list[[name]]) || 
           any(grepl("symbol", colnames(annotation_list[[name]]), ignore.case = TRUE)))) {
        annotation_list[[name]] <- add_symbol_column(annotation_list[[name]], "gene")
      }
    }
  }
  
  message("Loaded ", length(annotation_list), " annotation datasets")
  return(annotation_list)
}

# =============================================================================
# Data Standardization Functions
# =============================================================================

#' Standardize differential expression column names
standardize_de_columns <- function(de_data) {
  if (is.null(de_data)) return(NULL)
  
  # Common column name mappings
  col_mappings <- list(
    symbol = c("symbol", "gene", "gene_name", "Gene", "GeneSymbol", "gene_symbol", 
               "ensembl", "ensembl_id", "ensembl_gene_id", "ENSEMBL", "ENSEMBL_ID"),
    log2FoldChange = c("log2FoldChange", "logFC", "log2FC", "fold_change", "log2FoldChange", "lfc"),
    padj = c("padj", "adj.P.Val", "FDR", "qvalue", "adj_pval", "fdr"),
    pvalue = c("pvalue", "P.Value", "p_val", "pvalue", "PValue", "p_value"),
    baseMean = c("baseMean", "base_mean", "avg_expr", "baseMean")
  )
  
  # Apply mappings
  for (standard_name in names(col_mappings)) {
    found <- FALSE
    for (possible_name in col_mappings[[standard_name]]) {
      if (possible_name %in% colnames(de_data)) {
        if (possible_name != standard_name) {
          de_data[[standard_name]] <- de_data[[possible_name]]
        }
        found <- TRUE
        break
      }
    }
  }
  
  # Ensure symbol column exists
  if (!"symbol" %in% colnames(de_data)) {
    # Try to use row names if they are ENSEMBL IDs
    if (!is.null(rownames(de_data)) && all(grepl("^ENSG", rownames(de_data)))) {
      de_data$symbol <- convert_ensembl_to_symbol(rownames(de_data))
      message("Converted row names from ENSEMBL IDs to gene symbols")
    } else {
      # Create placeholder gene symbols
      de_data$symbol <- paste0("GENE_", 1:nrow(de_data))
      message("Created placeholder gene symbols")
    }
  }
  
  return(de_data)
}

#' Standardize ordinal regression column names
standardize_ordinal_columns <- function(ord_data) {
  if (is.null(ord_data)) return(NULL)
  
  # Column mappings for ordinal regression
  ord_mappings <- list(
    coefficient = c("coefficient", "Estimate", "estimate", "coef", "beta"),
    p_value = c("p_value", "P_value", "pvalue", "PValue", "p.val"),
    fdr = c("fdr", "FDR", "adj_pval", "padj", "adj.P.Val")
  )
  
  # Apply mappings
  for (standard_name in names(ord_mappings)) {
    found <- FALSE
    for (possible_name in ord_mappings[[standard_name]]) {
      if (possible_name %in% colnames(ord_data)) {
        if (possible_name != standard_name) {
          ord_data[[standard_name]] <- ord_data[[possible_name]]
        }
        found <- TRUE
        break
      }
    }
  }
  
  return(ord_data)
}

#' Validate differential expression data structure
validate_de_data <- function(de_data) {
  if (is.null(de_data)) {
    return(FALSE)
  }
  
  required <- c("symbol")
  recommended <- c("log2FoldChange", "padj")
  
  # Check required columns
  missing_required <- setdiff(required, colnames(de_data))
  if (length(missing_required) > 0) {
    warning("Missing required columns: ", paste(missing_required, collapse = ", "))
    return(FALSE)
  }
  
  # Check recommended columns
  missing_recommended <- setdiff(recommended, colnames(de_data))
  if (length(missing_recommended) > 0) {
    warning("Missing recommended columns: ", paste(missing_recommended, collapse = ", "))
  }
  
  # Check data quality
  if (nrow(de_data) == 0) {
    warning("Data has 0 rows")
    return(FALSE)
  }
  
  # Check if we have ENSEMBL IDs and warn
  if ("symbol" %in% colnames(de_data) && any(grepl("^ENSG", de_data$symbol))) {
    ensembl_count <- sum(grepl("^ENSG", de_data$symbol))
    message("Note: ", ensembl_count, " ENSEMBL IDs found in symbol column")
    message("Some features like gene set enrichment may work better with gene symbols")
  }
  
  return(TRUE)
}

# =============================================================================
# Data Filtering Functions
# =============================================================================

#' Filter data based on fold change and p-value thresholds
filter_de_data <- function(de_data, fc_threshold = c(-2, 2), p_threshold = 0.05) {
  if (is.null(de_data)) return(NULL)
  
  filtered_data <- de_data
  
  # Apply fold change filter if column exists
  if ("log2FoldChange" %in% colnames(filtered_data)) {
    filtered_data <- filtered_data[
      filtered_data$log2FoldChange >= fc_threshold[1] & 
      filtered_data$log2FoldChange <= fc_threshold[2], 
    ]
  }
  
  # Apply p-value filter if column exists
  if ("padj" %in% colnames(filtered_data)) {
    filtered_data <- filtered_data[filtered_data$padj <= p_threshold, ]
  } else if ("pvalue" %in% colnames(filtered_data)) {
    filtered_data <- filtered_data[filtered_data$pvalue <= p_threshold, ]
  }
  
  message("Filtered data: ", nrow(filtered_data), " rows remaining")
  return(filtered_data)
}

#' Get top N differentially expressed genes
get_top_genes <- function(de_data, n = 20, by = "padj") {
  if (is.null(de_data)) return(NULL)
  
  if (!by %in% colnames(de_data)) {
    warning("Column not found for ranking: ", by)
    return(de_data[1:min(n, nrow(de_data)), ])
  }
  
  # Remove NA values
  valid_data <- de_data[!is.na(de_data[[by]]), ]
  
  # Sort and get top N
  if (by %in% c("padj", "pvalue")) {
    # For p-values, lower is better
    top_genes <- valid_data[order(valid_data[[by]]), ]
  } else if (by %in% c("log2FoldChange", "minus_log10_padj")) {
    # For fold change, absolute value is better
    top_genes <- valid_data[order(abs(valid_data[[by]]), decreasing = TRUE), ]
  } else {
    # Default sorting
    top_genes <- valid_data[order(valid_data[[by]], decreasing = TRUE), ]
  }
  
  return(top_genes[1:min(n, nrow(top_genes)), ])
}

# =============================================================================
# Data Integration Functions
# =============================================================================

#' Merge differential expression data with annotation data
merge_with_annotations <- function(de_data, annotation_data) {
  if (is.null(de_data) || is.null(annotation_data)) {
    return(de_data)
  }
  
  merged_data <- de_data
  
  # Merge with HGNC mapping if available
  if (!is.null(annotation_data$hgnc_mapping) && "symbol" %in% colnames(de_data)) {
    # Simple merge on symbol
    common_cols <- setdiff(colnames(annotation_data$hgnc_mapping), colnames(merged_data))
    if (length(common_cols) > 0) {
      merged_data <- merge(merged_data, annotation_data$hgnc_mapping, 
                          by = "symbol", all.x = TRUE)
    }
  }
  
  return(merged_data)
}

#' Get gene symbols from various data sources
extract_gene_symbols <- function(data, symbol_col = "symbol") {
  if (is.null(data)) return(character(0))
  
  if (symbol_col %in% colnames(data)) {
    symbols <- unique(na.omit(data[[symbol_col]]))
    return(as.character(symbols))
  } else {
    warning("Symbol column not found: ", symbol_col)
    return(character(0))
  }
}

# =============================================================================
# Data Access Helper Functions
# =============================================================================

#' Get plasma protein information for genes
get_plasma_protein_info <- function(gene_symbols, plasma_data) {
  if (is.null(plasma_data) || is.null(gene_symbols)) {
    return(NULL)
  }
  
  # Find matching genes in plasma data
  if ("symbol" %in% colnames(plasma_data)) {
    protein_info <- plasma_data[plasma_data$symbol %in% gene_symbols, ]
  } else if ("gene" %in% colnames(plasma_data)) {
    protein_info <- plasma_data[plasma_data$gene %in% gene_symbols, ]
  } else {
    # Try first column
    protein_info <- plasma_data[plasma_data[[1]] %in% gene_symbols, ]
  }
  
  return(protein_info)
}

#' Get ordinal regression results for genes
get_ordinal_regression_info <- function(gene_symbols, ord_data) {
  if (is.null(ord_data) || is.null(gene_symbols)) {
    return(NULL)
  }
  
  # Find matching genes in ordinal regression data
  if ("symbol" %in% colnames(ord_data)) {
    ord_info <- ord_data[ord_data$symbol %in% gene_symbols, ]
  } else if ("gene" %in% colnames(ord_data)) {
    ord_info <- ord_data[ord_data$gene %in% gene_symbols, ]
  } else {
    # Try first column
    ord_info <- ord_data[ord_data[[1]] %in% gene_symbols, ]
  }
  
  return(ord_info)
}

#' Get count data for specific genes
get_gene_count_data <- function(gene_symbols, count_data, metadata = NULL) {
  if (is.null(count_data) || is.null(gene_symbols)) {
    return(NULL)
  }
  
  # Extract count data for specific genes
  # This depends on the structure of your count data
  if (is.matrix(count_data) || is.data.frame(count_data)) {
    # Assuming rows are genes and columns are samples
    if (all(gene_symbols %in% rownames(count_data))) {
      gene_counts <- count_data[gene_symbols, , drop = FALSE]
    } else if ("symbol" %in% colnames(count_data)) {
      gene_counts <- count_data[count_data$symbol %in% gene_symbols, ]
    } else {
      # Try to match by row names
      matched_genes <- intersect(gene_symbols, rownames(count_data))
      if (length(matched_genes) > 0) {
        gene_counts <- count_data[matched_genes, , drop = FALSE]
      } else {
        warning("No matching genes found in count data")
        return(NULL)
      }
    }
  } else {
    warning("Unsupported count data format")
    return(NULL)
  }
  
  return(gene_counts)
}

# =============================================================================
# Data Summary Functions
# =============================================================================

#' Generate summary statistics for differential expression data
summarize_de_data <- function(de_data) {
  if (is.null(de_data)) {
    return(list(
      total_genes = 0,
      message = "No data available"
    ))
  }
  
  summary_list <- list(
    total_genes = nrow(de_data),
    columns = colnames(de_data)
  )
  
  # Add statistical summaries if available
  if ("log2FoldChange" %in% colnames(de_data)) {
    summary_list$fc_summary <- summary(de_data$log2FoldChange)
  }
  
  if ("padj" %in% colnames(de_data)) {
    summary_list$significant_genes <- sum(de_data$padj < 0.05, na.rm = TRUE)
  }
  
  if ("pvalue" %in% colnames(de_data)) {
    summary_list$pvalue_summary <- summary(de_data$pvalue)
  }
  
  return(summary_list)
}

#' Print data loading summary
print_data_summary <- function(data_objects) {
  message("\n=== DATA LOADING SUMMARY ===")
  
  summary_items <- list(
    "Differential Expression" = data_objects$diffgenes,
    "Count Data" = data_objects$count_data,
    "Metadata" = data_objects$metadata,
    "Plasma Proteins" = data_objects$plasma_proteins,
    "Ordinal Regression" = data_objects$ordinal_regression,
    "STRING Database" = data_objects$string_db,
    "Annotation Data" = data_objects$annotation_data
  )
  
  for (name in names(summary_items)) {
    obj <- summary_items[[name]]
    if (is.data.frame(obj)) {
      message(name, ": ", nrow(obj), " rows, ", ncol(obj), " columns")
    } else if (is.matrix(obj)) {
      message(name, ": matrix ", paste(dim(obj), collapse = " x "))
    } else if (is.list(obj) && name == "annotation_data") {
      message(name, ": list with ", length(obj), " annotation datasets")
    } else if (!is.null(obj)) {
      message(name, ": loaded (type: ", class(obj), ")")
    } else {
      message(name, ": NOT LOADED")
    }
  }
  
  # Literature trends configuration
  if (!is.null(data_objects$literature_config)) {
    message("Literature Trends: configured")
    if (!is.null(data_objects$literature_config$disease_terms)) {
      message("  Disease terms: ", paste(data_objects$literature_config$disease_terms, collapse = ", "))
    }
    if (!is.null(data_objects$literature_config$years_back)) {
      message("  Years back: ", data_objects$literature_config$years_back)
    }
    if (!is.null(data_objects$literature_config$email)) {
      message("  Email: ", data_objects$literature_config$email)
    }
  } else {
    message("Literature Trends: NOT CONFIGURED")
  }
  
  message("Data status: ", data_objects$data_status)
  message("Load time: ", data_objects$load_time)
  message("============================\n")
}