# ENABL Biomarker/Target Analysis Tool - Global Configuration
# Updated for specific file formats: metadata with "Library" and ordinal regression with "Estimate", "P_value", "FDR"

# =============================================================================
# Package Availability Check
# =============================================================================

# Function to check package availability without installation
check_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    message("Package not available: ", package)
    return(FALSE)
  }

  # Try to load the package
  tryCatch({
    library(package, character.only = TRUE)
    message("Successfully loaded: ", package)
    return(TRUE)
  }, error = function(e) {
    warning("Failed to load package: ", package, " - ", e$message)
    return(FALSE)
  })
}

# List of required CRAN packages
cran_packages <- c(
  "shiny", "ggplot2", "magrittr", "dplyr", "ggrepel", "shinycssloaders",
  "readxl", "DT", "RCurl", "plot.matrix", "psych", "reshape2",
  "ggVennDiagram", "ggvenn", "readr", "ggupset", "UpSetR",
  "limma", "pheatmap", "hash", "gridExtra", "shinythemes",
  "tidyr", "purrr", "stringr", "tibble", "forcats", "scales", "yaml",
  "rentrez", "xml2", "jsonlite", "httr", "biomaRt"
)

# List of required Bioconductor packages
bioc_packages <- c(
  "clusterProfiler", "enrichplot", "DOSE", "org.Hs.eg.db"
)

# Check CRAN packages
message("Checking CRAN packages...")
cran_loaded <- sapply(cran_packages, check_package)

# Check Bioconductor packages
message("Checking Bioconductor packages...")
bioc_loaded <- sapply(bioc_packages, check_package)

# Check if all essential packages are available
essential_packages <- c("shiny", "ggplot2", "dplyr", "DT", "biomaRt")
missing_essential <- !sapply(essential_packages, function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
})

if (any(missing_essential)) {
  stop("Missing essential packages: ",
       paste(names(missing_essential)[missing_essential], collapse = ", "),
       "\nPlease install missing packages before running the application.")
} else {
  message("All essential packages available!")
}

# Summary of package availability
message("\n=== Package Availability Summary ===")
message("CRAN packages available: ", sum(cran_loaded), "/", length(cran_packages))
message("Bioconductor packages available: ", sum(bioc_loaded), "/", length(bioc_packages))

# Exit if any packages are missing
if (sum(cran_loaded) < length(cran_packages) || sum(bioc_loaded) < length(bioc_packages)) {
  missing_cran <- cran_packages[!cran_loaded]
  missing_bioc <- bioc_packages[!bioc_loaded]

  message("\nMissing CRAN packages: ", paste(missing_cran, collapse = ", "))
  message("Missing Bioconductor packages: ", paste(missing_bioc, collapse = ", "))
  message("\nPlease install missing packages using:")
  if (length(missing_cran) > 0) {
    message("install.packages(c(\"", paste(missing_cran, collapse = "\", \""), "\"))")
  }
  if (length(missing_bioc) > 0) {
    message("BiocManager::install(c(\"", paste(missing_bioc, collapse = "\", \""), "\"))")
  }
  stop("Application cannot start due to missing packages.")
}

# =============================================================================
# Utility Functions Loading
# =============================================================================

# Function to safely source utility files
safe_source <- function(file) {
  if (file.exists(file)) {
    tryCatch({
      source(file)
      message("Successfully sourced: ", file)
      TRUE
    }, error = function(e) {
      warning("Error sourcing ", file, ": ", e$message)
      FALSE
    })
  } else {
    warning("Utility file not found: ", file)
    FALSE
  }
}

# Source utility functions FIRST - before they are used
message("Loading utility functions...")
utility_files <- c(
  "utils/helpers.R",
  "utils/data_loading.R",
  "utils/plotting.R"
)

# Create utils directory if it doesn't exist
if (!dir.exists("utils")) {
  dir.create("utils", recursive = TRUE)
  message("Created utils directory")
}

# Source available utility files
utils_loaded <- sapply(utility_files, safe_source)

# =============================================================================
# Data Paths Configuration
# =============================================================================

# Load data paths configuration
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

# Get data paths
data_paths <- load_data_config()

# =============================================================================
# Enhanced ENSEMBL to Symbol Conversion with Local Files
# =============================================================================

# Load local biomart and HGNC mapping files
local_mappings <- load_local_mappings(
  biomart_path = if (!is.null(data_paths)) data_paths$biomart_path else NULL,
  hgnc_path = if (!is.null(data_paths)) data_paths$hgnc_path else NULL
)

# Override the convert_ensembl_to_symbol function in global environment
.GlobalEnv$convert_ensembl_to_symbol <- function(ensembl_ids) {
  convert_ensembl_to_symbol_enhanced(ensembl_ids, local_mappings)
}

# =============================================================================
# Consistent Minimal Data Creation
# =============================================================================

# Common gene set for all minimal data (1000 genes)
common_genes <- c(
  "TP53", "EGFR", "BRCA1", "BRCA2", "MYC", "KRAS", "PTEN", "VEGFA", "CDKN2A", "RB1",
  "AKT1", "PIK3CA", "ERBB2", "MET", "BRAF", "ALK", "ROS1", "FGFR1", "FGFR2", "FGFR3",
  "IDH1", "IDH2", "FLT3", "KIT", "PDGFRA", "JAK2", "STAT3", "NF1", "NF2", "APC",
  "SMAD4", "TGFBR2", "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", "DLL4", "JAG1", "JAG2",
  "WNT1", "WNT2", "WNT3", "WNT4", "WNT5A", "WNT5B", "WNT6", "WNT7A", "WNT7B", "WNT8A",
  "WNT8B", "WNT9A", "WNT9B", "WNT10A", "WNT10B", "WNT11", "WNT16", "FZD1", "FZD2", "FZD3",
  "FZD4", "FZD5", "FZD6", "FZD7", "FZD8", "FZD9", "FZD10", "LRP5", "LRP6", "DVL1",
  "DVL2", "DVL3", "AXIN1", "AXIN2", "GSK3B", "CTNNB1", "TCF7", "TCF7L1", "TCF7L2", "LEF1",
  "BCL2", "BCL2L1", "MCL1", "BAX", "BAK1", "BAD", "BID", "BIM", "PUMA", "NOXA",
  "CASP3", "CASP8", "CASP9", "APAF1", "CYCS", "DIABLO", "XIAP", "SURVIVIN", "C-FLIP", "TRAIL",
  "DR4", "DR5", "DCR1", "DCR2", "FADD", "TRADD", "RIPK1", "RIPK3", "MLKL", "FAS",
  "FASLG", "TNF", "TNFR1", "TNFR2", "IL1B", "IL6", "IL8", "IL10", "TGFB1", "TGFB2",
  "TGFB3", "IFNG", "IL2", "IL4", "IL12A", "IL12B", "IL17A", "IL17F", "IL23A", "IL23R",
  "STAT1", "STAT2", "STAT4", "STAT5A", "STAT5B", "STAT6", "IRF1", "IRF2", "IRF3", "IRF4",
  "IRF5", "IRF6", "IRF7", "IRF8", "IRF9", "NFKB1", "NFKB2", "RELA", "RELB", "CREL",
  "IKBKB", "IKBKG", "CHUK", "MAP3K1", "MAP3K2", "MAP3K3", "MAP3K4", "MAP3K5", "MAP3K6", "MAP3K7",
  "MAP3K8", "MAP3K9", "MAP3K10", "MAP3K11", "MAP3K12", "MAP3K13", "MAP3K14", "MAP3K15", "MAP2K1", "MAP2K2",
  "MAP2K3", "MAP2K4", "MAP2K5", "MAP2K6", "MAP2K7", "MAPK1", "MAPK3", "MAPK4", "MAPK6", "MAPK7",
  "MAPK8", "MAPK9", "MAPK10", "MAPK11", "MAPK12", "MAPK13", "MAPK14", "MAPK15", "JNK1", "JNK2",
  "JNK3", "P38", "ERK1", "ERK2", "ERK3", "ERK4", "ERK5", "ERK7", "ERK8", "RSK1",
  "RSK2", "RSK3", "RSK4", "MSK1", "MSK2", "MNK1", "MNK2", "MK2", "MK3", "MK5",
  "AKT2", "AKT3", "PDK1", "PDK2", "PDK3", "PDK4", "SGK1", "SGK2", "SGK3", "PKC",
  "PKCA", "PKCB", "PKCG", "PKCD", "PKCE", "PKCH", "PKCI", "PKCJ", "PKCL", "PKCM",
  "PKCN", "PKCO", "PKCP", "PKCQ", "PKCR", "PKCS", "PKCT", "PKCU", "PKCV", "PKCW",
  "PKCX", "PKCY", "PKCZ", "PKA", "PKB", "PKC", "PKD", "PKE", "PKF", "PKG", "PKH", "PKI", "PKJ"
)

# Extend to 1000 genes if needed
if (length(common_genes) < 1000) {
  common_genes <- c(common_genes, paste0("GENE", (length(common_genes)+1):1000))
} else {
  common_genes <- common_genes[1:1000]
}

# Create minimal example data with consistent genes
create_minimal_data <- function() {
  message("Creating minimal differential expression data with ", length(common_genes), " genes...")
  
  set.seed(123)
  example_data <- data.frame(
    symbol = common_genes,
    log2FoldChange = rnorm(length(common_genes), 0, 2),
    padj = runif(length(common_genes), 0, 1),
    pvalue = runif(length(common_genes), 0, 1),
    baseMean = runif(length(common_genes), 10, 1000),
    stringsAsFactors = FALSE
  )
  example_data$minus_log10_padj <- -log10(example_data$padj)
  
  return(example_data)
}

# Create minimal count data with consistent genes
create_minimal_count_data <- function() {
  message("Creating minimal count data with ", length(common_genes), " genes...")
  
  set.seed(456)
  count_data <- matrix(
    rnbinom(length(common_genes) * 20, mu = 100, size = 1),
    nrow = length(common_genes),
    ncol = 20,
    dimnames = list(
      common_genes,
      paste0("Sample", 1:20)
    )
  )
  
  return(count_data)
}

# Create minimal metadata with "Library" as sample_id
create_minimal_metadata <- function() {
  message("Creating minimal metadata structure with Library column...")
  
  metadata <- data.frame(
    Library = paste0("Sample", 1:20),
    sample_id = paste0("Sample", 1:20),
    condition = rep(c("Control", "Treatment"), each = 10),
    batch = rep(c("A", "B"), 10),
    Gender = sample(c("M", "F"), 20, replace = TRUE),
    Age_as_at_enrolment = sample(30:70, 20, replace = TRUE),
    Race = sample(c("I", "II", "III"), 20, replace = TRUE),
    Body_weight = runif(20, 50, 100),
    Body_height = runif(20, 1.5, 1.9),
    BMI = runif(20, 18, 35),
    Histology_Steatosis_grade = sample(0:3, 20, replace = TRUE),
    Histology_Lobular_inflammation = sample(0:3, 20, replace = TRUE),
    Histology_Ballooning = sample(0:2, 20, replace = TRUE),
    Histology_Fibrosis = sample(0:4, 20, replace = TRUE),
    RINe = runif(20, 6, 10),
    DISEASE_STAGES = sample(c("Early", "Advanced", "Severe"), 20, replace = TRUE),
    Cohort = sample(c("Cohort_A", "Cohort_B"), 20, replace = TRUE),
    Batch = sample(c("Batch1", "Batch2"), 20, replace = TRUE),
    Sex = sample(c("M", "F"), 20, replace = TRUE),
    stringsAsFactors = FALSE
  )
  
  return(metadata)
}

# Create minimal ordinal regression data with "Estimate", "P_value", "FDR" format
create_minimal_ordinal_regression <- function() {
  message("Creating minimal ordinal regression data with ", length(common_genes), " genes...")
  
  set.seed(789)
  ord_data <- data.frame(
    symbol = common_genes,
    Estimate = rnorm(length(common_genes), 0, 0.5),
    P_value = runif(length(common_genes), 0, 0.1),
    FDR = runif(length(common_genes), 0, 0.05),
    stringsAsFactors = FALSE
  )
  
  # Add coefficient and p_value columns for compatibility
  ord_data$coefficient <- ord_data$Estimate
  ord_data$p_value <- ord_data$P_value
  ord_data$fdr <- ord_data$FDR
  
  return(ord_data)
}

# =============================================================================
# Main Data Loading with Enhanced Conversion
# =============================================================================

message("Loading application data...")

# Load all data directly 
app_data <- list()

# Load core datasets using individual functions
app_data$diffgenes <- load_differential_expression(data_paths)
app_data$count_data <- load_count_data(data_paths)
app_data$metadata <- load_metadata(data_paths)
app_data$plasma_proteins <- load_plasma_proteins(data_paths)
app_data$ordinal_regression <- load_ordinal_regression(data_paths)
app_data$string_db <- load_string_db(data_paths)
app_data$annotation_data <- load_annotation_data(data_paths)

# Add literature trends configuration
app_data$literature_config <- if (!is.null(data_paths)) data_paths$literature_trends else NULL

# =============================================================================
# Ensure All Required Data Objects Exist
# =============================================================================

# Ensure we have data - with graceful fallbacks for all data types
if (is.null(app_data$diffgenes)) {
  message("Using minimal example data for differential expression")
  app_data$diffgenes <- create_minimal_data()
}

if (is.null(app_data$count_data)) {
  message("Using minimal example data for count data")
  app_data$count_data <- create_minimal_count_data()
}

if (is.null(app_data$metadata)) {
  message("Using minimal example data for metadata")
  app_data$metadata <- create_minimal_metadata()
}

if (is.null(app_data$ordinal_regression)) {
  message("Using minimal example data for ordinal regression")
  app_data$ordinal_regression <- create_minimal_ordinal_regression()
}

# Update data status
if (!is.null(app_data$diffgenes)) {
  app_data$data_status <- ifelse(exists("data_paths") && !is.null(data_paths$diffgenes_path) && file.exists(data_paths$diffgenes_path), 
                                "loaded", "minimal_example")
} else {
  app_data$data_status <- "no_data"
}

app_data$load_time <- Sys.time()

# =============================================================================
# Assign to Global Environment
# =============================================================================

# Assign ALL data objects to global environment
list2env(app_data, envir = .GlobalEnv)

# Double-check that all essential objects exist in global environment
required_objects <- c("diffgenes", "count_data", "metadata", "ordinal_regression")
for (obj in required_objects) {
  if (!exists(obj, envir = .GlobalEnv) || is.null(get(obj, envir = .GlobalEnv))) {
    message("WARNING: ", obj, " is not properly set in global environment")
    # Create emergency fallback
    switch(obj,
           "diffgenes" = assign("diffgenes", create_minimal_data(), envir = .GlobalEnv),
           "count_data" = assign("count_data", create_minimal_count_data(), envir = .GlobalEnv),
           "metadata" = assign("metadata", create_minimal_metadata(), envir = .GlobalEnv),
           "ordinal_regression" = assign("ordinal_regression", create_minimal_ordinal_regression(), envir = .GlobalEnv)
    )
    message("Created emergency fallback for: ", obj)
  }
}

# =============================================================================
# Debug Data Loading
# =============================================================================

# Check if data loaded properly
message("=== DATA LOADING VERIFICATION ===")
message("Differential expression data type: ", class(diffgenes))
message("Differential expression dimensions: ", if(!is.null(diffgenes)) paste(dim(diffgenes), collapse = " x ") else "NULL")
message("Count data type: ", class(count_data))
message("Count data dimensions: ", if(!is.null(count_data)) paste(dim(count_data), collapse = " x ") else "NULL")
message("Metadata type: ", class(metadata))
message("Metadata dimensions: ", if(!is.null(metadata)) paste(dim(metadata), collapse = " x ") else "NULL")
message("Ordinal regression type: ", class(ordinal_regression))
message("Ordinal regression dimensions: ", if(!is.null(ordinal_regression)) paste(dim(ordinal_regression), collapse = " x ") else "NULL")
message("Common genes count: ", length(common_genes))
message("==================================")

# Print loading summary
if (exists("print_data_summary")) {
  print_data_summary(app_data)
} else {
  # Fallback summary
  message("DATA LOADING SUMMARY:")
  message("Differential expression: ", ifelse(!is.null(diffgenes), nrow(diffgenes), "NOT LOADED"))
  message("Count data: ", ifelse(!is.null(count_data), "LOADED", "NOT LOADED"))
  message("Metadata: ", ifelse(!is.null(metadata), "LOADED", "NOT LOADED"))
  message("Plasma proteins: ", ifelse(!is.null(plasma_proteins), "LOADED", "NOT LOADED"))
  message("Ordinal regression: ", ifelse(!is.null(ordinal_regression), nrow(ordinal_regression), "NOT LOADED"))
  message("STRING database: ", ifelse(!is.null(string_db), "LOADED", "NOT LOADED"))
  message("Annotation data: ", length(annotation_data), " datasets")
  message("Local biomart mapping: ", ifelse(!is.null(local_mappings$biomart), "LOADED", "NOT AVAILABLE"))
  message("Local HGNC mapping: ", ifelse(!is.null(local_mappings$hgnc), "LOADED", "NOT AVAILABLE"))
}

# =============================================================================
# PubMed Literature Trends Functions
# =============================================================================

#' Search PubMed for literature trends
search_pubmed_trends <- function(gene_symbol, disease_terms, years_back = 10, email = NULL) {
  if (is.null(email)) {
    warning("Email is required for PubMed API. Please configure in data_paths.yml")
    return(NULL)
  }
  
  # Set email for PubMed API
  rentrez::set_entrez_email(email)
  
  current_year <- as.numeric(format(Sys.Date(), "%Y"))
  years <- (current_year - years_back):current_year
  
  trends_data <- data.frame(
    year = years,
    count = 0,
    gene = gene_symbol,
    stringsAsFactors = FALSE
  )
  
  # Build search query
  gene_query <- paste0(gene_symbol, "[Title/Abstract]")
  
  if (length(disease_terms) == 1) {
    disease_query <- paste0(disease_terms, "[Title/Abstract]")
  } else {
    # Combine multiple terms with AND
    disease_query <- paste0("(", paste(disease_terms, collapse = "[Title/Abstract] AND "), "[Title/Abstract])")
  }
  
  full_query <- paste(gene_query, "AND", disease_query)
  
  message("Searching PubMed for: ", full_query)
  
  for (i in seq_along(years)) {
    year <- years[i]
    
    tryCatch({
      # Search PubMed for this year
      search_result <- rentrez::entrez_search(
        db = "pubmed",
        term = paste(full_query, "AND (", year, "[PDAT])"),
        retmax = 0  # We only need the count
      )
      
      trends_data$count[i] <- search_result$count
      
      # Be polite to PubMed API
      Sys.sleep(0.3)
      
    }, error = function(e) {
      warning("Error searching PubMed for year ", year, ": ", e$message)
      trends_data$count[i] <- NA
    })
  }
  
  return(trends_data)
}

#' Get literature trends for multiple genes
get_literature_trends <- function(gene_symbols, disease_terms, years_back = 10, email = NULL) {
  if (is.null(gene_symbols) || length(gene_symbols) == 0) {
    return(NULL)
  }
  
  all_trends <- list()
  
  for (gene in gene_symbols) {
    message("Getting literature trends for: ", gene)
    gene_trends <- search_pubmed_trends(gene, disease_terms, years_back, email)
    
    if (!is.null(gene_trends)) {
      all_trends[[gene]] <- gene_trends
    }
    
    # Add delay between gene searches
    Sys.sleep(1)
  }
  
  if (length(all_trends) == 0) {
    return(NULL)
  }
  
  # Combine all trends
  combined_trends <- do.call(rbind, all_trends)
  return(combined_trends)
}

#' Create literature trends plot
create_literature_trends_plot <- function(trends_data, title = "Literature Trends") {
  if (is.null(trends_data) || nrow(trends_data) == 0) {
    return(create_empty_plot("No literature trends data available"))
  }
  
  p <- ggplot(trends_data, aes(x = year, y = count, color = gene, group = gene)) +
    geom_line(size = 1.2) +
    geom_point(size = 2) +
    theme_ebtat() +
    labs(
      title = title,
      x = "Year",
      y = "Number of Publications",
      color = "Gene"
    ) +
    scale_x_continuous(breaks = unique(trends_data$year)) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
  
  return(p)
}

# =============================================================================
# Application Options and Configuration
# =============================================================================

# Set application options
options(
  shiny.maxRequestSize = 10 * 1024^2,      # 10MB file upload limit
  shiny.sanitize.errors = FALSE,           # Show actual errors for debugging
  stringsAsFactors = FALSE,                # Don't automatically convert strings to factors
  warn = 1,                                # Show warnings immediately
  scipen = 999                             # Prefer fixed notation over scientific
)

# =============================================================================
# Color Palettes and Theme Configuration
# =============================================================================

# Color palettes
Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

# Application-specific color schemes
app_colors <- list(
  primary = "#2196F3",
  secondary = "#FF5722",
  success = "#4CAF50",
  warning = "#FFC107",
  danger = "#F44336",
  info = "#00BCD4",
  light = "#F5F5F5",
  dark = "#212121"
)

# Custom ggplot themes
theme_ebtat <- function(base_size = 14, base_family = "") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = base_size * 1.2),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.major = element_line(color = "grey80", size = 0.2),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey50", fill = NA, size = 0.5)
    )
}

# =============================================================================
# Application Constants and Defaults
# =============================================================================

# Application metadata
app_metadata <- list(
  name = "ENABL Biomarker / Target Analysis Tool (EB-TAT)",
  version = "2.0.0",
  release_date = "2024",
  author = "Vijay Singh (GIS)",
  description = "A Shiny app for analysis of transcriptomics data",
  maintainer = "ENABL Bioinformatics Team"
)

# =============================================================================
# Validation Functions
# =============================================================================

# Function to validate data structure
validate_data_structure <- function(data) {
  if (is.null(data)) return(FALSE)
  
  # Check for essential columns in differential expression data
  if ("log2FoldChange" %in% colnames(data) && "padj" %in% colnames(data)) {
    return(TRUE)
  }
  
  # Alternative column names
  essential_cols <- c("symbol", "gene", "logFC", "pvalue", "adj.P.Val")
  if (any(essential_cols %in% colnames(data))) {
    return(TRUE)
  }
  
  warning("Data missing required columns for analysis")
  return(FALSE)
}

# =============================================================================
# Final Initialization
# =============================================================================

message("=== EB-TAT Initialization Complete ===")
message("Time: ", Sys.time())
message("Data status: ", app_data$data_status)
message("Differential expression rows: ", ifelse(!is.null(diffgenes), nrow(diffgenes), "NOT LOADED"))
message("Count data loaded: ", !is.null(count_data))
message("Metadata loaded: ", !is.null(metadata))
message("Plasma proteins loaded: ", !is.null(plasma_proteins))
message("Ordinal regression loaded: ", !is.null(ordinal_regression))
message("STRING database loaded: ", !is.null(string_db))
message("Annotation data loaded: ", length(annotation_data))
message("PubMed integration: ", !is.null(literature_config))
message("Common genes defined: ", length(common_genes))
message("======================================")

# Cleanup temporary variables
rm(cran_packages, bioc_packages, cran_loaded, bioc_loaded, utility_files, utils_loaded)

# Create required directories if they don't exist
required_dirs <- c("www", "modules", "data", "config")
for (dir in required_dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    message("Created directory: ", dir)
  }
}

# =============================================================================
# Emergency Fallback Functions
# =============================================================================

# Emergency data provider if everything else fails
.get_safe_data <- function() {
  if (exists("diffgenes") && validate_data_structure(diffgenes)) {
    return(diffgenes)
  } else {
    return(create_minimal_data())
  }
}

# Emergency plot function
.safe_plot <- function() {
  data <- .get_safe_data()
  ggplot(data, aes(x = log2FoldChange, y = minus_log10_padj)) +
    geom_point(alpha = 0.6, color = app_colors$primary) +
    theme_ebtat() +
    labs(title = "Volcano Plot",
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-value")
}

# =============================================================================
# Global Error Handler
# =============================================================================

# Set up global error handling
options(shiny.error = function() {
  message("Global error caught at: ", Sys.time())
  # Don't break the app, just log the error
})

# =============================================================================
# Export Global Objects
# =============================================================================

# Make essential objects available globally
.GlobalEnv$Okabe_Ito <- Okabe_Ito
.GlobalEnv$app_colors <- app_colors
.GlobalEnv$theme_ebtat <- theme_ebtat
.GlobalEnv$app_metadata <- app_metadata
.GlobalEnv$common_genes <- common_genes

message("EB-TAT Global Configuration Complete!")