#!/bin/bash

# ENABL Biomarker/Target Analysis Tool (EB-TAT) Installation Script
# This script installs required R packages and sets up the application

echo "=============================================="
echo "ENABL Biomarker/Target Analysis Tool (EB-TAT)"
echo "Installation Script"
echo "=============================================="

# Check if R is installed
if ! command -v R &> /dev/null; then
    echo "ERROR: R is not installed. Please install R first."
    echo "Visit https://cran.r-project.org/"
    exit 1
fi

# Check if Rscript is available
if ! command -v Rscript &> /dev/null; then
    echo "ERROR: Rscript is not available. Please install R properly."
    exit 1
fi

echo "R version:"
R --version | head -1
echo ""

# Function to install packages
install_packages() {
    echo "Installing R packages..."
    
    # Install CRAN packages
    echo "Installing CRAN packages..."
    Rscript -e 'install.packages(c("shiny", "ggplot2", "dplyr", "ggrepel", "shinycssloaders", "readxl", "DT", "RCurl", "plot.matrix", "psych", "reshape2", "STRINGdb", "ggVennDiagram", "ggvenn", "readr", "clusterProfiler", "enrichplot", "DOSE", "ggupset", "UpSetR", "DESeq2", "limma", "pheatmap", "hash", "gridExtra", "shinythemes"), repos="https://cloud.r-project.org/")'
    
    # Install Bioconductor packages
    echo "Installing Bioconductor packages..."
    Rscript -e 'if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install(c("org.Hs.eg.db", "HPAanalyze"))'
    
    echo "Package installation completed."
}

# Function to check package installation
check_packages() {
    echo "Checking package installation..."
    
    # List of required packages
    packages=("shiny" "ggplot2" "dplyr" "ggrepel" "shinycssloaders" "readxl" "DT" "RCurl" 
              "plot.matrix" "psych" "reshape2" "STRINGdb" "ggVennDiagram" "ggvenn" "readr" 
              "clusterProfiler" "enrichplot" "DOSE" "ggupset" "UpSetR" "DESeq2" "limma" 
              "pheatmap" "hash" "gridExtra" "shinythemes" "org.Hs.eg.db" "HPAanalyze")
    
    missing_packages=()
    
    for pkg in "${packages[@]}"; do
        if Rscript -e "if(!require($pkg, quietly=TRUE, character.only=TRUE)) quit(status=1)" &> /dev/null; then
            echo "✓ $pkg"
        else
            echo "✗ $pkg"
            missing_packages+=("$pkg")
        fi
    done
    
    if [ ${#missing_packages[@]} -eq 0 ]; then
        echo "All required packages are installed."
        return 0
    else
        echo "Missing packages: ${missing_packages[*]}"
        return 1
    fi
}

# Function to create directory structure
create_directories() {
    echo "Creating directory structure..."
    
    # Create modules directory
    mkdir -p modules
    
    # Create utils directory
    mkdir -p utils
    
    # Create www directory
    mkdir -p www
    
    echo "Directory structure created."
}

# Function to create data loading configuration
create_data_config() {
    echo "Creating data loading configuration..."
    
    cat > utils/data_loading.R << 'EOL'
# Function to load all application data
load_app_data <- function() {
  # This function should be modified to point to your actual data files
  # Replace the paths below with the actual paths to your data files
  
  # Example structure:
  # diffgenes <- read.csv("path/to/your/differential_expression_data.csv")
  # count.data <- readRDS("path/to/your/count_data.RDS")
  # meta.data <- readRDS("path/to/your/metadata.RDS")
  
  # Return a list with all data objects
  return(list(
    # diffgenes = diffgenes,
    # count.data = count.data,
    # meta.data = meta.data,
    # Add other data objects as needed
  ))
}
EOL

    echo "Data configuration created. Please edit utils/data_loading.R with your actual data paths."
}

# Main installation process
main() {
    echo "Starting installation..."
    
    # Install packages
    install_packages
    
    # Check if packages were installed correctly
    if ! check_packages; then
        echo "Some packages failed to install. Please install them manually."
        read -p "Do you want to try installing missing packages manually? (y/n): " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            install_packages
        fi
    fi
    
    # Create directory structure
    create_directories
    
    # Create data loading configuration
    create_data_config
    
    echo ""
    echo "=============================================="
    echo "Installation completed!"
    echo "=============================================="
    echo ""
    echo "Next steps:"
    echo "1. Edit utils/data_loading.R to point to your data files"
    echo "2. Run the application with: R -e \"shiny::runApp()\""
    echo "3. Or open app.R in RStudio and click 'Run App'"
    echo ""
    echo "For detailed instructions, see README.md"
}

# Run main function
main