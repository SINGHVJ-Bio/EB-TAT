# ENABL Biomarker/Target Analysis Tool (EB-TAT)

A comprehensive Shiny application for analyzing transcriptomics data, specifically designed for biomarker and target identification in liver disease research.

![EB-TAT](https://img.shields.io/badge/Version-2.0.0-blue.svg)
![Shiny](https://img.shields.io/badge/Built%20with-Shiny-green.svg)
![R](https://img.shields.io/badge/R-%3E%3D%204.0-blue.svg)

## 📋 Table of Contents

- [Features](#-features)
- [Installation](#-installation)
- [Configuration](#-configuration)
- [Usage](#-usage)
- [Modules](#-modules)
- [Data Requirements](#-data-requirements)
- [Troubleshooting](#-troubleshooting)
- [Support](#-support)
- [License](#-license)

## ✨ Features

- **Interactive Volcano Plots**: Visualize differential expression with customizable thresholds and annotations
- **Comprehensive Analytics**: Gene set enrichment analysis, tissue expression, and cellular compartment visualization
- **Protein-Protein Interaction Networks**: Analyze interactions using STRING database with adjustable confidence scores
- **Publication Trend Analysis**: Track research interest in genes over time
- **Advanced Filtering**: Filter genes based on characteristics like secretion, transmembrane domains, tissue specificity
- **Multiple Data Input Options**: Load data from files or use built-in example data
- **Export Functionality**: Download plots and results in multiple formats (PDF, PNG, SVG, CSV, Excel)
- **Configurable Data Paths**: All data sources configured through YAML configuration file

## 🚀 Installation

### Prerequisites

- **R** (version 4.0 or higher)
- **RStudio** (recommended for development)
- **System dependencies** for R packages

### Automated Installation

1. **Clone or download** the application files
2. **Run the installation script**:
   ```bash
   chmod +x install.sh
   ./install.sh

3. **Edit the configuration file located at config data_paths.yml**
# Example configuration
diffgenes_path: "/path/to/your/differential_expression_data.tsv"
count_data_path: "/path/to/your/count_data.RDS"
meta_data_path: "/path/to/your/metadata.RDS"
# ... other paths

Configuration Options
The application uses a YAML configuration file (config/data_paths.yml) with the following options:

Configuration Key	Description	Required
diffgenes_path	Differential expression data (TSV/CSV)	✅
count_data_path	Count data (RDS)	❌
meta_data_path	Metadata (RDS)	❌
ord_reg_fib_path	Ordinal regression - fibrosis	❌
ord_reg_inf_path	Ordinal regression - inflammation	❌
lit_trends_f_path	Literature trends - fibrosis	❌
lit_trends_i_path	Literature trends - inflammation	❌
feature_data_path	Feature data (Excel)	❌
hgnc_data_path	HGNC gene mapping	❌
string_db_base_path	STRING database directory	❌
...	[See full list]	

4. **Run the application**:
# Method 1: From RStudio
# Open app.R and click "Run App"

# Method 2: From R console
shiny::runApp()

# Method 3: From command line
R -e "shiny::runApp()"


### File Structure ###
#### ############# ###
(base) singhvj@vijays-MacBook-Pro:EB-TAT % tree
.
├── app.R
├── check_encoding.R
├── config
│   └── data_paths.yml
├── data
│   ├── annotation_data
│   │   ├── hepatocytes.tsv
│   │   ├── hsc_cells.csv
│   │   ├── human_tf.csv
│   │   ├── innate_db_genes.xls
│   │   ├── kupffer_cells.tsv
│   │   ├── liver_expression.tsv
│   │   ├── receptor_ligand_pair.txt
│   │   ├── secreted_proteins.tsv
│   │   └── sp_tm_data.csv
│   ├── config
│   ├── count_data
│   │   ├── count_matrix.csv
│   │   └── count_matrix.RDS
│   ├── differential_expression
│   │   └── diff_expression_data.tsv
│   ├── hgnc_data
│   │   └── hgnc_mapping.txt
│   ├── literature_trends
│   │   ├── fibrosis_literature_trends.rds
│   │   └── inflammation_literature_trends.rds
│   ├── metadata
│   │   ├── metadata.csv
│   │   └── metadata.RDS
│   ├── ordinal_regression
│   │   └── fibrosis_ordinal_regression.csv
│   ├── plasma_proteins
│   │   └── blood_plasma_proteins.xlsx
│   └── string_db
│       ├── string_db200.rds
│       ├── string_db300.rds
│       ├── string_db400.rds
│       ├── string_db500.rds
│       ├── string_db600.rds
│       ├── string_db700.rds
│       ├── string_db800.rds
│       └── string_db900.rds
├── global.R
├── install.sh
├── Module_Description.pdf
├── modules
│   ├── about.R
│   ├── analytics.R
│   ├── data_input.R
│   ├── enrichment.R
│   ├── expression_plot.R
│   ├── literature_trends.R
│   └── volcano_plot.R
├── project_structure.txt
├── README.md
├── themes.R
├── utils
│   ├── data_loading.R
│   ├── helpers.R
│   └── plotting.R
└── www

16 directories, 46 files
(base) singhvj@vijays-MacBook-Pro:EB-TAT % 
