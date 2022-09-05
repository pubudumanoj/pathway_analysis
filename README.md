[![](https://img.shields.io/badge/Shiny-shinyapps.io-blue?style=flat&labelColor=white&logo=RStudio&logoColor=blue)](https://pubudumanoj.shinyapps.io/pathway_analysis/)
[![Github All Releases](https://img.shields.io/github/downloads/pubudumanoj/pathway_analysis/total.svg)]()

# Pathway Enrichment Analysis
This shiny app performs a pathway enrichment analysis using KEGG and Reactom databases utilizing the geneset provided by you. Currently it only supports human genes.

## Usage

The application is hosted on [Shiny.io](https://pubudumanoj.shinyapps.io/pathway_analysis/)

### Instructions to use

Insert a gene set file (text file with gene names seperated by a newline or paste gene set in the text area. Then change the options for pathway analysis or final figure.

## Installation

You can also deploy the application in your local machine using R Studio and/or R.

### Prerequisites

Following R packages are required to run this application
```
install.packages("shiny")
install.packages("data.table")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("DT")
install.packages("ggrepel")
install.packages("shinycssloaders")
if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("MAGeCKFlute")
BiocManager::install("ReactomePA")
BiocManager::install("clusterProfiler")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")

```

Clone the repo and run the application.


