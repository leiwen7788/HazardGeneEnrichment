#HazardGeneEnrichment

#Introduction

HazardGeneEnrichment is an R package designed to perform Hazard Ratio-based Gene Set Enrichment Analysis (GSEA) by integrating single-cell RNA sequencing (scRNA-seq) data with Cancer Genome Atlas (TCGA) survival data. This package provides a suite of functions to identify single-cell cluster-specific marker genes and assess the association of these gene sets with prognostic risk in TCGA cancer cohorts.
Objectives
To identify cell cluster-specific marker genes from single-cell RNA sequencing data.
To rank genes in TCGA cancer cohorts based on their Hazard Ratios.
To perform GSEA on single-cell marker gene sets to explore their association with cancer patient survival prognosis.
To visualize GSEA results, including enrichment plots and key statistical data.
#Installation
You can install the HazardGeneEnrichment package from GitHub using devtools:

If devtools is not already installed
install.packages("devtools")

devtools::install_github("leiwen7788/HazardGeneEnrichment") 

#Dependencies

The HazardGeneEnrichment package relies on the following CRAN and Bioconductor packages. devtools should automatically install these dependencies when you install HazardGeneEnrichment.
tidyverse (for data manipulation, including dplyr and ggplot2)
clusterProfiler (for core GSEA functionalities)
gridExtra (for plot arrangement)
Package Functions Overview
The HazardGeneEnrichment package primarily offers the following core functionalities:

display_genesets(): Displays a list of TCGA cancer types supported for analysis by this package.
HazardEnrichment(): The main function to execute Hazard Ratio-based GSEA.
gseaplot_hazard(): Plots GSEA enrichment curves, with an option to include a p-value table.

#Detailed Usage Workflow
This section will detail how to use the HazardGeneEnrichment package and reproduce the analysis workflow you've demonstrated in your code.
1. Load Necessary Libraries

# Load your developed package
library(HazardGeneEnrichment)

# Load other dependency packages
library(tidyverse)
library(clusterProfiler)
library(gridExtra) # For combining plots
2. Prepare Input Data
You will need to prepare a marker gene table from single-cell RNA sequencing data, along with gene expression and clinical survival data from a TCGA cancer cohort.
2.1 Load Marker Gene Table
The marker table contains differentially expressed genes identified from single-cell data. We assume it has been pre-processed to include columns like cluster and avg_log2FC.
code
R
# Load the marker gene table from the specified path
# This file is assumed to contain marker gene information for various single-cell clusters
marker <- read.table('all.diffgene_onlypos_padj005.txt',
                     header = TRUE, # Adjust based on your file's actual header status
                     sep = "\t")    # Adjust based on your file's actual separator
                     
2.2 Load Single-Cell Gene List
This is a list of all genes present in your single-cell dataset, which will be used as the background gene set for GSEA.

# Extract all gene names from the Seurat object
seuratobj <- readRDS('seuratObject.rds')
singlecell_gene <- rownames(seuratobj)

# Ensure singlecell_gene is a character vector
if (!is.character(singlecell_gene)) {
  singlecell_gene <- as.character(singlecell_gene)
}
2.3 Filter Top 50 Marker Genes per Cluster
This step prepares the gene sets for GSEA enrichment analysis. We select the top 50 genes with the highest avg_log2FC for each cell cluster.
code
R
top50_markers <- marker %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(n = 50, order_by = avg_log2FC)

# View a portion of the filtered marker genes
head(top50_markers)
2.4 Format Gene Set Data (TERM2GENE)
GSEA tools like clusterProfiler typically require gene set data in a specific format: a data frame with two columns, term and gene. The term column represents the gene set name (here, the cell cluster name), and the gene column represents the gene ID.
code
R
TERM2GENE <- top50_markers[, c('cluster', 'gene')]
colnames(TERM2GENE) <- c('term', 'gene')
TERM2GENE <- as.data.frame(TERM2GENE)

# Ensure the gene column is of character type to prevent subsequent matching issues
TERM2GENE$gene <- as.character(TERM2GENE$gene)

# View a portion of the formatted gene set data
head(TERM2GENE)
3. View Available TCGA Cancer Cohorts
The display_genesets() function shows a list of TCGA cancer types that HazardGeneEnrichment package supports for analysis.
code
R
# List all supported TCGA cancer cohorts
display_genesets()

# Example output:
# [1] "TCGA-ACC"     "TCGA-BLCA"    "TCGA-BRCA"    "TCGA-CESC"    ... "TCGA-UVM"     "NB(GSE85047)"
These are the values you can pass to the TCGA_cancer_type argument of the HazardEnrichment function.
4. Perform Hazard Ratio-based GSEA
Now, we can use the HazardEnrichment() function to perform GSEA for a specific TCGA cancer cohort. In this example, we select TCGA-HNSC (Head and Neck Squamous Cell Carcinoma).
code
R
# Execute Hazard Ratio-based GSEA
# TCGA_cancer_type: Specifies the TCGA cancer type to analyze, e.g., "TCGA-HNSC"
# TERM2GENE: The formatted gene set data frame
# pvalueCutoff: The p-value threshold for GSEA
# singlecell_gene: The background list of all single-cell genes, used for filtering
my_gsea_results <- HazardEnrichment(TCGA_cancer_type = 'TCGA-HNSC',
                                    TERM2GENE = TERM2GENE,
                                    pvalueCutoff = 0.05,
                                    singlecell_gene = singlecell_gene)

# View summary of GSEA results
summary(my_gsea_results)

# View the gene set IDs (cluster names) from the GSEA results
unique(my_gsea_results@result$ID)
5. Visualize GSEA Results
Use the gseaplot_hazard() function to draw the GSEA enrichment curves. You can choose to plot all significantly enriched gene sets or specify particular gene set IDs.
code
R
# Plot all significantly enriched gene sets
# geneSetID = unique(my_gsea_results@result$ID) will plot all gene sets present in the results
# pvalue_table = TRUE will include a table of p-values and q-values in the plot
gseaplot_hazard(my_gsea_results,
                geneSetID = unique(my_gsea_results@result$ID),
                pvalue_table = TRUE)

# If you want to plot a single gene set, e.g., the first significantly enriched gene set
# gseaplot_hazard(my_gsea_results,
#                 geneSetID = my_gsea_results@result$ID[1],
#                 pvalue_table = TRUE)
Interpretation of Results
The GSEA enrichment plot (output of gseaplot_hazard) will show the distribution of each single-cell cluster (as a gene set) within the gene list ranked by Hazard Ratio.
Peak of the Enrichment Curve:
If the peak appears at the front of the ranked gene list (Running Enrichment Score > 0), it indicates that the marker genes of that cell cluster are associated with a high Hazard Ratio (i.e., worse prognosis).
If the peak appears at the back of the ranked gene list (Running Enrichment Score < 0), it indicates that the marker genes of that cell cluster are associated with a low Hazard Ratio (i.e., better prognosis).
P-value and q-value (FDR): These values (typically shown in the legend or a table) are used to determine the statistical significance of the enrichment. A qvalue < 0.25 or qvalue < 0.05 is usually considered significant enrichment.
Through this analysis, you can identify which single-cell subpopulations' molecular characteristics (represented by their marker genes) are significantly associated with the survival prognosis of patients in specific cancer types.
Contributions
Feel free to submit Issues and Pull Requests!
License
GPL-3
