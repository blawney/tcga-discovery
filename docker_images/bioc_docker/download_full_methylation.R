library(TCGAbiolinks)
library(SummarizedExperiment)
library(biomaRt)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop("Usage: Rscript download_cohort.R <TCGA cohort> (e.g. TCGA-PAAD)")
}

# The TCGA project (e.g. TCGA-BRCA)
tcga_project <- args[1]

# Query the GDC
query <- GDCquery(
  project = tcga_project,
  data.category = "DNA Methylation",
  data.type = "Methylation Beta Value",
  platform = "Illumina Human Methylation 450",
  workflow.type = "SeSAMe Methylation Beta Estimation"
)

# Download the data
GDCdownload(query, files.per.chunk=10)

# Prepare the data to create a RangedSummarizedExperiment object
data <- GDCprepare(query)

save(data, file=sprintf('%s.methylation_data.rds', tcga_project))