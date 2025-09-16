library(TCGAbiolinks)
library(SummarizedExperiment)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop("Usage: Rscript download_cohort.R <TCGA cohort> (e.g. TCGA-PAAD)")
}

# The TCGA project (e.g. TCGA-BRCA)
tcga_project <- args[1]

# Query the GDC
query <- GDCquery(
  project = tcga_project,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# Download the data
GDCdownload(query)

# Prepare the data to create a SummarizedExperiment object
data <- GDCprepare(query)

# Access raw counts
raw_counts <- assay(data)        # genes x samples matrix
metadata <- colData(data)        # sample metadata

# remove the suffix, e.g. ENSG00000123456.11 -> ENSG00000123456
# so we don't have to handle it later when mapping to gene symbols.
rn <- sub("\\..*$", "", rownames(raw_counts))
not_duplicated <- !duplicated(rn)
raw_counts <- raw_counts[not_duplicated,]
rownames(raw_counts) <- rn[not_duplicated]

write.table(raw_counts, sprintf('%s.raw_counts.tsv', tcga_project), sep='\t', quote=F)
save(metadata, file=sprintf('%s.rnaseq_metadata.rds', tcga_project))
