library(TCGAbiolinks)
library(SummarizedExperiment)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop("Usage: Rscript download_cohort.R <TCGA cohort> (e.g. TCGA-PAAD)")
}

# The TCGA project (e.g. TCGA-BRCA)
tcga_project <- args[1]

# Query the GDC - all sample types (which can include normals)
query <- GDCquery(
  project = tcga_project,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# A query for normal samples. Note that this assumes we are making
# queries on cancer types where this makes sense. TCGA has some other
# normal types (typically blood normals for somatic mutation work), but 
# here we are focusing on solid tumors
normal_query <- GDCquery(
  project = tcga_project,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type='Solid Tissue Normal'
)

# Download the data
GDCdownload(query)
GDCdownload(normal_query)

# Prepare the data to create a SummarizedExperiment object
data <- GDCprepare(query)
normal_data <- GDCprepare(normal_query)

# Access raw counts for all samples
raw_counts <- assay(data)        # genes x samples matrix
metadata <- colData(data)        # sample metadata

# get the metadata for the normal samples
normal_metadata <- colData(normal_data)
selected_samples <- setdiff(rownames(metadata), rownames(normal_metadata))

# now keep only those non-normal:
raw_counts <- raw_counts[, selected_samples]
metadata <- metadata[selected_samples,]

# remove the suffix, e.g. ENSG00000123456.11 -> ENSG00000123456
# so we don't have to handle it later when mapping to gene symbols.
rn <- sub("\\..*$", "", rownames(raw_counts))
not_duplicated <- !duplicated(rn)
raw_counts <- raw_counts[not_duplicated,]
rownames(raw_counts) <- rn[not_duplicated]

write.table(raw_counts, sprintf('%s.raw_counts.tsv', tcga_project), sep='\t', quote=F)
save(metadata, file=sprintf('%s.rnaseq_metadata.rds', tcga_project))
