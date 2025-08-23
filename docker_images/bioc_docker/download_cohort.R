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

keep_cols <- c("barcode", "patient", "sample", "sample_submitter_id", 
               "sample_type_id", "tumor_descriptor", "sample_id", "sample_type", 
               "state", "specimen_type", "is_ffpe", "tissue_type", 
               "synchronous_malignancy", "ajcc_pathologic_stage", 
               "tissue_or_organ_of_origin", "days_to_last_follow_up", "age_at_diagnosis", 
               "primary_diagnosis", "prior_treatment", "diagnosis_is_primary_disease", 
               "ajcc_staging_system_edition", "ajcc_pathologic_t", "morphology", "ajcc_pathologic_n", 
               "ajcc_pathologic_m", "residual_disease", "classification_of_tumor", "diagnosis_id", 
               "icd_10_code", "site_of_resection_or_biopsy", "tumor_grade", "tumor_of_origin", 
               "race", "gender", "ethnicity", "vital_status", "age_at_index", 
               "days_to_birth", "demographic_id", "days_to_death", 
               "bcr_patient_barcode", "project_id")
write.table(metadata[, keep_cols], sprintf('%s.metadata.tsv', tcga_project), sep='\t', quote=F)
