library(TCGAbiolinks)
library(SummarizedExperiment)
library(GenomicRanges)
library(biomaRt)
library(dplyr)

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
  workflow.type = "SeSAMe Methylation Beta Estimation"
)

# Download the data
GDCdownload(query, files.per.chunk=10)

# Prepare the data to create a RangedSummarizedExperiment object
data <- GDCprepare(query)

# Extract out the parts we want
beta_vals <- assay(data)        # probes x samples matrix
metadata <- colData(data)        # sample metadata
probe_annot <- rowData(data)

# Use Ensembl Biomart to get the gene locations
ensembl <- useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl")
genes <- getBM(attributes = c("chromosome_name","start_position","end_position","strand","external_gene_name"),
               mart = ensembl)
genes <- genes %>% filter(external_gene_name != "")

# Create a GRanges object for the gene starts - below, we apply the `GRanges::promoters` function to get
# the upstream region
genes_gr <- GRanges(
    # The probe annotations have a 'chr' prefix
    seqnames = paste0("chr", genes$chromosome_name),
    ranges = IRanges(
        start = ifelse(genes$strand == 1, genes$start_position, genes$end_position),
        width=1
    ),
    strand = ifelse(genes$strand==1, "+", "-"),
    gene_name = genes$external_gene_name
)
genes_gr <- promoters(genes_gr, upstream=2000, downstream=100)


na_chroms <- is.na(probe_annot$chrm_A)
probe_annot <- probe_annot[!na_chroms,]
probe_ranges <- GRanges(
  seqnames = probe_annot$chrm_A,
  ranges = IRanges(start=probe_annot$probeBeg, end=probe_annot$probeEnd),
  probe_id = rownames(probe_annot)
)

# This looks like:
# Hits object with 212592 hits and 0 metadata columns:
#            queryHits subjectHits
#            <integer>   <integer>
#        [1]         2       55840
#        [2]         3       55841
#        [3]         4       55841
#        [4]         5       55841
#        [5]         8       55845
overlap_hits <- findOverlaps(probe_ranges, genes_gr)
probes_in_promoters <- as.data.frame(probe_ranges[queryHits(overlap_hits)]) %>% rename(probe_chrom=seqnames, probe_start=start, probe_end=end)
overlapping_promoters <- as.data.frame(genes_gr[subjectHits(overlap_hits)]) %>% rename(promoter_chrom=seqnames, promoter_start=start, promoter_end=end)

# this is now a dataframe which contains probes that overlap putative promoters (by our upstream/downstream definition)
# Note that the probe IDs are not unique - i.e. the same probe can fall into multiple promoter regions. Conversely, we
# can obviously have multiple probes that fall into a single gene's promoter region.
promoter_hits <- cbind(probes_in_promoters, overlapping_promoters)
write.table(promoter_hits, sprintf('%s.promoter_probe_meta.tsv', tcga_project), sep='\t', quote=F)

unique_promoter_probes <- unique(promoter_hits$probe_id)
beta_vals <- beta_vals[unique_promoter_probes,]
write.table(beta_vals, sprintf('%s.promoter_betas.tsv', tcga_project), sep='\t', quote=F, na='')

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
