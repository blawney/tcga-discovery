# This script does the following:
# - extracts the beta matrix from the methylation RDS input
# - re-maps the sample IDs for the RNA-seq and beta matrix so that we have matching samples
# - filters the beta matrix for probes that are in the passed 'probe file' (e.g. distal enhancers, promoters, etc.)

library(SummarizedExperiment)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop("Usage: Rscript align_and_prep.R <TCGA cohort> <meth RDS> <rna-seq counts> <probe file>")
}

# The TCGA project (e.g. TCGA-BRCA)
tcga_project <- args[1]

# the methylation RDS file which is the output of the TCGAbiolinks::GDCprepare function
methylation_data <- args[2]

# a normalized RNA-seq count matrix.
rnaseq_data <- args[3]

# A BED6-format file giving the probe locations.
probe_file <- args[4]

load(methylation_data)
meth_mtx <- assay(data)

counts <- read.table(rnaseq_data, sep='\t', row.names=1, header=T, check.names=F)

meth_cols <- colnames(meth_mtx)
rna_cols <- colnames(counts)

meth_cols <- sapply(strsplit(meth_cols, "-"), function(parts) paste(parts[1:4], collapse = "-"))
rna_cols <- sapply(strsplit(rna_cols, "-"), function(parts) paste(parts[1:4], collapse = "-"))

# there can be multiples after we transform the column names. 
meth_dupes <- duplicated(meth_cols)
rna_dupes <- duplicated(rna_cols)


if(length(unique(meth_cols)) != dim(meth_mtx)[2]){
    print(sprintf('In %s, methylation columns had repeated samples', tcga_project))
    meth_dupes <- duplicated(meth_cols) # array of bools
    meth_cols <- meth_cols[!meth_dupes]
    meth_mtx <- meth_mtx[,!meth_dupes]
}

if(length(unique(rna_cols)) != dim(counts)[2]){
    print(sprintf('In %s, rnaseq columns had repeated samples', tcga_project))
    rna_dupes <- duplicated(rna_cols)
    rna_cols <- rna_cols[!rna_dupes]
    counts <- counts[,!rna_dupes]
}

colnames(counts) <- rna_cols
colnames(meth_mtx) <- meth_cols

intersecting_samples <- intersect(meth_cols, rna_cols)

counts <- counts[, intersecting_samples]
meth_mtx <- meth_mtx[, intersecting_samples]

probes_df <- read.table(probe_file, sep='\t', header=F, col.names=c('chrom','start','end','probe_id','score', 'strand'))

meth_mtx <- meth_mtx[probes_df[, 'probe_id'], ]

# Use the probe file to name the output files
probe_prefix <- strsplit(basename(probe_file), '\\.')[[1]][1]

# write out the subset version of the normalized counts and the betas
write.table(counts, sprintf('%s.norm_counts.reformatted_sample_id.tsv', tcga_project), sep='\t', quote=F)
write.table(meth_mtx, file=sprintf('%s.betas.%s.reformatted_sample_id.tsv', tcga_project, probe_prefix), sep='\t', quote=FALSE, na='')
