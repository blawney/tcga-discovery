suppressMessages(suppressWarnings(
    library("DESeq2", character.only=T, warn.conflicts = F, quietly = T)))


###########################################################################

# args from command line:
args<-commandArgs(TRUE)

# the raw counts for the TCGA cohort, in TSV format
RAW_COUNT_MATRIX <- args[1]

# the gene we are using to segment the samples
SELECTED_GENE <- args[2]

# the quantile for the lower and upper cutoffs
LOW_Q <- as.numeric(args[3])
HIGH_Q <- as.numeric(args[4])

# the output file for differential expression results
DGE_OUTPUT <- args[5]

# normalized counts for the selected samples
# (those at the top and bottom of the selected gene's dist.)
NORM_COUNTS_OUTPUT <- args[6]

# An annotation file giving the samples that are high or low
# on the distribution, just for ease
ANN_OUTPUT <- args[7]

###########################################################################

# Don't mangle the column names quite yet. Keep them as-is and then we can convert them after we
# create a column name mapping
counts <- read.table(RAW_COUNT_MATRIX, 
                    sep='\t', 
                    stringsAsFactors=F, 
                    row.names=1, 
                    header=T, 
                    check.names=F)
orig_colnames = colnames(counts)

# cast the column names to "proper" names for R. This way we don't run into trouble with
# any called functions
new_colnames = make.names(orig_colnames)
colnames(counts) = new_colnames

# Create a map of the column names so we can re-assign them at the end
colname_mapping = setNames(orig_colnames, new_colnames)

# create a "dummy" annotation matrix which we will later fill in
ann_df <- data.frame(expression_state=rep(NA, length(new_colnames)))
rownames(ann_df) <- new_colnames

# Now run the normalization so we can find the "extremes"
# for the gene of interest:
deseq_sf = estimateSizeFactorsForMatrix(counts)
norm_mtx = sweep(counts,2,deseq_sf,'/')

# extract the expressions for the gene of interest and find
# the cutoffs for the desired quantiles
gene_array <- as.matrix(norm_mtx[SELECTED_GENE,])
cutoffs <- unname(
                  quantile(
                           gene_array, 
                           na.rm=T, 
                           probs=c(LOW_Q, HIGH_Q)
                  )
           )
low_cutoff <- cutoffs[1]
high_cutoff <- cutoffs[2]

# now get the samples corresponding to the low and high-expression groups
low_samples <- names(gene_array[,gene_array < low_cutoff])
high_samples <- names(gene_array[,gene_array > high_cutoff])

# complete the annotation dataframe so we can have for later
# in case we need it:
ann_df[low_samples, 'expression_state'] <- 'low'
ann_df[high_samples, 'expression_state'] <- 'high'

# the "un-interesting" samples (those with NAs) can be dropped
ann_df <- na.omit(ann_df)

# we will need to have the "expression_state" as a factor
# for use with DESeq2
ann_df$expression_state <- factor(ann_df$expression_state, levels=c('low', 'high'))

# Subset the raw count matrix to only include the low
# and high-expressed samples:
counts = counts[,rownames(ann_df)]


dds <- DESeqDataSetFromMatrix(countData = counts,
							  colData = ann_df,
							  design = ~expression_state)
dds <- DESeq(dds)
res <- results(dds, name='expression_state_high_vs_low')
resOrdered <- res[order(res$padj),]

# Write all our outputs:

# map the "R-mangled" names back to the original:
rownames(ann_df) <- colname_mapping[rownames(ann_df)]
write.table(ann_df, ANN_OUTPUT, sep='\t', quote=F)

write.table(resOrdered, DGE_OUTPUT, sep='\t', quote=F)

# map back to the original column/sample names for the normalized counts
# Note that this normalized matrix is ALL samples:
colnames(norm_mtx)<-colname_mapping[colnames(norm_mtx)]
write.table(norm_mtx, NORM_COUNTS_OUTPUT, sep='\t', quote=F)