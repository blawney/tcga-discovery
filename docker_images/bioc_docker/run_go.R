library(org.Hs.eg.db)
library(clusterProfiler)
library(rhdf5)

args <- commandArgs(trailingOnly = TRUE)

gmt_file <- args[1]
output_hdf5 <- args[2]

gmt_lines <- readLines(gmt_file)
gene_sets <- lapply(gmt_lines, function(line) {
  parts <- strsplit(line, "\t")[[1]]
  genes <- parts[-c(1,2)]  # Skip set name and description
  genes
})
names(gene_sets) <- sapply(gmt_lines, function(line) strsplit(line, "\t")[[1]][1])

ontologies = c('BP','CC','MF')
sapply(ontologies, function(ont) {h5createGroup(output_hdf5, ont)})

lapply(names(gene_sets), function(gs_name){
    gs <- gene_sets[[gs_name]]
    for(ont in ontologies){
        ego <- enrichGO(gene          = gs,
                        keyType       = "SYMBOL",
                        OrgDb         = org.Hs.eg.db,
                        ont           = ont,
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.01,
                        qvalueCutoff  = 0.05)
        df <- as.data.frame(ego)
        if(dim(df)[1] > 0){
            h5write(as.data.frame(ego), output_hdf5, paste0(ont, '/', gs_name))
        }
    }
})


