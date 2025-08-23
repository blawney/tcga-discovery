# TCGA transcriptome discovery and visualization

This repository contains a Nextflow script for processing RNA-seq data obtained from the The Cancer Genome Atlas (TCGA).

The original intent was to study how differing MUC1 levels might affect the transcriptome landscape, but the process is completely general and may be used to study your target of interest.

After selecting a candidate gene and cancer cohort (among several other parameters), the pipeline will stratify the cohort based on the expression of the target gene. It will then perform standard differential expression and pathway analyses to see how relatively low- and high-expression of the target gene affects the transcriptome.

Following the standard analyses, we use protein-protein interaction from stringDB (https://string-db.org/) to overlay high-confidence interactions on the differentially expressed transcripts. These networks are then clustered to discovery potentially related gene modules which can then be visualized with a companion dashboard application. 

## To run:

**Prerequisites:**
- Nextflow is installed (https://www.nextflow.io/docs/latest/getstarted.html)
- Docker is installed (software is containerized in Docker images)

Next we set the pipeline parameters. Copy `params.json.tmpl` and fill in variables as required. For a detailed explanation, see [parameters](#input-parameters)

Now you can run with:
```
nextflow run main.nf -c <NEXTFLOW CONFIG> -params-file <PARAMS JSON>
```

#### Output files

An example of the output folder:
```
results/250822.172517
└── TCGA-CHOL
    ├── annotations
    │   └── TCGA-CHOL.annotations.ENSG00000185499_0.3_0.7.tsv
    ├── dge_results
    │   ├── TCGA-CHOL.deseq2_results.ENSG00000185499_0.3_0.7.high_vs_low.symbol_remapped.tsv
    │   └── TCGA-CHOL.deseq2_results.ENSG00000185499_0.3_0.7.high_vs_low.tsv
    ├── gsea
    │   ├── TCGA-CHOL.ENSG00000185499_0.3_0.7.gsea_results.zip
    │   ├── TCGA-CHOL.ENSG00000185499_0.3_0.7.high_vs_low.cls
    │   └── TCGA-CHOL.ENSG00000185499_0.3_0.7.high_vs_low.gct
    ├── gsea_preranked
    │   └── TCGA-CHOL.ENSG00000185499_0.3_0.7.gsea_preranked_results.zip
    ├── normalized_counts
    │   ├── TCGA-CHOL.deseq2_norm_counts.all.tsv
    │   └── TCGA-CHOL.deseq2_norm_counts.symbol_remapped.all.tsv
    ├── ppi_networks
    │   ├── dge_only.pkl
    │   ├── full.pkl
    │   ├── louvain_communities.dge_only.gmt
    │   ├── louvain_communities.dge_only.h5
    │   ├── louvain_communities.full.gmt
    │   ├── louvain_communities.full.h5
    │   ├── network_stats.dge_only.tsv
    │   └── network_stats.full.tsv
    └── raw_data
        ├── TCGA-CHOL.metadata.tsv
        └── TCGA-CHOL.raw_counts.tsv
```
**Output file descriptions:**

This has annotations about which samples were in the lowest and highest quantiles (as set by your input paramters `low_q`, `high_q`):
```
    ├── annotations
    │   └── TCGA-CHOL.annotations.ENSG00000185499_0.3_0.7.tsv
```

These are both tables of differentially expressed genes. One is expressed in ENSG/Ensembl IDs, and the other has been mapped to gene symbols:
```
    ├── dge_results
    │   ├── TCGA-CHOL.deseq2_results.ENSG00000185499_0.3_0.7.high_vs_low.symbol_remapped.tsv
    │   └── TCGA-CHOL.deseq2_results.ENSG00000185499_0.3_0.7.high_vs_low.tsv
```

This has the traditional gene-set enrichment analysis. The zip file contains all the usual HTML files and associated figures. The `cls` and `gct` are the input files used to run GSEA:
```
    ├── gsea
    │   ├── TCGA-CHOL.ENSG00000185499_0.3_0.7.gsea_results.zip
    │   ├── TCGA-CHOL.ENSG00000185499_0.3_0.7.high_vs_low.cls
    │   └── TCGA-CHOL.ENSG00000185499_0.3_0.7.high_vs_low.gct
```
This has the results from the pre-ranked GSEA analysis. The `rnk` file has the ranking of the genes, which is determined by a combination of adjusted p-value and the direction of the log-fold change:
```
    ├── gsea_preranked
    │   ├── TCGA-CHOL.ENSG00000185499_0.3_0.7.gsea_preranked_results.zip
    │   └── TCGA-CHOL.ENSG00000185499_0.3_0.7.gsea_preranked_results.rnk
```

This has normalized counts following DESeq2's median-based normalization. One version has the original ENSG/Ensembl IDs, the other is mapped to the gene symbol:
```
    ├── normalized_counts
    │   ├── TCGA-CHOL.deseq2_norm_counts.all.tsv
    │   └── TCGA-CHOL.deseq2_norm_counts.symbol_remapped.all.tsv
```

This has the results of the stringDB network + clustering analysis. There are two versions of each file. The files marked "dge_only" are networks created where both genes/proteins associated with each edge are differentially expressed. The "full" files indicate that only one of the two proteins associated with each network edge are differentially expressed. 

- The `pkl` files are pickled networkx graphs
- The `gmt` files are GMT-format (e.g. like msigDB pathway files) which give the genes for each community.
- The `h5` files are HDF5-format files which contain the results of gene ontology analyses on the communities. They contain groups for the three ontologies (BP, MF, CC) at the top level, followed by datasets for each of the communities. For instance, to access the result table from the BP ontology analysis on gene set `gs3`, you would look for the key `/BP/gs3`.
- The network stats files have various network measurements such as centrality and "betweeness"
```
    ├── ppi_networks
    │   ├── dge_only.pkl
    │   ├── full.pkl
    │   ├── louvain_communities.dge_only.gmt
    │   ├── louvain_communities.dge_only.h5
    │   ├── louvain_communities.full.gmt
    │   ├── louvain_communities.full.h5
    │   ├── network_stats.dge_only.tsv
    │   └── network_stats.full.tsv
```
This has the raw data downloaded from the GDC. The `metadata` file has clinical/demographic metadata while the raw counts is simply the integer gene quantifications.
```
    └── raw_data
        ├── TCGA-CHOL.metadata.tsv
        └── TCGA-CHOL.raw_counts.tsv
```

#### Input parameters

- `gene`: The Ensembl/ENSG identifier for your target of interest
- `tcga_type`: The TCGA identifier/code for the cancer of interest. For example, `"TCGA-OV"` for ovarian cancer. See the [GDC portal](https://portal.gdc.cancer.gov/) for a full list of the TCGA cohort identifiers.
- `low_q`, `high_q`: The low and high cutoffs, expressed as a fraction. For example, values of `low_q=0.25` and `high_q=0.75` will compare the lowest and highest quartiles of expression for your target gene of interest.
- `padj_threshold`: A float on [0,1.0] which sets the adjusted p-value for significance. We use this to determine genes that are differentially expressed prior to performing clustering on the interaction network.
- `ppi_score_threshold`: A number on [0, 1000] which sets the threshold for a stringDB protein-protein interaction to be deemed high-confidence. This number is a heuristic which reflects the weight of evidence for determination of protein-protein interaction. Numbers closer to 1000 are more stringent and will give only the most high-confidence associations. Accordingly, larger values will result in fewer edges in the network graph.
- `output_dir`: The name of the output directory where results will be placed. Each run of the pipeline will create a time-stamped subdirectory underneath this.

As an example (for MUC1/ENSG00000185499) in cholangiocarcinoma with cutoffs in the lowest and highest quartile:
```
{
    "gene":"ENSG00000185499",
    "low_q": 0.25,
    "high_q": 0.75,
    "output_dir": "results",
    "tcga_type": "TCGA-CHOL",
    "padj_threshold": 0.05,
    "ppi_score_threshold": 900
}
```

## Visualizing the communities

To help with browsing the results of the community analysis, this repository contains a companion visualization plotly/Dash tool which is distributed as a Docker image.

To run this, you need to locate the results of your analysis. Below, we assume that our current working directory contains the timestamped results directory (`results/250822.172517`). Then we run: 
```
docker run -d \
    -v $PWD/results:/results \
    -p 8050:8050 ghcr.io/blawney/tcga-discovery/tcga-discovery-dash \
    --results_dir /results/250822.172517 --network_type dge_only
```

- The `-v` mounts the `results/` directory in the container so that the application has access to the pipeline results. 
- `-p 8050:8050` maps port 8050 in the container to the same port on your host machine
- The `--results_dir` and `--network_type` arguments are passed to the Dash app and tell it which files to select for display. 
  - The `--network_type` is either "full" or "dge_only". Refer to the notes above for the meaning of these and what they imply for the resulting network graphs.

Once the container is started. you can open your browser to `http://127.0.0.1:8050/` to view the application: 

![](dash_app.png)

**Shutting down the visualization app**

When you are done with the application, you can simply run:
```
docker ps -a
```
which should print some information, namely the container ID. With that identifier, 
```
docker stop <CONTAINER ID>
docker rm <CONTAINER ID>
```