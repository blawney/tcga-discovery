nextflow.enable.dsl=2

import java.text.SimpleDateFormat

// Define an output directory based on a timestamp:
def date = new java.util.Date()
def sdf = new SimpleDateFormat("yyMMdd.HHmmss")
def timestamp =  sdf.format(date)
output_dir = params.output_dir + "/" + timestamp 

// remove genes with cohort-wise means lower than this threshold:
// (where cohort is the subset of samples AFTER subsetting for our
// high and low expressed gene of interest)
params.min_reads = 5


process download_cohort {

    tag "Download TCGA cohort"
    publishDir "${output_dir}/${params.tcga_type}/raw_data", mode:"copy"
    container "ghcr.io/blawney/tcga-discovery/tcga-discovery-bioc"
    cpus 2
    memory '8 GB'

    output:
        path "${params.tcga_type}.raw_counts.tsv"
        path "${params.tcga_type}.metadata.tsv"

    script:
        """
        /opt/conda/envs/r-env/bin/Rscript /opt/software/scripts/download_cohort.R ${params.tcga_type}
        """
}


process run_dge {

    tag "Run differential expression on $raw_counts"
    publishDir "${output_dir}/${params.tcga_type}/normalized_counts", mode:"copy", pattern: "*.deseq2_norm_counts.all.tsv"
    publishDir "${output_dir}/${params.tcga_type}/dge_results", mode:"copy", pattern: "*.deseq2_results*"
    publishDir "${output_dir}/${params.tcga_type}/annotations", mode:"copy", pattern: "*.annotations*"
    container "ghcr.io/blawney/tcga-discovery/tcga-discovery-bioc"
    cpus 2
    memory '8 GB'

    input:
        path(raw_counts)

    output:
        path("${params.tcga_type}.deseq2_results.${params.gene}_${params.low_q}_${params.high_q}.high_vs_low.tsv")
        path("${params.tcga_type}.deseq2_norm_counts.all.tsv")
        path("${params.tcga_type}.annotations.${params.gene}_${params.low_q}_${params.high_q}.tsv")

    script:
        """
        /opt/conda/envs/r-env/bin/Rscript /opt/software/scripts/deseq2.R \
            ${raw_counts} \
            ${params.gene} \
            ${params.low_q} \
            ${params.high_q} \
            ${params.tcga_type}.deseq2_results.${params.gene}_${params.low_q}_${params.high_q}.high_vs_low.tsv \
            ${params.tcga_type}.deseq2_norm_counts.all.tsv \
            ${params.tcga_type}.annotations.${params.gene}_${params.low_q}_${params.high_q}.tsv
        """
}


process run_gsea {

    tag "Run GSEA on tcga type: $params.tcga_type"
    publishDir "${output_dir}/${params.tcga_type}/gsea", mode:"copy"
    container "ghcr.io/blawney/tcga-discovery/tcga-discovery-gsea"
    cpus 2
    memory '8 GB'

    input:
        path(norm_counts)
        path(annotations)

    output:
        path("${gct_file}")
        path("${cls_file}")
        path("${params.tcga_type}.${params.gene}_${params.low_q}_${params.high_q}.gsea_results.zip")

    script:
        def gct_file_template = "%s.%s_%s_%s.high_vs_low.gct"
        def cls_file_template = "%s.%s_%s_%s.high_vs_low.cls"
        gct_file = String.format(gct_file_template, params.tcga_type, params.gene, params.low_q, params.high_q)
        cls_file = String.format(cls_file_template, params.tcga_type, params.gene, params.low_q, params.high_q)
        """
        /py_venv/bin/python3 /opt/software/scripts/prep_files.py \
            -f ${norm_counts} \
            -a ${annotations} \
            -g ${gct_file} \
            -c ${cls_file} \
            -t ${params.min_reads}

        /opt/software/gsea/GSEA_4.3.2/gsea-cli.sh GSEA \
            -res "${gct_file}" \
            -cls "${cls_file}#high_versus_low" \
            -gmx ftp.broadinstitute.org://pub/gsea/msigdb/human/gene_sets/h.all.v2023.1.Hs.symbols.gmt \
            -chip ftp.broadinstitute.org://pub/gsea/msigdb/human/annotations/Human_Ensembl_Gene_ID_MSigDB.v2023.1.Hs.chip \
            -out /gsea/ \
            -rpt_label "${params.tcga_type}" \
            -zip_report true \
            -collapse Collapse \
            -mode Max_probe \
            -norm meandiv \
            -nperm 1000 \
            -permute phenotype \
            -rnd_seed timestamp \
            -rnd_type no_balance \
            -scoring_scheme weighted \
            -metric Signal2Noise \
            -sort real \
            -order descending \
            -create_gcts false \
            -create_svgs false \
            -include_only_symbols true \
            -make_sets true \
            -median false \
            -num 100 \
            -plot_top_x 20 \
            -save_rnd_lists false \
            -set_max 500 \
            -set_min 15

        /py_venv/bin/python3 /opt/software/scripts/move_final_files.py \
            -p "/gsea/${params.tcga_type}*/*.zip" \
            -o ${params.tcga_type}.${params.gene}_${params.low_q}_${params.high_q}.gsea_results.zip
        """
}


process run_gsea_preranked {

    tag "Run GSEA Preranked on tcga type: $params.tcga_type"
    publishDir "${output_dir}/${params.tcga_type}/gsea_preranked", mode:"copy"
    container "ghcr.io/blawney/tcga-discovery/tcga-discovery-gsea"
    cpus 2
    memory '8 GB'

    input:
        path(dge_results)

    output:
        path("${params.tcga_type}.${params.gene}_${params.low_q}_${params.high_q}.gsea_preranked_results.zip")
        path("${params.tcga_type}.${params.gene}_${params.low_q}_${params.high_q}.rnk")

    script:
        """
        /py_venv/bin/python3 /opt/software/scripts/create_rnk_file.py \
            -f ${dge_results} \
            -o ${params.tcga_type}.${params.gene}_${params.low_q}_${params.high_q}.rnk

        /opt/software/gsea/GSEA_4.3.2/gsea-cli.sh GSEAPreranked \
            -gmx ftp.broadinstitute.org://pub/gsea/msigdb/human/gene_sets/h.all.v2023.1.Hs.symbols.gmt \
            -chip ftp.broadinstitute.org://pub/gsea/msigdb/human/annotations/Human_Ensembl_Gene_ID_MSigDB.v2023.1.Hs.chip \
            -rnk ${params.tcga_type}.${params.gene}_${params.low_q}_${params.high_q}.rnk \
            -collapse Collapse \
            -mode Abs_max_of_probes \
            -norm meandiv \
            -nperm 1000 \
            -rnd_seed timestamp \
            -scoring_scheme weighted \
            -rpt_label "${params.tcga_type}.preranked" \
            -create_svgs false \
            -include_only_symbols true \
            -make_sets true \
            -plot_top_x 20 \
            -set_max 500 \
            -set_min 15 \
            -zip_report true \
            -out /gsea_preranked/

        ls /gsea_preranked/*/*.zip
        echo "******************"
        /py_venv/bin/python3 /opt/software/scripts/move_final_files.py \
            -p "/gsea_preranked/*/*.zip" \
            -o ${params.tcga_type}.${params.gene}_${params.low_q}_${params.high_q}.gsea_preranked_results.zip
        """
}


process map_ensg_to_symbol {

    tag "Run ENSG to symbol gene mapping on $exp_mtx and $dge_results"
    publishDir "${output_dir}/${params.tcga_type}/dge_results", mode:"copy", pattern: "*.deseq2_results*"
    publishDir "${output_dir}/${params.tcga_type}/normalized_counts", mode:"copy", pattern: "*.deseq2_norm_counts.symbol_remapped.*"
    container "ghcr.io/blawney/tcga-discovery/tcga-discovery-python"
    cpus 1
    memory '4 GB'

    input:
        path(exp_mtx)
        path(dge_results)

    output:
        path("${params.tcga_type}.deseq2_norm_counts.symbol_remapped.all.tsv")
        path("${params.tcga_type}.deseq2_results.${params.gene}_${params.low_q}_${params.high_q}.high_vs_low.symbol_remapped.tsv")

    script:
        """
        /py_venv/bin/python3 /opt/software/scripts/map_ensg_to_symbol.py \
            -i ${exp_mtx} \
            -m /opt/software/resources/Human_Ensembl_Gene_ID_MSigDB.v2023.1.Hs.chip \
            -o ${params.tcga_type}.deseq2_norm_counts.symbol_remapped.all.tsv

        /py_venv/bin/python3 /opt/software/scripts/map_ensg_to_symbol.py \
            -i ${dge_results} \
            -m /opt/software/resources/Human_Ensembl_Gene_ID_MSigDB.v2023.1.Hs.chip \
            -o ${params.tcga_type}.deseq2_results.${params.gene}_${params.low_q}_${params.high_q}.high_vs_low.symbol_remapped.tsv
        """
}


process overlay_ppi {

    tag "Overlay stringDB data"
    publishDir "${output_dir}/${params.tcga_type}/ppi_networks", mode:"copy", pattern: "*.gmt"
    publishDir "${output_dir}/${params.tcga_type}/ppi_networks", mode:"copy", pattern: "*.tsv"
    publishDir "${output_dir}/${params.tcga_type}/ppi_networks", mode:"copy", pattern: "*.pkl"
    container "ghcr.io/blawney/tcga-discovery/tcga-discovery-python"
    cpus 2
    memory '12 GB'

    input:
        path dge_results

    output:
        path "*.pkl", emit: graph_pkl
        path "*.gmt", emit: louvain_gmt
        path "network_stats.*.tsv", emit: network_stats

    script:
        """
        curl https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz -o links.txt.gz
        curl https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz -o info.txt.gz
        gunzip links.txt.gz
        gunzip info.txt.gz

        /py_venv/bin/python3 /opt/software/scripts/overlay_ppi.py --interaction_matrix links.txt \
                                            --info_matrix info.txt \
                                            --dge_matrix ${dge_results} \
                                            --score_threshold ${params.ppi_score_threshold} \
                                            --padj_threshold ${params.padj_threshold}
        """
}


process run_go_on_communities {

    container "ghcr.io/blawney/tcga-discovery/tcga-discovery-bioc"
    cpus 2
    memory '12 GB'

    publishDir "${output_dir}/${params.tcga_type}/ppi_networks", mode:"copy", pattern: "*.h5"


    input:
        path gmt_file

    output:
        path "${file_base}.h5"

    script:
        file_base = gmt_file.baseName
        """
        touch "${file_base}.h5"
        /opt/conda/envs/r-env/bin/Rscript /opt/software/scripts/run_go.R \
            ${gmt_file} \
            "${file_base}.h5"
        """
}

workflow {
    (raw_count_ch, metadata_ch) = download_cohort()
    (dge_results_ch, norm_counts_ch, ann_ch) = run_dge(raw_count_ch)
    run_gsea(norm_counts_ch, ann_ch)
    run_gsea_preranked(dge_results_ch)
    (norm_counts_symbol_remapped_ch, dge_results_remapped_ch) = map_ensg_to_symbol(norm_counts_ch, dge_results_ch)
    results = overlay_ppi(dge_results_remapped_ch)
    f_ch = results.louvain_gmt.flatten()
    h5_ch = run_go_on_communities(f_ch)
}