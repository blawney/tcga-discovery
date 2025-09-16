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


process pull_and_process_gtf {

    tag "Download GTF from GDC and create helper files"
    publishDir "${output_dir}/etc", mode:"copy"
    container "ghcr.io/blawney/tcga-discovery/tcga-discovery-python"
    cpus 2
    memory '12 GB'

    output:
        path "all_genes.sorted.bed"
        path "gene_mapping.tsv"

    script:
        """
        curl ${params.gtf} -o gene_annotations.gtf.gz
        gunzip gene_annotations.gtf.gz

        /py_venv/bin/python3 /opt/software/scripts/process_gtf.py \
            -i gene_annotations.gtf \
            -b all_genes.bed \
            -m gene_mapping.tsv

        /usr/bin/sort -k1,1 -k2,2n all_genes.bed > all_genes.sorted.bed
        """
}

process download_rnaseq_cohort {

    tag "Download TCGA cohort"
    publishDir "${output_dir}/${tcga_type}/raw_data", mode:"copy"
    container "ghcr.io/blawney/tcga-discovery/tcga-discovery-bioc"
    cpus 2
    memory '8 GB'

    input:
        val tcga_type

    output:
        tuple val(tcga_type), path("${tcga_type}.raw_counts.tsv"), path("${tcga_type}.rnaseq_metadata.rds")

    script:
        """
        /opt/conda/envs/r-env/bin/Rscript /opt/software/scripts/download_rnaseq_cohort.R "${tcga_type}"
        """
}


process run_dge {

    tag "Run differential expression on $raw_counts"
    publishDir "${output_dir}/${tcga_type}/normalized_counts", mode:"copy", pattern: "*.deseq2_norm_counts.all.tsv"
    publishDir "${output_dir}/${tcga_type}/dge_results", mode:"copy", pattern: "*.deseq2_results*"
    publishDir "${output_dir}/${tcga_type}/annotations", mode:"copy", pattern: "*.annotations*"
    container "ghcr.io/blawney/tcga-discovery/tcga-discovery-bioc"
    cpus 2
    memory '8 GB'

    input:
        tuple val(tcga_type), path(raw_counts)

    output:
        tuple val(tcga_type), \
            path("${tcga_type}.deseq2_results.${params.gene}_${params.low_q}_${params.high_q}.high_vs_low.tsv"), \
            path("${tcga_type}.deseq2_norm_counts.all.tsv"), \
            path("${tcga_type}.annotations.${params.gene}_${params.low_q}_${params.high_q}.tsv")

    script:
        """
        /opt/conda/envs/r-env/bin/Rscript /opt/software/scripts/deseq2.R \
            ${raw_counts} \
            ${params.gene} \
            ${params.low_q} \
            ${params.high_q} \
            ${tcga_type}.deseq2_results.${params.gene}_${params.low_q}_${params.high_q}.high_vs_low.tsv \
            ${tcga_type}.deseq2_norm_counts.all.tsv \
            ${tcga_type}.annotations.${params.gene}_${params.low_q}_${params.high_q}.tsv
        """
}


process run_gsea {

    tag "Run GSEA on tcga type: $tcga_type"
    publishDir "${output_dir}/${tcga_type}/gsea", mode:"copy"
    container "ghcr.io/blawney/tcga-discovery/tcga-discovery-gsea"
    cpus 2
    memory '8 GB'

    input:
        tuple val(tcga_type), \
            path(dge_results), \
            path(norm_counts), \
            path(annotations)

    output:
        tuple val(tcga_type), \
            path("${gct_file}")
            path("${cls_file}")
            path("${tcga_type}.${params.gene}_${params.low_q}_${params.high_q}.gsea_results.zip")

    script:
        def gct_file_template = "%s.%s_%s_%s.high_vs_low.gct"
        def cls_file_template = "%s.%s_%s_%s.high_vs_low.cls"
        gct_file = String.format(gct_file_template, tcga_type, params.gene, params.low_q, params.high_q)
        cls_file = String.format(cls_file_template, tcga_type, params.gene, params.low_q, params.high_q)
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
            -rpt_label "${tcga_type}" \
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
            -p "/gsea/${tcga_type}*/*.zip" \
            -o ${tcga_type}.${params.gene}_${params.low_q}_${params.high_q}.gsea_results.zip
        """
}


process run_gsea_preranked {

    tag "Run GSEA Preranked on tcga type: $tcga_type"
    publishDir "${output_dir}/${tcga_type}/gsea_preranked", mode:"copy"
    container "ghcr.io/blawney/tcga-discovery/tcga-discovery-gsea"
    cpus 2
    memory '8 GB'

    input:
        tuple val(tcga_type), \
            path(dge_results), \
            path(norm_counts), \
            path(annotations)
    output:
        path("${tcga_type}.${params.gene}_${params.low_q}_${params.high_q}.gsea_preranked_results.zip")
        path("${tcga_type}.${params.gene}_${params.low_q}_${params.high_q}.rnk")

    script:
        """
        /py_venv/bin/python3 /opt/software/scripts/create_rnk_file.py \
            -f ${dge_results} \
            -o ${tcga_type}.${params.gene}_${params.low_q}_${params.high_q}.rnk

        /opt/software/gsea/GSEA_4.3.2/gsea-cli.sh GSEAPreranked \
            -gmx ftp.broadinstitute.org://pub/gsea/msigdb/human/gene_sets/h.all.v2023.1.Hs.symbols.gmt \
            -chip ftp.broadinstitute.org://pub/gsea/msigdb/human/annotations/Human_Ensembl_Gene_ID_MSigDB.v2023.1.Hs.chip \
            -rnk ${tcga_type}.${params.gene}_${params.low_q}_${params.high_q}.rnk \
            -collapse Collapse \
            -mode Abs_max_of_probes \
            -norm meandiv \
            -nperm 1000 \
            -rnd_seed timestamp \
            -scoring_scheme weighted \
            -rpt_label "${tcga_type}.preranked" \
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
            -o ${tcga_type}.${params.gene}_${params.low_q}_${params.high_q}.gsea_preranked_results.zip
        """
}


process map_ensg_to_symbol {

    tag "Run ENSG to symbol gene mapping on $exp_mtx and $dge_results"
    publishDir "${output_dir}/${tcga_type}/dge_results", mode:"copy", pattern: "*.deseq2_results*"
    publishDir "${output_dir}/${tcga_type}/normalized_counts", mode:"copy", pattern: "*.deseq2_norm_counts.symbol_remapped.*"
    container "ghcr.io/blawney/tcga-discovery/tcga-discovery-python"
    cpus 1
    memory '4 GB'

    input:
        tuple val(tcga_type), \
            path(dge_results), \
            path(norm_counts), \
            path(annotations)
        path(gene_mapping)

    output:
        tuple val(tcga_type), \
            path("${tcga_type}.deseq2_norm_counts.symbol_remapped.all.tsv"), \
            path("${tcga_type}.deseq2_results.${params.gene}_${params.low_q}_${params.high_q}.high_vs_low.symbol_remapped.tsv")

    script:
        """
        /py_venv/bin/python3 /opt/software/scripts/map_ensg_to_symbol.py \
            -i ${norm_counts} \
            -m ${gene_mapping} \
            -o ${tcga_type}.deseq2_norm_counts.symbol_remapped.all.tsv

        /py_venv/bin/python3 /opt/software/scripts/map_ensg_to_symbol.py \
            -i ${dge_results} \
            -m ${gene_mapping} \
            -o ${tcga_type}.deseq2_results.${params.gene}_${params.low_q}_${params.high_q}.high_vs_low.symbol_remapped.tsv
        """
}


process overlay_ppi {

    tag "Overlay stringDB data"
    publishDir "${output_dir}/${tcga_type}/ppi_networks", mode:"copy", pattern: "*.gmt"
    publishDir "${output_dir}/${tcga_type}/ppi_networks", mode:"copy", pattern: "*.tsv"
    publishDir "${output_dir}/${tcga_type}/ppi_networks", mode:"copy", pattern: "*.pkl"
    container "ghcr.io/blawney/tcga-discovery/tcga-discovery-python"
    cpus 2
    memory '12 GB'

    input:
        tuple val(tcga_type), path(dge_results)

    output:
        tuple val(tcga_type), \
            path("${tcga_type}.*pkl", arity:2, glob:true), \
            path("${tcga_type}.*gmt", arity:2, glob:true), \
            path("${tcga_type}.*tsv", arity:2, glob:true)

    script:
        """
        curl https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz -o links.txt.gz
        curl https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz -o info.txt.gz
        gunzip links.txt.gz
        gunzip info.txt.gz

        /py_venv/bin/python3 /opt/software/scripts/overlay_ppi.py --interaction_matrix links.txt \
                                            --tcga_type ${tcga_type} \
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

    publishDir "${output_dir}/${tcga_type}/ppi_networks", mode:"copy", pattern: "*.h5"


    input:
        path gmt_file

    output:
        path "${file_base}.h5"

    script:
        file_base = gmt_file.baseName
        tcga_type = file_base.tokenize(".")[0]
        """
        touch "${file_base}.h5"
        /opt/conda/envs/r-env/bin/Rscript /opt/software/scripts/run_go.R \
            ${gmt_file} \
            "${file_base}.h5"
        """
}


process download_methylation_data {
    container "ghcr.io/blawney/tcga-discovery/tcga-discovery-bioc"
    cpus 2
    memory '12 GB'

    publishDir "${output_dir}/${tcga_type}/methylation", mode:"copy"

    input:
        val tcga_type

    output:
        tuple val(tcga_type), path("${tcga_type}.methylation_data.rds")

    script:
        """
        /opt/conda/envs/r-env/bin/Rscript /opt/software/scripts/download_full_methylation.R "${tcga_type}"
        """
}


process prep_methylation_data {
    container "ghcr.io/blawney/tcga-discovery/tcga-discovery-bioc"
    cpus 2
    memory '16 GB'

    publishDir "${output_dir}/${tcga_type}/methylation", mode:"copy"

    input:
        tuple val(tcga_type), \
            path(norm_counts), \
            path(meth_rds), \
            path(probe_file)

    output:
        tuple val(output_key), \
            path("${tcga_type}.norm_counts.reformatted_sample_id.tsv"), \
            path("${tcga_type}.betas.${probe_class}.reformatted_sample_id.tsv")

    script:
        probe_class = probe_file.baseName.tokenize(".")[0]
        output_key = tcga_type + "+" + probe_class
        """
        /opt/conda/envs/r-env/bin/Rscript /opt/software/scripts/prep_methylation_data.R ${tcga_type} ${meth_rds} ${norm_counts} ${probe_file}
        """
}


process get_differential_probes {
    container "ghcr.io/blawney/tcga-discovery/tcga-discovery-python"
    cpus 2
    memory '16 GB'

    publishDir "${output_dir}/${tcga_type}/methylation", mode:"copy", pattern: "*.supervised.*.tsv"

    input:
        tuple val(tcga_type), \
            val(tcga_and_probe_class_key), \
            path(betas), \
            path(annotations), \
            val(style)

    output:
        
        tuple val(tcga_and_probe_class_key), path("${tcga_type}.${probe_class}.supervised.${style}_na${params.na_fraction}.tsv")
    
    script:
        probe_class = tcga_and_probe_class_key.tokenize('+')[1]
        //out_key = tcga_and_probe_class_key + "+" style
        """
        /py_venv/bin/python3 /opt/software/scripts/supervised_test_for_diff_meth.py \
            -m ${betas} \
            -a ${annotations} \
            --style ${style} \
            --na_threshold ${params.na_fraction}
        """
}

process get_probe_neighbors {
    container "docker.io/biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1"
    cpus 2
    memory '8 GB'

    publishDir "${output_dir}/etc", mode:"copy"

    input:
        tuple path(probe_file), val(num_neighbors)
        path(all_genes_bed)

    output:
        tuple val(probe_class), path("${probe_class}.closest_genes.${num_neighbors}.sorted.bed")
    
    script:
        probe_class = probe_file.baseName.tokenize(".")[0]
        """
        bedtools closest -k ${num_neighbors} -d -D ref -io -iu -a ${probe_file} -b ${all_genes_bed} > downstream.bed`
        bedtools closest -k ${num_neighbors} -d -D ref -io -id -a ${probe_file} -b ${all_genes_bed} > upstream.bed`
        cat downstream.bed upstream.bed | awk -F"\t" '\$10 != "."' | /usr/bin/sort -k1,1 -k2,2n > ${probe_class}.closest_genes.${num_neighbors}.sorted.bed
        """
}


process test_probes_vs_genes {
    container "ghcr.io/blawney/tcga-discovery/tcga-discovery-python"
    cpus 2
    memory '8 GB'

    publishDir "${output_dir}/${tcga_type}/methylation", mode:"copy", pattern: "*tsv"

    input:
        tuple val(probe_class), val(tcga_type), path(normalized_counts), path(betas), path(diff_probe_test_results), path(annotations), path(gene_neighbors)

    output:
        path("*.probe_to_gene_results.tsv")

    script:
        """
        /py_venv/bin/python3 /opt/software/scripts/test_for_meth_regulation.py \
            -m ${diff_probe_test_results} \
            -b ${betas} \
            -n ${normalized_counts} \
            -t ${gene_neighbors} \
            -a ${annotations}
        """
}

workflow {

    tcga_types_ch = Channel.fromList(params.tcga_types)

    (all_genes_bed_ch, gene_mapping_ch) = pull_and_process_gtf()

    // has a tuple of tcga type, raw counts, and a clinical metadata RDS file:
    rnaseq_ch = download_rnaseq_cohort(tcga_types_ch)

    // don't need the metadata RDS in any downstream processes. We just have it in case
    // we do some custom processing after
    counts_ch = rnaseq_ch.map {
        it -> it[0..1]
    }

    // has a tuple of tcga type, dge results, norm counts, and high/low annotations
    dge_results_ch = run_dge(counts_ch)

    run_gsea(dge_results_ch)
    run_gsea_preranked(dge_results_ch)

    // has a tuple of tcga type, remapped norm counts, and remapped dge results
    symbol_remapped_ch = map_ensg_to_symbol(dge_results_ch, gene_mapping_ch)

    ppi_input_ch = symbol_remapped_ch.map {
        it -> [it[0], it[2]] // the tcga type and dge results only
    }

    // this channel has a tuple of tcga type followed by 3 lists 
    // corresponding to:
    // - pkl files of the graphs (networkx)
    // - GMT-format files of the communities
    // - TSV-format files giving the network stats
    // Each of those lists has length 2
    ppi_results_ch = overlay_ppi(ppi_input_ch)

    // This flatMap operation ultimately creates a single channel with
    // 2*n GMT files, where n is the number of TCGA types we are running.
    // Those all get sent to the same GO-process in parallel. File naming
    // is picked up by the names of the input files, so we don't need to
    // pass around the TCGA type.
    community_gmt_ch = ppi_results_ch.flatMap {
        it -> it[2]
    }
    run_go_on_communities(community_gmt_ch)

    // Start working on methylation data- 
    probe_file_ch = Channel.fromPath(params.probe_files)

    // a tuple of TCGA type to an RDS file with methylation data
    meth_rds_ch = download_methylation_data(tcga_types_ch)

    // join the normalized counts to the methylation RDS so we can
    // "align" the corresponding samples
    norm_count_ch = dge_results_ch.map {
        it -> [it[0], it[2]] // the tcga type and norm counts (with ENSG IDs) only
    }
    prep_meth_input_ch = norm_count_ch.join(meth_rds_ch)

    // This creates a cartesian product so we are able to prepare different beta
    // matrices corresponding to different probes (e.g. one for promoters, another for enhancers)
    prep_meth_input_ch = prep_meth_input_ch.combine(probe_file_ch)

    // a tuple of a unique key, subset normalized counts, and subset beta matrix. 
    // The unique key is a combination of the TCGA type and the probe style
    // (e.g. "TCGA-OV+enhancers"). Since the beta matrices are tied to both the cancer
    // type *and* probe style, we need both to keep all the files aligned properly
    // when merging/joining channels downstream. Note that this process creates a normalized count
    // matrix which contains only samples that have matching methylation data; it does not
    // subset that matrix in any other way
    prepped_meth_ch = prep_methylation_data(prep_meth_input_ch)

    // need to merge the betas with the 'high/low' annotations from the differential expression process
    diff_probe_input_ch = prepped_meth_ch.map {
        it -> 
            // to join with the annotations from above, we need only the TCGA type
            def tcga_type = it[0].tokenize("+")[0]

            // return a tuple of the tcga type, the tcga + probe class key, 
            // and the beta matrix. We need the lone tcga_type to permit
            // merging with channels created during the differential expression process
            [tcga_type, it[0], it[2]]
    }.combine(dge_results_ch.map{
        it -> [it[0], it[3]] // the tcga type and the high/low annotation file
    }, by:0)

    // also add in a string indicating whether we are testing for hypo/hypermethylation
    diff_probe_input_ch = diff_probe_input_ch.combine(Channel.fromList(["hypo", "hyper"]))

    // a tuple of the tcga type to differential probe table. There will be one for each 
    // probe style as well as hypo/hyper
    diff_meth_files_ch = get_differential_probes(diff_probe_input_ch)

    // This channel has size of: (# of TCGA typles) * (# probe files) * 2 (for hypo/hyper)
    meth_reg_test_input_ch = prepped_meth_ch.combine(diff_meth_files_ch, by: 0).map {
        it -> 
            // to join with the annotations from above, we need only the TCGA type
            def tcga_type = it[0].tokenize("+")[0]
            [tcga_type] + it
    }.combine(dge_results_ch.map{
        it -> [it[0], it[3]] // the tcga type and the high/low annotation file
    }, by: 0)

    // For each probe style, we can define a number of neighbors (e.g. 10 up/downstream).
    // The output of this process is a tuple of the probe class and the neighbor file
    num_gene_neighbors_ch = Channel.fromList(params.number_gene_neighbors)
    probe_neighbors_input_ch = probe_file_ch.merge(num_gene_neighbors_ch)
    neighbors_file_ch = get_probe_neighbors(probe_neighbors_input_ch, all_genes_bed_ch)

    // Finally, assign the neighbors file to the `meth_reg_test_input_ch`. Note that we merge by
    // the 'probe style' - for both the hypo and hypermethylated files, we still have the same
    // neighbors, so the `combine` operator with the `by` arg gives us this.
    meth_reg_test_input_ch = meth_reg_test_input_ch.map {
        it ->
        def probe_style = it[1].tokenize("+")[1]
        [probe_style] + [it[0]] + it[2..-1]
    }.combine(neighbors_file_ch, by:0)
    

    test_probes_vs_genes(meth_reg_test_input_ch)

}