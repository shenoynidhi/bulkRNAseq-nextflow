#!/usr/bin/env nextflow

process DEA{

    conda './envs/rnaseq_deseq2.yaml'

    publishDir "results/DEA", mode:'copy'

    input:
    val script
    path count_file
    path annotation_file
    path pathway_file
    
    output:
    path "filtered_count_matrix.csv", emit: filtered_counts
    path "normalized_count_matrix.csv", emit: normalized_counts
    path "DEGs_lncap.csv", emit: LNCAP_DEgenes
    path "DEGs_pc3.csv", emit: PC3_DEgenes
    path "*.png", optional: true
    path "*.csv", optional: true

    script:
    """
    Rscript ${script} ${count_file} ${annotation_file} filtered_count_matrix.csv normalized_count_matrix.csv ${pathway_file}
    """
}