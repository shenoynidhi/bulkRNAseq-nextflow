#!/usr/bin/env nextflow

process CONCAT{

    publishDir 'results/processed_data', mode: 'symlink'
    input:
    tuple val(sample), path(trimmed_reads)

    output:
    path "${sample}.fastq.gz", emit: concatenated_fastq

    script:
    """
    cat ${trimmed_reads.join(' ')} > ${sample}.fastq.gz
    """
}