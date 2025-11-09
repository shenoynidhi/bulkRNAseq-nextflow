#!/usr/bin/env nextflow

process FASTQC{
    
    container 'community.wave.seqera.io/library/fastqc:0.12.1--af7a5314d5015c29'
    publishDir "results/fastqc", mode: 'symlink'

    input:
    path fastq_files

    output:
    path "${fastq_files.simpleName}_fastqc.{zip,html}", emit: fastqc_reports

    script:
    """
    fastqc ${fastq_files}
    """

}