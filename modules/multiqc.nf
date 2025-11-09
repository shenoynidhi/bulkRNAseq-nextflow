#!/usr/bin/env nextflow 

process MULTIQC{
    
    container 'community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c'
    publishDir "results/multiqc", mode:'symlink'

    input:
    path fastqc_reports

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data", emit: data_dir

    script:
    """
    multiqc ${fastqc_reports} -o .
    """
}