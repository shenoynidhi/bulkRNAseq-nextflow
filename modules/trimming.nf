#!/usr/bin/env nextflow

process TRIMMING {
    
    conda 'bioconda::trimmomatic'
    publishDir "results/trim", mode: "symlink"

    input:
    path SRA_files 

    output:
    path "${SRA_files.simpleName}_trimmed.fastq.gz", emit: trimmed_reads


    script:
    """
    trimmomatic SE -threads 4 ${SRA_files} ${SRA_files.simpleName}_trimmed.fastq.gz TRAILING:10 -phred33
    """
}