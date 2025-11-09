#!/usr/bin/env nextflow

process FASTQ{
    
    conda 'bioconda::sra-tools'
    publishDir "results/", mode: 'symlink'

    input:
    path SRA_files

    output:
    path "fastq/${SRA_files.simpleName}_pass.fastq.gz", emit: fastq_files
    script:
    """
    fastq-dump --outdir fastq --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip ${SRA_files}
    """
}