#!/usr/bin/env nextflow

process ALIGN{
    
    conda 'bioconda::hisat2 bioconda::samtools'
    publishDir "results/align", mode: "symlink"

    input:
    path reads
    path genome_dir
    val genome_prefix

    output:
    path "${reads.simpleName}.bam", emit: bam

    script:
    """
    hisat2 -q -x ${genome_dir}/${genome_prefix} -U ${reads} | samtools sort -o ${reads.simpleName}.bam
    samtools index ${reads.simpleName}.bam
    """
}