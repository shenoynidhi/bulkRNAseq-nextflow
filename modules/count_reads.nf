#!/usr/bin/env nextflow

process COUNTING{

    container 'community.wave.seqera.io/library/subread:2.1.1--0ac4d7e46cd0c5d7'
    publishDir "results/read_counts", mode:'symlink'

    input:
    path bam
    path annotation_file

    output:
    path "${bam.simpleName}_featurecounts.txt", emit: count_files

    script:
    """
    featureCounts -S 2 -a ${annotation_file} -o ${bam.simpleName}_featurecounts.txt ${bam}
    """
}