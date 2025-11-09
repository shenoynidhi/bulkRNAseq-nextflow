#!/usr/bin/env nextflow

process COUNT_MATRIX{

    conda 'anaconda::python anaconda::pandas'
    publishDir "results/read_counts", mode: 'symlink'

    input:
    val script
    path count_files

    output:
    path "merged_matrix.csv", emit: count_matrix

    script:
    """
    python3 ${script} ${count_files.join('')} merged_matrix.csv
    """
}