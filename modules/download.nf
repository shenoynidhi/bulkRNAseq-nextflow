#!/ur/bin/env nextflow

process DOWNLOAD{
    
    conda 'bioconda::sra-tools'
    publishDir "results/prefetch", mode: 'symlink'
    
    input:
    val SRA_files

    output:
    path "${SRA_files}/*.sra", emit: sra_files

    script:
    """
    prefetch ${SRA_files}
    """
}