#!/usr/bin/env nextflow

include { DOWNLOAD } from './modules/download.nf'
include { FASTQ } from './modules/convert_fastq.nf'
include { FASTQC } from './modules/fastqc.nf'
include { MULTIQC } from './modules/multiqc.nf'
include { TRIMMING } from './modules/trimming.nf'
include { CONCAT } from './modules/concatenate.nf'
include { ALIGN } from './modules/align.nf'
include { COUNTING } from './modules/count_reads.nf'
include { COUNT_MATRIX } from './modules/count_matrix.nf'
include { DEA } from './modules/DEA_analysis.nf'

params.SRA_input='data/SRA_runs.csv'
params.sample_map='data/sample_map.csv'
params.genome_dir="/mnt/d/bulkRNAseq-nextflow/data/grch38"
params.genome_prefix = 'genome'
params.annotation_file="data/gencode.v48.primary_assembly.annotation.gtf"
params.count_matrix_script="scripts/generating_matrix.py"
params.count_dir="results/read_counts"
params.DEA_script="scripts/DEA_analysis.R"
params.DEA_annotation="data/GRCh38annotation.csv"
params.pathway_file="data/h.all.v2025.1.Hs.symbols.gmt"
params.do_trimming=true

workflow{
    read_ch=Channel.fromPath(params.SRA_input).splitCsv(header:false).map{ item -> item[0] }

    DOWNLOAD(read_ch)
    FASTQ(DOWNLOAD.out.sra_files)
    FASTQC(FASTQ.out.fastq_files)
    MULTIQC(FASTQC.out.fastqc_reports.collect())
    if (params.do_trimming){
        TRIMMING(FASTQ.out.fastq_files)
        trimmed_ch = TRIMMING.out.trimmed_reads.map { file ->
        def run = file.simpleName.replace('_pass_trimmed', '')
        tuple(run, file)
    }
    }else{
        trimmed_ch = FASTQ.out.fastq_files.map { file ->
        def run = file.simpleName.replace('_pass', '')
        tuple(run, file)
    }
    }    
    map_ch=Channel.fromPath(params.sample_map).splitCsv(header:true).map{ row -> tuple(row.run, row.sample) }

    group_ch=trimmed_ch.join(map_ch).map{run, file, sample -> tuple(sample, file) }.groupTuple()
    CONCAT(group_ch)
    ALIGN(CONCAT.out.concatenated_fastq, params.genome_dir, params.genome_prefix)
    COUNTING(ALIGN.out.bam, file(params.annotation_file))
    COUNT_MATRIX(file(params.count_matrix_script), COUNTING.out.count_files.collect())
    DEA(file(params.DEA_script),COUNT_MATRIX.out.count_matrix, file(params.DEA_annotation), file(params.pathway_file))
}