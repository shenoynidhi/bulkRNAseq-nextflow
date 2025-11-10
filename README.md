### Bulk RNAseq Analysis Pipeline with Nextflow automation

This repository contains a fully automated, reproducible Nextflow workflow for analyzing bulk RNA-Seq data, from raw reads to differential expression analysis and visualization.

**Note:** It builds upon my previously built Bulk RNA-Seq pipeline and integrates automation, reproducibility, and modular execution through Nextflow, Docker and Conda.

---
## Overview

This pipeline performs the complete RNA-Seq analysis workflow, including:
- Dataset download (from GEO)
- Quality control (FastQC and MultiQC)
- Read trimming (Trimmomatic)
- Alignment (HISAT2)
- Quantification (featureCounts)
- Differential expression analysis (DEA) (DESeq2)
- Gene set enrichment analysis (Fgsea, Reactome)
- Visualization of DEA results (volcano plots, PCA, heatmaps, and more)

---
## Main Repository Structure

- `data/` -> input data for the pipeline
- `envs/` -> conda environment for DEA analysis (rnaseq_deseq2.yaml)
- `modules/` -> .nf files required for each step of the processing and analysis
- `results/DEA` -> output files from the DEA analysis (output from other steps have been excluded due to space issues)
- `scripts` -> scripts for generation of count matrix and DEA analysis used in the respective module in the modules/ dir
- `README.md` -> readme file
- `nextflow.config` -> nextflow configuration file
- `rnaseq.nf` -> main nextflow file that is executed

---
## Environment Setup before execution
Depending on the process, Nextflow automatically chooses Docker or Conda environment for the download of tools/packages. Hence setting up both is a prerequisite

1. **Conda Setup**<br>
Some modules (`download`, `convert_fastq`, `trimming`, `align`, `count_matrix` ) run with Conda. [Refer documention for installation](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)<br>
To recreate the environment for DEA module, do as follows:
```
conda env create -f envs/rnaseq_deseq2.yaml
conda activate rnaseq_deseq2
```
2. **Docker Setup**<br>
Other modules (`fastqc`, `multiqc`, `count_reads`) are executed inside Docker containers. For installation, can refer- [Docker Installation Guide](https://dev.to/abhay_yt_52a8e72b213be229/how-to-install-docker-on-windows-macos-and-linux-a-step-by-step-guide-3a2i)

**Note:** 
- Ensure docker.enabled and conda.enabled is set to true in the `nextflow.config` file
---
## Reference Files Required
Before running the workflow, download the following into the `data/` directory:
- Prebuilt human genome index for alignment
```
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
tar -xvzf grch38_genome.tar.gz
```

- GTF file for read count estimation 
```
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.primary_assembly.annotation.gtf.gz
gunzip gencode.v48.primary_assembly.annotation.gtf.gz
```
**Note:** Other required input files have already been included in the `data/` directory

---
## Running the Workflow
To execute the pipeline, run the following in the terminal:
```
nextflow run rnaseq.nf
```
**Note:** Add `- resume` to resume after interruption 

---
## Output Summary
After successful execution, outputs get stored in the corresponding subdirectories of `results/` directory:

- `results/prefetch` -> SRA runs fetched from GEO
- `results/fastq` -> FASTQ files generated from respective SRA runs
- `results/fastqc` -> FASTQC reports
- `results/multiqc` -> MULTIQC report
- `results/trim` -> Trimmed fastq files
- `results/processed_data` -> Concatenated sample FASTQ files
- `results/align` -> Aligned Bam files
- `results/read_counts` -> Read count files for each sample and merged count matrix
- `results/DEA` -> Results from DEA analysis including filtered and normalized count matrices, DE genes for LNCAP and PC3 samples, sample variability visualizations, GSEA figures, and heatmaps

---
## Bonus
This workflow can be easily customized to different datasets and experimental setups.

- **Input samples:**
    1. `data/SRA_runs.csv` **(contains SRA run IDs for download)** -> Edit this file to include your own SRA run IDs

    2. `data/sample_map.csv`**(maps individual SRA runs to biological samples or conditions)** -> Update this file to define your sample types or experimental groups (e.g., “Control”, “Treatment”). Also modify `condition` in the `DEA_analysis.R` script

- **Trimming parameter:**<br>
By default, trimming is performed. To not perform trimming, set as follows in `rnaseq.nf`:
```
params.do_trimming= false
```

- **DEA outputs and visualization:**
    - The R script used for DEA and visualization `scripts/DEA_analysis.R` can be easily modified to fit your analysis goals
    - Change thresholds for log₂ fold change or adjusted p-values to redefine significant genes
    - Add new visualizations (e.g., MA plots, clustering heatmaps, pathway enrichment charts)
    - These modifications require only basic R editing and do not affect the Nextflow workflow itself. However, if additional analyses or visualizations are added, the corresponding R packages must be installed in the Conda environment.

---
**Visit my previous manual Bulk RNASeq Analysis Pipeline here**:[Bulk-RNASeq-Analysis](https://github.com/shenoynidhi/Bulk-RNASeq-Analysis)