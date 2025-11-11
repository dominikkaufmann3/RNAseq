#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --cpus-per-task=1
#SBATCH --mem=1000M
#SBATCH --time=1:00:00
#SBATCH --partition=pibu_el8

module load FastQC/0.11.9-Java-11  # load the FastQC module

# Make sure the output directory exists
mkdir -p ./QC_reports

# Run FastQC using 1 thread, saving output in QC_reports
fastqc -t 1 reads_Lung/*.fastq.gz -o ./QC_reports
