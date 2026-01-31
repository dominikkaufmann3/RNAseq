#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --cpus-per-task=1
#SBATCH --mem=1000M
#SBATCH --time=2:00:00
#SBATCH --partition=pibu_el8
#SBATCH --output=logs/slurm/fastqc_%j.out
#SBATCH --error=logs/slurm/fastqc_%j.err

# Make sure the output directory exists
mkdir -p ./QC/fastqc

# Run FastQC inside the Apptainer container
apptainer exec /containers/apptainer/fastqc-0.12.1.sif \
    fastqc -t 1 QC/reads_Lung/*.fastq.gz -o ./QC/fastqc
