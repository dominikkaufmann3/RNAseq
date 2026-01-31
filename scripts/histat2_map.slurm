#!/bin/bash
#SBATCH --job-name=hisat2_map
#SBATCH --output=logs/slurm/hisat2_map_%j.out
#SBATCH --error=logs/slurm/hisat2_map_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=8000
#SBATCH --time=04:00:00
#SBATCH --partition=pibu_el8

CONTAINER="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"

# UPDATED PATHS ONLY
READS_DIR="/data/users/dkaufmann3/RNAseq/QC/reads_Lung"
MAP_DIR="/data/users/dkaufmann3/RNAseq/map/bam"
INDEX_BASE="/data/users/dkaufmann3/RNAseq/data/reference/genome_index/genome"

mkdir -p ${MAP_DIR}

for R1 in ${READS_DIR}/*_1.fastq.gz
do
    SAMPLE=$(basename ${R1} _1.fastq.gz)
    R2=${READS_DIR}/${SAMPLE}_2.fastq.gz

    # Map reads to genome index
    apptainer exec --bind /data ${CONTAINER} hisat2 -x ${INDEX_BASE} \
        -1 ${R1} -2 ${R2} \
        --rna-strandness RF \
        -p ${SLURM_CPUS_PER_TASK} \
        -S ${MAP_DIR}/${SAMPLE}.sam
done