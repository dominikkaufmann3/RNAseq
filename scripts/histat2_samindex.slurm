#!/bin/bash
#SBATCH --job-name=samtools_index
#SBATCH --output=samtools_index_%j.out
#SBATCH --error=samtools_index_%j.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=4000
#SBATCH --time=01:00:00
#SBATCH --partition=pibu_el8

CONTAINER="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"
MAP_DIR="/data/users/dkaufmann3/RNAseq/map"

for BAM in ${MAP_DIR}/*.sorted.bam
do
    SAMPLE=$(basename ${BAM} .sorted.bam)

    apptainer exec --bind /data ${CONTAINER} samtools index ${BAM}
done

