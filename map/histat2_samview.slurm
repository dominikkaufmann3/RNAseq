#!/bin/bash
#SBATCH --job-name=samtools_view
#SBATCH --output=samtools_view_%j.out
#SBATCH --error=samtools_view_%j.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=4000
#SBATCH --time=02:00:00
#SBATCH --partition=pibu_el8

CONTAINER="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"
MAP_DIR="/data/users/dkaufmann3/RNAseq/map"

for SAM in ${MAP_DIR}/*.sam
do
    SAMPLE=$(basename ${SAM} .sam)

    apptainer exec --bind /data ${CONTAINER} samtools view -hbS ${SAM} > ${MAP_DIR}/${SAMPLE}.bam

    

done

