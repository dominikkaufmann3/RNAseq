#!/bin/bash
#SBATCH --job-name=featureCounts
#SBATCH --output=featureCounts_%j.out
#SBATCH --error=featureCounts_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=1000
#SBATCH --time=02:00:00
#SBATCH --partition=pibu_el8

CONTAINER="/containers/apptainer/subread_2.0.6.sif"
MAP_DIR="/data/users/dkaufmann3/RNAseq/map"
ANNOTATION="/data/users/dkaufmann3/RNAseq/map/Mus_musculus.GRCm39.115.gtf"
OUT_DIR="/data/users/dkaufmann3/RNAseq/counts"

mkdir -p ${OUT_DIR}

# Count reads per gene for all sorted BAMs
apptainer exec --bind /data ${CONTAINER} featureCounts \
    -T ${SLURM_CPUS_PER_TASK} \
    -p \
    -s 2 \
    -a ${ANNOTATION} \
    -o ${OUT_DIR}/gene_counts.txt \
    ${MAP_DIR}/*.sorted.bam
